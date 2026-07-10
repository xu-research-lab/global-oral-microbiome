#!/usr/bin/env Rscript
# ============================================================================
# MaAsLin3 v2 — Site Pairwise Differential Analysis
# For each pair of (Saliva, Tooth, Tongue, Oral_swab, Oral_Mucosa):
#   MaAsLin3: site_cleaned + age_cleaned_group + sex_cleaned +
#             sequencing_region_cleaned + continent_group + (1|bioproject_id)
#   abundance model only (evaluate_only = "abundance")
# ============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3_v2/site_pairwise")

MAASLIN3_CORES <- 4
Q_THRESHOLD <- 0.05
MIN_SAMPLES_PER_SITE <- 50

# ============================================================================
# 1. Load data
# ============================================================================
cat("=== Loading data ===\n")
genus_data <- read.table("../genus_table_relabund.tsv", header = TRUE, sep = "\t",
                          row.names = 1, check.names = FALSE,
                          comment.char = "", stringsAsFactors = FALSE)
metadata <- read.table("../../metadata_oral_qc_dedup.tsv", header = TRUE, sep = "\t",
                        fill = TRUE, comment.char = "", stringsAsFactors = FALSE, quote = "")

# Apply mappings
site_map <- c("saliva" = "Saliva", "tongue" = "Tongue", "tooth" = "Tooth",
              "oral swab" = "Oral_swab", "oral mucosa" = "Oral_Mucosa",
              "gingiva" = "Gingiva", "other" = "Others")
metadata$site_cleaned <- ifelse(metadata$site_cleaned %in% names(site_map),
                                site_map[metadata$site_cleaned], metadata$site_cleaned)
metadata$age_cleaned_group <- gsub("＞", ">", metadata$age_cleaned_group)

# ============================================================================
# 2. Match samples
# ============================================================================
sample_ids <- trimws(as.character(metadata$sample_id))
keep_ids <- !is.na(sample_ids) & sample_ids != "" & !duplicated(sample_ids)
metadata <- metadata[keep_ids, , drop = FALSE]
sample_ids <- sample_ids[keep_ids]

common_samples <- intersect(sample_ids, colnames(genus_data))
meta <- metadata[match(common_samples, sample_ids), , drop = FALSE]
rownames(meta) <- common_samples
abund_all <- genus_data[, common_samples, drop = FALSE]

# Clean required columns
needed_cols <- c("age_cleaned_group", "sex_cleaned", "sequencing_region_cleaned",
                 "bioproject_id", "site_cleaned", "continent_group")
for (col in needed_cols) {
  meta[[col]] <- trimws(as.character(meta[[col]]))
  meta[[col]][is.na(meta[[col]]) | meta[[col]] == "" | meta[[col]] == "NA"] <- NA_character_
}

# Keep only complete cases
valid <- complete.cases(meta[, needed_cols])
meta <- meta[valid, , drop = FALSE]
abund_all <- abund_all[, rownames(meta), drop = FALSE]
cat(sprintf("Complete cases: %d samples\n", nrow(meta)))

# Factorize
meta$age_cleaned_group <- factor(meta$age_cleaned_group)
meta$sex_cleaned <- factor(meta$sex_cleaned)
meta$sequencing_region_cleaned <- factor(meta$sequencing_region_cleaned)
meta$continent_group <- factor(meta$continent_group)
meta$bioproject_id <- factor(meta$bioproject_id)

# Reference levels
if ("18-60 years" %in% levels(meta$age_cleaned_group))
  meta$age_cleaned_group <- relevel(meta$age_cleaned_group, ref = "18-60 years")
if ("F" %in% levels(meta$sex_cleaned))
  meta$sex_cleaned <- relevel(meta$sex_cleaned, ref = "F")
if ("V4" %in% levels(meta$sequencing_region_cleaned))
  meta$sequencing_region_cleaned <- relevel(meta$sequencing_region_cleaned, ref = "V4")

# ============================================================================
# 3. Define sites and all pairs
# ============================================================================
TARGET_SITES <- c("Saliva", "Tooth", "Tongue", "Oral_swab", "Oral_Mucosa")

# Filter to target sites only
meta <- meta[meta$site_cleaned %in% TARGET_SITES, , drop = FALSE]
abund_all <- abund_all[, rownames(meta), drop = FALSE]
meta$site_cleaned <- droplevels(factor(meta$site_cleaned))

site_counts <- table(meta$site_cleaned)
cat("\nSite distribution:\n")
print(site_counts[sapply(TARGET_SITES, function(s) s)])  # ordered

# Generate all pairs
site_pairs <- combn(TARGET_SITES, 2, simplify = FALSE)
n_pairs <- length(site_pairs)

cat(sprintf("\nTotal pairs to run: %d\n", n_pairs))
for (i in seq_along(site_pairs)) {
  cat(sprintf("  %d: %s vs %s\n", i, site_pairs[[i]][1], site_pairs[[i]][2]))
}

# ============================================================================
# 4. Run pairwise analyses
# ============================================================================
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_along(site_pairs)) {
  s1 <- site_pairs[[i]][1]
  s2 <- site_pairs[[i]][2]
  pair_name <- paste0(gsub(" ", "_", s1), "_vs_", gsub(" ", "_", s2))
  out_dir <- file.path(results_dir, pair_name)

  cat(sprintf("\n==========================================\n"))
  cat(sprintf("[%d/%d] %s vs %s\n", i, n_pairs, s1, s2))
  cat(sprintf("==========================================\n"))

  # Skip if done
  done_file <- file.path(out_dir, "analysis.done")
  if (file.exists(done_file)) {
    cat("  Already done, skipping.\n")
    next
  }

  # Subset to this pair
  pair_meta <- meta[meta$site_cleaned %in% c(s1, s2), , drop = FALSE]
  pair_abund <- abund_all[, rownames(pair_meta), drop = FALSE]
  pair_meta$site_cleaned <- droplevels(factor(pair_meta$site_cleaned, levels = c(s1, s2)))
  pair_meta$bioproject_id <- droplevels(pair_meta$bioproject_id)

  n_s1 <- sum(pair_meta$site_cleaned == s1)
  n_s2 <- sum(pair_meta$site_cleaned == s2)
  cat(sprintf("  %s: %d, %s: %d, BioProjects: %d\n", s1, n_s1, s2, n_s2, nlevels(pair_meta$bioproject_id)))

  if (n_s1 < MIN_SAMPLES_PER_SITE || n_s2 < MIN_SAMPLES_PER_SITE) {
    cat("  SKIP: too few samples\n"); next
  }

  # Dynamic covariate filtering
  fixed_candidates <- c("age_cleaned_group", "sex_cleaned",
                        "sequencing_region_cleaned", "continent_group")
  valid_covars <- c()
  for (cov in fixed_candidates) {
    if (nlevels(pair_meta[[cov]]) >= 2) {
      valid_covars <- c(valid_covars, cov)
    } else {
      cat(sprintf("  WARNING: %s has only 1 level, dropping\n", cov))
    }
  }

  fixed_effects <- c("site_cleaned", valid_covars)
  n_bp <- nlevels(pair_meta$bioproject_id)

  if (n_bp >= 2) {
    random_effects <- c("bioproject_id")
  } else {
    random_effects <- NULL
    cat("  WARNING: only 1 BioProject, using fixed effects only\n")
  }

  cat(sprintf("  Fixed: %s\n", paste(fixed_effects, collapse=", ")))
  cat(sprintf("  Random: %s\n", ifelse(is.null(random_effects), "none", random_effects[[1]])))

  # Run MaAsLin3
  result <- tryCatch({
    maaslin3(
      input_data = t(pair_abund),
      input_metadata = pair_meta,
      output = out_dir,
      fixed_effects = fixed_effects,
      random_effects = random_effects,
      evaluate_only = "abundance",
      min_abundance = 0, min_prevalence = 0.1, max_significance = 0.25,
      normalization = "TSS", transform = "LOG", correction = "BH",
      warn_prevalence = FALSE,
      cores = MAASLIN3_CORES,
      verbosity = "INFO"
    )

    res_file <- file.path(out_dir, "all_results.tsv")
    res <- read.delim(res_file, stringsAsFactors = FALSE)
    res$Site1 <- s1
    res$Site2 <- s2
    res$Comparison <- pair_name

    # Count significant
    sig_rows <- res %>%
      filter(metadata == "site_cleaned", !is.na(qval_individual),
             qval_individual < Q_THRESHOLD)
    cat(sprintf("  Significant taxa (q<%.2f): %d\n", Q_THRESHOLD, nrow(sig_rows)))

    # Save summary
    summary <- data.frame(
      Site1 = s1, Site2 = s2, Comparison = pair_name,
      N_Significant = nrow(sig_rows),
      N_S1 = n_s1, N_S2 = n_s2,
      N_BioProjects = n_bp,
      Fixed_Effects = paste(fixed_effects, collapse=","),
      stringsAsFactors = FALSE
    )
    write.table(summary, file.path(out_dir, "summary.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    writeLines("done", done_file)
    cat("  SUCCESS!\n")
    TRUE
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    FALSE
  })

  if (!isTRUE(result)) next
  gc(verbose = FALSE)
}

cat("\n=== All site pairwise analyses complete ===\n")
