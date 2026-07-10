#!/usr/bin/env Rscript
# ============================================================================
# MaAsLin3 v2 — Analysis 4: Disease-associated microbiota (ALL sites)
# For each disease: Disease vs Healthy Control, matched by bioproject_id.
# Fixed effects: Disease_Status + age_cleaned_group + sex_cleaned +
#                site_cleaned + sequencing_region_cleaned + continent_group
# Random effect: bioproject_id
# ============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3_v2")

MAASLIN3_CORES <- 1
Q_THRESHOLD <- 0.1

PREFERRED_DISEASES <- c(
  "Dental Caries", "Head and Neck Cancer", "HIV", "Oral Cancer", "OSCC",
  "Periodontal Disease", "RA", "Benign Tumor", "Breast Cancer", "CRC", "IBD", "Lung Cancer"
)

OUT_DIR <- "disease_results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. Load data
# ============================================================================
cat("=== Loading data ===\n")
genus_data <- read.table("genus_table_relabund.tsv", header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE,
                         comment.char = "", stringsAsFactors = FALSE)
metadata <- read.table("../metadata_oral_qc_dedup.tsv", header = TRUE, sep = "\t",
                       fill = TRUE, comment.char = "", stringsAsFactors = FALSE, quote = "")

# Column mappings
site_map <- c("saliva" = "Saliva", "tongue" = "Tongue", "tooth" = "Tooth",
              "oral swab" = "Oral_swab", "oral mucosa" = "Oral_Mucosa",
              "gingiva" = "Gingiva", "other" = "Others")
metadata$site_cleaned <- ifelse(metadata$site_cleaned %in% names(site_map),
                                site_map[metadata$site_cleaned], metadata$site_cleaned)
metadata$status_cleaned[metadata$status_cleaned %in% c("HC", "Healthy")] <- "Healthy Control"
metadata$age_cleaned_group <- gsub("＞", ">", metadata$age_cleaned_group)

# ============================================================================
# 2. Match samples — ALL sites
# ============================================================================
sample_ids <- trimws(as.character(metadata$sample_id))
keep_ids <- !is.na(sample_ids) & sample_ids != "" & !duplicated(sample_ids)
metadata <- metadata[keep_ids, , drop = FALSE]
sample_ids <- sample_ids[keep_ids]

common_samples <- intersect(sample_ids, colnames(genus_data))
cat(sprintf("Common samples: %d\n", length(common_samples)))

abund <- genus_data[, common_samples, drop = FALSE]
meta <- metadata[match(common_samples, sample_ids), , drop = FALSE]
rownames(meta) <- common_samples

key_cols <- c("age_cleaned_group", "sex_cleaned", "sequencing_region_cleaned",
              "bioproject_id", "site_cleaned", "continent_group", "status_cleaned")
for (col in key_cols) {
  meta[[col]] <- trimws(as.character(meta[[col]]))
  meta[[col]][is.na(meta[[col]]) | meta[[col]] == "" | meta[[col]] == "NA"] <- NA_character_
}

valid <- !is.na(meta$bioproject_id) & meta$bioproject_id != "" &
  !is.na(meta$status_cleaned) & meta$status_cleaned != "" &
  meta$status_cleaned != "Other"
meta <- meta[valid, , drop = FALSE]
abund <- abund[, rownames(meta), drop = FALSE]

disease_counts <- sort(table(meta$status_cleaned), decreasing = TRUE)
cat("\n=== Disease/Status counts ===\n")
print(disease_counts)
write.csv(as.data.frame(disease_counts), "disease_status_counts.csv", row.names = FALSE)

available <- PREFERRED_DISEASES[PREFERRED_DISEASES %in% names(disease_counts)]
cat(sprintf("\n%d / %d preferred diseases available\n", length(available), length(PREFERRED_DISEASES)))
missing <- setdiff(PREFERRED_DISEASES, names(disease_counts))
if (length(missing) > 0)
  cat(sprintf("Missing: %s\n", paste(missing, collapse = ", ")))

# ============================================================================
# 3. Run per-disease analysis
# ============================================================================
all_results <- list()
run_summary <- list()

for (disease in available) {
  disease_bp <- unique(meta$bioproject_id[meta$status_cleaned == disease])
  keep <- meta$status_cleaned %in% c("Healthy Control", disease) &
    meta$bioproject_id %in% disease_bp
  sub_meta <- meta[keep, , drop = FALSE]
  sub_abund <- abund[, rownames(sub_meta), drop = FALSE]

  bp_tbl <- table(sub_meta$bioproject_id, sub_meta$status_cleaned)
  if (!("Healthy Control" %in% colnames(bp_tbl)) || !(disease %in% colnames(bp_tbl))) {
    cat(sprintf("\nSKIP %s: no matched HC/Case BioProjects\n", disease)); next
  }
  valid_bp <- rownames(bp_tbl)[bp_tbl[, "Healthy Control"] > 0 & bp_tbl[, disease] > 0]
  if (length(valid_bp) == 0) {
    cat(sprintf("\nSKIP %s: no BioProject with both Case and HC\n", disease)); next
  }
  sub_meta <- sub_meta[sub_meta$bioproject_id %in% valid_bp, , drop = FALSE]
  sub_abund <- sub_abund[, rownames(sub_meta), drop = FALSE]

  n_case <- sum(sub_meta$status_cleaned == disease)
  n_ctrl <- sum(sub_meta$status_cleaned == "Healthy Control")
  n_bp   <- length(unique(sub_meta$bioproject_id))
  cat(sprintf("\n=== %s: %d case / %d control / %d matched BioProjects ===\n",
              disease, n_case, n_ctrl, n_bp))

  if (n_case < 30 || n_ctrl < 30 || n_bp < 2) {
    cat("  SKIP: need >=30 each and >=2 BioProjects\n"); next
  }

  sub_meta$Disease_Status <- factor(
    ifelse(sub_meta$status_cleaned == disease, disease, "Healthy Control"),
    levels = c("Healthy Control", disease))
  sub_meta$bioproject_id <- droplevels(factor(sub_meta$bioproject_id))

  optional_effects <- c("age_cleaned_group", "sex_cleaned", "site_cleaned",
                        "sequencing_region_cleaned", "continent_group")
  for (eff in optional_effects) {
    sub_meta[[eff]][is.na(sub_meta[[eff]]) | sub_meta[[eff]] == ""] <- "Unknown"
    sub_meta[[eff]] <- droplevels(factor(sub_meta[[eff]]))
  }

  fixed_effects <- c("Disease_Status")
  for (eff in optional_effects) {
    if (nlevels(sub_meta[[eff]]) >= 2) fixed_effects <- c(fixed_effects, eff)
  }
  cat(sprintf("  Fixed effects: %s\n", paste(fixed_effects, collapse = ", ")))

  disease_dir <- file.path(OUT_DIR, gsub(" ", "_", disease))
  dir.create(disease_dir, showWarnings = FALSE, recursive = TRUE)

  res_file <- file.path(disease_dir, "all_results.tsv")
  if (!file.exists(res_file)) {
    result <- tryCatch({
      maaslin3(
        input_data = t(sub_abund), input_metadata = sub_meta,
        output = disease_dir, fixed_effects = fixed_effects,
        random_effects = c("bioproject_id"),
        evaluate_only = NULL,
        min_abundance = 0, min_prevalence = 0.1, max_significance = 0.25,
        normalization = "TSS", transform = "LOG", correction = "BH",
        warn_prevalence = FALSE, cores = MAASLIN3_CORES, verbosity = "INFO"
      )
      read.delim(file.path(disease_dir, "all_results.tsv"), stringsAsFactors = FALSE)
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message)); NULL
    })
    if (is.null(result)) next
    result$Disease <- disease
  } else {
    cat("  Using existing results\n")
    result <- read.delim(res_file, stringsAsFactors = FALSE)
    result$Disease <- disease
  }

  sig_count <- sum(result$metadata == "Disease_Status" &
                   result$qval_individual < Q_THRESHOLD & !is.na(result$qval_individual))
  run_summary[[disease]] <- data.frame(
    Disease = disease, Cases = n_case, Controls = n_ctrl,
    Matched_BioProjects = n_bp,
    Fixed_Effects = paste(fixed_effects, collapse = ","),
    Significant_Disease_Associations = sig_count)
  all_results[[disease]] <- result
  cat(sprintf("  Done: %d disease-sig associations\n", sig_count))
  gc(verbose = FALSE)
}

if (length(all_results) == 0) {
  cat("No disease analyses completed.\n")
  quit(save = "no", status = 0)
}

combined <- bind_rows(all_results)
write.table(combined, "disease_results_all.tsv", sep = "\t", quote = TRUE, row.names = FALSE)
write.csv(bind_rows(run_summary), "disease_run_summary.csv", row.names = FALSE)

cat("\n=== Disease run summary ===\n")
print(bind_rows(run_summary))
cat("\n=== Disease analysis done ===\n")
