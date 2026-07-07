#!/usr/bin/env Rscript
# ============================================================================
# MaAsLin3 Site-stratified Analysis
# Runs one MaAsLin3 model per Site_Cleaned.
# Fixed effects: Age_Group, Sex, Region_custom, Geography
# Random effect: BioProject
# ============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3")
source("input_helpers.R")

MAASLIN3_CORES <- 4

cat("=== Loading data ===\n")
genus_data <- load_genus_relabund()
metadata <- load_standard_metadata(dedup = TRUE)

sample_ids <- trimws(as.character(metadata$SampleID))
keep_ids <- !is.na(sample_ids) & sample_ids != "" & !duplicated(sample_ids)
metadata <- metadata[keep_ids, , drop = FALSE]
sample_ids <- sample_ids[keep_ids]

common_samples <- intersect(sample_ids, colnames(genus_data))
cat(sprintf("Common samples: %d\n", length(common_samples)))

abund <- genus_data[, common_samples, drop = FALSE]
meta <- metadata[match(common_samples, sample_ids), , drop = FALSE]
rownames(meta) <- common_samples

for (col in c("Age_Group", "Sex", "Region_custom", "BioProject", "Site_Cleaned", "Geography")) {
  if (col %in% colnames(meta)) meta[[col]] <- trimws(as.character(meta[[col]]))
}

valid <- complete.cases(meta[, c("Age_Group", "Sex", "Region_custom", "BioProject", "Site_Cleaned", "Geography")]) &
  meta$Age_Group != "" & meta$Sex != "" & meta$Region_custom != "" &
  meta$BioProject != "" & meta$Site_Cleaned != "" & meta$Geography != "" &
  meta$Age_Group != "NA" & meta$Sex != "NA" & meta$Geography != "Unknown"
meta <- meta[valid, , drop = FALSE]
abund <- abund[, rownames(meta), drop = FALSE]
cat(sprintf("Samples with complete metadata: %d\n", nrow(meta)))

meta$Age_Group <- factor(meta$Age_Group)
meta$Sex <- factor(meta$Sex)
meta$Region_custom <- factor(meta$Region_custom)
meta$Geography <- factor(meta$Geography)
meta$BioProject <- factor(meta$BioProject)

MIN_N <- 50
site_table <- meta %>%
  group_by(Site_Cleaned) %>%
  summarise(n = n(), n_bioproject = n_distinct(BioProject), .groups = "drop") %>%
  arrange(desc(n)) %>%
  filter(n >= MIN_N)

cat("\n=== Site sizes ===\n")
print(site_table)
write.csv(site_table, "site_stratified_group_sizes.csv", row.names = FALSE)

out_root <- "site_stratified_results"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

all_results <- list()
run_summary <- list()

for (i in seq_len(nrow(site_table))) {
  site <- site_table$Site_Cleaned[i]
  n_samp <- site_table$n[i]
  site_key <- gsub(" ", "_", site)
  out_dir <- file.path(out_root, site_key)

  cat(sprintf("\n[%d/%d] %s (n=%d)\n", i, nrow(site_table), site_key, n_samp))

  site_idx <- meta$Site_Cleaned == site
  sub_meta <- meta[site_idx, , drop = FALSE]
  sub_abund <- abund[, rownames(sub_meta), drop = FALSE]

  sub_meta$Age_Group <- droplevels(sub_meta$Age_Group)
  sub_meta$Sex <- droplevels(sub_meta$Sex)
  sub_meta$Region_custom <- droplevels(sub_meta$Region_custom)
  sub_meta$Geography <- droplevels(sub_meta$Geography)
  sub_meta$BioProject <- droplevels(sub_meta$BioProject)

  if ("18-60 years" %in% levels(sub_meta$Age_Group)) sub_meta$Age_Group <- relevel(sub_meta$Age_Group, ref = "18-60 years")
  if ("F" %in% levels(sub_meta$Sex)) sub_meta$Sex <- relevel(sub_meta$Sex, ref = "F")
  if ("V4" %in% levels(sub_meta$Region_custom)) sub_meta$Region_custom <- relevel(sub_meta$Region_custom, ref = "V4")
  if ("North America" %in% levels(sub_meta$Geography)) sub_meta$Geography <- relevel(sub_meta$Geography, ref = "North America")

  fixed_effects <- c("Age_Group", "Sex", "Region_custom", "Geography")
  fixed_effects <- fixed_effects[sapply(fixed_effects, function(x) nlevels(sub_meta[[x]]) >= 2)]
  n_bp <- nlevels(sub_meta$BioProject)

  cat(sprintf("  Levels: Age=%d, Sex=%d, Region=%d, Geography=%d, BioProject=%d\n",
              nlevels(sub_meta$Age_Group), nlevels(sub_meta$Sex),
              nlevels(sub_meta$Region_custom), nlevels(sub_meta$Geography), n_bp))
  cat(sprintf("  Fixed effects: %s\n", paste(fixed_effects, collapse = ", ")))

  if (length(fixed_effects) == 0 || n_bp < 2) {
    cat("  SKIP: insufficient factor levels\n")
    next
  }

  res_file_existing <- file.path(out_dir, "all_results.tsv")
  if (file.exists(res_file_existing)) {
    cat("  SKIP: existing completed result found\n")
    res <- read.delim(res_file_existing, stringsAsFactors = FALSE)
    res$Site_Cleaned <- site
    all_results[[site_key]] <- res
    run_summary[[site_key]] <- data.frame(
      Site_Cleaned = site, Samples = nrow(sub_meta), BioProjects = n_bp,
      Fixed_Effects = paste(fixed_effects, collapse = ","),
      Associations = nrow(res), Significant_q0.1 = sum(res$qval_individual < 0.1, na.rm = TRUE)
    )
    next
  }

  result <- tryCatch({
    maaslin3(
      input_data = t(sub_abund),
      input_metadata = sub_meta,
      output = out_dir,
      fixed_effects = fixed_effects,
      random_effects = c("BioProject"),
      min_abundance = 0,
      min_prevalence = 0.1,
      max_significance = 0.25,
      normalization = "TSS",
      transform = "LOG",
      correction = "BH",
      warn_prevalence = FALSE,
      cores = MAASLIN3_CORES,
      verbosity = "INFO"
    )
    res <- read.delim(file.path(out_dir, "all_results.tsv"), stringsAsFactors = FALSE)
    res$Site_Cleaned <- site
    res
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    NULL
  })

  if (!is.null(result)) {
    all_results[[site_key]] <- result
    run_summary[[site_key]] <- data.frame(
      Site_Cleaned = site, Samples = nrow(sub_meta), BioProjects = n_bp,
      Fixed_Effects = paste(fixed_effects, collapse = ","),
      Associations = nrow(result), Significant_q0.1 = sum(result$qval_individual < 0.1, na.rm = TRUE)
    )
    cat(sprintf("  Done: %d associations, %d significant (q<0.1)\n",
                nrow(result), sum(result$qval_individual < 0.1, na.rm = TRUE)))
  }
  gc(verbose = FALSE)
}

if (length(all_results) > 0) {
  combined <- bind_rows(all_results, .id = "Site_Group")
  write.table(combined, "site_stratified_results_all.tsv", sep = "\t", quote = TRUE, row.names = FALSE)
  saveRDS(all_results, file = "site_stratified_results.rds")
  write.csv(bind_rows(run_summary), "site_stratified_run_summary.csv", row.names = FALSE)
}

cat("\n=== Done ===\n")
