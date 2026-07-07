#!/usr/bin/env Rscript
# ============================================================================
# MaAsLin3 Unified Mixed Model for Sequencing Region Effects
# Fixed effects: Region_custom + available covariates
# Random effect: BioProject
# ============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3")
source("input_helpers.R")

OUT_ROOT <- "unified_region_model"
MAASLIN3_CORES <- 10
Q_THRESHOLDS <- c(0.05, 0.1)
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

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

needed_cols <- c("Region_custom", "BioProject", "Geography", "Site_Cleaned", "Age_Group", "Sex", "Status_Cleaned")
for (col in needed_cols) {
  if (!col %in% colnames(meta)) meta[[col]] <- NA_character_
  meta[[col]] <- trimws(as.character(meta[[col]]))
  meta[[col]][is.na(meta[[col]]) | meta[[col]] == "" | meta[[col]] == "NA"] <- NA_character_
}

valid <- !is.na(meta$Region_custom) & !is.na(meta$BioProject)
meta <- meta[valid, , drop = FALSE]
abund <- abund[, rownames(meta), drop = FALSE]
cat(sprintf("Samples with Region_custom and BioProject: %d\n", nrow(meta)))

optional_effects <- c("Geography", "Site_Cleaned", "Age_Group", "Sex", "Status_Cleaned")
for (eff in optional_effects) {
  meta[[eff]][is.na(meta[[eff]]) | meta[[eff]] == ""] <- "Unknown"
}

meta$Region_custom <- droplevels(factor(meta$Region_custom))
meta$BioProject <- droplevels(factor(meta$BioProject))
for (eff in optional_effects) meta[[eff]] <- droplevels(factor(meta[[eff]]))

if ("V4" %in% levels(meta$Region_custom)) meta$Region_custom <- relevel(meta$Region_custom, ref = "V4")
if ("North America" %in% levels(meta$Geography)) meta$Geography <- relevel(meta$Geography, ref = "North America")
if ("Saliva" %in% levels(meta$Site_Cleaned)) meta$Site_Cleaned <- relevel(meta$Site_Cleaned, ref = "Saliva")
if ("18-60 years" %in% levels(meta$Age_Group)) meta$Age_Group <- relevel(meta$Age_Group, ref = "18-60 years")
if ("F" %in% levels(meta$Sex)) meta$Sex <- relevel(meta$Sex, ref = "F")
if ("Healthy Control" %in% levels(meta$Status_Cleaned)) meta$Status_Cleaned <- relevel(meta$Status_Cleaned, ref = "Healthy Control")

fixed_effects <- c("Region_custom")
for (eff in optional_effects) {
  if (nlevels(meta[[eff]]) >= 2) fixed_effects <- c(fixed_effects, eff)
}

cat("\n=== Model factors ===\n")
cat(sprintf("Region_custom levels: %s\n", paste(levels(meta$Region_custom), collapse = ", ")))
cat(sprintf("BioProject levels: %d\n", nlevels(meta$BioProject)))
cat(sprintf("Fixed effects: %s\n", paste(fixed_effects, collapse = ", ")))

factor_summary <- bind_rows(lapply(c("Region_custom", optional_effects, "BioProject"), function(col) {
  data.frame(Variable = col, Level = names(table(meta[[col]])), N = as.integer(table(meta[[col]])), stringsAsFactors = FALSE)
}))
write.csv(factor_summary, file.path(OUT_ROOT, "unified_model_factor_counts.csv"), row.names = FALSE)
write.csv(data.frame(
  Samples = nrow(meta),
  Features = nrow(abund),
  BioProjects = nlevels(meta$BioProject),
  Fixed_Effects = paste(fixed_effects, collapse = ","),
  Cores = MAASLIN3_CORES,
  stringsAsFactors = FALSE
), file.path(OUT_ROOT, "unified_model_run_summary.csv"), row.names = FALSE)

if (nlevels(meta$Region_custom) < 2 || nlevels(meta$BioProject) < 2) {
  stop("Insufficient Region_custom or BioProject levels for unified model.")
}

res_file <- file.path(OUT_ROOT, "all_results.tsv")
if (!file.exists(res_file)) {
  maaslin3(
    input_data = t(abund),
    input_metadata = meta,
    output = OUT_ROOT,
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
} else {
  cat("Existing all_results.tsv found; skipping MaAsLin3 fit.\n")
}

results <- read.delim(res_file, stringsAsFactors = FALSE)
region_results <- results %>%
  filter(metadata == "Region_custom") %>%
  arrange(qval_individual, desc(abs(coef)))
write.table(region_results, file.path(OUT_ROOT, "region_custom_results_all.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

for (q in Q_THRESHOLDS) {
  sig <- region_results %>% filter(!is.na(qval_individual), qval_individual <= q)
  q_label <- gsub("\\.", "p", as.character(q))
  write.table(sig, file.path(OUT_ROOT, paste0("region_custom_significant_q", q_label, ".tsv")), sep = "\t", quote = TRUE, row.names = FALSE)
}

sig_summary <- bind_rows(lapply(Q_THRESHOLDS, function(q) {
  data.frame(
    q_threshold = q,
    significant_associations = sum(!is.na(region_results$qval_individual) & region_results$qval_individual <= q),
    significant_features = length(unique(region_results$feature[!is.na(region_results$qval_individual) & region_results$qval_individual <= q])),
    stringsAsFactors = FALSE
  )
}))
write.csv(sig_summary, file.path(OUT_ROOT, "region_custom_significant_summary.csv"), row.names = FALSE)

cat("\n=== Region_custom significant summary ===\n")
print(sig_summary)
cat("\n=== Done ===\n")
