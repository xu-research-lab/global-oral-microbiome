#!/usr/bin/env Rscript
# ============================================================================
# Run MaAsLin3 for North America saliva samples ONLY, using stratified
# downsampled subset and full default models (abundance + prevalence).
# ============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
})

BASE <- "/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3"
setwd(BASE)

OUT_DIR <- file.path(BASE, "saliva_sex_downsampled_na", "results", "Saliva_North_America")
META_FILE <- file.path(BASE, "saliva_sex_downsampled_na", "downsampled_metadata.tsv")
ABUND_FILE <- file.path(BASE, "saliva_sex_downsampled_na", "downsampled_genus_relabund.tsv")
Q_THRESHOLD <- 0.1
CORES <- max(1, min(8, parallel::detectCores() - 1))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Loading downsampled data ===\n")
meta <- read.delim(META_FILE, stringsAsFactors = FALSE, check.names = FALSE)
abund <- read.delim(ABUND_FILE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

# Filter to North America saliva samples only
meta <- meta[meta$Geography == "North America" & meta$Site_Cleaned == "Saliva", , drop = FALSE]
rownames(meta) <- meta$SampleID; abund <- abund[, rownames(meta), drop = FALSE]
cat(sprintf("North America saliva samples: %d\n", nrow(meta)))

meta$Age_Group <- droplevels(factor(meta$Age_Group))
meta$Sex <- droplevels(factor(meta$Sex))
meta$Region_custom <- droplevels(factor(meta$Region_custom))
meta$BioProject <- droplevels(factor(meta$BioProject))

if ("18-60 years" %in% levels(meta$Age_Group)) meta$Age_Group <- relevel(meta$Age_Group, ref = "18-60 years")
if ("F" %in% levels(meta$Sex)) meta$Sex <- relevel(meta$Sex, ref = "F")
if ("V4" %in% levels(meta$Region_custom)) meta$Region_custom <- relevel(meta$Region_custom, ref = "V4")

cat(sprintf("Levels: Age=%d, Sex=%d, Region=%d, BioProject=%d\n",
            nlevels(meta$Age_Group), nlevels(meta$Sex), nlevels(meta$Region_custom), nlevels(meta$BioProject)))

res_file <- file.path(OUT_DIR, "all_results.tsv")
if (!file.exists(res_file)) {
  cat("Running MaAsLin3 (default abundance + prevalence models)...\n")
  maaslin3(
    input_data = t(abund),
    input_metadata = meta,
    output = OUT_DIR,
    fixed_effects = c("Age_Group", "Sex", "Region_custom"),
    random_effects = c("BioProject"),
    min_abundance = 0,
    min_prevalence = 0.1,
    max_significance = 0.25,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    warn_prevalence = FALSE,
    cores = CORES,
    verbosity = "INFO"
  )
} else {
  cat("Existing result found; skipping.\n")
}

cat("Done.\n")
