#!/usr/bin/env Rscript
# ============================================================================
# MaAsLin3 v2 — Analysis 2&3: Sex & Region in Saliva, by Geography
# - Saliva only, per continent
# - North America: stratified downsampling to ~2500
# - Model: age_cleaned_group + sex_cleaned + sequencing_region_cleaned + (1|bioproject_id)
# - Output: combined results for downstream plotting
# ============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

BASE <- "/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3_v2"
setwd(BASE)

OUT_DIR <- "by_geography_saliva"
TARGET_NA_N <- 2500
Q_THRESHOLD <- 0.1
MAASLIN3_CORES <- 1

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "plots"), showWarnings = FALSE, recursive = TRUE)

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
# 2. Match samples — Saliva only
# ============================================================================
sample_ids <- trimws(as.character(metadata$sample_id))
keep_ids <- !is.na(sample_ids) & sample_ids != "" & !duplicated(sample_ids)
metadata <- metadata[keep_ids, , drop = FALSE]
sample_ids <- sample_ids[keep_ids]

common_samples <- intersect(sample_ids, colnames(genus_data))
abund <- genus_data[, common_samples, drop = FALSE]
meta <- metadata[match(common_samples, sample_ids), , drop = FALSE]
rownames(meta) <- common_samples

key_cols <- c("age_cleaned_group", "sex_cleaned", "sequencing_region_cleaned",
              "bioproject_id", "site_cleaned", "continent_group")
for (col in key_cols) {
  meta[[col]] <- trimws(as.character(meta[[col]]))
  meta[[col]][is.na(meta[[col]]) | meta[[col]] == "" | meta[[col]] == "NA"] <- NA_character_
}

# Saliva + complete cases
meta <- meta[meta$site_cleaned == "Saliva" & !is.na(meta$site_cleaned), , drop = FALSE]
abund <- abund[, rownames(meta), drop = FALSE]
valid <- complete.cases(meta[, key_cols])
meta <- meta[valid, , drop = FALSE]
abund <- abund[, rownames(meta), drop = FALSE]
cat(sprintf("Saliva samples with complete metadata: %d\n", nrow(meta)))

# Factorize
meta$age_cleaned_group <- factor(meta$age_cleaned_group)
meta$sex_cleaned      <- factor(meta$sex_cleaned)
meta$sequencing_region_cleaned <- factor(meta$sequencing_region_cleaned)
meta$bioproject_id    <- factor(meta$bioproject_id)
meta$continent_group  <- factor(meta$continent_group)

# ============================================================================
# 3. Stratified downsampling — North America
# ============================================================================
na_samples <- rownames(meta)[meta$continent_group == "North America"]
na_meta <- meta[na_samples, , drop = FALSE]
cat(sprintf("\nOriginal North America saliva samples: %d\n", length(na_samples)))

na_meta$stratum <- interaction(
  na_meta$bioproject_id, na_meta$sex_cleaned,
  na_meta$age_cleaned_group, na_meta$sequencing_region_cleaned,
  drop = TRUE, sep = "__"
)

# Collapse tiny strata (<3)
stratum_counts <- table(na_meta$stratum)
small_strata <- names(stratum_counts[stratum_counts < 3])
if (length(small_strata) > 0) {
  cat(sprintf("Collapsing %d small strata (<3 samples)\n", length(small_strata)))
  for (s in small_strata) {
    parts <- strsplit(s, "__", fixed = TRUE)[[1]]
    na_meta$stratum[na_meta$stratum == s] <- paste(parts[1], parts[2],
                                                   "__pooled__", parts[4], sep = "__")
  }
  na_meta$stratum <- factor(na_meta$stratum)
}

stratum_n <- table(na_meta$stratum)
prop_n <- floor(TARGET_NA_N * stratum_n / sum(stratum_n))
prop_n[prop_n == 0] <- 1
if (sum(prop_n) > TARGET_NA_N) {
  extra <- sum(prop_n) - TARGET_NA_N
  sorted <- order(prop_n, decreasing = TRUE)
  for (i in seq_len(extra)) {
    k <- sorted[((i - 1) %% length(sorted)) + 1]
    if (prop_n[k] > 1) prop_n[k] <- prop_n[k] - 1
  }
}
cat(sprintf("Target NA allocation sum: %d\n", sum(prop_n)))

selected_in_na <- unlist(lapply(names(stratum_n), function(s) {
  rows <- which(na_meta$stratum == s)
  n_take <- min(prop_n[s], length(rows))
  sample(rows, n_take)
}))
na_selected <- na_samples[selected_in_na]
cat(sprintf("Selected NA samples after downsampling: %d\n", length(na_selected)))

# Combine: non-NA + downsampled NA
non_na <- rownames(meta)[meta$continent_group != "North America"]
keep_samples <- c(non_na, na_selected)
meta_ds <- meta[keep_samples, , drop = FALSE]
abund_ds <- abund[, keep_samples, drop = FALSE]

meta_ds$age_cleaned_group <- droplevels(meta_ds$age_cleaned_group)
meta_ds$sex_cleaned      <- droplevels(meta_ds$sex_cleaned)
meta_ds$sequencing_region_cleaned <- droplevels(meta_ds$sequencing_region_cleaned)
meta_ds$bioproject_id    <- droplevels(meta_ds$bioproject_id)
meta_ds$continent_group  <- droplevels(meta_ds$continent_group)

if ("18-60 years" %in% levels(meta_ds$age_cleaned_group))
  meta_ds$age_cleaned_group <- relevel(meta_ds$age_cleaned_group, ref = "18-60 years")
if ("F" %in% levels(meta_ds$sex_cleaned))
  meta_ds$sex_cleaned <- relevel(meta_ds$sex_cleaned, ref = "F")
if ("V4" %in% levels(meta_ds$sequencing_region_cleaned))
  meta_ds$sequencing_region_cleaned <- relevel(meta_ds$sequencing_region_cleaned, ref = "V4")

write.table(meta_ds, file.path(OUT_DIR, "downsampled_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

geo_counts <- meta_ds %>% count(continent_group)
write.csv(geo_counts, file.path(OUT_DIR, "geography_counts_downsampled.csv"), row.names = FALSE)
cat("\nGeography counts after downsampling:\n")
print(geo_counts)

# ============================================================================
# 4. Run MaAsLin3 per geography
# ============================================================================
geographies <- as.character(unique(meta_ds$continent_group))
geographies <- geographies[geographies != "Africa"]

all_results <- list()

for (geo in geographies) {
  cat(sprintf("\n[%s] Running MaAsLin3...\n", geo))
  geo_idx <- meta_ds$continent_group == geo
  sub_meta <- meta_ds[geo_idx, , drop = FALSE]
  sub_abund <- abund_ds[, rownames(sub_meta), drop = FALSE]

  sub_meta$age_cleaned_group <- droplevels(sub_meta$age_cleaned_group)
  sub_meta$sex_cleaned      <- droplevels(sub_meta$sex_cleaned)
  sub_meta$sequencing_region_cleaned <- droplevels(sub_meta$sequencing_region_cleaned)
  sub_meta$bioproject_id    <- droplevels(sub_meta$bioproject_id)

  n_age <- nlevels(sub_meta$age_cleaned_group)
  n_sex <- nlevels(sub_meta$sex_cleaned)
  n_region <- nlevels(sub_meta$sequencing_region_cleaned)
  n_bp <- nlevels(sub_meta$bioproject_id)
  cat(sprintf("  n=%d; Age=%d; Sex=%d; Region=%d; BioProject=%d\n",
              nrow(sub_meta), n_age, n_sex, n_region, n_bp))

  if (n_age < 2 || n_sex < 2 || n_region < 2 || n_bp < 2) {
    cat("  SKIP: insufficient factor levels\n"); next
  }

  geo_dir <- file.path(OUT_DIR, "results", paste0("Saliva_", gsub(" ", "_", geo)))
  dir.create(geo_dir, showWarnings = FALSE, recursive = TRUE)

  res_file <- file.path(geo_dir, "all_results.tsv")
  if (!file.exists(res_file)) {
    maaslin3(
      input_data = t(sub_abund), input_metadata = sub_meta,
      output = geo_dir,
      fixed_effects = c("age_cleaned_group", "sex_cleaned",
                        "sequencing_region_cleaned"),
      random_effects = c("bioproject_id"),
      evaluate_only = NULL,
      min_abundance = 0, min_prevalence = 0.1, max_significance = 0.25,
      normalization = "TSS", transform = "LOG", correction = "BH",
      warn_prevalence = FALSE, cores = MAASLIN3_CORES, verbosity = "INFO"
    )
  }

  res <- read.delim(res_file, stringsAsFactors = FALSE)
  res$continent_group <- geo
  all_results[[geo]] <- res
  cat(sprintf("  %d associations, %d sig (q<%.2f)\n", nrow(res),
              sum(res$qval_individual < Q_THRESHOLD & !is.na(res$qval_individual)), Q_THRESHOLD))
}

combined <- bind_rows(all_results, .id = "Geography")
write.table(combined, file.path(OUT_DIR, "all_results_by_geography.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)

# Quick summary
clean_genus <- function(x) {
  x <- trimws(as.character(x))
  has_genus <- grepl("g__", x)
  extracted <- regmatches(x, regexpr("g__[^;]*", x))
  x[has_genus] <- sub("^g__", "", extracted)
  x <- sub("_.+$", "", x)
  ifelse(is.na(x) | x == "" | x == "g__", NA_character_, x)
}
combined$Genus <- clean_genus(combined$feature)

sig_summary <- combined %>%
  filter(!is.na(qval_individual), qval_individual < Q_THRESHOLD) %>%
  group_by(metadata, continent_group) %>%
  summarise(n_associations = n(), n_genera = n_distinct(Genus), .groups = "drop") %>%
  arrange(metadata, desc(n_associations))

cat("\n=== Significant summary ===\n")
print(sig_summary)
write.csv(sig_summary, file.path(OUT_DIR, "significant_summary_by_geography.csv"), row.names = FALSE)
cat("\n=== By-geography saliva analysis done ===\n")
