#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3")
source("input_helpers.R")
out_root <- "disease_results"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

cat("=== Loading data ===\n")
genus_data <- load_genus_relabund()
metadata <- load_standard_metadata(dedup = TRUE)

sample_ids <- trimws(as.character(metadata$SampleID))
keep_ids <- !is.na(sample_ids) & sample_ids != "" & !duplicated(sample_ids)
metadata <- metadata[keep_ids, , drop = FALSE]
sample_ids <- sample_ids[keep_ids]

common_samples <- intersect(sample_ids, colnames(genus_data))
meta <- metadata[match(common_samples, sample_ids), , drop = FALSE]
rownames(meta) <- common_samples
abund <- genus_data[, common_samples, drop = FALSE]

for (col in c("Age_Group", "Sex", "Region_custom", "BioProject",
              "Site_Cleaned", "Geography", "Status_Cleaned")) {
  meta[[col]] <- trimws(as.character(meta[[col]]))
}

valid <- !is.na(meta$BioProject) & !is.na(meta$Status_Cleaned) &
  meta$BioProject != "" & meta$Status_Cleaned != "" &
  meta$Status_Cleaned != "Other" & meta$Status_Cleaned != "Healthy"
meta <- meta[valid, , drop = FALSE]
abund <- abund[, rownames(meta), drop = FALSE]

disease_counts <- sort(table(meta$Status_Cleaned), decreasing = TRUE)
write.csv(as.data.frame(disease_counts), "disease_status_counts.csv", row.names = FALSE)
print(disease_counts)

preferred_diseases <- c(
  "Dental Caries",
  "Head and Neck Cancer",
  "HIV",
  "Oral Cancer",
  "OSCC",
  "Periodontal Disease",
  "RA",
  "Benign Tumor",
  "Breast Cancer",
  "CRC",
  "IBD",
  "Lung Cancer"
)
missing_diseases <- setdiff(preferred_diseases, names(disease_counts))
if (length(missing_diseases) > 0) {
  cat(sprintf("Missing diseases in metadata, will skip: %s\n",
              paste(missing_diseases, collapse = ", ")))
}
available <- preferred_diseases[preferred_diseases %in% names(disease_counts)]

extract_genus <- function(x) {
  g <- regmatches(x, regexpr("g__[^;]*", x))
  if (length(g) == 0 || identical(g, "g__")) return(x)
  gsub("g__", "", g)
}

all_results <- list()
run_summary <- list()

for (disease in available) {
  disease_bp <- unique(meta$BioProject[meta$Status_Cleaned == disease])
  keep <- meta$Status_Cleaned %in% c("Healthy Control", disease) &
    meta$BioProject %in% disease_bp
  sub_meta <- meta[keep, , drop = FALSE]
  sub_abund <- abund[, rownames(sub_meta), drop = FALSE]

  bp_tbl <- table(sub_meta$BioProject, sub_meta$Status_Cleaned)
  if (!("Healthy Control" %in% colnames(bp_tbl)) || !(disease %in% colnames(bp_tbl))) {
    cat("  SKIP: no matched Healthy Control/case BioProjects after filtering\n")
    next
  }
  valid_bp <- rownames(bp_tbl)[bp_tbl[, "Healthy Control"] > 0 & bp_tbl[, disease] > 0]
  if (length(valid_bp) == 0) {
    cat("  SKIP: no BioProject contains both case and Healthy Control\n")
    next
  }
  sub_meta <- sub_meta[sub_meta$BioProject %in% valid_bp, , drop = FALSE]
  sub_abund <- sub_abund[, rownames(sub_meta), drop = FALSE]

  n_case <- sum(sub_meta$Status_Cleaned == disease)
  n_ctrl <- sum(sub_meta$Status_Cleaned == "Healthy Control")
  n_bp <- length(unique(sub_meta$BioProject))
  cat(sprintf("\n=== %s: %d case / %d control / %d matched BioProjects ===\n",
              disease, n_case, n_ctrl, n_bp))

  if (n_case < 30 || n_ctrl < 30 || n_bp < 2) {
    cat("  SKIP: insufficient matched support\n")
    next
  }

  sub_meta$Disease_Status <- factor(
    ifelse(sub_meta$Status_Cleaned == disease, disease, "Healthy Control"),
    levels = c("Healthy Control", disease)
  )
  sub_meta$BioProject <- droplevels(factor(sub_meta$BioProject))
  optional_effects <- c("Age_Group", "Sex", "Site_Cleaned", "Region_custom", "Geography")
  for (eff in optional_effects) {
    sub_meta[[eff]][is.na(sub_meta[[eff]]) | sub_meta[[eff]] == ""] <- "Unknown"
    sub_meta[[eff]] <- droplevels(factor(sub_meta[[eff]]))
  }

  fixed_effects <- c("Disease_Status")
  for (eff in optional_effects) {
    if (nlevels(sub_meta[[eff]]) >= 2) fixed_effects <- c(fixed_effects, eff)
  }

  out_dir <- file.path(out_root, gsub(" ", "_", disease))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  res_file_existing <- file.path(out_dir, "all_results.tsv")
  if (file.exists(res_file_existing)) {
    cat("  SKIP: existing completed result found\n")
    res <- read.delim(res_file_existing, stringsAsFactors = FALSE)
    res$Disease <- disease
    sig_count <- sum(res$metadata == "Disease_Status" & res$qval_individual < 0.1, na.rm = TRUE)
    run_summary[[disease]] <- data.frame(
      Disease = disease,
      Cases = n_case,
      Controls = n_ctrl,
      Matched_BioProjects = n_bp,
      Fixed_Effects = paste(fixed_effects, collapse = ","),
      Significant_Disease_Associations = sig_count
    )
    all_results[[disease]] <- res
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
      cores = max(1, min(8, parallel::detectCores() - 1)),
      verbosity = "INFO"
    )

    res <- read.delim(file.path(out_dir, "all_results.tsv"), stringsAsFactors = FALSE)
    res$Disease <- disease
    res
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    NULL
  })

  if (is.null(result)) next

  sig_count <- sum(result$metadata == "Disease_Status" & result$qval_individual < 0.1, na.rm = TRUE)
  run_summary[[disease]] <- data.frame(
    Disease = disease,
    Cases = n_case,
    Controls = n_ctrl,
    Matched_BioProjects = n_bp,
    Fixed_Effects = paste(fixed_effects, collapse = ","),
    Significant_Disease_Associations = sig_count
  )
  all_results[[disease]] <- result
}

if (length(all_results) == 0) {
  cat("No disease analyses completed.\n")
  quit(save = "no", status = 0)
}

combined_results <- bind_rows(all_results)
write.table(combined_results, "disease_results_all.tsv",
            sep = "\t", quote = TRUE, row.names = FALSE)
write.csv(bind_rows(run_summary), "disease_run_summary.csv", row.names = FALSE)

disease_sig <- combined_results %>%
  filter(metadata == "Disease_Status", !is.na(qval_individual), qval_individual < 0.1) %>%
  mutate(Genus = vapply(feature, extract_genus, character(1)),
         Direction = ifelse(coef > 0, 1, -1)) %>%
  select(Disease, Genus, Direction) %>%
  distinct()

if (nrow(disease_sig) > 0) {
  hm <- disease_sig %>%
    pivot_wider(names_from = Disease, values_from = Direction, values_fill = 0)
  hm_mat <- as.matrix(hm[, -1, drop = FALSE])
  rownames(hm_mat) <- hm$Genus
  pdf("disease_association_heatmap.pdf", width = 7, height = max(6, nrow(hm_mat) * 0.15))
  pheatmap(hm_mat,
           color = c("#4575b4", "white", "#d73027"),
           breaks = c(-1.5, -0.5, 0.5, 1.5),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           main = "Disease-associated genera (q < 0.1)")
  dev.off()
}

cat("\n=== Done ===\n")

