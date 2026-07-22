#!/usr/bin/env Rscript
# ============================================================================
# da_v2 MaAsLin3: Disease vs HC, all sites, bioproject_id random effect.
# Filter: matched BP (>0 both), Case/HC >=30, matched BP >=2. q <= 0.05.
# ============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
})

WORK_DIR <- "/home/zhangcaiyi/project/global_oral_diseae_da/da_v2"
META_PATH <- paste0(
  "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/",
  "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv"
)
OUT_DIR <- file.path(WORK_DIR, "maaslin3")
RESULTS_DIR <- file.path(OUT_DIR, "disease_results")
Q_THRESHOLD <- 0.05
MAASLIN3_CORES <- 1

PREFERRED_DISEASES <- c(
  "Dental Caries", "Head and Neck Cancer", "HIV", "Oral Cancer", "OSCC",
  "Periodontal Disease", "RA", "Benign Tumor", "Breast Cancer", "CRC",
  "IBD", "Lung Cancer"
)
OPTIONAL_EFFECTS <- c(
  "age_cleaned_group", "sex_cleaned", "site_cleaned",
  "sequencing_region_cleaned", "continent_group"
)

args <- commandArgs(trailingOnly = TRUE)
MERGE_ONLY <- "--merge-only" %in% args
args <- setdiff(args, "--merge-only")
TARGET_DISEASE <- if (length(args) >= 1) args[[1]] else NULL

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(WORK_DIR)

clean_genus <- function(x) {
  x <- trimws(as.character(x))
  has_genus <- grepl("g__", x)
  extracted <- regmatches(x, regexpr("g__[^;]*", x))
  x[has_genus] <- sub("^g__", "", extracted)
  x <- sub("_.+$", "", x)
  ifelse(is.na(x) | x == "" | x == "g__", NA_character_, x)
}

disease_dir_name <- function(disease) gsub(" ", "_", disease)

merge_all <- function() {
  cat("=== Merging MaAsLin3 results ===\n")
  res_files <- list.files(RESULTS_DIR, pattern = "^all_results\\.tsv$",
                          recursive = TRUE, full.names = TRUE)
  if (length(res_files) == 0) stop("No all_results.tsv under ", RESULTS_DIR)
  all_res <- bind_rows(lapply(res_files, function(f) {
    df <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE,
                     colClasses = "character")
    num_cols <- c("coef", "null_hypothesis", "stderr", "pval_individual",
                  "qval_individual", "pval_joint", "qval_joint", "N", "N_not_zero")
    for (nc in intersect(num_cols, names(df))) {
      df[[nc]] <- suppressWarnings(as.numeric(df[[nc]]))
    }
    if (!("Disease" %in% names(df))) {
      dname <- basename(dirname(f))
      df$Disease <- gsub("_", " ", dname)
    }
    df
  }))
  write.table(all_res, file.path(OUT_DIR, "maaslin3_all_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  summ_files <- list.files(RESULTS_DIR, pattern = "^run_summary\\.csv$",
                           recursive = TRUE, full.names = TRUE)
  if (length(summ_files) > 0) {
    summ <- bind_rows(lapply(summ_files, read.csv, stringsAsFactors = FALSE))
    write.csv(summ, file.path(OUT_DIR, "maaslin3_run_summary.csv"), row.names = FALSE)
  }

  sig <- all_res %>%
    filter(metadata == "Disease_Status",
           !is.na(qval_individual), qval_individual <= Q_THRESHOLD) %>%
    mutate(
      coef = as.numeric(coef),
      Genus = clean_genus(feature),
      Direction = ifelse(coef > 0, "+", "-"),
      qval = qval_individual,
      lfc = coef
    ) %>%
    filter(!is.na(Genus), !is.na(coef)) %>%
    group_by(Disease, Genus) %>%
    arrange(qval, desc(abs(coef)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  write.table(sig, file.path(OUT_DIR, "maaslin3_significant_genera_q0.05.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Significant genera rows: %d\n", nrow(sig)))
}

if (MERGE_ONLY) {
  merge_all()
  quit(save = "no", status = 0)
}

cat("=== Loading data ===\n")
genus_data <- read.table(file.path(WORK_DIR, "genus_table_relabund.tsv"),
                         header = TRUE, sep = "\t", row.names = 1,
                         check.names = FALSE, comment.char = "",
                         stringsAsFactors = FALSE)
metadata <- read.table(META_PATH, header = TRUE, sep = "\t", fill = TRUE,
                       comment.char = "", stringsAsFactors = FALSE, quote = "")

site_map <- c("saliva" = "Saliva", "tongue" = "Tongue", "tooth" = "Tooth",
              "oral swab" = "Oral_swab", "oral mucosa" = "Oral_Mucosa",
              "gingiva" = "Gingiva", "other" = "Others")
metadata$site_cleaned <- ifelse(metadata$site_cleaned %in% names(site_map),
                                site_map[metadata$site_cleaned], metadata$site_cleaned)
metadata$status_cleaned[metadata$status_cleaned %in% c("HC", "Healthy")] <- "Healthy Control"
metadata$age_cleaned_group <- gsub("＞", ">", metadata$age_cleaned_group)

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
available <- PREFERRED_DISEASES[PREFERRED_DISEASES %in% names(disease_counts)]
if (!is.null(TARGET_DISEASE)) {
  if (!(TARGET_DISEASE %in% available)) {
    stop("Disease not available: ", TARGET_DISEASE)
  }
  available <- TARGET_DISEASE
}
cat(sprintf("Diseases to run: %s\n", paste(available, collapse = ", ")))

run_one <- function(disease) {
  cat(sprintf("\n########## %s ##########\n", disease))
  disease_bp <- unique(meta$bioproject_id[meta$status_cleaned == disease])
  keep <- meta$status_cleaned %in% c("Healthy Control", disease) &
    meta$bioproject_id %in% disease_bp
  sub_meta <- meta[keep, , drop = FALSE]
  sub_abund <- abund[, rownames(sub_meta), drop = FALSE]

  bp_tbl <- table(sub_meta$bioproject_id, sub_meta$status_cleaned)
  if (!("Healthy Control" %in% colnames(bp_tbl)) || !(disease %in% colnames(bp_tbl))) {
    cat("SKIP: no matched HC/Case BioProjects\n"); return(invisible(NULL))
  }
  valid_bp <- rownames(bp_tbl)[bp_tbl[, "Healthy Control"] > 0 & bp_tbl[, disease] > 0]
  if (length(valid_bp) == 0) {
    cat("SKIP: no BioProject with both Case and HC\n"); return(invisible(NULL))
  }
  sub_meta <- sub_meta[sub_meta$bioproject_id %in% valid_bp, , drop = FALSE]
  sub_abund <- sub_abund[, rownames(sub_meta), drop = FALSE]

  n_case <- sum(sub_meta$status_cleaned == disease)
  n_ctrl <- sum(sub_meta$status_cleaned == "Healthy Control")
  n_bp <- length(unique(sub_meta$bioproject_id))
  cat(sprintf("Matched: %d case / %d HC / %d BP\n", n_case, n_ctrl, n_bp))
  if (n_case < 30 || n_ctrl < 30 || n_bp < 2) {
    cat("SKIP: need >=30 each and >=2 BioProjects\n"); return(invisible(NULL))
  }

  sub_meta$Disease_Status <- factor(
    ifelse(sub_meta$status_cleaned == disease, disease, "Healthy Control"),
    levels = c("Healthy Control", disease)
  )
  sub_meta$bioproject_id <- droplevels(factor(sub_meta$bioproject_id))

  for (eff in OPTIONAL_EFFECTS) {
    sub_meta[[eff]][is.na(sub_meta[[eff]]) | sub_meta[[eff]] == ""] <- "Unknown"
    sub_meta[[eff]] <- droplevels(factor(sub_meta[[eff]]))
  }
  fixed_effects <- c("Disease_Status")
  for (eff in OPTIONAL_EFFECTS) {
    if (nlevels(sub_meta[[eff]]) >= 2) fixed_effects <- c(fixed_effects, eff)
  }
  cat(sprintf("Fixed effects: %s\n", paste(fixed_effects, collapse = " + ")))

  disease_dir <- file.path(RESULTS_DIR, disease_dir_name(disease))
  dir.create(disease_dir, showWarnings = FALSE, recursive = TRUE)
  res_file <- file.path(disease_dir, "all_results.tsv")

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
    read.delim(res_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message)); NULL
  })
  if (is.null(result)) return(invisible(NULL))
  result$Disease <- disease
  write.table(result, res_file, sep = "\t", quote = FALSE, row.names = FALSE)

  n_sig <- sum(result$metadata == "Disease_Status" &
                 !is.na(result$qval_individual) &
                 result$qval_individual <= Q_THRESHOLD)
  summ <- data.frame(
    Disease = disease, Cases = n_case, Controls = n_ctrl,
    Matched_BioProjects = n_bp,
    Fixed_Effects = paste(fixed_effects, collapse = ","),
    Significant_q0.05 = n_sig, stringsAsFactors = FALSE
  )
  write.csv(summ, file.path(disease_dir, "run_summary.csv"), row.names = FALSE)
  cat(sprintf("Done %s: sig q<=0.05 = %d\n", disease, n_sig))
  invisible(summ)
}

for (disease in available) {
  tryCatch(run_one(disease), error = function(e) {
    cat(sprintf("FATAL in %s: %s\n", disease, e$message))
  })
  gc(verbose = FALSE)
}

tryCatch(merge_all(), error = function(e) cat("Merge failed:", e$message, "\n"))
cat("\n=== MaAsLin3 da_v2 done ===\n")
