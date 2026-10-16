#!/usr/bin/env Rscript
# ============================================================================
# da_v4 ANCOM-BC2: Disease vs HC, mixed model with (1|bioproject_id).
# Diseases come from eligible_diseases*.txt (matched BP, Case/HC >=30, BP >=2).
# q <= 0.05.
# ============================================================================

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(phyloseq)
  library(dplyr)
})

# Repository root: run from it, or set DA_WORK_DIR. Raw metadata defaults to
# data/ under the root (not distributed; see README), or set DA_META_PATH.
WORK_DIR <- normalizePath(Sys.getenv("DA_WORK_DIR", "."), mustWork = TRUE)
META_PATH <- Sys.getenv(
  "DA_META_PATH", file.path(WORK_DIR, "data", "metadata_oral_qc_dedup.tsv")
)
# Track selector: "" (raw) / "_rare1k" / "_rare5k"
RARE_TAG <- Sys.getenv("DA_RARE_TAG", "")
COUNTS_TABLE <- Sys.getenv(
  "DA_COUNTS_TABLE", paste0("genus_table_counts", RARE_TAG, ".tsv")
)
DISEASE_LIST <- Sys.getenv("DA_DISEASE_LIST", "eligible_diseases.txt")
OUT_DIR <- file.path(WORK_DIR, paste0("ancombc2", RARE_TAG))
RESULTS_DIR <- file.path(OUT_DIR, "disease_results")
Q_THRESHOLD <- 0.05
PRV_CUT <- 0.01

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

read_disease_list <- function(fname) {
  path <- if (file.exists(fname)) fname else file.path(WORK_DIR, fname)
  if (!file.exists(path)) {
    stop("Disease list not found: ", path,
         " (run compute_eligible_diseases.py first)")
  }
  d <- trimws(readLines(path, warn = FALSE))
  d[nzchar(d)]
}

clean_genus <- function(x) {
  x <- trimws(as.character(x))
  has_genus <- grepl("g__", x)
  extracted <- regmatches(x, regexpr("g__[^;]*", x))
  x[has_genus] <- sub("^g__", "", extracted)
  ifelse(is.na(x) | x == "" | x == "g__", NA_character_, x)
}

disease_dir_name <- function(disease) gsub(" ", "_", disease)

pick_col <- function(cn, prefixes) {
  for (prefix in prefixes) {
    hits <- grep(paste0("^", prefix), cn, value = TRUE)
    hits <- hits[!grepl("Intercept", hits, ignore.case = TRUE)]
    hits <- hits[grepl("Disease_Status", hits)]
    if (length(hits) >= 1) return(hits[[1]])
  }
  NA_character_
}

extract_disease_stats <- function(res_df) {
  cn <- colnames(res_df)
  lfc_col <- pick_col(cn, c("lfc_"))
  se_col <- pick_col(cn, c("se_"))
  p_col <- pick_col(cn, c("p_"))
  q_col <- pick_col(cn, c("q_", "qval_"))
  if (is.na(lfc_col) || is.na(se_col)) {
    stop("Cannot find Disease_Status lfc/se: ", paste(cn, collapse = ", "))
  }
  taxon_col <- if ("taxon" %in% cn) "taxon" else cn[[1]]
  out <- data.frame(
    feature = as.character(res_df[[taxon_col]]),
    lfc = as.numeric(res_df[[lfc_col]]),
    se = as.numeric(res_df[[se_col]]),
    pval = if (!is.na(p_col)) as.numeric(res_df[[p_col]]) else NA_real_,
    qval = if (!is.na(q_col)) as.numeric(res_df[[q_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )
  if (all(is.na(out$qval)) && !all(is.na(out$pval))) {
    out$qval <- p.adjust(out$pval, method = "BH")
  }
  out <- out[!is.na(out$lfc) & is.finite(out$lfc), , drop = FALSE]
  out
}

merge_all <- function() {
  cat("=== Merging ANCOM-BC2 results ===\n")
  res_files <- list.files(RESULTS_DIR, pattern = "^all_results\\.tsv$",
                          recursive = TRUE, full.names = TRUE)
  if (length(res_files) == 0) stop("No all_results.tsv under ", RESULTS_DIR)
  combined <- bind_rows(lapply(res_files, read.delim, stringsAsFactors = FALSE,
                               check.names = FALSE))
  # Recompute so the merged table cannot carry a stale Genus column written by
  # an earlier run under a different naming rule.
  combined$Genus <- clean_genus(combined$feature)
  combined$Direction <- ifelse(combined$lfc > 0, "+", "-")
  write.table(combined, file.path(OUT_DIR, "ancombc2_all_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  summ_files <- list.files(RESULTS_DIR, pattern = "^run_summary\\.csv$",
                           recursive = TRUE, full.names = TRUE)
  if (length(summ_files) > 0) {
    summ <- bind_rows(lapply(summ_files, read.csv, stringsAsFactors = FALSE))
    write.csv(summ, file.path(OUT_DIR, "ancombc2_run_summary.csv"), row.names = FALSE)
  }

  sig <- combined %>%
    filter(!is.na(qval), qval <= Q_THRESHOLD) %>%
    mutate(
      Genus = clean_genus(feature),
      Direction = ifelse(lfc > 0, "+", "-")
    ) %>%
    filter(!is.na(Genus)) %>%
    group_by(Disease, Genus) %>%
    arrange(qval, desc(abs(lfc)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  write.table(sig, file.path(OUT_DIR, "ancombc2_significant_genera_q0.05.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Significant genera rows: %d\n", nrow(sig)))
}

if (MERGE_ONLY) {
  merge_all()
  quit(save = "no", status = 0)
}

cat("=== Loading counts & metadata ===\n")
cat(sprintf("OUT_DIR=%s  COUNTS_TABLE=%s  DISEASE_LIST=%s\n",
            OUT_DIR, COUNTS_TABLE, DISEASE_LIST))
ELIGIBLE_DISEASES <- read_disease_list(DISEASE_LIST)
cat(sprintf("Eligible diseases (%d): %s\n", length(ELIGIBLE_DISEASES),
            paste(ELIGIBLE_DISEASES, collapse = ", ")))

counts_path <- file.path(WORK_DIR, COUNTS_TABLE)
genus_counts <- read.table(counts_path, header = TRUE, sep = "\t", row.names = 1,
                           check.names = FALSE, comment.char = "",
                           stringsAsFactors = FALSE)
genus_counts <- as.matrix(genus_counts)
storage.mode(genus_counts) <- "double"

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
common_samples <- intersect(sample_ids, colnames(genus_counts))
cat(sprintf("Common samples: %d\n", length(common_samples)))

counts <- genus_counts[, common_samples, drop = FALSE]
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
counts <- counts[, rownames(meta), drop = FALSE]

disease_counts <- sort(table(meta$status_cleaned), decreasing = TRUE)
available <- ELIGIBLE_DISEASES[ELIGIBLE_DISEASES %in% names(disease_counts)]
if (!is.null(TARGET_DISEASE)) {
  if (!(TARGET_DISEASE %in% available)) stop("Disease not available: ", TARGET_DISEASE)
  available <- TARGET_DISEASE
}
cat(sprintf("Diseases to run: %s\n", paste(available, collapse = ", ")))

run_one <- function(disease) {
  cat(sprintf("\n########## %s ##########\n", disease))
  disease_bp <- unique(meta$bioproject_id[meta$status_cleaned == disease])
  keep <- meta$status_cleaned %in% c("Healthy Control", disease) &
    meta$bioproject_id %in% disease_bp
  sub_meta <- meta[keep, , drop = FALSE]
  sub_counts <- counts[, rownames(sub_meta), drop = FALSE]

  bp_tbl <- table(sub_meta$bioproject_id, sub_meta$status_cleaned)
  if (!("Healthy Control" %in% colnames(bp_tbl)) || !(disease %in% colnames(bp_tbl))) {
    cat("SKIP: no matched HC/Case\n"); return(invisible(NULL))
  }
  valid_bp <- rownames(bp_tbl)[bp_tbl[, "Healthy Control"] > 0 & bp_tbl[, disease] > 0]
  if (length(valid_bp) == 0) {
    cat("SKIP: no matched BP\n"); return(invisible(NULL))
  }
  sub_meta <- sub_meta[sub_meta$bioproject_id %in% valid_bp, , drop = FALSE]
  sub_counts <- sub_counts[, rownames(sub_meta), drop = FALSE]

  n_case <- sum(sub_meta$status_cleaned == disease)
  n_ctrl <- sum(sub_meta$status_cleaned == "Healthy Control")
  n_bp <- length(unique(sub_meta$bioproject_id))
  cat(sprintf("Matched: %d case / %d HC / %d BP\n", n_case, n_ctrl, n_bp))
  if (n_case < 30 || n_ctrl < 30 || n_bp < 2) {
    cat("SKIP: need >=30 each and >=2 BP\n"); return(invisible(NULL))
  }

  sub_meta$Disease_Status <- factor(
    ifelse(sub_meta$status_cleaned == disease, disease, "Healthy Control"),
    levels = c("Healthy Control", disease)
  )
  sub_meta$bioproject_id <- droplevels(factor(as.character(sub_meta$bioproject_id)))

  # Fill missing covariates with "Unknown" (no complete-case drop)
  for (eff in OPTIONAL_EFFECTS) {
    sub_meta[[eff]][is.na(sub_meta[[eff]]) | sub_meta[[eff]] == ""] <- "Unknown"
    sub_meta[[eff]] <- droplevels(factor(sub_meta[[eff]]))
  }
  fixed_effects <- c("Disease_Status")
  for (eff in OPTIONAL_EFFECTS) {
    if (nlevels(sub_meta[[eff]]) >= 2) fixed_effects <- c(fixed_effects, eff)
  }
  fix_formula <- paste(fixed_effects, collapse = " + ")

  lib <- colSums(sub_counts)
  keep_s <- lib > 0
  sub_counts <- sub_counts[, keep_s, drop = FALSE]
  sub_meta <- sub_meta[keep_s, , drop = FALSE]
  keep_f <- rowSums(sub_counts) > 0
  sub_counts <- sub_counts[keep_f, , drop = FALSE]

  # ANCOM-BC2 works on log abundance with zeros as missing, so each genus is
  # fitted only on the samples where it is detected. A genus detected in a
  # single BioProject leaves (1|bioproject_id) with one sampled level, and one
  # detected in a single group leaves Disease_Status without a contrast; either
  # aborts the whole disease because ancombc2 smoke-tests only the first genus.
  detected <- sub_counts > 0
  bp_ind <- model.matrix(~ 0 + sub_meta$bioproject_id)
  grp_ind <- model.matrix(~ 0 + sub_meta$Disease_Status)
  keep_fit <- rowSums((detected %*% bp_ind) > 0) >= 2 &
    rowSums((detected %*% grp_ind) > 0) >= 2
  cat(sprintf("Unfittable genera dropped: %d / %d\n",
              sum(!keep_fit), length(keep_fit)))
  sub_counts <- sub_counts[keep_fit, , drop = FALSE]

  n_case <- sum(sub_meta$Disease_Status == disease)
  n_ctrl <- sum(sub_meta$Disease_Status == "Healthy Control")
  n_bp <- nlevels(droplevels(sub_meta$bioproject_id))
  cat(sprintf("Final n=%d (case=%d HC=%d BP=%d); formula: %s + (1|bioproject_id)\n",
              nrow(sub_meta), n_case, n_ctrl, n_bp, fix_formula))

  disease_dir <- file.path(RESULTS_DIR, disease_dir_name(disease))
  dir.create(disease_dir, showWarnings = FALSE, recursive = TRUE)

  ps <- phyloseq(otu_table(sub_counts, taxa_are_rows = TRUE),
                 sample_data(sub_meta))

  run_once <- function(ff) {
    tryCatch({
      ancombc2(
        data = ps,
        fix_formula = ff,
        rand_formula = "(1 | bioproject_id)",
        p_adj_method = "BH",
        prv_cut = PRV_CUT,
        lib_cut = 0,
        group = NULL,
        struc_zero = FALSE,
        neg_lb = FALSE,
        alpha = Q_THRESHOLD,
        n_cl = 1,
        verbose = TRUE,
        global = FALSE,
        pairwise = FALSE,
        dunnet = FALSE,
        trend = FALSE
      )
    }, error = function(e) {
      cat(sprintf("ERROR ancombc2 (%s): %s\n", ff, e$message)); NULL
    })
  }

  ancom_res <- run_once(fix_formula)
  used_formula <- fix_formula
  if ((is.null(ancom_res) || is.null(ancom_res$res)) && fix_formula != "Disease_Status") {
    cat("Retry Disease_Status only\n")
    ancom_res <- run_once("Disease_Status")
    used_formula <- "Disease_Status"
  }
  if (is.null(ancom_res) || is.null(ancom_res$res)) {
    cat("No result\n"); return(invisible(NULL))
  }

  stats <- tryCatch(extract_disease_stats(ancom_res$res), error = function(e) {
    cat(sprintf("ERROR extract: %s\n", e$message)); NULL
  })
  if (is.null(stats) || nrow(stats) == 0) {
    cat("Empty stats\n"); return(invisible(NULL))
  }
  stats$Disease <- disease
  stats$n_case <- n_case
  stats$n_ctrl <- n_ctrl
  stats$n_bp <- n_bp
  stats$fix_formula <- used_formula
  stats$rand_formula <- "(1 | bioproject_id)"
  stats$Direction <- ifelse(stats$lfc > 0, "+", "-")
  stats$Genus <- clean_genus(stats$feature)
  write.table(stats, file.path(disease_dir, "all_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  n_sig <- sum(!is.na(stats$qval) & stats$qval <= Q_THRESHOLD)
  summ <- data.frame(
    Disease = disease, Cases = n_case, Controls = n_ctrl,
    Matched_BioProjects = n_bp, Fixed_Effects = used_formula,
    Significant_q0.05 = n_sig, stringsAsFactors = FALSE
  )
  write.csv(summ, file.path(disease_dir, "run_summary.csv"), row.names = FALSE)
  cat(sprintf("Done %s: features=%d sig=%d\n", disease, nrow(stats), n_sig))
  invisible(summ)
}

for (disease in available) {
  tryCatch(run_one(disease), error = function(e) {
    cat(sprintf("FATAL in %s: %s\n", disease, e$message))
  })
  gc(verbose = FALSE)
}

if (is.null(TARGET_DISEASE)) {
  tryCatch(merge_all(), error = function(e) cat("Merge failed:", e$message, "\n"))
}
cat("\n=== ANCOM-BC2 da_v4 done ===\n")
