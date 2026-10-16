#!/usr/bin/env Rscript
# ============================================================================
# da_v4 LinDA: Disease vs HC, mixed model with (1|bioproject_id).
# Diseases come from eligible_diseases*.txt (matched BP, Case/HC >=30, BP >=2).
# is.winsor = FALSE (no winsorization). q <= 0.05.
# ============================================================================

suppressPackageStartupMessages({
  library(MicrobiomeStat)
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
OUT_DIR <- file.path(WORK_DIR, paste0("linda", RARE_TAG))
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

pick_linda_disease_term <- function(output_names, disease) {
  hits <- output_names[grepl("Disease_Status", output_names)]
  hits <- hits[!grepl("Intercept", hits, ignore.case = TRUE)]
  if (length(hits) == 0) return(NA_character_)
  prefer <- hits[grepl(disease, hits, fixed = TRUE)]
  if (length(prefer) >= 1) return(prefer[[1]])
  hits[[1]]
}

extract_linda_stats <- function(linda_obj, disease) {
  term <- pick_linda_disease_term(names(linda_obj$output), disease)
  if (is.na(term)) {
    stop("No Disease_Status term in: ", paste(names(linda_obj$output), collapse = ", "))
  }
  res_df <- linda_obj$output[[term]]
  out <- data.frame(
    feature = rownames(res_df),
    lfc = as.numeric(res_df$log2FoldChange),
    se = as.numeric(res_df$lfcSE),
    pval = as.numeric(res_df$pvalue),
    qval = as.numeric(res_df$padj),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$lfc) & is.finite(out$lfc), , drop = FALSE]
  attr(out, "term") <- term
  out
}

merge_all <- function() {
  cat("=== Merging LinDA results ===\n")
  res_files <- list.files(RESULTS_DIR, pattern = "^all_results\\.tsv$",
                          recursive = TRUE, full.names = TRUE)
  if (length(res_files) == 0) stop("No all_results.tsv under ", RESULTS_DIR)
  combined <- bind_rows(lapply(res_files, read.delim, stringsAsFactors = FALSE,
                               check.names = FALSE))
  # Recompute so the merged table cannot carry a stale Genus column written by
  # an earlier run under a different naming rule.
  combined$Genus <- clean_genus(combined$feature)
  combined$Direction <- ifelse(combined$lfc > 0, "+", "-")
  write.table(combined, file.path(OUT_DIR, "linda_all_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  summ_files <- list.files(RESULTS_DIR, pattern = "^run_summary\\.csv$",
                           recursive = TRUE, full.names = TRUE)
  if (length(summ_files) > 0) {
    summ <- bind_rows(lapply(summ_files, read.csv, stringsAsFactors = FALSE))
    write.csv(summ, file.path(OUT_DIR, "linda_run_summary.csv"), row.names = FALSE)
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
  write.table(sig, file.path(OUT_DIR, "linda_significant_genera_q0.05.tsv"),
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

genus_counts <- read.table(file.path(WORK_DIR, COUNTS_TABLE),
                           header = TRUE, sep = "\t", row.names = 1,
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

run_linda_once <- function(feature_mat, meta_df, formula_str) {
  tryCatch({
    linda(
      feature.dat = feature_mat,
      meta.dat = meta_df,
      formula = paste0("~", formula_str),
      feature.dat.type = "count",
      prev.filter = PRV_CUT,
      mean.abund.filter = 0,
      max.abund.filter = 0,
      is.winsor = FALSE,
      outlier.pct = 0.03,
      p.adj.method = "BH",
      alpha = Q_THRESHOLD,
      n.cores = 1,
      verbose = TRUE
    )
  }, error = function(e) {
    cat(sprintf("ERROR linda (%s): %s\n", formula_str, e$message)); NULL
  })
}

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

  for (eff in OPTIONAL_EFFECTS) {
    sub_meta[[eff]][is.na(sub_meta[[eff]]) | sub_meta[[eff]] == ""] <- "Unknown"
    sub_meta[[eff]] <- droplevels(factor(sub_meta[[eff]]))
  }
  fixed_effects <- c("Disease_Status")
  for (eff in OPTIONAL_EFFECTS) {
    if (nlevels(sub_meta[[eff]]) >= 2) fixed_effects <- c(fixed_effects, eff)
  }
  fix_formula <- paste(c(fixed_effects, "(1|bioproject_id)"), collapse = " + ")
  fix_only <- paste(fixed_effects, collapse = " + ")

  lib <- colSums(sub_counts)
  keep_s <- lib > 0
  sub_counts <- sub_counts[, keep_s, drop = FALSE]
  sub_meta <- sub_meta[keep_s, , drop = FALSE]
  keep_f <- rowSums(sub_counts) > 0
  sub_counts <- sub_counts[keep_f, , drop = FALSE]

  n_case <- sum(sub_meta$Disease_Status == disease)
  n_ctrl <- sum(sub_meta$Disease_Status == "Healthy Control")
  n_bp <- nlevels(droplevels(sub_meta$bioproject_id))
  cat(sprintf("Final n=%d (case=%d HC=%d BP=%d); formula: ~%s\n",
              nrow(sub_meta), n_case, n_ctrl, n_bp, fix_formula))

  disease_dir <- file.path(RESULTS_DIR, disease_dir_name(disease))
  dir.create(disease_dir, showWarnings = FALSE, recursive = TRUE)

  linda_res <- run_linda_once(sub_counts, sub_meta, fix_formula)
  used_formula <- fix_formula
  if (is.null(linda_res)) {
    cat("Retry Disease_Status + (1|bioproject_id) only\n")
    used_formula <- "Disease_Status + (1|bioproject_id)"
    linda_res <- run_linda_once(sub_counts, sub_meta, used_formula)
  }
  if (is.null(linda_res)) {
    cat("No result\n"); return(invisible(NULL))
  }

  stats <- tryCatch(extract_linda_stats(linda_res, disease), error = function(e) {
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
  stats$Direction <- ifelse(stats$lfc > 0, "+", "-")
  stats$Genus <- clean_genus(stats$feature)
  write.table(stats, file.path(disease_dir, "all_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  n_sig <- sum(!is.na(stats$qval) & stats$qval <= Q_THRESHOLD)
  summ <- data.frame(
    Disease = disease, Cases = n_case, Controls = n_ctrl,
    Matched_BioProjects = n_bp, Fixed_Effects = fix_only,
    Significant_q0.05 = n_sig, stringsAsFactors = FALSE
  )
  write.csv(summ, file.path(disease_dir, "run_summary.csv"), row.names = FALSE)
  cat(sprintf("Done %s: features=%d sig=%d (term=%s)\n",
              disease, nrow(stats), n_sig, attr(stats, "term")))
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
cat("\n=== LinDA da_v4 done (is.winsor=FALSE) ===\n")
