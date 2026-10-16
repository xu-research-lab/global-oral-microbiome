#!/usr/bin/env Rscript
# ============================================================================
# da_v4 MaAsLin2: Disease vs HC, all sites, bioproject_id random effect.
# Diseases come from eligible_diseases*.txt (matched BP, Case/HC >=30, BP >=2).
# q <= 0.05. Runs on R4 (Maaslin2 package).
# ============================================================================

suppressPackageStartupMessages({
  library(Maaslin2)
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
RELAB_TABLE <- Sys.getenv(
  "DA_RELAB_TABLE", paste0("genus_table_relabund", RARE_TAG, ".tsv")
)
DISEASE_LIST <- Sys.getenv("DA_DISEASE_LIST", "eligible_diseases.txt")
OUT_DIR <- file.path(WORK_DIR, paste0("maaslin2", RARE_TAG))
RESULTS_DIR <- file.path(OUT_DIR, "disease_results")
Q_THRESHOLD <- 0.05
PRV_CUT <- 0.01
MAASLIN2_CORES <- 1

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

# MaAsLin2 pushes feature names through make.names(), so GTDB labels such as
# "CAG-177" or "1XD8-76" come back as "CAG.177" / "X1XD8.76" and stop matching
# the names ANCOM-BC2 and LinDA report. Map them back to the abundance table.
canonical_features <- function() {
  path <- if (file.exists(RELAB_TABLE)) RELAB_TABLE else file.path(WORK_DIR, RELAB_TABLE)
  if (!file.exists(path)) return(character())
  con <- pipe(sprintf("cut -f1 %s", shQuote(path)), "r")
  feats <- readLines(con)
  close(con)
  if (length(feats) <= 1) character() else feats[-1]
}

unmangle_features <- function(x, canon) {
  if (!length(canon)) return(x)
  mangled <- make.names(canon)
  ambiguous <- unique(mangled[duplicated(mangled)])
  usable <- !(mangled %in% ambiguous) & mangled != canon
  lut <- setNames(canon[usable], mangled[usable])
  hit <- x %in% names(lut)
  x[hit] <- unname(lut[x[hit]])
  x
}

# MaAsLin2 uses qval/pval; tolerate qval_individual if present
pick_qval <- function(df) {
  if ("qval" %in% names(df)) return(as.numeric(df$qval))
  if ("qval_individual" %in% names(df)) return(as.numeric(df$qval_individual))
  rep(NA_real_, nrow(df))
}

merge_all <- function() {
  cat("=== Merging MaAsLin2 results ===\n")
  res_files <- list.files(RESULTS_DIR, pattern = "^all_results\\.tsv$",
                          recursive = TRUE, full.names = TRUE)
  if (length(res_files) == 0) stop("No all_results.tsv under ", RESULTS_DIR)
  all_res <- bind_rows(lapply(res_files, function(f) {
    df <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE,
                     colClasses = "character")
    num_cols <- c("coef", "stderr", "pval", "qval", "pval_individual",
                  "qval_individual", "N", "N.not.0", "N_not_zero")
    for (nc in intersect(num_cols, names(df))) {
      df[[nc]] <- suppressWarnings(as.numeric(df[[nc]]))
    }
    if (!("Disease" %in% names(df))) {
      dname <- basename(dirname(f))
      df$Disease <- gsub("_", " ", dname)
    }
    df
  }))
  canon <- canonical_features()
  n_before <- sum(!(all_res$feature %in% canon))
  all_res$feature <- unmangle_features(as.character(all_res$feature), canon)
  cat(sprintf("Feature names restored from make.names(): %d -> %d unmatched\n",
              n_before, sum(!(all_res$feature %in% canon))))
  all_res$Genus <- clean_genus(all_res$feature)
  write.table(all_res, file.path(OUT_DIR, "maaslin2_all_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  summ_files <- list.files(RESULTS_DIR, pattern = "^run_summary\\.csv$",
                           recursive = TRUE, full.names = TRUE)
  if (length(summ_files) > 0) {
    summ <- bind_rows(lapply(summ_files, read.csv, stringsAsFactors = FALSE))
    write.csv(summ, file.path(OUT_DIR, "maaslin2_run_summary.csv"), row.names = FALSE)
  }

  all_res$qval_use <- pick_qval(all_res)
  sig <- all_res %>%
    filter(metadata == "Disease_Status",
           !is.na(qval_use), qval_use <= Q_THRESHOLD) %>%
    mutate(
      coef = as.numeric(coef),
      Genus = clean_genus(feature),
      Direction = ifelse(coef > 0, "+", "-"),
      qval = qval_use,
      lfc = coef
    ) %>%
    filter(!is.na(Genus), !is.na(coef)) %>%
    group_by(Disease, Genus) %>%
    arrange(qval, desc(abs(coef)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(-qval_use)
  write.table(sig, file.path(OUT_DIR, "maaslin2_significant_genera_q0.05.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Significant genera rows: %d\n", nrow(sig)))
}

if (MERGE_ONLY) {
  merge_all()
  quit(save = "no", status = 0)
}

cat("=== Loading data ===\n")
cat(sprintf("OUT_DIR=%s  RELAB_TABLE=%s  DISEASE_LIST=%s\n",
            OUT_DIR, RELAB_TABLE, DISEASE_LIST))
ELIGIBLE_DISEASES <- read_disease_list(DISEASE_LIST)
cat(sprintf("Eligible diseases (%d): %s\n", length(ELIGIBLE_DISEASES),
            paste(ELIGIBLE_DISEASES, collapse = ", ")))

genus_data <- read.table(file.path(WORK_DIR, RELAB_TABLE),
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
available <- ELIGIBLE_DISEASES[ELIGIBLE_DISEASES %in% names(disease_counts)]
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
    Maaslin2(
      input_data = t(sub_abund),
      input_metadata = sub_meta,
      output = disease_dir,
      fixed_effects = fixed_effects,
      random_effects = c("bioproject_id"),
      reference = "Disease_Status,Healthy Control",
      min_abundance = 0,
      min_prevalence = PRV_CUT,
      max_significance = 0.25,
      normalization = "TSS",
      transform = "LOG",
      analysis_method = "LM",
      correction = "BH",
      standardize = TRUE,
      plot_heatmap = FALSE,
      plot_scatter = FALSE,
      save_scatter = FALSE,
      cores = MAASLIN2_CORES
    )
    if (!file.exists(res_file)) {
      cand <- list.files(disease_dir, pattern = "all_results\\.tsv$",
                         recursive = TRUE, full.names = TRUE)
      if (length(cand) == 0) stop("MaAsLin2 did not write all_results.tsv")
      file.copy(cand[[1]], res_file, overwrite = TRUE)
    }
    read.delim(res_file, stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message)); NULL
  })
  if (is.null(result)) return(invisible(NULL))
  result$Disease <- disease
  write.table(result, res_file, sep = "\t", quote = FALSE, row.names = FALSE)

  qv <- pick_qval(result)
  n_sig <- sum(result$metadata == "Disease_Status" & !is.na(qv) & qv <= Q_THRESHOLD)
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

if (is.null(TARGET_DISEASE)) {
  tryCatch(merge_all(), error = function(e) cat("Merge failed:", e$message, "\n"))
}
cat("\n=== MaAsLin2 da_v4 done ===\n")
