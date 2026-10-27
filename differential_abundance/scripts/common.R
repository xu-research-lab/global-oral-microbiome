# Shared input schema, model settings and result harmonization.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

DA_SEED <- 20260615

DA_ALPHA <- 0.05

PREV_CUT <- 0.01

MIN_ABUNDANCE <- 0

LIB_CUT <- 0

DROP_BIOPROJECTS <- NULL

DROP_SITES <- NULL

COL_SAMPLE <- "sample_id"

COL_CONTINENT <- "continent_group"

COL_SITE <- "site_cleaned"

COL_SITE_GROUP <- "site_group_cleaned"

COL_AGE <- "age_cleaned_group"

COL_SEX <- "sex_cleaned"

COL_STATUS <- "status_cleaned"

COL_BMI <- "bmi_cleaned_group"

COL_PROJECT <- "bioproject_id"

COL_PLATFORM <- "platform_cleaned"

COL_SEQ_REGION <- "sequencing_region_cleaned"

COMMON_COLS <- c(
  "feature", "genus_short", "variable", "baseMean", "log2FoldChange", "lfcSE", "stat",
  "pvalue", "padj", "reject", "df", "direction", "enriched_in"
)

sanitize <- function(x) gsub(" ", "_", x)

load_genus_table <- function(which = c("counts", "relabund")) {
  which <- match.arg(which)
  path <- if (which == "counts") {
    COUNTS_TSV
  } else {
    RELABUND_TSV
  }
  if (!file.exists(path)) {
    stop("Genus tsv not found: ", path, " (run prepare_genus_tables.py / sbatch first)")
  }
  df <- read.table(path,
    header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, comment.char = "",
    stringsAsFactors = FALSE
  )
  m <- as.matrix(df)
  storage.mode(m) <- "double"
  m
}

align_samples <- function(meta, genus_mat) {
  if (!is.null(DROP_BIOPROJECTS)) {
    n0 <- nrow(meta)
    meta <- meta[!(as.character(meta[[COL_PROJECT]]) %in% DROP_BIOPROJECTS), , drop = FALSE]
    if (nrow(meta) < n0) {
      cat("    Dropped ", n0 - nrow(meta), " samples by DROP_BIOPROJECTS\n", sep = "")
    }
  }
  if (COL_SITE %in% names(meta) && !is.null(DROP_SITES)) {
    n0 <- nrow(meta)
    meta <- meta[!(as.character(meta[[COL_SITE]]) %in% DROP_SITES), , drop = FALSE]
    if (nrow(meta) < n0) {
      cat("    Dropped ", n0 - nrow(meta), " samples by DROP_SITES\n", sep = "")
    }
  }
  common <- intersect(rownames(meta), colnames(genus_mat))
  if (length(common) == 0) {
    stop("No common samples between metadata and genus table!")
  }
  meta <- meta[common, , drop = FALSE]
  mat <- genus_mat[, common, drop = FALSE]
  list(meta = meta, mat = mat)
}

make_genus_short <- function(feature) {
  g <- sub(".*g__", "", feature)
  sub(";.*", "", g)
}

harmonize_results <- function(
    df, group1, group2, comparison, target_site = NA_character_, group1_col = "Group1",
    group2_col = "Group2") {
  if (nrow(df) == 0) {
    return(data.frame(
      feature = character(0), genus_short = character(0), variable = character(0),
      baseMean = numeric(0), log2FoldChange = numeric(0), lfcSE = numeric(0), stat = numeric(0),
      pvalue = numeric(0), padj = numeric(0), reject = logical(0), df = numeric(0), direction = character(0),
      enriched_in = character(0), stringsAsFactors = FALSE
    ))
  }
  out <- data.frame(feature = df$feature, genus_short = make_genus_short(df$feature), variable = if (!is.null(df$variable)) {
    df$variable
  } else {
    NA_character_
  }, baseMean = if (!is.null(df$baseMean)) {
    df$baseMean
  } else {
    NA_real_
  }, log2FoldChange = df$log2FoldChange, lfcSE = if (!is.null(df$lfcSE)) {
    df$lfcSE
  } else {
    NA_real_
  }, stat = if (!is.null(df$stat)) {
    df$stat
  } else {
    NA_real_
  }, pvalue = df$pvalue, padj = df$padj, reject = if (!is.null(df$reject)) {
    df$reject
  } else {
    (df$padj < DA_ALPHA)
  }, df = if (!is.null(df$df)) {
    df$df
  } else {
    NA_real_
  }, stringsAsFactors = FALSE)
  out$direction <- ifelse(is.na(out$log2FoldChange), NA_character_, ifelse(out$log2FoldChange > 0,
    "up_in_target", ifelse(out$log2FoldChange < 0, "up_in_reference", "zero")
  ))
  out$enriched_in <- ifelse(is.na(out$log2FoldChange), NA_character_, ifelse(out$log2FoldChange > 0,
    group2, ifelse(out$log2FoldChange < 0, group1, "none")
  ))
  out[[group1_col]] <- group1
  out[[group2_col]] <- group2
  out$Comparison <- comparison
  if (!is.na(target_site)) {
    out$Target_Site <- target_site
  }
  keep <- c(COMMON_COLS, group1_col, group2_col, "Comparison", if (!is.na(target_site)) "Target_Site" else NULL)
  out[, keep, drop = FALSE]
}

linda_to_common <- function(linda_fit, group_coef) {
  if (is.null(linda_fit$output) || !group_coef %in% names(linda_fit$output)) {
    stop("LinDA coefficient not in output: ", group_coef, " (have: ", paste(names(linda_fit$output),
      collapse = ", "
    ), ")")
  }
  res <- linda_fit$output[[group_coef]]
  data.frame(
    feature = rownames(res), variable = group_coef, baseMean = res$baseMean, log2FoldChange = res$log2FoldChange,
    lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj, reject = res$reject,
    df = res$df, stringsAsFactors = FALSE
  )
}

.maalind2_rev_map_cache <- NULL

restore_feature_names <- function(feature_vec) {
  if (is.null(.maalind2_rev_map_cache)) {
    feats <- fread(COUNTS_TSV, sep = "\t", header = TRUE, select = "V1", showProgress = FALSE)$V1
    mk <- make.names(feats)
    if (anyDuplicated(mk)) {
      stop("make.names not 1:1 on feature ids; cannot restore safely.")
    }
    .maalind2_rev_map_cache <<- setNames(feats, mk)
  }
  m <- .maalind2_rev_map_cache
  known <- feature_vec %in% names(m)
  ifelse(known, m[feature_vec], feature_vec)
}

maaslin2_to_common <- function(maaslin_all, group_col, group2) {
  grp_rows <- maaslin_all[maaslin_all$metadata == group_col & as.character(maaslin_all$value) == group2, ,
    drop = FALSE
  ]
  if (nrow(grp_rows) == 0) {
    grp_rows <- maaslin_all[maaslin_all$metadata == group_col, , drop = FALSE]
  }
  feat_fixed <- restore_feature_names(as.character(grp_rows$feature))
  data.frame(
    feature = feat_fixed, variable = paste0(group_col, group2), baseMean = NA_real_, log2FoldChange = grp_rows$coef,
    lfcSE = grp_rows$stderr, stat = grp_rows$coef / pmax(grp_rows$stderr, .Machine$double.xmin), pvalue = grp_rows$pval,
    padj = grp_rows$qval, reject = !is.na(grp_rows$qval) & grp_rows$qval < DA_ALPHA, df = NA_real_,
    stringsAsFactors = FALSE
  )
}

ancombc2_to_common <- function(ancombc_res, group_coef) {
  lfc_col <- paste0("lfc_", group_coef)
  se_col <- paste0("se_", group_coef)
  W_col <- paste0("W_", group_coef)
  p_col <- paste0("p_", group_coef)
  q_col <- paste0("q_", group_coef)
  data.frame(
    feature = ancombc_res$taxon, variable = group_coef, baseMean = NA_real_, log2FoldChange = ancombc_res[[lfc_col]],
    lfcSE = ancombc_res[[se_col]], stat = ancombc_res[[W_col]], pvalue = ancombc_res[[p_col]], padj = ancombc_res[[q_col]],
    reject = !is.na(ancombc_res[[q_col]]) & ancombc_res[[q_col]] < DA_ALPHA, df = NA_real_, stringsAsFactors = FALSE
  )
}

ancombc_group_coef <- function(group_col, group2) {
  paste0(group_col, group2)
}

collinear_covars_from_error <- function(error_msg, valid_covars) {
  if (length(valid_covars) == 0) {
    return(character(0))
  }
  lines <- strsplit(error_msg, "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines)
  hdr <- grep("Estimation failed", lines)
  fin <- grep("Please ensure|collinearity|missing values", lines)
  if (length(hdr) == 0) {
    return(character(0))
  }
  end <- if (length(fin)) {
    min(fin[fin > hdr[1]]) - 1L
  } else {
    length(lines)
  }
  tokens <- lines[(hdr[1] + 1L):end]
  tokens <- tokens[nzchar(tokens)]
  drop <- character(0)
  for (tk in tokens) {
    cand <- valid_covars[vapply(valid_covars, function(v) substr(tk, 1L, nchar(v)) == v, logical(1))]
    if (length(cand) >= 1L) {
      drop <- unique(c(drop, cand[which.max(nchar(cand))]))
    }
  }
  drop
}

# Persist one complete main-contrast table and its fit configuration.
write_comparison_outputs <- function(
    all_res, output_subdir, comparison, group1, group2,
    group1_col, group2_col, main_results_var_pattern, summary_row, run_config_row,
    main_results_filename, sanitize_names = TRUE) {
  dir.create(output_subdir, recursive = TRUE, showWarnings = FALSE)
  main <- all_res[grepl(main_results_var_pattern, all_res$variable), , drop = FALSE]
  main <- main[order(main$padj, -abs(main$log2FoldChange)), , drop = FALSE]
  sig <- main[!is.na(main$padj) & main$padj < DA_ALPHA, , drop = FALSE]
  fwrite(main, file.path(output_subdir, main_results_filename), sep = "\t", quote = FALSE)
  fwrite(summary_row, file.path(output_subdir, "summary.tsv"), sep = "\t", quote = FALSE)
  fwrite(run_config_row, file.path(output_subdir, "run_config.tsv"), sep = "\t", quote = FALSE)
  writeLines("done", file.path(output_subdir, "analysis.done"))
  invisible(list(
    n_sig = nrow(sig), n_up_g1 = sum(sig$enriched_in == group1, na.rm = TRUE),
    n_up_g2 = sum(sig$enriched_in == group2, na.rm = TRUE)
  ))
}
