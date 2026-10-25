# ==============================================================================
# 02-merge-and-consensus.R  (continent OVA da_v2 — tsv based)
# Merged post-processing for the three-method continent DA pipeline
# (LinDA + Maaslin2 + ANCOMBC2, OVA one-vs-ALL comparisons: each valid continent
# vs the pooled "Other" of the other valid continents).
#
# Stages (select with --stages=merge,consensus,summary,enhanced,fulltable;
#         default = all, in this order):
#   merge     - bind per-comparison results of each method
#   consensus - 3-way intersection flags per comparison. ROBUST = significant
#               (padj<0.05) in >=2 methods (n2 stratification — the primary
#               output). A method missing for a comparison is handled
#               automatically: >=2 of the available methods define robust.
#   summary   - per-comparison venn-region counts + grand summary
#   enhanced  - robust features joined with per-method LFC/padj/enriched_in
#   fulltable - ALL features tested by any method, per-method stats side by side
#
# NOTE: the OVO pipeline's symmetric summary_diff_counts.tsv (reversed
# contrasts) has no OVA analogue — the reverse contrast is the same
# comparison with flipped signs — so that stage is dropped here.
#
# Outputs (under <DA_RESULTS_ROOT>/<continent_ova | continent_hc_ova>/continent_da_consensus/):
#   consensus_per_comparison/ova__<comp>.tsv
#   consensus_all.tsv / consensus_robust.tsv (n_methods>=2) / consensus_all3.tsv (n==3)
#   consensus_summary.tsv / consensus_grand_summary.tsv
#   consensus_per_comparison/ova__<comp>_enhanced.tsv
#   consensus_all_enhanced.tsv / consensus_robust_enhanced.tsv / consensus_all3_enhanced.tsv
#   full_table.tsv / full_table_per_comparison/ova__<comp>.tsv
#
# Usage (from project root):
#   conda activate R4
#   Rscript differential_abundance/continent_differential_ova/02-merge-and-consensus.R
#   Rscript differential_abundance/continent_differential_ova/02-merge-and-consensus.R --scope=hc
#   Rscript differential_abundance/continent_differential_ova/02-merge-and-consensus.R --stages=merge,consensus
# ==============================================================================

suppressPackageStartupMessages({ library(data.table); library(dplyr) })
source("differential_abundance/continent_differential_ova/da_helpers.R")

# ==========================================================================
# 0. Configuration (continent scope, hardcoded for this directory)
# ==========================================================================
args <- commandArgs(trailingOnly = TRUE)
stages_arg <- "all"
scope <- parse_scope_arg(args)
for (a in args) if (grepl("^--stages=", a)) stages_arg <- sub("^--stages=", "", a)
ALL_STAGES <- c("merge", "consensus", "summary", "enhanced", "fulltable")
stages <- if (stages_arg == "all") ALL_STAGES else strsplit(stages_arg, ",", fixed = TRUE)[[1]]
bad <- setdiff(stages, ALL_STAGES)
if (length(bad) > 0) stop("Unknown --stages value(s): ", paste(bad, collapse = ", "))

alpha     <- DA_ALPHA
scope_tag <- if (scope == "hc") "continent OVA (HC only)" else "continent OVA (all)"
group_col <- COL_CONTINENT
g1_col    <- "Continent1"; g2_col <- "Continent2"
main_res_file <- "continent_results.tsv"
alt_res_files <- c("main_results.tsv", "significant_results.tsv")

results_root <- continent_root(scope)
method_dirs <- list(
  linda    = file.path(results_root, "continent_linda_ova"),
  maaslin2 = file.path(results_root, "continent_maaslin2_ova"),
  ancombc2 = file.path(results_root, "continent_ancombc2_ova")
)
out_dir <- file.path(results_root, "continent_da_consensus")

# ==========================================================================
# Shared helpers
# ==========================================================================
read_full <- function(dir_path, main_file, alt_files) {
  for (f in c(main_file, alt_files)) {
    p <- file.path(dir_path, f)
    if (file.exists(p)) {
      df <- fread(p, data.table = FALSE)
      if ("padj" %in% names(df) && "feature" %in% names(df)) return(df)
    }
  }
  NULL
}

read_main_sig <- function(dir_path, main_file, alt_files) {
  df <- read_full(dir_path, main_file, alt_files)
  if (is.null(df)) return(NULL)
  df <- df[!is.na(df$padj) & df$padj < alpha, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df
}

enrich_by_method <- function(df, feats) {
  if (is.null(df)) return(setNames(rep(NA_character_, length(feats)), feats))
  sub <- df[match(feats, df$feature), , drop = FALSE]
  enriched <- if ("enriched_in" %in% names(sub)) as.character(sub$enriched_in) else NA_character_
  setNames(enriched, feats)
}
direction_consensus_of <- function(e_list, feats) {
  vapply(feats, function(f) {
    vals <- unlist(lapply(e_list, function(e) e[f]))
    vals <- vals[!is.na(vals) & vals != "" & vals != "none"]
    if (length(vals) == 0) return("none")
    tab <- sort(table(vals), decreasing = TRUE)
    if (tab[1] >= 2) paste0("agreed:", names(tab)[1]) else "discordant"
  }, character(1))
}

cat("==========================================\n")
cat("Merge + consensus (3-method continent DA) —", scope_tag, "\n")
cat("Stages: ", paste(stages, collapse = ", "), "\n", sep = "")
cat("Robust definition: significant in >=2 methods (n2)\n")
cat("Start: ", format(Sys.time()), "\n")
cat("==========================================\n\n")

# ==========================================================================
# Stage 1: merge per-method results
# ==========================================================================
if ("merge" %in% stages) {
  cat("==== [merge] binding per-comparison results per method ====\n")

  merge_one_method <- function(results_dir) {
    if (!dir.exists(results_dir)) {
      cat("    [", results_dir, "] not found, skipping\n", sep = "")
      return(NULL)
    }
    run_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
    run_dirs <- run_dirs[file.exists(file.path(run_dirs, "analysis.done"))]
    if (length(run_dirs) == 0) {
      cat("    [", results_dir, "] no completed runs\n", sep = "")
      return(NULL)
    }
    all_list <- list(); main_list <- list(); sum_list <- list()
    for (d in run_dirs) {
      comp <- basename(d)
      all_file <- file.path(d, "all_results.tsv")
      if (file.exists(all_file)) {
        m <- fread(all_file, data.table = FALSE)
        if (!"Comparison" %in% names(m)) m$Comparison <- comp
        all_list[[comp]] <- m
      }
      main <- read_full(d, main_res_file, alt_res_files)
      if (!is.null(main)) main_list[[comp]] <- main
      sum_file <- file.path(d, "summary.tsv")
      if (file.exists(sum_file)) sum_list[[comp]] <- fread(sum_file, data.table = FALSE)
    }
    if (length(all_list) == 0 && length(main_list) == 0) return(NULL)
    list(all    = if (length(all_list))  bind_rows(all_list)  else NULL,
         main   = if (length(main_list)) bind_rows(main_list) else NULL,
         summary = if (length(sum_list)) bind_rows(sum_list)  else NULL)
  }

  for (mname in names(method_dirs)) {
    m <- merge_one_method(method_dirs[[mname]])
    if (is.null(m)) { cat("[", mname, "] nothing to merge\n", sep = ""); next }
    rdir <- method_dirs[[mname]]
    if (!is.null(m$all)) {
      fwrite(m$all, file.path(rdir, "all_ova_results.tsv"), sep = "\t", quote = FALSE)
      cat("[", mname, "] all_ova_results.tsv: ", nrow(m$all), " rows\n", sep = "")
    }
    if (!is.null(m$main)) {
      fwrite(m$main, file.path(rdir, "all_continent_ova_results.tsv"), sep = "\t", quote = FALSE)
      sig <- m$main[!is.na(m$main$padj) & m$main$padj < alpha, , drop = FALSE]
      fwrite(sig, file.path(rdir, "significant_ova_results.tsv"), sep = "\t", quote = FALSE)
      cat("[", mname, "] significant_ova_results.tsv: ", nrow(sig), " rows\n", sep = "")

      if (mname == "linda" && nrow(sig) > 0 && "enriched_in" %in% names(sig)) {
        sig_sorted <- sig[order(sig$Comparison, sig$padj, -abs(sig$log2FoldChange)), ]
        fwrite(sig_sorted, file.path(rdir, "significant_ova_results_with_direction.tsv"),
               sep = "\t", quote = FALSE)
      }

      # summary counts incl. enrichment direction, every completed comparison present
      if (!is.null(m$summary)) {
        all_comparisons <- m$summary %>%
          select(any_of(c(g1_col, g2_col, "Site", "Comparison", "N1", "N2",
                          paste0("N_Enriched_", g1_col), paste0("N_Enriched_", g2_col)))) %>%
          distinct()
      } else if (nrow(sig) > 0) {
        all_comparisons <- unique(sig[, c(g1_col, g2_col, "Comparison"), drop = FALSE])
      } else next

      if (nrow(sig) > 0 && "enriched_in" %in% names(sig)) {
        dir_counts <- sig %>%
          group_by(.data[[g1_col]], .data[[g2_col]], Comparison) %>%
          summarise(
            N_Significant = n_distinct(feature),
            !!paste0("N_Enriched_", g1_col) :=
              sum(enriched_in == .data[[g1_col]][1], na.rm = TRUE),
            !!paste0("N_Enriched_", g2_col) :=
              sum(enriched_in == .data[[g2_col]][1], na.rm = TRUE),
            .groups = "drop"
          ) %>% as.data.frame()
      } else {
        dir_counts <- all_comparisons
        dir_counts$N_Significant <- 0L
        dir_counts[[paste0("N_Enriched_", g1_col)]] <- 0L
        dir_counts[[paste0("N_Enriched_", g2_col)]] <- 0L
      }
      summary_counts <- merge(all_comparisons, dir_counts,
                              by = c(g1_col, g2_col, "Comparison"),
                              all.x = TRUE, suffixes = c(".cfg", ""))
      if ("N_Significant" %in% names(summary_counts))
        summary_counts$N_Significant[is.na(summary_counts$N_Significant)] <- 0
      for (cc in c(paste0("N_Enriched_", g1_col), paste0("N_Enriched_", g2_col))) {
        if (cc %in% names(summary_counts))
          summary_counts[[cc]][is.na(summary_counts[[cc]])] <- 0
      }

      # (no symmetric reversal table for OVA: the reverse contrast is the
      # same comparison with flipped signs, so it carries no extra info)
      fwrite(summary_counts, file.path(rdir, "summary_diff_counts.tsv"),
             sep = "\t", quote = FALSE)
      cat("[", mname, "] summary_diff_counts.tsv: ", nrow(summary_counts),
          " rows\n", sep = "")
    }
  }
  cat("\n")
}

# ==========================================================================
# Stage 2: build consensus (robust = >=2 methods significant)
# ==========================================================================
if ("consensus" %in% stages) {
  cat("==== [consensus] 3-way intersection flags (robust = n_methods>=2) ====\n")
  per_cmp_dir <- file.path(out_dir, "consensus_per_comparison")
  if (!dir.exists(per_cmp_dir)) dir.create(per_cmp_dir, recursive = TRUE)

  results_by_comp <- list()
  for (mname in names(method_dirs)) {
    ddir <- method_dirs[[mname]]
    if (!dir.exists(ddir)) next
    for (d in list.dirs(ddir, recursive = FALSE, full.names = TRUE)) {
      if (!file.exists(file.path(d, "analysis.done"))) next
      comp <- basename(d)
      sig <- read_main_sig(d, main_res_file, alt_res_files)
      if (!is.null(sig)) {
        if (mname == "maaslin2" && exists("restore_feature_names", mode = "function"))
          sig$feature <- restore_feature_names(as.character(sig$feature))
        results_by_comp[[paste0("ova::", comp)]][[mname]] <- sig
      }
    }
  }

  cat("    Comparisons found:", length(results_by_comp), "\n")
  if (length(results_by_comp) == 0) stop("No completed comparisons found across methods.")

  all_consensus <- list()
  for (key in names(results_by_comp)) {
    comp <- strsplit(key, "::", fixed = TRUE)[[1]][2]
    msets <- results_by_comp[[key]]

    feats <- unique(unlist(lapply(msets, function(x) x$feature)))
    if (length(feats) == 0) next

    flag <- data.frame(feature = feats, stringsAsFactors = FALSE)
    flag$in_linda    <- feats %in% if (!is.null(msets$linda))    msets$linda$feature    else character(0)
    flag$in_maaslin2 <- feats %in% if (!is.null(msets$maaslin2)) msets$maaslin2$feature else character(0)
    flag$in_ancombc2 <- feats %in% if (!is.null(msets$ancombc2)) msets$ancombc2$feature else character(0)
    flag$n_methods <- rowSums(flag[, c("in_linda", "in_maaslin2", "in_ancombc2")])
    flag$all_three <- flag$n_methods == 3
    flag$robust <- flag$n_methods >= 2   # n2 stratification (primary)

    e_list <- list(
      linda    = enrich_by_method(msets$linda, feats),
      maaslin2 = enrich_by_method(msets$maaslin2, feats),
      ancombc2 = enrich_by_method(msets$ancombc2, feats)
    )
    flag$direction_consensus <- direction_consensus_of(e_list, feats)

    flag$design <- "ova"
    flag$comparison <- comp
    fwrite(flag, file.path(per_cmp_dir, paste0("ova__", comp, ".tsv")),
           sep = "\t", quote = FALSE)
    all_consensus[[key]] <- flag

    cat(sprintf("    %s: %d union features, %d robust (>=2/3), %d all-three\n",
                comp, length(feats), sum(flag$robust), sum(flag$all_three)))
  }

  if (length(all_consensus) == 0) stop("No consensus built.")
  all_df <- bind_rows(all_consensus)
  fwrite(all_df, file.path(out_dir, "consensus_all.tsv"), sep = "\t", quote = FALSE)
  robust_df <- all_df[all_df$robust == TRUE, , drop = FALSE]
  fwrite(robust_df, file.path(out_dir, "consensus_robust.tsv"), sep = "\t", quote = FALSE)
  all3_df <- all_df[all_df$all_three == TRUE, , drop = FALSE]
  fwrite(all3_df, file.path(out_dir, "consensus_all3.tsv"), sep = "\t", quote = FALSE)

  cat("\n    Consensus built\n")
  cat("    Total union features (all comparisons): ", nrow(all_df), "\n", sep = "")
  cat("    Robust (>=2/3 methods): ", nrow(robust_df), "\n", sep = "")
  cat("    All three (3/3 methods): ", nrow(all3_df), "\n", sep = "")
  cat("\n")
}

# ==========================================================================
# Stage 3: consensus summary (venn regions per comparison + grand)
# ==========================================================================
if ("summary" %in% stages) {
  cat("==== [summary] venn-region counts per comparison ====\n")
  cons_file <- file.path(out_dir, "consensus_all.tsv")
  if (!file.exists(cons_file))
    stop("Run the consensus stage first; missing ", cons_file)
  df <- fread(cons_file, data.table = FALSE)

  sum_df <- df %>%
    group_by(design, comparison) %>%
    summarise(
      N_union          = dplyr::n(),
      N_linda_only     = sum(in_linda & !in_maaslin2 & !in_ancombc2),
      N_maaslin2_only  = sum(!in_linda & in_maaslin2 & !in_ancombc2),
      N_ancombc2_only  = sum(!in_linda & !in_maaslin2 & in_ancombc2),
      N_linda_maaslin2 = sum(in_linda & in_maaslin2 & !in_ancombc2),
      N_linda_ancombc2 = sum(in_linda & !in_maaslin2 & in_ancombc2),
      N_maaslin2_ancombc2 = sum(!in_linda & in_maaslin2 & in_ancombc2),
      N_all_three      = sum(all_three),
      N_robust_any2    = sum(robust),
      .groups = "drop"
    ) %>% as.data.frame()
  fwrite(sum_df, file.path(out_dir, "consensus_summary.tsv"), sep = "\t", quote = FALSE)
  cat("Per-comparison consensus counts:\n")
  print(sum_df)

  grand <- data.frame(
    n_comparisons = nrow(sum_df),
    total_union_features = nrow(df),
    total_robust_features = sum(df$robust),        # n2 (primary)
    total_all3_features = sum(df$all_three),       # 3/3 (reference)
    mean_robust_per_comparison = round(mean(sum_df$N_robust_any2), 1)
  )
  fwrite(grand, file.path(out_dir, "consensus_grand_summary.tsv"), sep = "\t", quote = FALSE)
  cat("\nGrand summary:\n")
  print(grand)
  cat("\n")
}

# ==========================================================================
# Stage 4: enhanced consensus (robust features + per-method LFC/padj/dir)
# ==========================================================================
if ("enhanced" %in% stages) {
  cat("==== [enhanced] robust consensus with per-method statistics ====\n")
  per_cmp_dir <- file.path(out_dir, "consensus_per_comparison")
  if (!dir.exists(per_cmp_dir)) dir.create(per_cmp_dir, recursive = TRUE)

  keep_cols <- c("feature", "log2FoldChange", "padj", "enriched_in")
  rename_method_cols <- function(df, mname) {
    nm <- names(df)
    nm[nm == "log2FoldChange"] <- paste0(mname, "_lfc")
    nm[nm == "padj"]           <- paste0(mname, "_padj")
    nm[nm == "enriched_in"]    <- paste0(mname, "_enriched_in")
    setNames(df, nm)
  }

  full_by_comp <- list(); sig_by_comp <- list()
  for (mname in names(method_dirs)) {
    ddir <- method_dirs[[mname]]
    if (!dir.exists(ddir)) next
    for (d in list.dirs(ddir, recursive = FALSE, full.names = TRUE)) {
      if (!file.exists(file.path(d, "analysis.done"))) next
      comp <- basename(d)
      key <- paste0("ova::", comp)
      full <- read_full(d, main_res_file, alt_res_files)
      if (is.null(full)) next
      if (mname == "maaslin2" && exists("restore_feature_names", mode = "function"))
        full$feature <- restore_feature_names(as.character(full$feature))
      full_by_comp[[key]][[mname]] <- full
      sig <- full[!is.na(full$padj) & full$padj < alpha, , drop = FALSE]
      if (nrow(sig) > 0) sig_by_comp[[key]][[mname]] <- sig
    }
  }
  cat("    Comparisons found:", length(full_by_comp), "\n")
  if (length(full_by_comp) == 0) stop("No completed comparisons found.")

  all_consensus <- list()
  for (key in names(full_by_comp)) {
    comp <- strsplit(key, "::", fixed = TRUE)[[1]][2]
    fsets <- full_by_comp[[key]]
    ssets <- sig_by_comp[[key]]

    feats <- unique(unlist(lapply(ssets, function(x) x$feature)))
    if (length(feats) == 0) next

    flag <- data.frame(feature = feats, stringsAsFactors = FALSE)
    flag$in_linda    <- feats %in% if (!is.null(ssets$linda))    ssets$linda$feature    else character(0)
    flag$in_maaslin2 <- feats %in% if (!is.null(ssets$maaslin2)) ssets$maaslin2$feature else character(0)
    flag$in_ancombc2 <- feats %in% if (!is.null(ssets$ancombc2)) ssets$ancombc2$feature else character(0)
    flag$n_methods <- rowSums(flag[, c("in_linda", "in_maaslin2", "in_ancombc2")])
    flag$all_three <- flag$n_methods == 3
    flag$robust <- flag$n_methods >= 2   # n2 stratification (primary)

    e_list <- list(
      linda    = enrich_by_method(fsets$linda, feats),
      maaslin2 = enrich_by_method(fsets$maaslin2, feats),
      ancombc2 = enrich_by_method(fsets$ancombc2, feats)
    )
    flag$direction_consensus <- direction_consensus_of(e_list, feats)

    for (mname in c("linda", "maaslin2", "ancombc2")) {
      mfull <- fsets[[mname]]
      if (is.null(mfull)) {
        flag[[paste0(mname, "_lfc")]]          <- NA_real_
        flag[[paste0(mname, "_padj")]]         <- NA_real_
        flag[[paste0(mname, "_enriched_in")]]  <- NA_character_
        next
      }
      j <- mfull[match(feats, mfull$feature),
                 intersect(keep_cols, names(mfull)), drop = FALSE]
      j <- rename_method_cols(j, mname)
      flag <- cbind(flag, j[, setdiff(names(j), "feature"), drop = FALSE])
    }

    flag$design <- "ova"
    flag$comparison <- comp
    fwrite(flag, file.path(per_cmp_dir, paste0("ova__", comp, "_enhanced.tsv")),
           sep = "\t", quote = FALSE)
    all_consensus[[key]] <- flag
    cat(sprintf("    %s: %d union, %d robust (>=2/3), %d all-three\n",
                comp, length(feats), sum(flag$robust), sum(flag$all_three)))
  }

  if (length(all_consensus) == 0) stop("No consensus built.")
  all_df <- bind_rows(all_consensus)
  fwrite(all_df, file.path(out_dir, "consensus_all_enhanced.tsv"), sep = "\t", quote = FALSE)
  robust_df <- all_df[all_df$robust == TRUE, , drop = FALSE]
  fwrite(robust_df, file.path(out_dir, "consensus_robust_enhanced.tsv"), sep = "\t", quote = FALSE)
  all3_df <- all_df[all_df$all_three == TRUE, , drop = FALSE]
  fwrite(all3_df, file.path(out_dir, "consensus_all3_enhanced.tsv"), sep = "\t", quote = FALSE)

  cat("\n    Enhanced consensus built\n")
  cat("    Total union features: ", nrow(all_df), "\n", sep = "")
  cat("    Robust (>=2/3):       ", nrow(robust_df), "\n", sep = "")
  cat("    All three (3/3):      ", nrow(all3_df), "\n", sep = "")
  cat("\n")
}

# ==========================================================================
# Stage 5: full feature table (no padj filter; per-method stats side by side)
# ==========================================================================
if ("fulltable" %in% stages) {
  cat("==== [fulltable] all features, 3-method padj side by side ====\n")
  per_cmp_dir <- file.path(out_dir, "full_table_per_comparison")
  if (!dir.exists(per_cmp_dir)) dir.create(per_cmp_dir, recursive = TRUE)

  keep_cols <- c("feature", "log2FoldChange", "padj", "enriched_in")
  rename_method_cols <- function(df, mname) {
    nm <- names(df)
    nm[nm == "log2FoldChange"] <- paste0(mname, "_lfc")
    nm[nm == "padj"]           <- paste0(mname, "_padj")
    nm[nm == "enriched_in"]    <- paste0(mname, "_enriched_in")
    setNames(df, nm)
  }

  full_by_comp <- list()
  for (mname in names(method_dirs)) {
    ddir <- method_dirs[[mname]]
    if (!dir.exists(ddir)) next
    for (d in list.dirs(ddir, recursive = FALSE, full.names = TRUE)) {
      if (!file.exists(file.path(d, "analysis.done"))) next
      comp <- basename(d)
      full <- read_full(d, main_res_file, alt_res_files)
      if (is.null(full)) next
      if (mname == "maaslin2" && exists("restore_feature_names", mode = "function"))
        full$feature <- restore_feature_names(as.character(full$feature))
      full_by_comp[[paste0("ova::", comp)]][[mname]] <- full
    }
  }
  cat("    Comparisons found:", length(full_by_comp), "\n")
  if (length(full_by_comp) == 0) stop("No completed comparisons found.")

  all_tables <- list()
  for (key in names(full_by_comp)) {
    comp <- strsplit(key, "::", fixed = TRUE)[[1]][2]
    fsets <- full_by_comp[[key]]

    feats <- unique(unlist(lapply(fsets, function(x) x$feature)))
    if (length(feats) == 0) next
    out <- data.frame(feature = feats, genus = make_genus_short(feats),
                      design = "ova", comparison = comp,
                      stringsAsFactors = FALSE)

    for (mname in c("linda", "maaslin2", "ancombc2")) {
      mfull <- fsets[[mname]]
      if (is.null(mfull)) {
        out[[paste0(mname, "_lfc")]]          <- NA_real_
        out[[paste0(mname, "_padj")]]         <- NA_real_
        out[[paste0(mname, "_enriched_in")]]  <- NA_character_
        next
      }
      j <- mfull[match(feats, mfull$feature),
                 intersect(keep_cols, names(mfull)), drop = FALSE]
      j <- rename_method_cols(j, mname)
      out <- cbind(out, j[, setdiff(names(j), "feature"), drop = FALSE])
    }

    out$linda_sig    <- as.integer(!is.na(out$linda_padj)    & out$linda_padj    < alpha)
    out$maaslin2_sig <- as.integer(!is.na(out$maaslin2_padj) & out$maaslin2_padj < alpha)
    out$ancombc2_sig <- as.integer(!is.na(out$ancombc2_padj) & out$ancombc2_padj < alpha)
    out$n_methods_sig <- out$linda_sig + out$maaslin2_sig + out$ancombc2_sig
    out$all_three <- as.integer(out$n_methods_sig == 3)
    out$robust <- as.integer(out$n_methods_sig >= 2)

    fwrite(out, file.path(per_cmp_dir, paste0("ova__", comp, ".tsv")),
           sep = "\t", quote = FALSE)
    all_tables[[key]] <- out
    cat(sprintf("    %s: %d features (%d robust >=2/3, %d all-three)\n",
                comp, nrow(out), sum(out$robust), sum(out$all_three)))
  }

  if (length(all_tables) == 0) stop("Nothing built.")
  all_df <- bind_rows(all_tables)
  fwrite(all_df, file.path(out_dir, "full_table.tsv"), sep = "\t", quote = FALSE)

  cat("\n    Full table built\n")
  cat("    Total rows (all features, all comparisons): ", nrow(all_df), "\n", sep = "")
  cat("    Robust (>=2/3 sig): ", sum(all_df$robust), "\n", sep = "")
  cat("    All three (3/3 sig): ", sum(all_df$all_three), "\n", sep = "")
  cat("\n")
}

cat("==========================================\n")
cat("Merge + consensus complete —", scope_tag, "\n")
cat("End: ", format(Sys.time()), "\n")
cat("==========================================\n")