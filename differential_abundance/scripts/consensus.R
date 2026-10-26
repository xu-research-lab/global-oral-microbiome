# Read each completed comparison once, then build the tables used by the figures.
read_main <- function(directory) {
  path <- file.path(directory, main_res_file)
  if (!file.exists(path)) stop("Completed fit is missing ", path)
  df <- fread(path, data.table = FALSE)
  if (!all(c("feature", "log2FoldChange", "padj", "enriched_in") %in% names(df))) {
    stop("Invalid result schema: ", path)
  }
  df
}

direction_consensus <- function(sets, feats) {
  enriched <- lapply(sets, function(df) as.character(df$enriched_in[match(feats, df$feature)]))
  vapply(seq_along(feats), function(i) {
    values <- vapply(enriched, `[`, "", i)
    values <- values[!is.na(values) & values != "" & values != "none"]
    if (!length(values)) {
      return("none")
    }
    tab <- sort(table(values), decreasing = TRUE)
    if (tab[1] >= 2) paste0("agreed:", names(tab)[1]) else "discordant"
  }, "")
}

method_stats <- function(feats, sets) {
  out <- data.frame(row.names = seq_along(feats))
  for (m in names(method_dirs)) {
    df <- sets[[m]]
    idx <- if (!is.null(df)) match(feats, df$feature) else rep(NA_integer_, length(feats))
    out[[paste0(m, "_lfc")]] <- if (is.null(df)) NA_real_ else df$log2FoldChange[idx]
    out[[paste0(m, "_padj")]] <- if (is.null(df)) NA_real_ else df$padj[idx]
    out[[paste0(m, "_enriched_in")]] <- if (is.null(df)) NA_character_ else as.character(df$enriched_in[idx])
  }
  out
}

consensus_flags <- function(feats, significant, directions, comparison) {
  flag <- data.frame(feature = feats)
  for (m in names(method_dirs)) flag[[paste0("in_", m)]] <- feats %in% significant[[m]]$feature
  flag$n_methods <- rowSums(flag[, paste0("in_", names(method_dirs))])
  flag$all_three <- flag$n_methods == 3L
  flag$robust <- flag$n_methods >= 2L
  flag$direction_consensus <- direction_consensus(directions, feats)
  flag
}

by_comparison <- list()
for (m in names(method_dirs)) {
  for (d in list.dirs(method_dirs[[m]], recursive = FALSE, full.names = TRUE)) {
    if (!file.exists(file.path(d, "analysis.done"))) next
    df <- read_main(d)
    if (m == "maaslin2") df$feature <- restore_feature_names(as.character(df$feature))
    by_comparison[[basename(d)]][[m]] <- df
  }
}
if (!length(by_comparison)) stop("No completed comparisons found")
design <- if (analysis == "sex") "pairwise" else "ova"
enhanced <- list()
full_tables <- list()
robust <- list()
for (comparison in names(by_comparison)) {
  sets <- by_comparison[[comparison]]
  significant <- lapply(sets, function(df) df[!is.na(df$padj) & df$padj < DA_ALPHA, , drop = FALSE])
  significant <- Filter(function(df) nrow(df) > 0L, significant)
  sig_features <- unique(unlist(lapply(significant, `[[`, "feature")))
  if (length(sig_features)) {
    flag <- consensus_flags(sig_features, significant, sets, comparison)
    enhanced[[comparison]] <- data.frame(flag, method_stats(sig_features, sets), design = design, comparison = comparison)
    if (analysis == "sex") {
      # Sex bar directions use votes from significant methods, as in the model consensus.
      flag <- consensus_flags(sig_features, significant, significant, comparison)
      flag$design <- design
      flag$comparison <- comparison
      robust[[comparison]] <- flag[flag$robust, , drop = FALSE]
    }
  }
  feats <- unique(unlist(lapply(sets, `[[`, "feature")))
  out <- data.frame(
    feature = feats, genus = make_genus_short(feats), design = design,
    comparison = comparison, method_stats(feats, sets)
  )
  for (m in names(method_dirs)) out[[paste0(m, "_sig")]] <- as.integer(!is.na(out[[paste0(m, "_padj")]]) & out[[paste0(m, "_padj")]] < DA_ALPHA)
  out$n_methods_sig <- rowSums(out[, paste0(names(method_dirs), "_sig")])
  out$all_three <- as.integer(out$n_methods_sig == 3L)
  out$robust <- as.integer(out$n_methods_sig >= 2L)
  full_tables[[comparison]] <- out
}
if (!length(enhanced)) stop("No significant features found")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(bind_rows(enhanced), file.path(out_dir, "consensus_all_enhanced.tsv"), sep = "\t", quote = FALSE)
fwrite(bind_rows(full_tables), file.path(out_dir, "full_table.tsv"), sep = "\t", quote = FALSE)
if (analysis == "sex") fwrite(bind_rows(robust), file.path(out_dir, "consensus_robust.tsv"), sep = "\t", quote = FALSE)
message("Saved consensus for ", length(by_comparison), " comparisons in ", out_dir)
