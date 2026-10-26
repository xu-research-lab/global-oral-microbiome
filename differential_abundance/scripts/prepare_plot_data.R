# Recompute the selected genera and all three panel tables from consensus + inputs.
# Plotting itself only needs the saved tables; this stage reads the source data.
if (analysis != "sex") {
  sys.source(file.path(script_dir, "plot_common.R"), envir = environment())
  sys.source(file.path(script_dir, "selection.R"), envir = environment())
  enh <- fread(file.path(out_dir, "consensus_all_enhanced.tsv"), data.table = FALSE)
  parsed <- lapply(enh$comparison, parse_ova_comparison)
  enh$target <- vapply(parsed, `[[`, "", "target")
  enh$site <- if (analysis == "age") vapply(parsed, `[[`, "", "site") else enh$target
  enh$pair_lab <- short_target(enh$target)
  enh$genus_lab <- make.unique(vapply(enh$feature, genus_label, character(1)))
  meta <- load_age_metadata(hc_only = (scope == "hc"))
  meta <- meta[meta[[main_col]] %in% group_order, , drop = FALSE]
  meta$group_lab <- if (analysis == "age") paste(meta[[COL_SITE]], meta[[COL_AGE]]) else meta[[main_col]]
  plotted_groups <- group_order[group_order %in% meta[[main_col]]]
  group_lab_order <- if (analysis == "age") paste("saliva", plotted_groups) else plotted_groups
  robust_feats <- unique(enh$feature[enh$n_methods >= 2 &
    if (analysis == "age") !is.na(enh$maaslin2_padj) else TRUE])
  rel <- read_relabund_rows(robust_feats)
  samples <- intersect(colnames(rel), rownames(meta))
  sdf <- as.data.frame(t(rel[, samples]))
  colnames(sdf) <- rel[[1]]
  means <- do.call(rbind, lapply(plotted_groups, function(tgt) {
    idx <- which(meta[rownames(sdf), main_col] == tgt)
    m <- colMeans(sdf[idx, , drop = FALSE])
    data.frame(feature = names(m), target = tgt, target_mean_abundance = as.numeric(m))
  }))
  feats <- select_target_topn(enh, 2L, means, 10L, group_levels = plotted_groups, alpha = DA_ALPHA)
  if (!length(feats)) stop("No eligible target-enriched genera")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  write_plot_table <- function(df, name) {
    fwrite(df, file.path(
      figures_dir,
      paste0("figureB_", name, "_2method.tsv")
    ), sep = "\t")
  }
  write_plot_table(attr(feats, "selection"), "selected")
  write_plot_table(attr(feats, "candidates"), "candidates")
  counts <- data.frame(target = group_order) %>%
    left_join(attr(feats, "candidates") %>% filter(passes_peak_filter) %>%
      count(target, name = "eligible_genera"), by = "target") %>%
    left_join(attr(feats, "selection") %>% count(target, name = "selected_genera"), by = "target") %>%
    mutate(across(c(eligible_genera, selected_genera), ~ replace_na(.x, 0L)))
  write_plot_table(counts, "selection_counts")
  sdf <- sdf[, feats, drop = FALSE]
  sdf$group_lab <- meta[rownames(sdf), "group_lab"]
  group_n <- sdf %>%
    group_by(group_lab) %>%
    summarise(n_total = n(), .groups = "drop")
  heat <- sdf %>%
    pivot_longer(-group_lab, names_to = "feature", values_to = "abundance") %>%
    group_by(group_lab, feature) %>%
    summarise(mean_abundance = mean(abundance), n_present = sum(abundance > 0), .groups = "drop") %>%
    left_join(group_n, by = "group_lab") %>%
    group_by(feature) %>%
    mutate(
      z_abun = {
        x <- mean_abundance
        s <- sd(x)
        if (is.finite(s) && s > 0) (x - mean(x)) / s else NA_real_
      },
      raw_abun = ifelse(mean_abundance > 0, log10(mean_abundance), NA_real_)
    ) %>%
    ungroup()
  lab_map <- setNames(enh$genus_lab[match(feats, enh$feature)], feats)
  row_levels <- rev(unname(lab_map))
  bottom_blocks <- rev(attr(feats, "block"))
  block_index <- cumsum(c(TRUE, tail(bottom_blocks, -1L) != head(bottom_blocks, -1L)))
  row_y <- setNames(seq_along(row_levels) + (block_index - 1L) * OVA_BLOCK_GAP, row_levels)
  heat$genus_lab <- factor(lab_map[heat$feature], levels = row_levels)
  heat$group_lab <- factor(heat$group_lab, levels = group_lab_order)
  heat$y_plot <- unname(row_y[as.character(heat$genus_lab)])
  if (analysis == "age") heat$age_lab <- factor(sub("^\\S+ ", "", heat$group_lab), levels = plotted_groups)
  heat$x_plot <- as.numeric(heat$group_lab)
  prev <- heat %>%
    filter(n_total > 0) %>%
    mutate(prevalence = 100 * n_present / n_total) %>%
    arrange(genus_lab, group_lab)
  set.seed(42)
  prev$y_jit <- prev$y_plot + runif(nrow(prev), -0.15, 0.15)
  bubble <- enh[enh$feature %in% feats & !is.na(enh$target), , drop = FALSE]
  bubble$maaslin2_coef <- bubble$maaslin2_lfc
  bubble$genus_lab <- factor(lab_map[bubble$feature], levels = row_levels)
  bubble$y_plot <- unname(row_y[as.character(bubble$genus_lab)])
  sig_mat <- as.matrix(bubble[, c("in_linda", "in_maaslin2", "in_ancombc2")])
  bubble$maa_q <- bubble$maaslin2_padj
  bubble <- bubble[is.finite(bubble$maa_q) & is.finite(bubble$maaslin2_coef) & rowSums(sig_mat) >= 2L, , drop = FALSE]
  bubble$neg_log10_maaq <- -log10(pmax(bubble$maa_q, .Machine$double.xmin))
  write_plot_table(heat, "heatmap_data")
  write_plot_table(bubble, "bubble_data")
  write_plot_table(prev, "prevalence_data")
}
