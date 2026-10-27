# Shared genus selection and font settings for OVA figures.
# Select a real italic font: the cluster's default sans family can resolve
# both regular and italic requests to the same upright font file.
ova_genus_font <- function() {
  fonts <- systemfonts::system_fonts()
  candidates <- c("Arial", "Liberation Sans", "DejaVu Sans", "Ubuntu", "Nimbus Sans L")
  for (family in candidates) {
    italic_fonts <- fonts[fonts$family == family & fonts$italic, , drop = FALSE]
    if (nrow(italic_fonts) == 0L) next
    matched <- systemfonts::match_fonts(family, italic = TRUE)
    if (matched$path %in% italic_fonts$path) {
      return(family)
    }
  }
  stop("OVA plots require a sans-serif family with an actual italic font installed")
}
OVA_GENUS_FONT <- ova_genus_font()
OVA_BLOCK_GAP <- 0.6

# A genus must peak uniquely in its target across all plotted groups, and
# have a significant positive MaAsLin2 coefficient plus support from at least
# min_methods - 1 other significant positive methods. Select by coefficient;
# display within each target by descending raw mean abundance (including zeros).
select_target_topn <- function(enh, min_methods, mean_abundance_df, n_each = 10L,
                               group_levels, alpha = 0.05) {
  stopifnot(length(n_each) == 1L, is.finite(n_each), n_each >= 1L)
  methods <- c("linda", "maaslin2", "ancombc2")
  positive_sig <- vapply(methods, function(method) {
    q <- enh[[paste0(method, "_padj")]]
    effect <- enh[[paste0(method, "_lfc")]]
    is.finite(q) & q < alpha & is.finite(effect) & effect > 0
  }, logical(nrow(enh)))
  # vapply simplifies a single comparison row to a vector; retain a matrix.
  dim(positive_sig) <- c(nrow(enh), length(methods))
  colnames(positive_sig) <- methods
  enh$n_positive_sig_methods <- rowSums(positive_sig)
  keep <- enh[
    positive_sig[, "maaslin2"] &
      enh$n_positive_sig_methods >= min_methods & !is.na(enh$target), ,
    drop = FALSE
  ]
  means <- mean_abundance_df %>%
    select(feature, target, target_mean_abundance)
  if (anyDuplicated(means[c("feature", "target")]) ||
    any(!is.finite(means$target_mean_abundance))) {
    stop("Global peak selection requires one finite mean per genus and group")
  }
  peaks <- means %>%
    group_by(feature) %>%
    mutate(
      is_global_peak = target_mean_abundance == max(target_mean_abundance, -Inf),
      n_peak_ties = sum(is_global_peak),
      complete_groups = setequal(target, group_levels),
      peak_target = if (sum(is_global_peak) == 1L) {
        target[which.max(target_mean_abundance)]
      } else {
        NA_character_
      }
    ) %>%
    ungroup()
  if (any(!peaks$complete_groups)) {
    stop("Global peak selection requires mean abundance for every plotted group")
  }
  keep <- keep %>% left_join(peaks, by = c("feature", "target"))
  if (any(!is.finite(keep$target_mean_abundance))) {
    stop("Missing mean abundance for an eligible target-enriched genus")
  }
  candidates <- keep %>%
    group_by(feature) %>%
    mutate(
      n_enriched_targets = n(),
      passes_peak_filter = is_global_peak & n_peak_ties == 1L
    ) %>%
    ungroup()
  assigned <- candidates %>%
    filter(passes_peak_filter) %>%
    group_by(target) %>%
    arrange(desc(maaslin2_lfc), maaslin2_padj, feature, .by_group = TRUE) %>%
    mutate(coef_rank_in_target = row_number(), selected = coef_rank_in_target <= n_each) %>%
    ungroup()
  selection <- assigned %>%
    filter(selected) %>%
    arrange(
      match(target, group_levels), desc(target_mean_abundance),
      desc(maaslin2_lfc), feature
    ) %>%
    group_by(target) %>%
    mutate(display_rank_in_target = row_number()) %>%
    ungroup()
  ordered_rows <- selection$feature
  attr(ordered_rows, "block") <- selection$target
  attr(ordered_rows, "selection") <- selection %>%
    select(feature, target, target_mean_abundance,
      maaslin2_coef = maaslin2_lfc,
      maaslin2_padj, n_methods, n_positive_sig_methods, peak_target,
      n_enriched_targets, coef_rank_in_target, display_rank_in_target
    )
  attr(ordered_rows, "candidates") <- candidates %>%
    select(feature, target, target_mean_abundance,
      maaslin2_coef = maaslin2_lfc,
      maaslin2_padj, n_methods, n_positive_sig_methods, peak_target,
      n_enriched_targets, n_peak_ties, passes_peak_filter
    ) %>%
    left_join(assigned %>% select(feature, target, coef_rank_in_target, selected),
      by = c("feature", "target")
    ) %>%
    mutate(selected = tidyr::replace_na(selected, FALSE))
  ordered_rows
}
