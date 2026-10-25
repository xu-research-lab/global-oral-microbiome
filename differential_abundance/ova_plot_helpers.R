# Shared only by the age, site and continent OVA visualization scripts.
# Select a real italic font: the cluster's default sans family can resolve
# both regular and italic requests to the same upright font file.
ova_genus_font <- function() {
  fonts <- systemfonts::system_fonts()
  candidates <- c("Arial", "Liberation Sans", "DejaVu Sans", "Ubuntu", "Nimbus Sans L")
  for (family in candidates) {
    italic_fonts <- fonts[fonts$family == family & fonts$italic, , drop = FALSE]
    if (nrow(italic_fonts) == 0L) next
    matched <- systemfonts::match_fonts(family, italic = TRUE)
    if (matched$path %in% italic_fonts$path) return(family)
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
  keep <- enh[positive_sig[, "maaslin2"] &
                enh$n_positive_sig_methods >= min_methods & !is.na(enh$target), ,
              drop = FALSE]
  means <- mean_abundance_df %>%
    select(feature, target, target_mean_abundance)
  if (anyDuplicated(means[c("feature", "target")]) ||
      any(!is.finite(means$target_mean_abundance)))
    stop("Global peak selection requires one finite mean per genus and group")
  peaks <- means %>%
    group_by(feature) %>%
    mutate(is_global_peak = target_mean_abundance == max(target_mean_abundance, -Inf),
           n_peak_ties = sum(is_global_peak),
           complete_groups = setequal(target, group_levels),
           peak_target = if (sum(is_global_peak) == 1L)
             target[which.max(target_mean_abundance)] else NA_character_) %>%
    ungroup()
  if (any(!peaks$complete_groups))
    stop("Global peak selection requires mean abundance for every plotted group")
  keep <- keep %>% left_join(peaks, by = c("feature", "target"))
  if (any(!is.finite(keep$target_mean_abundance)))
    stop("Missing mean abundance for an eligible target-enriched genus")
  candidates <- keep %>%
    group_by(feature) %>%
    mutate(n_enriched_targets = n(),
           passes_peak_filter = is_global_peak & n_peak_ties == 1L) %>%
    ungroup()
  assigned <- candidates %>%
    filter(passes_peak_filter) %>%
    group_by(target) %>%
    arrange(desc(maaslin2_lfc), maaslin2_padj, feature, .by_group = TRUE) %>%
    mutate(coef_rank_in_target = row_number(), selected = coef_rank_in_target <= n_each) %>%
    ungroup()
  selection <- assigned %>%
    filter(selected) %>%
    arrange(match(target, group_levels), desc(target_mean_abundance),
            desc(maaslin2_lfc), feature) %>%
    group_by(target) %>% mutate(display_rank_in_target = row_number()) %>% ungroup()
  ordered_rows <- selection$feature
  attr(ordered_rows, "block") <- selection$target
  attr(ordered_rows, "selection") <- selection %>%
    select(feature, target, target_mean_abundance, maaslin2_coef = maaslin2_lfc,
           maaslin2_padj, n_methods, n_positive_sig_methods, peak_target,
           n_enriched_targets, coef_rank_in_target, display_rank_in_target)
  attr(ordered_rows, "candidates") <- candidates %>%
    select(feature, target, target_mean_abundance, maaslin2_coef = maaslin2_lfc,
           maaslin2_padj, n_methods, n_positive_sig_methods, peak_target,
           n_enriched_targets, n_peak_ties, passes_peak_filter) %>%
    left_join(assigned %>% select(feature, target, coef_rank_in_target, selected),
              by = c("feature", "target")) %>%
    mutate(selected = replace_na(selected, FALSE))
  ordered_rows
}

# Count unique group-mean extrema with matching significant method directions.
ova_extreme_counts <- function(cons, means, group_levels, tested_targets,
                               min_methods = 2L, alpha = 0.05) {
  stopifnot(!anyDuplicated(means[c("feature", "target")]),
            all(is.finite(means$target_mean_abundance)))
  extrema <- means %>% group_by(feature) %>%
    mutate(complete_groups = setequal(target, group_levels),
           is_max = target_mean_abundance == max(target_mean_abundance),
           is_min = target_mean_abundance == min(target_mean_abundance),
           n_max = sum(is_max), n_min = sum(is_min),
           direction = case_when(is_max & n_max == 1L & !is_min ~ "Enriched",
                                 is_min & n_min == 1L & !is_max ~ "Depleted",
                                 TRUE ~ NA_character_)) %>% ungroup()
  stopifnot(all(extrema$complete_groups))
  detail <- cons %>% filter(n_methods >= min_methods, target %in% tested_targets) %>%
    distinct(feature, target, .keep_all = TRUE) %>%
    left_join(extrema, by = c("feature", "target"))
  stopifnot(all(is.finite(detail$target_mean_abundance)))
  methods <- c("linda", "maaslin2", "ancombc2")
  for (sign_name in c("positive", "negative")) {
    sig <- vapply(methods, function(method) {
      q <- detail[[paste0(method, "_padj")]]
      effect <- detail[[paste0(method, "_lfc")]]
      is.finite(q) & q < alpha & is.finite(effect) &
        if (sign_name == "positive") effect > 0 else effect < 0
    }, logical(nrow(detail)))
    dim(sig) <- c(nrow(detail), length(methods))
    detail[[paste0("n_", sign_name, "_sig_methods")]] <- rowSums(sig)
    detail[[paste0("maaslin2_", sign_name, "_sig")]] <- sig[, 2L]
  }
  detail <- detail %>% mutate(
    passes_direction_filter = case_when(
      direction == "Enriched" ~ maaslin2_positive_sig & n_positive_sig_methods >= min_methods,
      direction == "Depleted" ~ maaslin2_negative_sig & n_negative_sig_methods >= min_methods,
      TRUE ~ FALSE))
  counts <- detail %>% filter(passes_direction_filter) %>%
    count(target, direction, name = "n_genera") %>%
    right_join(expand_grid(target = tested_targets,
                           direction = c("Enriched", "Depleted")),
               by = c("target", "direction")) %>%
    mutate(n_genera = replace_na(n_genera, 0L))
  list(counts = counts, detail = detail)
}

ova_count_plot <- function(counts, group_order, font_size = 7.5) {
  counts$direction <- factor(counts$direction, levels = c("Depleted", "Enriched"))
  ggplot(counts, aes(n_genera, factor(target, levels = rev(group_order)),
                    fill = direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = n_genera), position = position_dodge(width = 0.8),
              hjust = -0.2, family = OVA_GENUS_FONT, size = font_size / .pt) +
    scale_fill_manual(values = c(Enriched = "#B2182B", Depleted = "#2166AC"),
                      breaks = c("Enriched", "Depleted"), drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = "Number of genera", y = NULL, fill = NULL) +
    theme_classic(base_family = OVA_GENUS_FONT, base_size = font_size) +
    theme(legend.position = "top", legend.key.size = unit(2.5, "mm"),
          legend.margin = margin(0, 0, 0, 0),
          axis.text = element_text(colour = "black"),
          axis.line.y = element_blank(), axis.ticks.y = element_blank())
}
