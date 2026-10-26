# Age three-panel figure with the selected genera in their saved row order.
render_age <- function(heat, bubble, prev) {
  mean_by_group <- heat
  row_table <- heat %>%
    distinct(genus_lab, y_plot) %>%
    arrange(y_plot)
  row_levels <- row_table$genus_lab
  row_y <- setNames(row_table$y_plot, row_levels)
  y_limits <- range(row_y) + c(-0.5, 0.5)
  fill_col <- "z_abun"
  fill_name <- "mean relative abundance\n(z-score)"
  fill_lims <- sym_limits(mean_by_group$z_abun)
  mean_by_group$genus_lab <- factor(mean_by_group$genus_lab, levels = row_levels)
  mean_heatmap <- ggplot(mean_by_group, aes(x = x_plot, y = y_plot, fill = .data[[fill_col]])) +
    geom_tile(
      color = NA,
      width = 0.9, height = 1
    ) +
    scale_fill_viridis_c(
      option = "mako", direction = -1, limits = fill_lims,
      na.value = "grey85", name = fill_name
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(
        size = 8,
        face = "italic", family = OVA_GENUS_FONT, hjust = 1
      ), axis.text.x = element_text(
        size = 8, angle = 60,
        hjust = 1, face = "bold"
      ), axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
      legend.position = "right", legend.direction = "vertical", legend.title = element_text(
        size = 9,
        face = "bold"
      ), legend.text = element_text(size = 8), legend.spacing = unit(0, "cm"), legend.box.spacing = unit(
        0,
        "cm"
      ), legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.25, "cm"), panel.border = element_rect(
        fill = NA,
        colour = "black", linewidth = 0.3
      ), plot.margin = margin(t = 5.5, r = 2, b = 5.5, l = 2)
    ) +
    scale_x_continuous(
      breaks = seq_along(age_levels), labels = unname(age_display[age_levels]),
      limits = c(0.55, length(age_levels) + 0.45), expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = unname(row_y),
      labels = row_levels, limits = y_limits, expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off")
  prev_df <- prev
  prev_df$genus_lab <- factor(prev_df$genus_lab, levels = row_levels)
  prev_df$age_lab <- factor(prev_df$age_lab, levels = age_levels)
  prev_ranges <- prev_df %>%
    group_by(genus_lab, y_plot) %>%
    summarise(
      prev_min = min(prevalence),
      prev_max = max(prevalence), .groups = "drop"
    )
  prevalence_plot <- ggplot(prev_df, aes(x = prevalence, y = y_jit)) +
    geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "grey45", linewidth = 0.4
    ) +
    geom_segment(data = prev_ranges, aes(
      x = prev_min, xend = prev_max,
      y = y_plot, yend = y_plot
    ), inherit.aes = FALSE, colour = "grey75", linewidth = 0.3) +
    geom_point(aes(fill = age_lab),
      shape = 21, size = 1.6, stroke = 0.15, alpha = 0.6
    ) +
    labs(x = "Prevalence (%)", y = NULL) +
    scale_fill_manual(values = age_colors, labels = age_short, name = "Age group", guide = guide_legend(override.aes = list(
      size = 2.4,
      alpha = 1
    ))) +
    scale_x_continuous(limits = c(-8, 108), breaks = c(0, 50, 100), expand = c(
      0,
      0
    )) +
    scale_y_continuous(limits = y_limits, breaks = unname(row_y), labels = row_levels, expand = c(
      0,
      0
    )) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(), axis.text.y = element_text(
        size = 8,
        face = "italic", family = OVA_GENUS_FONT, color = "black", hjust = 1
      ), axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 7, face = "bold", margin = margin(t = 3)), axis.text.x = element_text(size = 6),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3), legend.position = "right",
      legend.direction = "vertical", legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8),
      legend.spacing = unit(0, "cm"), legend.box.spacing = unit(0, "cm"), legend.key.width = unit(
        0.2,
        "cm"
      ), legend.key.height = unit(0.25, "cm"), plot.margin = margin(
        t = 5.5, r = 2, b = 5.5,
        l = 2
      )
    ) +
    coord_cartesian(clip = "off")
  sub <- bubble
  sub$pair_lab <- factor(sub$pair_lab, levels = short_target(age_levels))
  diff_bubble <- NULL
  if (nrow(sub) > 0) {
    diff_bubble <- ggplot(sub, aes(x = pair_lab, y = y_plot, size = neg_log10_maaq, fill = maaslin2_coef)) +
      geom_point(shape = 21, alpha = 0.9, color = "black", stroke = 0.4) +
      scale_size_continuous(range = c(
        0.8,
        3.2
      ), name = "-log10(FDR)") +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white",
        high = "#B2182B", midpoint = 0, na.value = "grey80", name = "MaAsLin2 coefficient\n(Target vs. Rest)"
      ) +
      scale_x_discrete(expand = expansion(add = 0.55)) +
      scale_y_continuous(
        breaks = unname(row_y),
        labels = row_levels, limits = y_limits, expand = c(0, 0)
      ) +
      labs(x = NULL, y = NULL) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(
          size = 8, face = "italic", family = OVA_GENUS_FONT,
          color = "black", hjust = 1
        ), axis.text.x = element_text(
          size = 8, angle = 60, hjust = 1,
          face = "bold"
        ), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(
          fill = NA, colour = "black",
          linewidth = 0.3
        ), legend.position = "right", legend.direction = "vertical", legend.title = element_text(
          size = 9,
          face = "bold"
        ), legend.text = element_text(size = 8), legend.spacing = unit(0, "cm"),
        legend.box.spacing = unit(0, "cm"), legend.key.width = unit(0.2, "cm"), legend.key.height = unit(
          0.25,
          "cm"
        ), legend.box = "vertical", plot.margin = margin(t = 5.5, r = 2, b = 5.5, l = 2)
      ) +
      coord_cartesian(clip = "on")
  }
  prev_panel <- prevalence_plot + theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_text(hjust = 0)
  )
  prev_panel <- patchwork::free(prev_panel, type = "label", side = "b")
  combined <- if (is.null(diff_bubble)) {
    (mean_heatmap | prev_panel) + patchwork::plot_layout(widths = c(6, 3), guides = "collect")
  } else {
    (diff_bubble | mean_heatmap + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) |
      prev_panel) + patchwork::plot_layout(widths = c(6, 6, 3), guides = "collect")
  }
  combined
}
