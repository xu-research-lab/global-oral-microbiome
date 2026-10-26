# Shared publication layout for the site and continent three-panel figures.
render_ova <- function(heat, bubble, prev) {
  rows <- heat %>%
    distinct(feature, genus_lab, y_plot) %>%
    arrange(y_plot)
  width_mm <- 120
  detail_height_mm <- 205
  old_theme <- theme_get()
  on.exit(theme_set(old_theme), add = TRUE)
  # Match regular and italic text within the same installed sans-serif family.
  font_family <- OVA_GENUS_FONT
  font_body <- 7.5
  font_axis <- 7
  font_legend <- 7
  theme_set(theme_minimal(base_family = font_family, base_size = font_body) +
    theme(
      legend.position = "right", legend.box = "vertical",
      legend.direction = "vertical", legend.title.position = "top",
      legend.box.just = "left", legend.spacing.y = unit(2, "mm")
    ))
  colors <- setNames(viridis(length(group_order)), group_order)
  group_labels <- group_order
  legend_labels <- group_order
  y_limits <- range(rows$y_plot) + c(-0.5, 0.5)
  bubble$target <- factor(bubble$target, levels = group_order)
  prev$group_lab <- factor(prev$group_lab, levels = group_order)
  fill_col <- "z_abun"
  m <- max(abs(heat$z_abun), na.rm = TRUE)
  fill_limits <- c(-m, m)
  fill_title <- "Mean relative\nabundance\n(z-score)"

  panel_theme <- theme_minimal(base_family = font_family, base_size = font_body) +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(colour = "black"),
      axis.text.y = element_text(
        size = font_body, face = "italic", hjust = 1,
        margin = margin(r = 2)
      ),
      axis.text.x = element_text(size = font_axis, angle = 60, hjust = 1, vjust = 1),
      axis.ticks = element_blank(), panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.25),
      legend.title = element_text(size = font_legend, face = "plain"),
      legend.position = "right", legend.box = "vertical",
      legend.direction = "vertical", legend.title.position = "top",
      legend.text = element_text(size = font_legend),
      legend.key.width = unit(2.4, "mm"),
      legend.key.height = unit(2.5, "mm"),
      legend.spacing.x = unit(1.2, "mm"),
      legend.spacing.y = unit(0.5, "mm"),
      legend.margin = margin(0, 0.5, 0, 0.5, unit = "mm"),
      legend.box.spacing = unit(1.5, "mm"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm")
    )
  row_scale <- function(labels = NULL) {
    scale_y_continuous(
      limits = y_limits, breaks = rows$y_plot,
      labels = labels, expand = c(0, 0)
    )
  }
  colorbar <- function(order) {
    guide_colorbar(
      order = order, direction = "vertical",
      theme = theme(
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(12, "mm")
      )
    )
  }

  bubble_plot <- ggplot(bubble, aes(target, y_plot,
    size = neg_log10_maaq,
    fill = maaslin2_coef
  )) +
    geom_point(shape = 21, alpha = 0.9, colour = "black", stroke = 0.22) +
    scale_size_continuous(
      range = c(0.45, 1.8), name = "-log10(FDR)",
      guide = guide_legend(order = 2, ncol = 1)
    ) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "MaAsLin2 coefficient\n(Target vs. Rest)",
      guide = colorbar(1)
    ) +
    scale_x_discrete(
      labels = paste0(group_labels, " vs rest"), drop = FALSE,
      expand = expansion(add = 0.55)
    ) +
    row_scale(rows$genus_lab) +
    panel_theme

  heat_plot <- ggplot(heat, aes(x_plot, y_plot, fill = .data[[fill_col]])) +
    geom_tile(width = 0.9, height = 1) +
    scale_fill_viridis_c(
      option = "mako", direction = -1,
      limits = fill_limits, na.value = "grey85",
      name = fill_title, guide = colorbar(3)
    ) +
    scale_x_continuous(
      breaks = seq_along(group_order), labels = group_labels,
      limits = c(0.55, length(group_order) + 0.45), expand = c(0, 0)
    ) +
    row_scale() +
    panel_theme +
    theme(axis.text.y = element_blank())

  prev_plot <- ggplot(prev, aes(prevalence, y_jit)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.25) +
    geom_point(aes(fill = group_lab), shape = 21, size = 1.05, stroke = 0.1, alpha = 0.6) +
    scale_fill_manual(
      values = colors, labels = legend_labels, name = group_name,
      guide = guide_legend(
        order = 4, ncol = 1,
        override.aes = list(size = 1.7, alpha = 1)
      )
    ) +
    scale_x_continuous(limits = c(-8, 108), breaks = c(0, 50, 100), expand = c(0, 0)) +
    row_scale() +
    labs(x = "Prevalence (%)") +
    panel_theme +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = font_axis, angle = 0, hjust = 0.5, vjust = 1),
      axis.title.x = element_text(size = font_axis, hjust = 0, margin = margin(t = 2))
    )
  prev_plot <- patchwork::free(prev_plot, type = "label", side = "b")

  detail_plot <- (bubble_plot | heat_plot | prev_plot) +
    plot_layout(widths = c(6, 6, 4), guides = "collect")
  # Use the same font-aware device for layout measurements and PDF export.
  grDevices::cairo_pdf(
    filename = "/dev/null", width = width_mm / 25.4,
    height = detail_height_mm / 25.4
  )
  detail_grob <- patchworkGrob(detail_plot)

  invisible(grDevices::dev.off())

  detail_grob
}
