# Save standalone OVA count and three-panel figures at final publication size.
# Requires the same plotting libraries and OVA_GENUS_FONT as the OVA drivers.
save_ova_publication_figures <- function(input_dir, group_order, group_name,
                                         versions = c("2method", "3method"),
                                         width_mm = 120, count_height_mm = 57,
                                         detail_height_mm = 205,
                                         count_width_mm = 70) {
  stopifnot(width_mm >= 85, count_width_mm >= 50,
            count_height_mm >= 30, detail_height_mm >= 170)
  old_theme <- theme_get()
  on.exit(theme_set(old_theme), add = TRUE)
  for (ver in versions) {
    if (!file.exists(file.path(input_dir, paste0("figureB_heatmap_data_", ver, ".tsv")))) next
    for (hm_mode in c("zscore", "raw")) {
      read_plot_table <- function(name) fread(file.path(input_dir, name), data.table = FALSE)
      counts <- read_plot_table(paste0("figureA_consensus_counts_", ver, ".tsv")) %>%
        arrange(desc(n_genera), match(target, group_order))
      heat <- read_plot_table(paste0("figureB_heatmap_data_", ver, ".tsv"))
      bubble <- read_plot_table(paste0("figureB_bubble_data_", ver, ".tsv"))
      prev <- read_plot_table(paste0("figureB_prevalence_data_", ver, ".tsv"))
      rows <- heat %>% distinct(feature, genus_lab, y_plot) %>% arrange(y_plot)
      stopifnot(!anyDuplicated(rows$feature), !anyDuplicated(rows$y_plot),
                all(diff(counts$n_genera) <= 0),
                all(bubble$y_plot == rows$y_plot[match(bubble$feature, rows$feature)]),
                all(prev$y_plot == rows$y_plot[match(prev$feature, rows$feature)]),
                isTRUE(all.equal(prev$prevalence, 100 * prev$n_present / prev$n_total)))
      
      # Match regular and italic text within the same installed sans-serif family.
      font_family <- OVA_GENUS_FONT
      font_body <- 7.5
      font_axis <- 7
      font_legend <- 7
      theme_set(theme_minimal(base_family = font_family, base_size = font_body) +
                  theme(legend.position = "right", legend.box = "vertical",
                        legend.direction = "vertical", legend.title.position = "top",
                        legend.box.just = "left", legend.spacing.y = unit(2, "mm")))
      colors <- setNames(viridis(length(group_order)), group_order)
      group_labels <- group_order
      legend_labels <- group_order
      count_order <- counts$target
      y_limits <- range(rows$y_plot) + c(-0.5, 0.5)
      bubble$target <- factor(bubble$target, levels = group_order)
      prev$group_lab <- factor(prev$group_lab, levels = group_order)
      fill_col <- if (hm_mode == "raw") "raw_abun" else "z_abun"
      fill_limits <- if (hm_mode == "raw") c(floor(min(heat$raw_abun, na.rm = TRUE)), ceiling(max(heat$raw_abun, na.rm = TRUE))) else {
        m <- max(abs(heat$z_abun), na.rm = TRUE); c(-m, m)
      }
      fill_title <- if (hm_mode == "raw") "Mean relative\nabundance\n(log10)" else "Mean relative\nabundance\n(z-score)"
      
      count_plot <- ova_count_plot(counts, group_order[group_order %in% counts$target], font_body)

      panel_theme <- theme_minimal(base_family = font_family, base_size = font_body) +
        theme(axis.title = element_blank(),
              axis.text = element_text(colour = "black"),
              axis.text.y = element_text(size = font_body, face = "italic", hjust = 1,
                                         margin = margin(r = 2)),
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
              plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm"))
      row_scale <- function(labels = NULL) {
        scale_y_continuous(limits = y_limits, breaks = rows$y_plot,
                           labels = labels, expand = c(0, 0))
      }
      colorbar <- function(order) guide_colorbar(
        order = order, direction = "vertical",
        theme = theme(legend.key.width = unit(2, "mm"),
                      legend.key.height = unit(12, "mm")))
      
      bubble_plot <- ggplot(bubble, aes(target, y_plot, size = neg_log10_maaq,
                                       fill = maaslin2_coef)) +
        geom_point(shape = 21, alpha = 0.9, colour = "black", stroke = 0.22) +
        scale_size_continuous(range = c(0.45, 1.8), name = "-log10(FDR)",
                              guide = guide_legend(order = 2, ncol = 1)) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                              midpoint = 0, name = "MaAsLin2 coefficient\n(Target vs. Rest)",
                              guide = colorbar(1)) +
        scale_x_discrete(labels = paste0(group_labels, " vs rest"), drop = FALSE,
                          expand = expansion(add = 0.55)) +
        row_scale(rows$genus_lab) + panel_theme
      
      heat_plot <- ggplot(heat, aes(x_plot, y_plot, fill = .data[[fill_col]])) +
        geom_tile(width = 0.9, height = 1) +
        scale_fill_viridis_c(option = "mako", direction = -1,
                              limits = fill_limits, na.value = "grey85",
                              name = fill_title, guide = colorbar(3)) +
        scale_x_continuous(breaks = seq_along(group_order), labels = group_labels,
                           limits = c(0.55, length(group_order) + 0.45), expand = c(0, 0)) +
        row_scale() + panel_theme + theme(axis.text.y = element_blank())
      
      prev_plot <- ggplot(prev, aes(prevalence, y_jit)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.25) +
        geom_point(aes(fill = group_lab), shape = 21, size = 1.05, stroke = 0.1, alpha = 0.6) +
        scale_fill_manual(values = colors, labels = legend_labels, name = group_name,
                           guide = guide_legend(order = 4, ncol = 1,
                                                override.aes = list(size = 1.7, alpha = 1))) +
        scale_x_continuous(limits = c(-8, 108), breaks = c(0, 50, 100), expand = c(0, 0)) +
        row_scale() + labs(x = "Prevalence (%)") + panel_theme +
        theme(axis.text.y = element_blank(),
              axis.text.x = element_text(size = font_axis, angle = 0, hjust = 0.5, vjust = 1),
              axis.title.x = element_text(size = font_axis, hjust = 0, margin = margin(t = 2)))
      prev_plot <- patchwork::free(prev_plot, type = "label", side = "b")
      
      detail_plot <- (bubble_plot | heat_plot | prev_plot) +
        plot_layout(widths = c(6, 6, 4), guides = "collect")
      # Use the same font-aware device for layout measurements and PDF export.
      grDevices::cairo_pdf(filename = "/dev/null", width = width_mm / 25.4,
                           height = detail_height_mm / 25.4)
      detail_grob <- patchworkGrob(detail_plot)
      count_grob <- ggplotGrob(count_plot)
      
      invisible(grDevices::dev.off())
      
      save_publication <- function(plot, stem, height, width = width_mm) {
        ggsave(file.path(input_dir, paste0(stem, ".pdf")), plot,
               device = grDevices::cairo_pdf, width = width, height = height,
               units = "mm", limitsize = FALSE, bg = "white")
        ggsave(file.path(input_dir, paste0(stem, ".png")), plot,
               device = ragg::agg_png, width = width, height = height,
               units = "mm", dpi = 300, limitsize = FALSE, bg = "white")
      }
      if (hm_mode == "zscore") {
        save_publication(count_grob, paste0("figureA_consensus_counts_by_target_", ver),
                         count_height_mm, width = count_width_mm)
      }
      if (ver == "2method" && hm_mode == "zscore") {
        save_publication(count_grob, "figureA_consensus_counts_by_target", count_height_mm,
                         width = count_width_mm)
      }
      save_publication(detail_grob, paste0("figureB_abundance_", ver, "_combined_", hm_mode), detail_height_mm)
      if (hm_mode == "zscore") {
        save_publication(detail_grob, paste0("figureB_abundance_", ver, "_combined"), detail_height_mm)
      }
    }
  }
  writeLines(c(
    paste0("Standalone Figure A: ", count_width_mm, " x ", count_height_mm, " mm."),
    paste0("Standalone Figure B: ", width_mm, " x ", detail_height_mm, " mm."),
    "Genus/count text: 7.5 pt; axes/legends: 7 pt. Figure B legends are on the right.",
    "PNG: 300 dpi. PDF: vector artwork with embedded fonts.",
    "Saved plotting values, genus order and group colors are preserved."
  ), file.path(input_dir, "figure_publication_layout.txt"))
  message("Saved standalone publication figures in ", input_dir)
}
