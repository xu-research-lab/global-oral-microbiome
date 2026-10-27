###############################################################################
# 05_pcoa_matrix_EC.R — EC PCoA 两两比较矩阵图（仅 PC1 vs PC2）
# 上三角=Bray-Curtis, 下三角=Aitchison, 对角线=大洲名+n
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrastr); library(patchwork)
})

proj   <- "/home/hexincheng/projects/hxc/oral_project/picrust2"
outdir <- file.path(proj, "downstream_EC_nsti2/02.global_structure")

continent_colors <- c(
  "North America" = "#E41A1C", "Asia" = "#FF7F00", "Europe" = "#377EB8",
  "Africa" = "#4DAF4A", "Oceania" = "#984EA3", "South America" = "#A65628"
)
continents <- names(continent_colors)
n <- length(continents)

ve_bray <- c(41.43, 13.95)
ve_ait  <- c(15.17, 14.20)

diag_plot <- function(g, n_samples, color) {
  ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.58, label = g, size = 6, fontface = "bold", color = color) +
    annotate("text", x = 0.5, y = 0.38, label = paste0("n = ", n_samples), size = 4, color = "grey40") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)
}

pair_plot <- function(df, g1, g2, xlabel, ylabel) {
  sub <- df[continent_group %in% c(g1, g2)]
  grp_order <- names(sort(table(sub$continent_group), decreasing = TRUE))
  sub$continent_group <- factor(sub$continent_group, levels = grp_order)
  ggplot(sub, aes(PC1, PC2, color = continent_group)) +
    rasterise(geom_point(size = 0.6, alpha = 0.3, shape = 16), dpi = 300) +
    stat_ellipse(type = "t", level = 0.95, linewidth = 0.7, show.legend = FALSE) +
    scale_color_manual(values = continent_colors) +
    labs(x = xlabel, y = ylabel) +
    theme_bw(base_size = 9) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
          axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          plot.margin = margin(2, 2, 2, 2))
}

df_bray <- fread(file.path(outdir, "pcoa_braycurtis_coords.tsv"))
df_ait  <- fread(file.path(outdir, "pcoa_aitchison_coords.tsv"))
n_samples <- table(df_bray$continent_group)[continents]

xl_bray <- sprintf("PC1 (%.1f%%)", ve_bray[1]); yl_bray <- sprintf("PC2 (%.1f%%)", ve_bray[2])
xl_ait  <- sprintf("PC1 (%.1f%%)", ve_ait[1]);  yl_ait  <- sprintf("PC2 (%.1f%%)", ve_ait[2])

plots <- list()
for (i in 1:n) {
  for (j in 1:n) {
    gi <- continents[i]; gj <- continents[j]
    if (i == j) {
      plots[[length(plots) + 1]] <- diag_plot(gi, n_samples[gi], continent_colors[gi])
    } else if (i < j) {
      plots[[length(plots) + 1]] <- pair_plot(df_bray, gi, gj, xl_bray, yl_bray)
    } else {
      plots[[length(plots) + 1]] <- pair_plot(df_ait, gi, gj, xl_ait, yl_ait)
    }
  }
}

matrix_plot <- wrap_plots(plots, ncol = n, nrow = n) +
  plot_annotation(
    title = "EC PCoA Matrix — Pairwise by Continent",
    subtitle = "Upper triangle: Bray-Curtis (TSS) | Lower triangle: Aitchison (CLR) | Ellipses: 95% CI",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 14, hjust = 0.5)))

plot_size <- n * 3.2
ggsave(file.path(outdir, "pcoa_matrix.pdf"), matrix_plot,
       width = plot_size, height = plot_size, limitsize = FALSE, device = cairo_pdf)
ggsave(file.path(outdir, "pcoa_matrix.png"), matrix_plot,
       width = plot_size, height = plot_size, dpi = 250, type = "cairo", bg = "white", limitsize = FALSE)
cat("Saved: EC pcoa_matrix.png/pdf\n")
cat("=== 05_pcoa_matrix_EC done ===\n")
