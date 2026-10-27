###############################################################################
# 05_pcoa_matrix_density.R — PCoA 两两比较矩阵图（密度等高线版）
# EC / KO / pathway 三套数据统一处理
# 上三角=Bray-Curtis, 下三角=Aitchison, 对角线=大洲名+n（加边框）
#
# 说明：
#   - 参考自 beta_diversity 的 04-plot-density-contour-pcoa.R（密度等高线画法）
#   - 数据不同：本脚本直接读各 outdir 下的 pcoa_*_coords.tsv，其中已含
#     continent_group 列，无需再合并元数据（参考脚本需按 SampleID 合并
#     metadata_oral_qc_dedup.tsv，此处不适用）
#   - EC/KO/pathway 距离矩阵没有成对按洲的 PERMANOVA 表，故不加 R²/P 值注释
# 输出：各 outdir 下 pcoa_matrix_density.pdf / .png（不覆盖原图）
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

proj <- "/home/hexincheng/projects/hxc/oral_project/picrust2"

continent_colors <- c(
  "North America" = "#E41A1C", "Asia" = "#FF7F00", "Europe" = "#377EB8",
  "Africa" = "#4DAF4A", "Oceania" = "#984EA3", "South America" = "#A65628"
)
continents <- names(continent_colors)
n <- length(continents)

# 每个panel只有两种颜色：行分组=红，列分组=蓝
color_row <- "#CC3311"
color_col <- "#0072B2"

configs <- list(
  list(type = "EC",       outdir = file.path(proj, "downstream_EC_nsti2/02.global_structure"),
       ve_bray = c(41.43, 13.95), ve_ait = c(15.17, 14.20)),
  list(type = "KO",       outdir = file.path(proj, "downstream_KO_nsti2/02.global_structure"),
       ve_bray = c(41.47, 15.82), ve_ait = c(14.42, 13.94)),
  list(type = "pathway",  outdir = file.path(proj, "downstream_pathway_nsti2/02.global_structure"),
       ve_bray = c(47.47, 10.09), ve_ait = c(18.05, 16.34))
)

# 对角线：大洲名 + 样本数，加黑色边框，文字统一黑色（与参考脚本一致）
# 坐标设置与密度panel保持一致（scale expansion + clip="off"），
# 否则 expand=FALSE 会让边框描边被裁一半，看起来比其他panel细
diag_plot <- function(g, n_samples) {
  ggplot() + theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          plot.margin = margin(2, 2, 2, 2)) +
    annotate("text", x = 0.5, y = 0.58, label = g, size = 6, fontface = "bold", color = "black") +
    annotate("text", x = 0.5, y = 0.38, label = paste0("n = ", n_samples), size = 4, color = "grey40") +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off")
}

# 密度等高线两两panel：每组各画一层半透明填充等高线（bins=6, alpha=0.15），
# 内层等高线通过重叠叠加不透明度，显示各组样本的实际聚集区域
# 不加轴标题（与参考脚本一致），保证所有panel方框大小一致；
# 方差解释比例改放在总subtitle
pair_plot_density <- function(df, g_row, g_col) {
  data_row <- df[continent_group == g_row, .(PC1, PC2)]
  data_col <- df[continent_group == g_col, .(PC1, PC2)]
  if (nrow(data_row) < 5 || nrow(data_col) < 5) return(plot_spacer())

  p <- ggplot() +
    theme_bw(base_size = 9) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          plot.margin = margin(2, 2, 2, 2))

  # 行分组（红色）——先画，背景层
  p <- p + stat_density_2d(
    data = data_row, aes(x = PC1, y = PC2),
    geom = "polygon", fill = color_row, color = color_row,
    alpha = 0.15, linewidth = 0.3, bins = 6, contour = TRUE
  )

  # 列分组（蓝色）——后画，前景层
  p <- p + stat_density_2d(
    data = data_col, aes(x = PC1, y = PC2),
    geom = "polygon", fill = color_col, color = color_col,
    alpha = 0.15, linewidth = 0.3, bins = 6, contour = TRUE
  )

  p <- p +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    coord_cartesian(clip = "off")
  p
}

for (cfg in configs) {
  cat("=== Processing", cfg$type, "===\n")
  df_bray <- fread(file.path(cfg$outdir, "pcoa_braycurtis_coords.tsv"))
  df_ait  <- fread(file.path(cfg$outdir, "pcoa_aitchison_coords.tsv"))
  n_samples <- table(df_bray$continent_group)[continents]

  xl_bray <- sprintf("PC1 (%.1f%%)", cfg$ve_bray[1]); yl_bray <- sprintf("PC2 (%.1f%%)", cfg$ve_bray[2])
  xl_ait  <- sprintf("PC1 (%.1f%%)", cfg$ve_ait[1]);  yl_ait  <- sprintf("PC2 (%.1f%%)", cfg$ve_ait[2])

  plots <- list()
  for (i in 1:n) {
    for (j in 1:n) {
      gi <- continents[i]; gj <- continents[j]
      if (i == j) {
        plots[[length(plots) + 1]] <- diag_plot(gi, n_samples[gi])
      } else if (i < j) {
        plots[[length(plots) + 1]] <- pair_plot_density(df_bray, gi, gj)
      } else {
        plots[[length(plots) + 1]] <- pair_plot_density(df_ait, gi, gj)
      }
    }
  }

  matrix_plot <- wrap_plots(plots, ncol = n, nrow = n) +
    plot_annotation(
      title = sprintf("%s PCoA Matrix — Pairwise by Continent", cfg$type),
      subtitle = paste0(
        "Upper triangle: Bray-Curtis (TSS) | Lower triangle: Aitchison (CLR)",
        "\nBray-Curtis: ", xl_bray, ", ", yl_bray,
        " | Aitchison: ", xl_ait, ", ", yl_ait,
        "\nDensity contours: two semi-transparent layers per group pair (red = row, blue = column)"
      ),
      theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 14, hjust = 0.5)))

  plot_size <- n * 3.2
  ggsave(file.path(cfg$outdir, sprintf("pcoa_matrix_density_%s.pdf", cfg$type)), matrix_plot,
         width = plot_size, height = plot_size, limitsize = FALSE, device = cairo_pdf)
  ggsave(file.path(cfg$outdir, sprintf("pcoa_matrix_density_%s.png", cfg$type)), matrix_plot,
         width = plot_size, height = plot_size, dpi = 250, bg = "white", limitsize = FALSE)
  cat("Saved:", cfg$type, sprintf("pcoa_matrix_density_%s.png/pdf", cfg$type), "\n")
}

cat("=== 05_pcoa_matrix_density done (EC/KO/pathway) ===\n")