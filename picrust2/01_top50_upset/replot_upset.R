###############################################################################
# replot_upset.R — 三层级 Top50 UpSet 图重绘（ComplexUpset 引擎）
# 仅改两处：
#   (1) 配色:
#         交集矩阵点   = #d95f5f (珊瑚红实心点, 非成员 = grey85)
#         交集大小柱   = #7bb5de (浅蓝)
#         集合大小柱   = #f7cf6b (淡黄)
#   (2) 大洲缩写: SA/OC/NA/EU/AS/AF (自上而下, 字母序自下而上)
# 其余全部保持之前设定:
#   - 交集按大小从大到小排序 (ComplexUpset 默认)
#   - height_ratio = 0.3, width_ratio = 0.2
#   - 版面 10x6 @ 300 dpi, PNG + PDF
# 直接读已保存的 membership matrix，不重跑 Top50 流程
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ComplexUpset); library(ggplot2)
})

proj <- "/home/hexincheng/projects/hxc/oral_project/picrust2"
TOP_N <- 50
types <- c("EC", "KO", "pathway")   # 与 01.top50_upset.R 的 short 名一致
W <- 10; H <- 6; DPI <- 300
HR <- 0.3   # Intersection size 面板高度比（默认 0.5）

RED  <- "#d95f5f"   # 矩阵成员点
BLUE <- "#7bb5de"   # Intersection size 柱
YEL  <- "#f7cf6b"   # Set size 柱
GREY <- "grey85"    # 矩阵非成员点

abbr <- c("Africa" = "AF", "Asia" = "AS", "Europe" = "EU",
          "North America" = "NA", "Oceania" = "OC", "South America" = "SA")

for (type in types) {
  outdir <- file.path(proj, paste0("downstream_", type, "_nsti2/01.top50_upset"))

  m <- fread(file.path(outdir, paste0("top", TOP_N, "_", type, "_membership_matrix.tsv")))
  cont_cols <- setdiff(names(m), "feature_id")
  for (cc in cont_cols) m[[cc]] <- as.logical(m[[cc]])

  ia <- intersection_size(counts = TRUE, fill = BLUE,
         text_colors = c(on_background = "black", on_bar = "black"))
  s1 <- upset_set_size(mapping = aes(fill = YEL)) +
        scale_fill_identity(guide = "none")
  mx <- intersection_matrix(geom = ggplot2::geom_point(size = 3),
         outline_color = list(active = RED, inactive = GREY)) +
        scale_color_manual(values = c("TRUE" = RED, "FALSE" = GREY), guide = "none")

  p <- upset(m, intersect = cont_cols,
             name = paste0("Top", TOP_N, " ", type),
             width_ratio = 0.2,
             height_ratio = HR,
             labeller = function(x) unname(abbr[x]),
             base_annotations = list("Intersection size" = ia),
             set_sizes = s1,
             matrix = mx,
             themes = upset_modify_themes(list(
               "overall_sizes" = theme(legend.position = "none")))) +
    ggtitle(paste0("Top", TOP_N, " ", type, " overlap by continent"))

  ggsave(file.path(outdir, paste0("top", TOP_N, "_", type, "_UpSet.pdf")), p,
         width = W, height = H)
  ggsave(file.path(outdir, paste0("top", TOP_N, "_", type, "_UpSet.png")), p,
         width = W, height = H, dpi = DPI, bg = "white")
  cat("replotted:", type, "\n")
}

cat("=== replot_upset done ===\n")