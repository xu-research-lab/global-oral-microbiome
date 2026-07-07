#!/usr/bin/env Rscript

# 仅重新绘制 dbRDA 双序图（优化版），不重新运行分析
# 散点：用形状区分组别（不同测序区间），颜色统一为灰色
# 箭头：用颜色区分不同因子，标签仅显示组别名
library(vegan)
library(ggplot2)
library(data.table)

setwd("/beegfs/huangxiaochang/global_oral/genus_correction/dbRDA")

# 加载已有结果
cat("Loading saved results...\n")
load("dbRDA_results.RData")

rda_eig <- eigenvals(rda_model)
rda_total_var <- sum(rda_eig[rda_eig > 0])
rda1_pct <- round(rda_eig[1] / rda_total_var * 100, 2)
rda2_pct <- round(rda_eig[2] / rda_total_var * 100, 2)

sample_scores <- scores(rda_model, display = "sites", scaling = 3)
sample_scores_df <- as.data.frame(sample_scores)
sample_scores_df$SampleID <- rownames(sample_scores_df)

species_scores <- scores(rda_model, display = "bp", scaling = 3)
if (!is.null(species_scores) && nrow(species_scores) > 0) {
  species_scores_df <- as.data.frame(species_scores)
  species_scores_df$Variable <- rownames(species_scores_df)
} else {
  species_scores_df <- data.frame(RDA1 = numeric(0), RDA2 = numeric(0), Variable = character(0))
}

# 重新计算 anova
anova_df <- as.data.frame(anova_result)
anova_for_plot <- anova_df[!rownames(anova_df) %in% c("Model", "Residual"), , drop = FALSE]
anova_for_plot$F_val <- as.numeric(anova_df$F[!rownames(anova_df) %in% c("Model", "Residual")])

most_significant_factor <- rownames(anova_for_plot)[which.max(anova_for_plot$F_val)]
cat("Most significant factor:", most_significant_factor, "\n")

# 提取最显著因子的分组信息
point_group_df <- data.frame(
  SampleID = metadata$`#SampleID`,
  PointGroup = metadata[[most_significant_factor]]
)

# 过滤 Unknown 和低样本组别（含空字符串）
if (most_significant_factor == "Sex") {
  point_group_df <- point_group_df[point_group_df$PointGroup != "Unknown", ]
} else if (most_significant_factor == "Region_custom") {
  point_group_df <- point_group_df[!point_group_df$PointGroup %in% c("V2", "V3", "V5-V7", ""), ]
}

sample_color_df <- merge(sample_scores_df, point_group_df, by = "SampleID", all.x = TRUE)
sample_color_df <- sample_color_df[!is.na(sample_color_df$PointGroup), ]
point_groups <- sort(unique(sample_color_df$PointGroup))
n_point_groups <- length(point_groups)
cat("Number of point groups:", n_point_groups, "\n")
cat("Point groups:", paste(point_groups, collapse = ", "), "\n")

# 箭头：匹配因子
rda_terms <- attr(terms(rda_model), "term.labels")
rda_terms <- rda_terms[!grepl("^Condition\\(", rda_terms)]

match_factor <- function(var_name, factors) {
  for (f in factors) {
    if (startsWith(var_name, f)) return(f)
  }
  return(NA_character_)
}
species_scores_df$Factor <- sapply(species_scores_df$Variable, match_factor, rda_terms)

# 从变量名中提取纯组别名
strip_factor_prefix <- function(var_name, factor) {
  sub(paste0("^", gsub(".", "\\\\.", factor, fixed = TRUE)), "", var_name)
}
species_scores_df$Label <- mapply(strip_factor_prefix, species_scores_df$Variable, species_scores_df$Factor)

# 过滤 Health_Status 和 Unknown 组别
species_scores_clean <- species_scores_df[
  species_scores_df$Factor != "Health_Status" &
    !is.na(species_scores_df$Factor) &
    !species_scores_df$Variable %in% c("SexUnknown", "GeographyUnknown"),
]
# 过滤低样本测序区间
species_scores_clean <- species_scores_clean[!species_scores_clean$Label %in% c("V2", "V3", "V5-V7"), ]

# 重叠标签处理
species_scores_clean$LabelMult <- ifelse(
  species_scores_clean$Factor %in% c("Age_Group", "Geography"),
  1.35,
  1.15
)

# 箭头因子颜色
factor_colors <- c(
  "Sex" = "#1B9E77",
  "Age_Group" = "#D95F02",
  "Region_custom" = "#7570B3",
  "Geography" = "#E7298A"
)
available_factors <- unique(species_scores_clean$Factor)
factor_colors_used <- factor_colors[available_factors]

# 散点形状分配（最多 8 个不同形状）
all_shapes <- c(16, 17, 15, 3, 4, 8, 9, 10)
shape_vals <- all_shapes[seq_len(min(n_point_groups, length(all_shapes)))]
if (n_point_groups > length(all_shapes)) {
  shape_vals <- c(shape_vals, rep(0, n_point_groups - length(all_shapes)))
}
shape_palette <- setNames(shape_vals, point_groups)

# 绘图
if (nrow(species_scores_clean) > 0) {
  rda_plot <- ggplot() +
    # 散点：统一灰色，用形状区分组别
    geom_point(data = sample_color_df,
               aes(x = RDA1, y = RDA2, shape = PointGroup),
               color = "gray50", alpha = 0.7, size = 2) +
    # 箭头：按因子着色
    geom_segment(data = species_scores_clean,
                 aes(x = 0, y = 0, xend = RDA1, yend = RDA2, color = Factor),
                 arrow = arrow(length = unit(0.2, "cm")),
                 linewidth = 1.2) +
    geom_text(data = species_scores_clean,
              aes(x = RDA1 * LabelMult,
                  y = RDA2 * LabelMult,
                  label = Label, color = Factor),
              size = 4, fontface = "bold") +
    # 颜色只给箭头
    scale_color_manual(
      name = "Arrows (Factor)",
      values = factor_colors_used
    ) +
    # 形状给组别
    scale_shape_manual(
      name = most_significant_factor,
      values = shape_palette
    ) +
    labs(
      x = paste0("db-RDA1 (", rda1_pct, "%)"),
      y = paste0("db-RDA2 (", rda2_pct, "%)"),
      title = "db-RDA Biplot: Global Oral Microbiome",
      subtitle = paste0("Points (", most_significant_factor, ") | Arrows: ",
                        paste(available_factors, collapse = ", "))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
      panel.grid = element_line(color = "gray85"),
      aspect.ratio = 1,
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing = unit(0.8, "cm")
    )
} else {
  rda_plot <- ggplot(sample_color_df, aes(x = RDA1, y = RDA2)) +
    geom_point(color = "gray50", alpha = 0.7, size = 2) +
    labs(
      x = paste0("db-RDA1 (", rda1_pct, "%)"),
      y = paste0("db-RDA2 (", rda2_pct, "%)"),
      title = "db-RDA: Global Oral Microbiome"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_line(color = "gray85"),
      aspect.ratio = 1
    )
}

ggsave("dbRDA_biplot.png", rda_plot, width = 12, height = 9, dpi = 300)
cat("Saved: dbRDA_biplot.png\n")
