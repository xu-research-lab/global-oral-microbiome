#!/usr/bin/env Rscript

# =============================================================================
# 全球口腔微生物组 db-RDA 分析（使用预计算的 PCoA 坐标）
# 分析 Sex, Age_Group, Region_custom, Health_Status, Geography 对微生物组成的影响
# 同时排除 Bioproject 因素的影响
#
# 修改说明：
# 1. 使用 qiime2_metadata_with_geography.tsv（已有 Geography 列，无需重新生成）
# 2. 使用 QIIME2 预计算的 Bray-Curtis PCoA 坐标（从 bray_curtis_pcoa_results.qza 提取）
# 3. 跳过距离矩阵计算和 PCoA 步骤，直接对 PCoA 坐标做 RDA
# 4. 显著节省内存和时间（避免 13GB+ 距离矩阵计算）
# =============================================================================

library(vegan)
library(ggplot2)
library(data.table)

# 输出到 log 文件
log_con <- file("dbRDA_analysis.log", open = "wt")
sink(log_con, append = TRUE, type = "output")
sink(log_con, append = TRUE, type = "message")

cat("=== db-RDA Analysis Log ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

mem_report <- function(label = "") {
  gc_result <- gc(verbose = FALSE)
  cat(sprintf("[%s] Memory: Used = %.1f MB, Max = %.1f MB\n",
              label, gc_result[2, 2], gc_result[2, 3]))
  return(invisible(NULL))
}

theme_set(theme_minimal(base_size = 14))

setwd("/beegfs/huangxiaochang/global_oral/genus_correction/dbRDA")

# =============================================================================
# 1. 加载 metadata
# =============================================================================
cat("\nLoading metadata...\n")

metadata <- fread("../../qiime2_metadata_with_geography.tsv", header = TRUE,
                  sep = "\t", check.names = FALSE, data.table = FALSE)
cat("Metadata loaded:", nrow(metadata), "samples\n")

if (!"Geography" %in% colnames(metadata)) {
  stop("Geography column not found in metadata file!")
}
cat("Geography column found in metadata.\n")
cat("Geography distribution:\n")
print(table(metadata$Geography))

# =============================================================================
# 2. 数据预处理
# =============================================================================
cat("\nPreprocessing data...\n")

metadata$Sex <- ifelse(metadata$Sex %in% c("M", "Male", "male"), "Male",
                       ifelse(metadata$Sex %in% c("F", "Female", "female"), "Female", "Unknown"))
metadata$Age_Group <- as.factor(metadata$Age_Group)
metadata$Region_custom <- as.factor(metadata$Region_custom)
metadata$Health_Status <- as.factor(metadata$Healthy_Status)
metadata$Geography <- as.factor(metadata$Geography)
metadata$Bioproject <- as.factor(metadata$BioProject)

# 过滤缺失值
metadata <- metadata[!is.na(metadata$Sex) & !is.na(metadata$Age_Group) &
                       !is.na(metadata$Region_custom) & !is.na(metadata$Health_Status) &
                       !is.na(metadata$Geography) & !is.na(metadata$Bioproject), ]
cat("Metadata after filtering:", nrow(metadata), "samples\n")

sample_ids <- metadata$`#SampleID`

# =============================================================================
# 3. 加载预提取的 PCoA 坐标和特征值
#    从 QIIME2 bray_curtis_pcoa_results.qza 提取
#    运行: python extract_pcoa_from_qza.py <qza_file> <output_dir> [max_axes]
# =============================================================================
cat("\nLoading pre-computed PCoA coordinates...\n")
mem_report("Before loading PCoA")

pcoa_coords <- fread("pcoa_coordinates.tsv", header = TRUE,
                     sep = "\t", check.names = FALSE, data.table = FALSE)
cat("PCoA coordinates loaded:", nrow(pcoa_coords), "samples x", ncol(pcoa_coords) - 1, "axes\n")

eig_df <- fread("pcoa_eigenvalues.tsv", header = TRUE,
                sep = "\t", check.names = FALSE, data.table = FALSE)
cat("Eigenvalues loaded:", nrow(eig_df), "axes\n")

# 确保样本顺序与 metadata 一致
common_samples <- intersect(sample_ids, pcoa_coords$SampleID)
cat("Common samples:", length(common_samples), "\n")

if (length(common_samples) < 100) {
  stop("ERROR: Too few common samples (<100)!")
}

# 子集化 PCoA 坐标
pcoa_coords <- pcoa_coords[pcoa_coords$SampleID %in% common_samples, ]
metadata <- metadata[metadata$`#SampleID` %in% common_samples, ]

# 按样本 ID 排序以确保一致
pcoa_coords <- pcoa_coords[order(pcoa_coords$SampleID), ]
metadata <- metadata[order(metadata$`#SampleID`), ]

# 验证顺序
if (!all(pcoa_coords$SampleID == metadata$`#SampleID`)) {
  stop("ERROR: Sample order mismatch between PCoA coords and metadata!")
}

cat("Samples aligned:", length(common_samples), "\n")

# 提取 PCoA 坐标矩阵（去除 SampleID 列）
pcoa_axes <- setdiff(colnames(pcoa_coords), "SampleID")
n_axes_total <- length(pcoa_axes)

# 特征值统计
positive_eig <- eig_df$Eigenvalue[eig_df$Eigenvalue > 0]
total_var <- sum(positive_eig)

cat("\n=== PCoA Summary ===\n")
cat("Total PCoA axes available:", n_axes_total, "\n")
cat("Positive eigenvalues:", length(positive_eig), "\n")
cat("Top 3 eigenvalues:", head(eig_df$Eigenvalue, 3), "\n")

n_axes_10 <- min(10, nrow(eig_df))
n_axes_50 <- min(50, nrow(eig_df))
cat("Cumulative variance - top 10 axes:", round(eig_df$CumulativeVariance[n_axes_10], 2), "%\n")
cat("Cumulative variance - top 50 axes:", round(eig_df$CumulativeVariance[n_axes_50], 2), "%\n")
cat("Cumulative variance - all", n_axes_total, "axes:", round(eig_df$CumulativeVariance[n_axes_total], 2), "%\n")

mem_report("After loading PCoA")

# =============================================================================
# 4. RDA 分析
#    对 PCoA 坐标做 RDA，等价于 db-RDA
#    PCoA 坐标已经由特征值加权，因此直接做 RDA 即可
# =============================================================================
cat("\nRunning RDA on PCoA scores...\n")
start_time <- Sys.time()

Y_pcoa <- as.matrix(pcoa_coords[, pcoa_axes, drop = FALSE])
rownames(Y_pcoa) <- pcoa_coords$SampleID

rda_model <- rda(
  Y_pcoa ~ Sex + Age_Group + Region_custom + Health_Status + Geography +
    Condition(Bioproject),
  data = metadata
)

end_time <- Sys.time()
cat("RDA completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

# =============================================================================
# 5. 统计检验
# =============================================================================
cat("\nCalculating statistics...\n")

# 方差分解（按项检验）
anova_result <- anova(rda_model, by = "terms", permutations = 999)
cat("\n=== ANOVA-like Permutation Test (by terms) ===\n")
print(anova_result)

# 整体模型显著性
overall_anova <- anova(rda_model, permutations = 999)
cat("\n=== Overall Model Significance ===\n")
print(overall_anova)

# 调整后 R²
r2_result <- RsquareAdj(rda_model)$adj.r.squared
cat("\n=== Adjusted R² ===\n")
cat("Adjusted R²:", round(r2_result, 4), "\n")

# 提取得分
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

# =============================================================================
# 6. 保存结果
# =============================================================================
cat("\nSaving results...\n")

save(rda_model, sample_scores_df, species_scores_df, anova_result,
     overall_anova, r2_result, metadata, pcoa_coords, eig_df,
     file = "dbRDA_results.RData")

anova_df <- as.data.frame(anova_result)
write.table(anova_df, "dbRDA_anova_results.tsv", sep = "\t",
            quote = FALSE, row.names = TRUE)

write.table(sample_scores_df, "dbRDA_sample_scores.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

if (nrow(species_scores_df) > 0) {
  write.table(species_scores_df, "dbRDA_variable_scores.tsv", sep = "\t",
              quote = FALSE, row.names = FALSE)
}

write.table(eig_df, "dbRDA_pcoa_eigenvalues.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("Results saved.\n")

# =============================================================================
# 7. 可视化
# =============================================================================
cat("\nCreating visualizations...\n")

rda_eig <- eigenvals(rda_model)
rda_total_var <- sum(rda_eig[rda_eig > 0])
rda1_pct <- round(rda_eig[1] / rda_total_var * 100, 2)
rda2_pct <- round(rda_eig[2] / rda_total_var * 100, 2)

# --- 图 1: db-RDA 双序图（优化版）---
# 策略：
# 1. 箭头按因子着色，同一因子内所有组别共用一色
# 2. Health_Status 组别过多，不显示其箭头
# 3. 散点按最显著的因子着色

if (nrow(species_scores_df) > 0) {
  # 从变量名中提取因子名（去掉水平名后缀）
  # RDA 的 bp 得分命名格式为 "因子名水平名"，如 "SexMale", "GeographyAsia"
  # 需要根据 RDA 模型的项名来匹配
  rda_terms <- attr(terms(rda_model), "term.labels")
  # 去掉 "Condition(...)" 项
  rda_terms <- rda_terms[!grepl("^Condition\\(", rda_terms)]
  cat("RDA constraint factors:", paste(rda_terms, collapse = ", "), "\n")

  # 为每个箭头变量匹配其所属因子
  match_factor <- function(var_name, factors) {
    for (f in factors) {
      if (startsWith(var_name, f)) return(f)
    }
    return(NA_character_)
  }
  species_scores_df$Factor <- sapply(species_scores_df$Variable, match_factor, rda_terms)

  # 过滤掉 Health_Status（组别过多，箭头太密集）
  species_scores_clean <- species_scores_df[species_scores_df$Factor != "Health_Status" &
                                              !is.na(species_scores_df$Factor), ]

  # 为每个因子分配颜色
  factor_colors <- c(
    "Sex" = "#1B9E77",
    "Age_Group" = "#D95F02",
    "Region_custom" = "#7570B3",
    "Geography" = "#E7298A",
    "Health_Status" = "#66A61E"  # 保留定义但不在图中显示
  )
  available_factors <- unique(species_scores_clean$Factor)
  factor_colors_used <- factor_colors[available_factors]

  # 确定最显著的因子（用于散点着色）
  anova_sorted <- anova_for_plot[order(anova_for_plot$F.Model, decreasing = TRUE), ]
  most_significant_factor <- anova_sorted$Variable[1]
  cat("Most significant factor:", most_significant_factor, "\n")

  # 提取最显著因子的分组信息
  if (most_significant_factor %in% colnames(metadata)) {
    point_color_df <- data.frame(
      SampleID = metadata$`#SampleID`,
      ColorGroup = metadata[[most_significant_factor]]
    )
    sample_color_df <- merge(sample_scores_df, point_color_df, by = "SampleID", all.x = TRUE)
  } else {
    sample_color_df <- sample_scores_df
    sample_color_df$ColorGroup <- "Unknown"
  }

  n_color_groups <- length(unique(sample_color_df$ColorGroup))
  cat("Number of color groups:", n_color_groups, "\n")

  # 根据组数选择合适的散点调色板
  if (n_color_groups <= 8) {
    point_palette <- RColorBrewer::brewer.pal(n_color_groups, "Set1")
  } else if (n_color_groups <= 12) {
    point_palette <- RColorBrewer::brewer.pal(n_color_groups, "Set3")
  } else {
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      point_palette <- viridisLite::viridis(n_color_groups, option = "turbo")
    } else {
      point_palette <- scales::hue_pal()(n_color_groups)
    }
  }
  names(point_palette) <- unique(sample_color_df$ColorGroup)

  # 合并所有颜色（散点组别 + 箭头因子），避免重复键
  all_colors <- c(factor_colors_used, point_palette)

  # 绘制双序图
  if (nrow(species_scores_clean) > 0) {
    rda_plot <- ggplot() +
      # 散点：按最显著因子着色
      geom_point(data = sample_color_df,
                 aes(x = RDA1, y = RDA2, color = ColorGroup),
                 alpha = 0.6, size = 1.5) +
      # 箭头：按因子着色
      geom_segment(data = species_scores_clean,
                   aes(x = 0, y = 0, xend = RDA1, yend = RDA2, color = Factor),
                   arrow = arrow(length = unit(0.2, "cm")),
                   linewidth = 1.2) +
      geom_text(data = species_scores_clean,
                aes(x = RDA1 * 1.15, y = RDA2 * 1.15, label = Variable, color = Factor),
                size = 4, fontface = "bold") +
      scale_color_manual(
        values = all_colors,
        guide = guide_legend(override.aes = list(
          shape = c(rep(NA, length(factor_colors_used)), rep(16, n_color_groups)),
          linetype = c(rep(1, length(factor_colors_used)), rep(0, n_color_groups)),
          size = c(rep(1.2, length(factor_colors_used)), rep(3, n_color_groups))
        ))
      ) +
      labs(
        x = paste0("db-RDA1 (", rda1_pct, "%)"),
        y = paste0("db-RDA2 (", rda2_pct, "%)"),
        title = "db-RDA Biplot: Global Oral Microbiome",
        color = "Factor / Group",
        subtitle = paste0("Points colored by: ", most_significant_factor,
                          " | Arrows: ", paste(available_factors, collapse = ", "))
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
        panel.grid = element_line(color = "gray85"),
        aspect.ratio = 1,
        legend.position = "right",
        legend.box = "vertical"
      )
  } else {
    # 如果没有可显示的箭头，只画散点
    rda_plot <- ggplot(sample_color_df,
                       aes(x = RDA1, y = RDA2, color = ColorGroup)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = all_colors) +
      labs(
        x = paste0("db-RDA1 (", rda1_pct, "%)"),
        y = paste0("db-RDA2 (", rda2_pct, "%)"),
        title = "db-RDA: Global Oral Microbiome",
        subtitle = paste0("Points colored by: ", most_significant_factor),
        color = most_significant_factor
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
        panel.grid = element_line(color = "gray85"),
        aspect.ratio = 1
      )
  }

  ggsave("dbRDA_biplot.png", rda_plot, width = 12, height = 9, dpi = 300)
  cat("  - Saved: dbRDA_biplot.png\n")
}

# --- 图 2: 按 Geography 着色 ---
geo_scores <- merge(sample_scores_df,
                    metadata[, c("#SampleID", "Geography")],
                    by.x = "SampleID", by.y = "#SampleID")

geo_plot <- ggplot(geo_scores, aes(x = RDA1, y = RDA2, color = Geography)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = paste0("db-RDA1 (", rda1_pct, "%)"),
    y = paste0("db-RDA2 (", rda2_pct, "%)"),
    title = "db-RDA by Geography",
    color = "Geography"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1
  )

ggsave("dbRDA_by_geography.png", geo_plot, width = 10, height = 8, dpi = 300)
cat("  - Saved: dbRDA_by_geography.png\n")

# --- 图 3: 按 Health_Status 着色 ---
health_scores <- merge(sample_scores_df,
                       metadata[, c("#SampleID", "Health_Status")],
                       by.x = "SampleID", by.y = "#SampleID")

health_plot <- ggplot(health_scores, aes(x = RDA1, y = RDA2, color = Health_Status)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = paste0("db-RDA1 (", rda1_pct, "%)"),
    y = paste0("db-RDA2 (", rda2_pct, "%)"),
    title = "db-RDA by Health Status",
    color = "Health Status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1
  )

ggsave("dbRDA_by_health.png", health_plot, width = 10, height = 8, dpi = 300)
cat("  - Saved: dbRDA_by_health.png\n")

# --- 图 4: 方差解释条形图 ---
anova_for_plot <- anova_df
anova_for_plot$Variable <- rownames(anova_for_plot)
anova_for_plot <- anova_for_plot[!anova_for_plot$Variable %in% c("Model", "Residual"), ]
anova_for_plot$F.Model <- as.numeric(anova_for_plot$F.Model)
anova_for_plot$Pr <- as.numeric(anova_for_plot$`Pr(>F)`)
anova_for_plot$Variable <- factor(anova_for_plot$Variable,
                                  levels = anova_for_plot$Variable[order(anova_for_plot$F.Model, decreasing = TRUE)])

variance_bar_plot <- ggplot(anova_for_plot, aes(x = Variable, y = F.Model, fill = Variable)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("p = %.3f", Pr)), vjust = -0.5, size = 4) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    x = "Variable",
    y = "F-statistic",
    title = "Variance Explained by Each Factor (db-RDA)",
    subtitle = paste("Adjusted R² =", round(r2_result, 4))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("dbRDA_variance_bar.png", variance_bar_plot, width = 10, height = 7, dpi = 300)
cat("  - Saved: dbRDA_variance_bar.png\n")

# --- 图 5: PCoA 特征值（前20轴）---
eig_plot <- ggplot(head(eig_df[eig_df$Eigenvalue > 0, ], 20),
                   aes(x = Axis, y = VarianceExplained)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "PCoA Axis",
    y = "Variance Explained (%)",
    title = "PCoA Eigenvalues (Top 20 Axes)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("dbRDA_pcoa_eigenvalues.png", eig_plot, width = 10, height = 6, dpi = 300)
cat("  - Saved: dbRDA_pcoa_eigenvalues.png\n")

# --- 图 6: 累积解释度曲线 ---
cumulative_plot <- ggplot(eig_df[eig_df$Eigenvalue > 0, ],
                          aes(x = Axis, y = CumulativeVariance)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 1) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red", alpha = 0.7) +
  labs(
    x = "PCoA Axis",
    y = "Cumulative Variance Explained (%)",
    title = "PCoA Cumulative Variance"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("dbRDA_cumulative_variance.png", cumulative_plot, width = 10, height = 6, dpi = 300)
cat("  - Saved: dbRDA_cumulative_variance.png\n")

# --- 组合图 ---
if (requireNamespace("gridExtra", quietly = TRUE) && nrow(species_scores_df) > 0) {
  library(gridExtra)
  combined_plot <- grid.arrange(
    rda_plot + theme(plot.title = element_text(size = 12)),
    variance_bar_plot + theme(plot.title = element_text(size = 12)),
    ncol = 2
  )
  ggsave("dbRDA_combined.png", combined_plot, width = 16, height = 8, dpi = 300)
  cat("  - Saved: dbRDA_combined.png\n")
}

cat("Visualizations completed.\n")

# =============================================================================
# 8. 结果摘要
# =============================================================================
cat("\n")
cat("========================================\n")
cat("     DB-RDA ANALYSIS SUMMARY\n")
cat("========================================\n\n")

cat("1. Sample size:", nrow(metadata), "samples\n")
cat("2. Number of Bioproject factors (conditioned):", length(unique(metadata$Bioproject)), "\n")
cat("3. PCoA axes available:", n_axes_total, "\n")
cat("4. PCoA axes used for RDA:", n_axes_total, "\n")
cat("5. Cumulative variance explained by PCoA axes:",
    round(eig_df$CumulativeVariance[n_axes_total], 2), "%\n")
cat("6. Adjusted R² (db-RDA):", round(r2_result, 4), "\n\n")

cat("7. Significance of each factor:\n")
print(anova_df)

cat("\n8. Output files:\n")
cat("   - dbRDA_results.RData (all results)\n")
cat("   - dbRDA_anova_results.tsv\n")
cat("   - dbRDA_sample_scores.tsv\n")
cat("   - dbRDA_variable_scores.tsv\n")
cat("   - dbRDA_pcoa_eigenvalues.tsv\n")
cat("   - dbRDA_biplot.png\n")
cat("   - dbRDA_by_geography.png\n")
cat("   - dbRDA_by_health.png\n")
cat("   - dbRDA_variance_bar.png\n")
cat("   - dbRDA_pcoa_eigenvalues.png\n")
cat("   - dbRDA_cumulative_variance.png\n")

cat("\n========================================\n")
cat("           ANALYSIS COMPLETE\n")
cat("========================================\n")

sink(type = "output")
sink(type = "message")
close(log_con)
