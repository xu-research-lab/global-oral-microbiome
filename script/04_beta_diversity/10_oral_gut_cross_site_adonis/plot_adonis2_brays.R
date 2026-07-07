library(ggplot2)
library(dplyr)
library(tidyr)

# 读取四个部位的 adonis 结果
saliva <- read.csv("adonis2_result_Saliva_brays.csv", row.names = 1)
mucosa <- read.csv("adonis2_result_Oral_Mucosa_brays.csv", row.names = 1)
tongue <- read.csv("adonis2_result_Tongue_brays.csv", row.names = 1)
tooth <- read.csv("adonis2_result_Tooth_brays.csv", row.names = 1)

# 将行名转换为 Factor 列
saliva$Factor <- rownames(saliva)
mucosa$Factor <- rownames(mucosa)
tongue$Factor <- rownames(tongue)
tooth$Factor <- rownames(tooth)

# 重命名 Region_custom 为 Sequencing_Region
saliva$Factor[saliva$Factor == "Region_custom"] <- "Sequencing_Region"
mucosa$Factor[mucosa$Factor == "Region_custom"] <- "Sequencing_Region"
tongue$Factor[tongue$Factor == "Region_custom"] <- "Sequencing_Region"
tooth$Factor[tooth$Factor == "Region_custom"] <- "Sequencing_Region"

# 添加部位信息
saliva$Site <- "Saliva"
mucosa$Site <- "Oral Mucosa"
tongue$Site <- "Tongue"
tooth$Site <- "Tooth"

# 合并数据
all_data <- bind_rows(saliva, mucosa, tongue, tooth)

# 移除 Residual 和 Total 行，只保留影响因素
all_data <- all_data[!all_data$Factor %in% c("Residual", "Total"), ]

# 将 Site 转换为因子，并设置顺序
all_data$Site <- factor(all_data$Site,
                        levels = c("Saliva", "Oral Mucosa", "Tongue", "Tooth"))

# 获取 Saliva 的 R2 值用于排序
saliva_order <- saliva[saliva$Factor %in% all_data$Factor, ]
saliva_order <- saliva_order[order(saliva_order$R2, decreasing = TRUE), ]
factor_order <- saliva_order$Factor

# 将因素转换为因子，按照 Saliva R2 从高到低排序
all_data$Factor <- factor(all_data$Factor, levels = rev(factor_order))

# 创建显著性标注函数
get_sig <- function(p) {
  ifelse(p <= 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

# 添加显著性标注（在 merge 之前）
all_data$Significance <- get_sig(all_data$Pr..F.)

# 设置颜色
colors <- c("Saliva" = "#E41A1C",
            "Oral Mucosa" = "#377EB8",
            "Tongue" = "#4DAF4A",
            "Tooth" = "#984EA3")

# 计算每个分面的最大 R2 用于设置 x 轴范围
max_r2_per_site <- all_data %>%
  group_by(Site) %>%
  summarise(max_r2 = max(R2, na.rm = TRUE))

# 使用 left_join 合并，保持原始数据顺序
all_data <- left_join(all_data, max_r2_per_site, by = "Site")

# 计算标注的 x 位置（基于各分面最大 R2 的固定比例）
all_data$label_x <- all_data$R2 + all_data$max_r2 * 0.03

# 创建分面柱状图，排成一行
p <- ggplot(all_data, aes(x = R2*100, y = Factor, fill = Site)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Significance, x = label_x, group = 1),
            vjust = 0.5, hjust = 0, size = 3.5, color = "black", fontface = "bold") +
  facet_grid(~Site, scales = "free_x") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = colors) +
  labs(x = paste(expression(R^2), "(%)"),
       y = NULL) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(5, 20, 5, 5)
  )

# 保存图片
png_width <- 7
png_height <- 2.5
ggsave("adonis_results_facet_brays.png", p, width = png_width, height = png_height, dpi = 300)
ggsave("adonis_results_facet_brays.pdf", p, width = png_width, height = png_height)

print(p)

# 输出排序后的因素
print("Factors ordered by Saliva R2 (high to low):")
print(factor_order)

# --- 额外：创建带图例的版本 ---
p_legend <- ggplot(all_data, aes(x = R2*100, y = Factor, fill = Site)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Significance, x = label_x, group = 1),
            vjust = 0.5, hjust = 0, size = 3.5, color = "black", fontface = "bold") +
  facet_grid(~Site, scales = "free_x") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = colors, name = "Site") +
  labs(x = paste(expression(R^2), "(%)"),
       y = NULL) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(5, 20, 5, 5)
  )

ggsave("adonis_results_facet_brays_with_legend.png", p_legend, width = 8, height = 3, dpi = 300)
ggsave("adonis_results_facet_brays_with_legend.pdf", p_legend, width = 8, height = 3)

print("Plots saved successfully!")
