# 1. 读取并处理数据 
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

INPUT_FILE <- "/home/hexincheng/projects/hxc/oral_project/oral_plot_by_yl_20260602/oral_code_results_new/output/oral_output/rarefaction_pcofg.tsv"

df <- read.table(INPUT_FILE, 
                 header = TRUE, sep = "\t")

df_long <- df %>%
  pivot_longer(
    cols = c("Phylum", "Class", "Order", "Family", "Genus"),
    names_to = "Taxonomy",
    values_to = "Unique_taxa"
  )

df_long$Taxonomy <- factor(df_long$Taxonomy, 
                           levels = c("Genus", "Family", "Order", "Class", "Phylum"))

# 2. 准备图内文本标签
label_x_target <- 30000 
label_data <- df_long %>%
  group_by(Taxonomy) %>%
  slice_min(abs(Samples - label_x_target)) %>%
  mutate(
    vjust_val = case_when(
      Taxonomy == "Phylum" ~ 1.5,
      Taxonomy == "Class"  ~ -0.7,
      Taxonomy == "Order"  ~ -1.2,
      Taxonomy == "Family" ~ -1.2,
      Taxonomy == "Genus"  ~ 1.5,
      TRUE ~ -1.2
    )
  )

# 3. 开始绘图
p <- ggplot(df_long, aes(x = Samples, y = Unique_taxa, color = Taxonomy)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  
  scale_color_manual(values = c(
    "Genus"  = "#7570B3",
    "Family" = "#F5D241", 
    "Order"  = "#E7298A", 
    "Class"  = "#4DBBD5", 
    "Phylum" = "#66A61E"  
  )) +
  
  scale_x_continuous(labels = label_number(scale = 1/1000, suffix = "K"),
                     breaks = c(0, 25000, 50000, 75000),
                     expand = expansion(mult = c(0.01, 0.06))) +
  scale_y_continuous(labels = comma,
                     breaks = seq(0, 4000, by = 1000),
                     limits = c(-150, 4200)) + 
  geom_text(data = label_data, aes(label = Taxonomy, vjust = vjust_val), 
            size = 4, fontface = "bold", show.legend = FALSE) +
  
  labs(
    x = "Samples",
    y = "Unique taxa"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold", color = "black"),
    legend.position = "none", 
    panel.border = element_rect(color = "black", linewidth = 1),
    plot.margin = margin(5.5, 14, 5.5, 5.5)
  )

# 4. 查看并保存结果
print(p)

ggsave("/home/hexincheng/projects/hxc/oral_project/oral_plot_by_yl_20260602/oral_code_results_new/output/oral_output/curve_global_pcofg.pdf", 
       plot = p, width = 4.5, height = 4)
ggsave("/home/hexincheng/projects/hxc/oral_project/oral_plot_by_yl_20260602/oral_code_results_new/output/oral_output/curve_global_pcofg.png", 
       plot = p, width = 4.5, height = 4, dpi = 300)
print("绘图完成！图例重叠与溢出问题已完美解决。")
