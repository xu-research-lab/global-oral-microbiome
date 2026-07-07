library(tidyverse)
library(viridis)
library(grid)

# 1. Use the current directory for both input and output.
current_dir <- "E:/oral/26/metadata-2604/Oral-Metadata-merge-YL/results/260519/260522"
input_csv <- file.path(current_dir, "oral_chim_perc_YL_260519.csv")
out_dir <- "E:/oral/26/metadata-2604/Oral-Metadata-merge-YL/results/260519/260522/results"

out_plot_pdf <- file.path(out_dir, "oral_chim_perc_boxplot_ppt_4rows_desc.pdf")
out_plot_png <- file.path(out_dir, "oral_chim_perc_boxplot_ppt_4rows_desc.png")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# 2. Read data.
cat("Reading CSV data...\n")
df <- read.csv(input_csv)

# 3. Rank projects by median chimera percentage and split them into four panels.
cat("Calculating project-level summary statistics...\n")

project_stats <- df %>%
  group_by(Project) %>%
  summarise(
    median_chim = median(chim_perc, na.rm = TRUE),
    sample_count = n(),
    
    # 统计该项目中，嵌合体 > 25% 的样本数量
    high_chim_count = sum(chim_perc > 25, na.rm = TRUE),
    
    # 计算该项目中，嵌合体 > 25% 的样本占比 (%)
    high_chim_perc = (high_chim_count / sample_count) * 100,
    
    .groups = "drop"
  ) %>%
  arrange(desc(median_chim), Project) %>%
  mutate(
    project_rank = row_number(),
    panel_row = paste("Row", ntile(project_rank, 4))
  )

# ==========================================
# 提取“超过30%样本的嵌合体>25%”的数据集
# ==========================================
flagged_projects <- project_stats %>%
  filter(high_chim_perc > 30) %>%
  select(Project, sample_count, median_chim, high_chim_count, high_chim_perc) %>%
  arrange(desc(high_chim_perc))

cat("\n=== 异常预警: 超过30%样本嵌合体占比>25%的项目 ===\n")
print(flagged_projects)
write.csv(flagged_projects, 
          file.path(out_dir, "flagged_projects_gt30perc_samples_gt25chim.csv"), 
          row.names = FALSE)

# ==========================================
# 提取图表最左侧（中位数最高）的前30个项目
# ==========================================
top_30_projects <- project_stats %>%
  slice_head(n = 30) %>%
  select(project_rank, Project, sample_count, median_chim, high_chim_perc)

cat("\n=== 图表最左侧 (Top 30) 异常数据集 ===\n")
print(top_30_projects)
write.csv(top_30_projects, 
          file.path(out_dir, "top_30_highest_median_projects.csv"), 
          row.names = FALSE)
# ==========================================

plot_data <- df %>%
  left_join(project_stats, by = "Project") %>%
  mutate(
    Project = factor(Project, levels = project_stats$Project),
    panel_row = factor(panel_row, levels = paste("Row", 1:4))
  )

anno_data <- project_stats %>%
  mutate(
    Project = factor(Project, levels = project_stats$Project),
    panel_row = factor(panel_row, levels = paste("Row", 1:4)),
    is_flagged = high_chim_perc > 30,
    label_text = ifelse(high_chim_perc == 0, "0%", sprintf("%.1f%%", high_chim_perc))
  )

# 4. Use four independent stacked panels so each row reflows its own projects.
cat("\nGenerating four-row PPT-friendly figure...\n")

p <- ggplot(plot_data, aes(x = Project, y = chim_perc)) +
  
  geom_jitter(
    width = 0.25,        
    size = 0.7,          
    alpha = 0.8,         
    color = "#4A90E2",   
    stroke = 0           
  ) +
  
  geom_boxplot(
    fill = "transparent", 
    outlier.shape = NA,   
    linewidth = 0.45,     
    color = "#2C3E50"    
  ) +
  
  geom_text(
    data = anno_data,
    aes(x = Project, y = 101, label = label_text, color = is_flagged), 
    angle = 45,
    hjust = 0,
    vjust = 0,
    size = 2.4,          
    fontface = "bold"    
  ) +
  
  scale_color_manual(
    values = c(
      "TRUE" = "#E74C3C",  
      "FALSE" = "grey50"   
    )
  ) +
  
  facet_wrap(~panel_row, ncol = 1, scales = "free_x") +
  theme_minimal() +
  labs(
    title = "Distribution of Chimera Percentages Across All Samples",
    subtitle = paste(
      "Total Projects:", length(unique(plot_data$Project)),
      "| Total Samples:", nrow(plot_data),
      "| Ranked by median chimera percentage (descending)\n",
      "🔴 Red text indicates projects where >30% of samples have >25% chimera" 
    ),
    x = "Bioprojects",
    y = "Chimera Percentage of Individual Samples (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 4.8, color = "black"),
    
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 11, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E0E0E0", linewidth = 0.4),
    panel.grid.minor.y = element_blank(),
    panel.spacing.y = unit(0.35, "lines"),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 9.5, color = "grey30", margin = margin(b = 8)),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none", 
    plot.margin = margin(t = 8, r = 10, b = 8, l = 8)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    limits = c(0, 122),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_discrete(drop = TRUE, expand = expansion(add = c(0.15, 0.15)))

# 5. Save outputs with a 16:9 PPT-friendly canvas.
cat("Saving outputs to current directory...\n")

ggsave(out_plot_pdf, plot = p, width = 13.333, height = 7.5, limitsize = FALSE)
ggsave(out_plot_png, plot = p, width = 13.333, height = 7.5, dpi = 300, limitsize = FALSE)

cat("Completed. Projects now reflow within each of the four stacked rows and outputs were saved in the current directory.\n")