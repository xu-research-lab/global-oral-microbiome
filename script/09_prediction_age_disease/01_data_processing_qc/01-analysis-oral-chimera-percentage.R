# ==============================================================================
# 批量提取嵌合体占比并复刻 Ranked Area Line Plot
# ==============================================================================

# 加载所需的 R 包
# 如果缺少包，请在 R 中运行: install.packages(c("tidyverse", "ggrepel", "scales"))
library(tidyverse)
library(ggrepel)
library(scales)

# ==============================================================================
# 1. 定义路径与文件变量
# ==============================================================================
project_list_file <- "/home/yanglun/YSD/global_oral/global_oral_v2/project_list.txt"
base_dir <- "/home/yanglun/YSD/global_oral/global_oral_v2/projects"

out_dir <- "/home/yanglun/YSD/global_oral/global_oral_v2/results/2605"
out_csv <- file.path(out_dir, "oral_chim_perc_YL_260519.csv")
out_log <- file.path(out_dir, "chim_perc_log_YL_260519.txt")
out_plot_pdf <- file.path(out_dir, "oral_chim_perc_rank_YL_260519.pdf")
out_plot_png <- file.path(out_dir, "oral_chim_perc_rank_YL_260519.png")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# 2. 读取项目列表并提取【每个样本】的数据
# ==============================================================================
cat("正在读取项目列表...\n")
projects <- readLines(project_list_file)
projects <- projects[projects != ""]

all_data <- list()
log_messages <- c()
success_proj_count <- 0
total_sample_count <- 0

for (proj in projects) {
  file_path <- file.path(base_dir, proj, "results", "summary.tsv")
  
  if (file.exists(file_path)) {
    df <- tryCatch({
      read.delim(file_path, row.names = 1, check.names = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(df) && "chim_perc" %in% colnames(df)) {
      temp_df <- data.frame(
        Project = proj,
        SampleID = rownames(df),
        chim_perc = as.numeric(df$chim_perc),
        stringsAsFactors = FALSE
      )
      all_data[[proj]] <- temp_df
      success_proj_count <- success_proj_count + 1
      total_sample_count <- total_sample_count + nrow(temp_df)
      log_messages <- c(log_messages, paste("[成功] 文件夹:", proj, "- 提取样本数:", nrow(temp_df)))
    } else {
      log_messages <- c(log_messages, paste("[警告] 文件夹:", proj, "- 存在，但缺少 chim_perc 列或读取失败。"))
    }
  } else {
    log_messages <- c(log_messages, paste("[失败] 文件夹:", proj, "- 找不到 summary.tsv 文件。"))
  }
}

# 合并并保存为 CSV 输入文件
final_df <- bind_rows(all_data)
write.csv(final_df, out_csv, row.names = FALSE, quote = FALSE)
cat("样本级数据已保存至:", out_csv, "\n")

# 输出日志文件
log_summary <- paste(
  "========================================\n",
  "数据提取总结报告\n",
  "========================================\n",
  "要求检查的文件夹总数 :", length(projects), "\n",
  "成功提取的文件夹数   :", success_proj_count, "\n",
  "失败/缺失的文件夹数  :", length(projects) - success_proj_count, "\n",
  "最终获取的总样本数量 :", total_sample_count, "\n\n",
  "--- 详细日志 ---\n", sep = ""
)
writeLines(c(log_summary, log_messages), out_log)
cat("日志记录已保存至:", out_log, "\n")

# ==============================================================================
# 3. 读取 CSV 并进行作图数据准备 (计算各项目嵌合体中位数并排序)
# ==============================================================================
cat("正在生成图片...\n")

# 从刚刚保存的 CSV 中读取数据（模拟完整的分析流）
plot_data <- read.csv(out_csv)

# 聚合到 Project 层面：计算中位数，并降序排列
project_stats <- plot_data %>%
  group_by(Project) %>%
  summarise(
    median_chim_perc = median(chim_perc, na.rm = TRUE),
    sample_count = n() # 记录每个项目的样本数
  ) %>%
  arrange(desc(median_chim_perc)) %>%
  mutate(Rank = row_number()) # 生成 X 轴的排序编号

# 提取 Top 10 和 Bottom 5 用于高亮显示
top10_data <- project_stats %>% slice_head(n = 10)
bottom5_data <- project_stats %>% slice_tail(n = 5)

# 动态生成图表副标题 (与参考图保持一致)
subtitle_text <- sprintf("Total: %d bioprojects, %d samples | Red = top 10, Green = bottom 5", 
                         success_proj_count, total_sample_count)

# ==============================================================================
# 4. 使用 ggplot2 完美复刻参考图形式
# ==============================================================================
p <- ggplot(project_stats, aes(x = Rank, y = median_chim_perc)) +
  # 1. 下方的浅蓝色阴影区域
  geom_area(fill = "#F4F7F9", alpha = 0.8) +
  # 2. 主折线 (蓝色)
  geom_line(color = "#1F77B4", linewidth = 1) +
  
  # 3. 底部 5 个项目打上绿点
  geom_point(data = bottom5_data, aes(x = Rank, y = median_chim_perc), 
             color = "#2CA02C", size = 2) +
  
  # 4. 顶部 10 个项目打上红点
  geom_point(data = top10_data, aes(x = Rank, y = median_chim_perc), 
             color = "#D62728", size = 2) +
  
  # 5. 为 Top 10 添加引线标签 (参考图中带有具体数值的形式)
  geom_text_repel(
    data = top10_data,
    aes(label = paste0(Project, " (", round(median_chim_perc, 2), "%)")),
    nudge_x = 20,          # 标签统一向右偏
    nudge_y = max(project_stats$median_chim_perc) * 0.05, # 标签微调向上偏
    segment.color = "grey50", # 引线颜色
    segment.size = 0.3,
    size = 3.5,
    hjust = 0,             # 文本左对齐
    direction = "y"        # 尽量在Y轴方向错开以防重叠
  ) +
  
  # 6. 设置主题和背景格式，复刻图片中的水平网格线
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 12, color = "black", margin = margin(b = 15)),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major.x = element_blank(), # 去除垂直网格线
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5"), # 保留水平灰色网格线
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  
  # 7. 坐标轴标签设置
  labs(
    title = "Median Chimera Percentage per Bioproject (descending)",
    subtitle = subtitle_text,
    x = "Bioproject Rank",
    y = "Median Chimera Percentage (%)"
  ) +
  
  # 8. Y 轴调整：采用 pseudo_log_trans 处理含 0 的百分比数据，呈现 log 视觉效果
  # 若你不想用对数轴（因为百分比通常用线性轴），可以将下面这行注释掉
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 20, 50, 100), # 自定义 Y 轴刻度
    labels = c("0", "1", "5", "10", "20", "50", "100")
  ) +
  scale_x_continuous(expand = c(0, 0)) # 让图形两端紧贴坐标轴

# ==============================================================================
# 5. 保存图表
# ==============================================================================
ggsave(out_plot_pdf, plot = p, width = 12, height = 7)
ggsave(out_plot_png, plot = p, width = 12, height = 7, dpi = 300)

cat("图表已成功保存为:\n", out_plot_pdf, "\n", out_plot_png, "\n")
cat("全流程运行完毕！\n")