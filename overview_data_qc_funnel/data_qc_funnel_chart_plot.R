#!/usr/bin/env Rscript
# Data QC Funnel Chart
# Dataset-level-only QC pipeline funnel plot (v4).
# Based on data_qc_pipeline.md: 3 filter steps + oral subset.
#
# Output: data_qc_funnel_chart.pdf / data_qc_funnel_chart.png

library(tidyverse)
library(scales)

# ---------------------------------------------------------------------------
# Global theme & font sizing
# ---------------------------------------------------------------------------
base_sz <- 17
text_sz <- base_sz / .pt

update_geom_defaults("text", list(family = "sans"))
update_geom_defaults("label", list(family = "sans"))

theme_set(theme_minimal(base_size = base_sz, base_family = "sans"))

# ---------------------------------------------------------------------------
# Funnel data: 4 steps (Raw + 3 filters)
#   Step 0: Raw merge (116,099 samples, 334 studies)
#   Step 1: Chimera filter (114,246 samples, 312 studies)
#   Step 2: CR mapping filter (109,257 samples, 288 studies)
#   Step 3: Oral subset (81,069 samples, 288 studies)
# ---------------------------------------------------------------------------
funnel_data <- data.frame(
  step     = 1:4,
  n_inside = c("116,099", "114,246", "109,257", "81,069"),
  s_inside = c("334 studies", "312 studies", "288 studies", "288 studies"),
  n_num    = c(116099, 114246, 109257, 81069),
  y        = c(12.0, 8.0, 4.0, 0.0)
)

right_anno <- data.frame(
  y_mid = c(10.0, 6.0, 2.0),
  line1 = c("Step 1: Chimera Filter",
            "Step 2: CR Mapping Filter",
            "Step 3: Oral Subset"),
  line2 = c("(>30% oral samples with chimera >25%)",
            "(oral aggregate GG2 CR mapping <70%)",
            "(remove non-oral samples)"),
  line3 = c("\u221222 studies, \u22121,853 samples",
            "\u221224 studies, \u22124,989 samples",
            "\u221228,188 non-oral samples")
)

max_w <- 13.5
min_w <- 8.5
bar_h <- 0.60

funnel_data$half_w <- min_w +
  (funnel_data$n_num - min(funnel_data$n_num)) /
  (max(funnel_data$n_num) - min(funnel_data$n_num)) * (max_w - min_w)

poly_blocks <- do.call(rbind, lapply(1:nrow(funnel_data), function(i) {
  b <- funnel_data[i, ]
  hw <- b$half_w
  yc <- b$y
  data.frame(
    x    = c(-hw, hw, hw, -hw),
    y    = c(yc - bar_h, yc - bar_h, yc + bar_h, yc + bar_h),
    step = b$step,
    grp  = paste0("bar_", i)
  )
}))

conn_indices <- c(1, 2, 3)
conn_blocks <- do.call(rbind, lapply(conn_indices, function(i) {
  top_hw      <- funnel_data$half_w[i]
  bot_hw      <- funnel_data$half_w[i + 1]
  top_y_bottom <- funnel_data$y[i] - bar_h
  bot_y_top   <- funnel_data$y[i + 1] + bar_h
  data.frame(
    x    = c(-top_hw, top_hw, bot_hw, -bot_hw),
    y    = c(top_y_bottom, top_y_bottom, bot_y_top, bot_y_top),
    step = i,
    grp  = paste0("conn_", i)
  )
}))

arrow_data <- data.frame(
  x       = 0,
  y_start = funnel_data$y[conn_indices] - bar_h,
  y_end   = funnel_data$y[conn_indices + 1] + bar_h
)

bar_fills      <- c("#263238", "#455A64", "#1976D2", "#43A047")
bar_text_colors <- rep("white", 4)

# Annotation box sizing — tight fit to text content
anno_x <- max_w + 2.0
anno_w <- 13.0
anno_h <- 1.5

right_anno$rect_xmin <- anno_x - 0.2
right_anno$rect_xmax <- anno_x + anno_w + 0.5
right_anno$rect_ymin <- right_anno$y_mid - anno_h / 2
right_anno$rect_ymax <- right_anno$y_mid + anno_h / 2

p_a <- ggplot() +
  # Dataset-Level Filters background
  annotate("rect",
           xmin = -max_w - 0.5, xmax = anno_x + anno_w + 0.5,
           ymin = -1.5, ymax = 14.0,
           fill = "#F8F9FA", color = "#B0BEC5",
           linetype = "dashed", alpha = 0.5) +
  annotate("text",
           x = anno_x + anno_w - 1.0, y = 13.5,
           label = "Dataset-Level Filters",
           fontface = "bold", hjust = 1,
           size = text_sz * 1.1, color = "#37474F") +
  # Connector trapezoids
  geom_polygon(data = conn_blocks,
               aes(x = x, y = y, group = grp),
               fill = "#CFD8DC", alpha = 0.45) +
  # Main funnel bars
  geom_polygon(data = poly_blocks,
               aes(x = x, y = y, group = grp),
               fill = rep(bar_fills, each = 4), color = NA) +
  # White dividers on bars
  geom_segment(data = funnel_data,
               aes(x = -half_w, xend = half_w, y = y + bar_h, yend = y + bar_h),
               color = "white", linewidth = 0.4) +
  geom_segment(data = funnel_data,
               aes(x = -half_w, xend = half_w, y = y - bar_h, yend = y - bar_h),
               color = "white", linewidth = 0.4) +
  # Inter-step arrows — open arrowhead, sleek style
  geom_segment(data = arrow_data,
               aes(x = x, y = y_start, xend = x, yend = y_end),
               arrow = arrow(length = unit(0.28, "cm"), type = "closed", angle = 18),
               linewidth = 0.6, color = "#455A64") +
  # Sample count labels inside bars
  geom_text(data = funnel_data,
            aes(x = 0, y = y + 0.25, label = n_inside),
            size = text_sz * 1.0, color = bar_text_colors,
            fontface = "bold", hjust = 0.5) +
  # Study count labels inside bars
  geom_text(data = funnel_data,
            aes(x = 0, y = y - 0.30, label = s_inside),
            size = text_sz * 0.75, color = bar_text_colors,
            fontface = "plain", hjust = 0.5) +
  # Right-side annotation boxes
  geom_rect(data = right_anno,
            aes(xmin = rect_xmin, xmax = rect_xmax,
                ymin = rect_ymin, ymax = rect_ymax),
            fill = "#FFFFFF", color = "#B0BEC5",
            linewidth = 0.4, alpha = 0.95) +
  geom_text(data = right_anno,
            aes(x = anno_x + 0.4, y = y_mid + 0.45, label = line1),
            size = text_sz * 0.85, fontface = "bold",
            color = "#212121", hjust = 0) +
  geom_text(data = right_anno,
            aes(x = anno_x + 0.4, y = y_mid, label = line2),
            size = text_sz * 0.65, fontface = "italic",
            color = "#616161", hjust = 0) +
  geom_text(data = right_anno,
            aes(x = anno_x + 0.4, y = y_mid - 0.45, label = line3),
            size = text_sz * 0.80, fontface = "bold",
            color = "#D32F2F", hjust = 0) +
  # Dashed connector lines
  geom_segment(data = right_anno,
               aes(x = max_w, xend = anno_x - 0.15,
                   y = y_mid, yend = y_mid),
               color = "#78909C", linewidth = 0.5,
               alpha = 0.8, linetype = "dashed") +
  # Coordinate and theme
  coord_cartesian(xlim = c(-14.5, 33.0), ylim = c(-1.5, 14.5), clip = "off") +
  theme_void(base_size = base_sz, base_family = "sans") +
  theme(plot.margin = margin(10, 5, 10, 5, unit = "pt"))

# ---------------------------------------------------------------------------
# Save output
# ---------------------------------------------------------------------------
output_dir <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-pipeline-map/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "data_qc_funnel_chart.pdf"),
       plot = p_a, width = 14, height = 8, device = "pdf")
ggsave(file.path(output_dir, "data_qc_funnel_chart.png"),
       plot = p_a, width = 14, height = 8, device = "png", dpi = 600, bg = "white")

message("Data QC funnel chart saved to: ", file.path(output_dir, "data_qc_funnel_chart.pdf"))
