#!/usr/bin/env Rscript
# ==============================================================================
# 01_composite_figure.R
# Generate composite figure: HOMD Coverage Analysis
# Target: Nature/Cell level, 5 panels (A-E), 183mm full width
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)
library(ggrepel)

BASE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"
OUTDIR <- file.path(BASE, "pipeline", "analysis_report")
RDS <- file.path(OUTDIR, "R", "validated_data.rds")

# ==============================================================================
# Load validated data
# ==============================================================================
cat("Loading validated data ...\n")
d <- readRDS(RDS)

# ==============================================================================
# Color palette — Nature-friendly, color-blind safe
# ==============================================================================
cols_capture <- c("HOMD captured" = "#4DAF4A", "HOMD uncaptured" = "#E41A1C")
cols_reason <- c("Low identity" = "#E41A1C", "No HOMD hit" = "#FF7F00", "Low query coverage" = "#984EA3")
cols_novelty <- c("novel_all_refs" = "#E41A1C", "known_not_in_homd" = "#377EB8")
cols_region <- c("V3-V4" = "#1B9E77", "V4" = "#D95F02", "V1-V3" = "#7570B3",
                 "V4-V5" = "#E7298A", "V6-V8" = "#66A61E", "Other" = "#999999")

theme_pub <- theme_classic(base_size = 8) +
  theme(
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, face = "bold"),
    strip.text = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 9, face = "bold"),
    plot.margin = margin(4, 4, 4, 4)
  )

# ==============================================================================
# Panel A: Data overview — horizontal stacked bar
# Show: Total derep → Captured vs Uncaptured
# ==============================================================================
cat("Panel A: Data overview ...\n")

panel_a_data <- tibble(
  category = factor(c("All dereplicated ASVs"), levels = "All dereplicated ASVs"),
  captured = d$n_homd_captured,
  uncaptured = d$n_homd_uncaptured
) %>%
  mutate(total = captured + uncaptured,
         pct_captured = captured / total * 100,
         pct_uncaptured = uncaptured / total * 100)

# Long format for stacking
panel_a_long <- panel_a_data %>%
  select(category, captured, uncaptured) %>%
  pivot_longer(cols = c(captured, uncaptured), names_to = "status", values_to = "count") %>%
  mutate(label = ifelse(status == "captured", "HOMD captured", "HOMD uncaptured"),
         pct = count / sum(panel_a_data$total) * 100)

p_a <- ggplot(panel_a_long, aes(x = category, y = count, fill = label)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(count, big.mark = ","), pct)),
            position = position_stack(vjust = 0.5), size = 2.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = cols_capture, name = NULL) +
  scale_y_continuous(labels = label_comma(scale = 1e-6, suffix = "M")) +
  coord_flip() +
  labs(x = NULL, y = "Number of dereplicated ASVs",
       title = "A  HOMD Capture Overview") +
  theme_pub +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# ==============================================================================
# Panel B: HOMD capture by 16S region — grouped bar
# ==============================================================================
cat("Panel B: HOMD capture by region ...\n")

# Select main regions (those with >50K ASVs)
main_regions <- d$homd_region_df %>%
  filter(total > 50000) %>%
  mutate(region = factor(region, levels = region)) %>%
  arrange(desc(total))

panel_b_data <- main_regions %>%
  select(region, captured, uncaptured) %>%
  pivot_longer(cols = c(captured, uncaptured), names_to = "status", values_to = "count") %>%
  mutate(label = ifelse(status == "captured", "HOMD captured", "HOMD uncaptured"),
         capture_rate = ifelse(status == "captured",
                               main_regions$capture_rate[match(region, main_regions$region)],
                               100 - main_regions$capture_rate[match(region, main_regions$region)]))

p_b <- ggplot(panel_b_data, aes(x = region, y = count, fill = label)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = cols_capture, name = NULL) +
  scale_y_continuous(labels = label_comma(scale = 1e-3, suffix = "K")) +
  labs(x = "16S rRNA amplicon region", y = "Number of ASVs",
       title = "B  Capture Rate by Amplicon Region") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        legend.position = "bottom")

# ==============================================================================
# Panel C: Uncaptured reasons — donut chart
# ==============================================================================
cat("Panel C: Uncaptured reasons ...\n")

panel_c_data <- d$uncaptured_reasons %>%
  mutate(reason = case_when(
    reason == "low_identity" ~ "Low identity\n(<97%)",
    reason == "no_hit" ~ "No HOMD hit",
    reason == "low_qcov" ~ "Low query\ncoverage (<99%)",
    TRUE ~ reason
  ),
  pct = count / sum(count) * 100,
  label = sprintf("%s\n%s (%.1f%%)", reason, format(count, big.mark = ","), pct))

p_c <- ggplot(panel_c_data, aes(x = 2, y = count, fill = reason)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  scale_fill_manual(values = cols_reason, name = NULL) +
  labs(title = "C  Why ASVs Escape HOMD") +
  geom_text(aes(label = ifelse(pct > 3, label, "")),
            position = position_stack(vjust = 0.5), size = 2.2, color = "white", fontface = "bold") +
  theme_pub +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 6.5))

# ==============================================================================
# Panel D: Candidate novel phylotype landscape — bubble scatter
# ==============================================================================
cat("Panel D: Phylotype landscape ...\n")

panel_d_data <- d$phylotypes %>%
  mutate(novelty_label = ifelse(novelty_class == "novel_all_refs",
                                "Novel to all references", "Known, absent from HOMD"))

p_d <- ggplot(panel_d_data, aes(x = max_n_bioprojects, y = max_subject_prevalence)) +
  geom_point(aes(size = sum_total_reads, color = novelty_label), alpha = 0.7) +
  scale_size_continuous(range = c(1.5, 8), labels = label_comma(scale = 1e-3, suffix = "K"),
                        name = "Total reads") +
  scale_color_manual(values = c("Novel to all references" = "#E41A1C",
                                 "Known, absent from HOMD" = "#377EB8"),
                     name = NULL) +
  geom_text_repel(
    data = . %>% filter(max_subject_prevalence > 1200 | max_n_bioprojects > 28),
    aes(label = phylotype_id), size = 2, max.overlaps = 10,
    segment.color = "grey50", segment.size = 0.3
  ) +
  labs(x = "Number of BioProjects", y = "Number of Subjects",
       title = "D  Novel Phylotype Landscape") +
  theme_pub +
  theme(legend.position = "bottom",
        legend.box = "vertical")

# ==============================================================================
# Panel E: Novelty classification — simple bar
# ==============================================================================
cat("Panel E: Novelty classification ...\n")

panel_e_data <- d$phylotypes %>%
  count(novelty_class) %>%
  mutate(label = ifelse(novelty_class == "novel_all_refs",
                        "Novel to all\nreferences", "Known, absent\nfrom HOMD"),
         pct = n / sum(n) * 100)

p_e <- ggplot(panel_e_data, aes(x = reorder(label, n), y = n, fill = novelty_class)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d (%.0f%%)", n, pct)),
            hjust = -0.2, size = 3, fontface = "bold") +
  scale_fill_manual(values = cols_novelty, guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Number of candidate\nnovel phylotypes",
       title = "E  Novelty Classification") +
  theme_pub +
  ylim(0, max(panel_e_data$n) * 1.25)

# ==============================================================================
# Assemble composite figure
# ==============================================================================
cat("Assembling composite figure ...\n")

composite <- (p_a | p_b) / (p_c | p_d | p_e) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "",
    tag_suffix = "",
    theme = theme(
      plot.tag = element_text(size = 10, face = "bold", hjust = 0)
    )
  )

# ==============================================================================
# Save outputs
# ==============================================================================
cat("Saving outputs ...\n")

ggsave(file.path(OUTDIR, "fig_homd_coverage.pdf"),
       composite, width = 183, height = 130, units = "mm", dpi = 300)

ggsave(file.path(OUTDIR, "fig_homd_coverage.png"),
       composite, width = 183, height = 130, units = "mm", dpi = 300)

cat(sprintf("Saved: fig_homd_coverage.pdf\n"))
cat(sprintf("Saved: fig_homd_coverage.png\n"))
cat("\n=== Composite Figure Complete ===\n")
