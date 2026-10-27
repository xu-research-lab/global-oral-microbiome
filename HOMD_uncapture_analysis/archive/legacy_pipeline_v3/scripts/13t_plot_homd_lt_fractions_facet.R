#!/usr/bin/env Rscript
# pipeline_v3: HOMD identity curves, axes swapped, "fraction with identity < x",
# x-axis reversed (100% on the left -> 0% on the right).
# One row of four panels (facet_grid, shared x-axis):
#   [All phyla combined] [ASV count] [Reads] [Subject prevalence]
# x = Fraction (0-100%) of the metric with identity < threshold (normal axis);
# y = Best identity to HOMD reference (%) 70-100.
# 97% identity dashed horizontal line in all panels.
# Two colour scales via ggnewscale: phylum (15) + metric (3).
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggnewscale)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

phy <- read.delim(file.path(PIPE_V3, "13r_homd_phylum_lt_fractions.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
tot <- read.delim(file.path(PIPE_V3, "13r_homd_total_lt_fractions.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cov <- read.delim(file.path(PIPE_V3, "13d_homd_phylum_coverage.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# ---------------- phylum metadata (legend order: frac_ge97 high -> low) ----------------
meta <- cov %>% transmute(
  phylum    = phylum,
  frac_ge97 = n_asvs_ge97 / n_asvs_total
) %>% arrange(desc(frac_ge97))
meta$label <- sprintf("%s (%.1f%%)", meta$phylum, 100 * meta$frac_ge97)
meta$label <- factor(meta$label, levels = meta$label)
PAL <- hue_pal()(15)

# ---------------- long data (fraction with identity < x) ----------------
phy <- phy %>% filter(identity >= 70 & identity <= 100) %>%
  left_join(meta, by = "phylum")

phy_long <- bind_rows(
  phy %>% transmute(panel = "ASV count", identity, frac = frac_n_lt, label),
  phy %>% transmute(panel = "Reads", identity, frac = frac_r_lt, label),
  phy %>% transmute(panel = "Subject prevalence", identity, frac = frac_s_lt, label)
)
phy_long$label <- factor(phy_long$label, levels = levels(meta$label))

tot_long <- tot %>% filter(identity >= 70 & identity <= 100) %>%
  pivot_longer(cols = c(frac_n_lt, frac_r_lt, frac_s_lt),
               names_to = "metric", values_to = "frac") %>%
  mutate(metric = recode(metric,
                         frac_n_lt = "ASV count",
                         frac_r_lt = "Reads",
                         frac_s_lt = "Subject prevalence"),
         metric = factor(metric, levels = c("ASV count", "Reads", "Subject prevalence")),
         panel = "All phyla combined")

panel_levels <- c("All phyla combined", "ASV count", "Reads", "Subject prevalence")
phy_long$panel <- factor(phy_long$panel, levels = panel_levels)
tot_long$panel <- factor(tot_long$panel, levels = panel_levels)

# ---------------- plot (x = fraction identity<x, reversed: 100% left -> 0% right) ----------------
p <- ggplot() +
  geom_hline(yintercept = 97, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_line(data = phy_long, aes(x = frac, y = identity, colour = label), linewidth = 0.6) +
  scale_colour_manual(values = PAL,
                      name = "Phylum (ASV counts fraction of identity>=97%)") +
  new_scale_colour() +
  geom_line(data = tot_long, aes(x = frac, y = identity, colour = metric), linewidth = 0.9) +
  scale_colour_manual(values = c("ASV count" = "grey40",
                                 "Reads" = "#1f77b4",
                                 "Subject prevalence" = "#d62728"),
                      name = "Metric") +
  facet_grid(. ~ panel) +
  scale_x_continuous(labels = percent_format(), limits = c(0, 1),
                  name = "Fraction of ASVs with identity < threshold") +
  scale_y_continuous(breaks = c(70, 75, 80, 85, 90, 95, 100), limits = c(70, 100),
                     name = "Best identity to HOMD reference (%)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 9))

fig <- file.path(FIGOUT, "fig_v3_homd_identity_lt_fraction_facet")
ggsave(paste0(fig, ".png"), p, width = 10, height = 3.5, dpi = 300)
ggsave(paste0(fig, ".pdf"), p, width = 10, height = 3.5)
message("[13t] faceted 4-panel figure (x = fraction identity<x, reversed axis) written")
