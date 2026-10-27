#!/usr/bin/env Rscript
# pipeline_v3: HOMD identity "fraction with identity < x" curves.
# Four panels in one row (aplot, shared legends):
#   p1: ASV count      (per phylum)   y = fraction of ASVs with identity < x
#   p2: Reads          (per phylum)   y = fraction of reads with identity < x
#   p3: Subject prevalence (per phylum) y = fraction of subject prevalence with identity < x
#   p4: Total (all phyla combined)    y = fraction of ASVs / reads / subject prevalence with identity < x
# x = 70-100% identity for all panels; y = 0-100%.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(aplot)
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

phy <- phy %>% filter(identity >= 70 & identity <= 100) %>%
  left_join(meta, by = "phylum")
phy$label <- factor(phy$label, levels = levels(meta$label))

tot_long <- tot %>% filter(identity >= 70 & identity <= 100) %>%
  pivot_longer(cols = c(frac_n_lt, frac_r_lt, frac_s_lt),
               names_to = "metric", values_to = "frac") %>%
  mutate(metric = recode(metric,
                         frac_n_lt = "ASV count",
                         frac_r_lt = "Reads",
                         frac_s_lt = "Subject prevalence"),
         metric = factor(metric, levels = c("ASV count", "Reads", "Subject prevalence")))

XBR <- c(70, 75, 80, 85, 90, 95, 100)
base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 10, hjust = 0.5))

# ---- p1: ASV count (per phylum) ----
p1 <- ggplot(phy, aes(x = identity, y = frac_n_lt, colour = label)) +
  geom_vline(xintercept = 97, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_line(linewidth = 0.7) +
  ggtitle("ASV count") +
  scale_x_continuous(breaks = XBR, limits = c(70, 100),
                     name = "Best identity to HOMD reference (%)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1),
                     name = "Fraction of ASVs (identity < x)") +
  scale_colour_manual(values = PAL, name = "Phylum") +
  base_theme + theme(legend.position = "right")

# ---- p2: reads (per phylum) ----
p2 <- ggplot(phy, aes(x = identity, y = frac_r_lt, colour = label)) +
  geom_vline(xintercept = 97, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_line(linewidth = 0.7) +
  ggtitle("Reads") +
  scale_x_continuous(breaks = XBR, limits = c(70, 100),
                     name = "Best identity to HOMD reference (%)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1),
                     name = "Fraction of reads (identity < x)") +
  scale_colour_manual(values = PAL, name = "Phylum") +
  base_theme + theme(legend.position = "right")

# ---- p3: subject prevalence (per phylum) ----
p3 <- ggplot(phy, aes(x = identity, y = frac_s_lt, colour = label)) +
  geom_vline(xintercept = 97, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_line(linewidth = 0.7) +
  ggtitle("Subject prevalence") +
  scale_x_continuous(breaks = XBR, limits = c(70, 100),
                     name = "Best identity to HOMD reference (%)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1),
                     name = "Fraction of subject prevalence (identity < x)") +
  scale_colour_manual(values = PAL, name = "Phylum") +
  base_theme + theme(legend.position = "right")

# ---- p4: total (all phyla combined, 3 metrics) ----
p4 <- ggplot(tot_long, aes(x = identity, y = frac, colour = metric)) +
  geom_vline(xintercept = 97, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_line(linewidth = 0.9) +
  ggtitle("All phyla combined") +
  scale_x_continuous(breaks = XBR, limits = c(70, 100),
                     name = "Best identity to HOMD reference (%)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1),
                     name = "Fraction (identity < x)") +
  scale_colour_manual(values = c("ASV count" = "grey40",
                                 "Reads" = "#1f77b4",
                                 "Subject prevalence" = "#d62728"),
                      name = "Metric") +
  base_theme + theme(legend.position = "right")

p_row <- aplot::plot_list(p1, p2, p3, p4, ncol = 4, guides = "collect")

fig <- file.path(FIGOUT, "fig_v3_homd_identity_lt_fraction_row")
ggsave(paste0(fig, ".png"), p_row, width = 14.5, height = 4.6, dpi = 300)
ggsave(paste0(fig, ".pdf"), p_row, width = 14.5, height = 4.6)
message("[13s] 4-panel row written")
