#!/usr/bin/env Rscript
# pipeline_v3: HOMD identity ratio curves.
#
# Fig 1 (total, placed first): all phyla combined, three ratio curves
#   y = log10( sum(metric of ASVs with identity >= x) / sum(metric of ASVs with identity < x) )
#   metrics: ASV count, total reads, subject prevalence.
#
# Fig 2 (phylum): one row via aplot, shared phylum legend
#   p1: original cumulative identity curve (2 facets, as before)
#   p2: reads-count ratio curve (70-100%)
#   p3: subject-prevalence ratio curve (70-100%)
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

tot <- read.delim(file.path(PIPE_V3, "13p_homd_total_ratio_curves.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
phy <- read.delim(file.path(PIPE_V3, "13p_homd_phylum_ratio_curves.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
d   <- read.delim(file.path(PIPE_V3, "13f_homd_phylum_ecdf.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cov <- read.delim(file.path(PIPE_V3, "13d_homd_phylum_coverage.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# ---------------- phylum metadata (legend order: frac_ge97 high -> low) ----------------
meta <- cov %>% transmute(
  phylum    = phylum,
  n_total   = n_asvs_total,
  n_ge97    = n_asvs_ge97,
  frac_ge97 = n_asvs_ge97 / n_asvs_total,
  low_ge97  = n_asvs_ge97 < 500
) %>% arrange(desc(frac_ge97))
meta$label <- sprintf("%s (%.1f%%)", meta$phylum, 100 * meta$frac_ge97)
meta$label <- factor(meta$label, levels = meta$label)
PAL <- hue_pal()(15)

# ================= Fig 1: total ratio curves =================
tot_long <- tot %>%
  filter(identity >= 70 & identity <= 100) %>%
  pivot_longer(cols = c(ratio_n_asvs, ratio_reads, ratio_prevalence),
               names_to = "metric", values_to = "ratio") %>%
  mutate(metric = recode(metric,
                         ratio_n_asvs    = "ASV count",
                         ratio_reads     = "Reads",
                         ratio_prevalence = "Subject prevalence"),
         metric = factor(metric, levels = c("ASV count", "Reads", "Subject prevalence"))) %>%
  filter(is.finite(log10(ratio)))

p_tot <- ggplot(tot_long, aes(x = identity, y = log10(ratio), colour = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(breaks = seq(70, 100, 5), limits = c(70, 100),
                     name = "Best identity to HOMD reference (%)") +
  scale_y_continuous(name = "log10 ratio (sum >= x / sum < x)") +
  scale_colour_manual(values = c("ASV count" = "grey40",
                                 "Reads" = "#1f77b4",
                                 "Subject prevalence" = "#d62728"),
                      name = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9), panel.grid.minor = element_blank())

fig1 <- file.path(FIGOUT, "fig_v3_homd_total_identity_ratio")
ggsave(paste0(fig1, ".png"), p_tot, width = 6.5, height = 4.2, dpi = 300)
ggsave(paste0(fig1, ".pdf"), p_tot, width = 6.5, height = 4.2)
message("[13q] total ratio figure written")

# ================= Fig 2: phylum row (aplot) =================
# ---- p1: original cumulative curve (2 facets) ----
d1 <- d %>% filter(identity >= 70 & identity <= 97) %>%
  left_join(meta, by = "phylum") %>%
  mutate(panel = "All ASVs (70-97%)", frac = frac_ge)
d2 <- d %>% filter(identity >= 97) %>%
  group_by(phylum) %>% mutate(frac2 = frac_ge / max(frac_ge)) %>% ungroup() %>%
  left_join(meta, by = "phylum") %>%
  mutate(panel = "ASVs >=97% (97-100%)", frac = frac2)
dd <- bind_rows(d1, d2)
dd$panel <- factor(dd$panel, levels = c("All ASVs (70-97%)", "ASVs >=97% (97-100%)"))
dd$lty <- ifelse(dd$panel == "ASVs >=97% (97-100%)" & dd$low_ge97, "dashed", "solid")

p1 <- ggplot(dd, aes(x = identity, y = frac, colour = label)) +
  geom_line(aes(linetype = lty), linewidth = 0.8) +
  scale_linetype_identity(guide = "none") +
  geom_vline(xintercept = 97, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  scale_x_continuous(name = "Best identity to HOMD reference (%)",
                     breaks = c(70, 75, 80, 85, 90, 95, 99, 100)) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.08),
                     name = "Cumulative fraction of ASVs") +
  scale_colour_manual(values = PAL, name = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 7.5), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 8.5))

# ---- p2 / p3: ratio curves ----
pr <- phy %>% filter(identity >= 70 & identity <= 100) %>%
  left_join(meta, by = "phylum") %>%
  filter(is.finite(log10(ratio_reads)), is.finite(log10(ratio_prevalence)))

p2 <- ggplot(pr, aes(x = identity, y = log10(ratio_reads), colour = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  ggtitle("Reads ratio") +
  scale_x_continuous(breaks = seq(70, 100, 5), limits = c(70, 100),
                     name = "Best identity to HOMD reference (%)") +
  scale_y_continuous(name = "log10 ratio (>=x / <x)") +
  scale_colour_manual(values = PAL, name = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 7.5), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5))

p3 <- ggplot(pr, aes(x = identity, y = log10(ratio_prevalence), colour = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  ggtitle("Subject prevalence ratio") +
  scale_x_continuous(breaks = seq(70, 100, 5), limits = c(70, 100),
                     name = "Best identity to HOMD reference (%)") +
  scale_y_continuous(name = "log10 ratio (>=x / <x)") +
  scale_colour_manual(values = PAL, name = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 7.5), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5))

p_row <- aplot::plot_list(p1, p2, p3, ncol = 3, widths = c(2, 1, 1), guides = "collect")

fig2 <- file.path(FIGOUT, "fig_v3_homd_phylum_identity_ratio_row")
ggsave(paste0(fig2, ".png"), p_row, width = 12.5, height = 4.6, dpi = 300)
ggsave(paste0(fig2, ".pdf"), p_row, width = 12.5, height = 4.6)
message("[13q] phylum ratio row written")
