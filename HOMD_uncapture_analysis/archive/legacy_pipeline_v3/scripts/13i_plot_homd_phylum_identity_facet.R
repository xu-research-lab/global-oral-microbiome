#!/usr/bin/env Rscript
# pipeline_v3: combined phylum-identity figure.
#   Left block (one faceted ggplot, 2 rows x 2 columns, shared axes):
#     All phyla combined | ASV count
#     Reads | Subject prevalence
#   Axes are swapped: x = fraction of the metric with identity >= x,
#                     y = best HOMD identity (%).
#   Cumulative direction: fraction of the metric with identity >= x.
#   Right column: 100% stacked bar of phyla among ASVs with
#     70 <= identity < 97 (best-hit species phylum), sorted by that share.
# Run 13u_homd_lt97_ge70_phylum.py first to create the subset count table.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(aplot)
  library(ggnewscale)
  library(ggh4x)     # strip_themed(): one strip background per facet
})

# ---- palette ---------------------------------------------------------------
# 15 phyla need 15 distinguishable colours, but only the handful that dominate
# the two stacked bars carries the message.  The eight phyla that are visible in
# the bars therefore get saturated, well-separated hues; the long tail of rare
# phyla gets deliberately muted tones, which also pushes their (already
# bunched) curves into the background.
PHYLUM_COL <- c(
  Fusobacteriota    = "#8E2F52",
  Bacteroidota      = "#E07B39",
  Spirochaetota     = "#3AA6A0",
  Actinomycetota    = "#7B5EA7",
  Patescibacteria   = "#C9A227",
  Bacillota         = "#2E6DA4",
  Pseudomonadota    = "#3F9C6D",
  Synergistota      = "#8C613C",
  Mycoplasmatota    = "#A47FB5",
  Verrucomicrobiota = "#7FB3D5",
  Methanobacteriota = "#D9A66B",
  Chloroflexota     = "#9BAE7A",
  Chlorobiota       = "#C97F7F",
  Chlamydiota       = "#6E7B85",
  Cyanobacteriota   = "#5BC8C8")

# the three metrics, used both for the "All phyla" curves and for the
# facet strip that carries each metric
MET_COL   <- c("ASV count" = "#4D4D4D", "Reads" = "#0072B2",
               "Prevalence" = "#D55E00")
MET_STRIP <- c("ASV count" = "#D2D2D2", "Reads" = "#D3E6F5",
               "Prevalence" = "#FAE0CE")
MET_TEXT  <- c("ASV count" = "#3A3A3A", "Reads" = "#00527F",
               "Prevalence" = "#A84400")

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

cov <- read.delim(file.path(PIPE_V3, "13d_homd_phylum_coverage.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)

phy_frac_file <- file.path(PIPE_V3, "13r_homd_phylum_lt_fractions.tsv")
tot_frac_file <- file.path(PIPE_V3, "13r_homd_total_lt_fractions.tsv")
sub_file      <- file.path(PIPE_V3, "13u_homd_lt97_ge70_phylum.tsv")
gt_file       <- file.path(PIPE_V3, "13u_homd_gt97_phylum.tsv")
stopifnot(file.exists(phy_frac_file), file.exists(tot_frac_file),
          file.exists(sub_file), file.exists(gt_file))

phy_frac <- read.delim(phy_frac_file, sep = "\t", header = TRUE,
                       check.names = FALSE, stringsAsFactors = FALSE)
tot_frac <- read.delim(tot_frac_file, sep = "\t", header = TRUE,
                       check.names = FALSE, stringsAsFactors = FALSE)
sub_cnt  <- read.delim(sub_file, sep = "\t", header = TRUE,
                       check.names = FALSE, stringsAsFactors = FALSE)
gt_cnt   <- read.delim(gt_file, sep = "\t", header = TRUE,
                       check.names = FALSE, stringsAsFactors = FALSE)

# shared phylum metadata; legend ordered by the first percentage (frac_ge97)
meta2 <- cov %>%
  left_join(sub_cnt, by = "phylum") %>%
  left_join(gt_cnt, by = "phylum") %>%
  transmute(phylum,
            n_total   = n_asvs_total,
            n_ge97    = n_asvs_ge97,
            n_70_97   = n_asvs_70_97,
            n_gt97    = n_asvs_gt97,
            frac_ge97 = n_asvs_ge97 / n_asvs_total,
            frac_70_97 = n_asvs_70_97 / sum(n_asvs_70_97)) %>%
  arrange(desc(frac_ge97))
meta2$label <- sprintf("%s (%.1f%%)", meta2$phylum, 100 * meta2$frac_ge97)
meta2$label <- factor(meta2$label, levels = meta2$label)
stopifnot(all(meta2$phylum %in% names(PHYLUM_COL)))
PAL <- setNames(unname(PHYLUM_COL[meta2$phylum]), as.character(meta2$label))
legend_levels <- levels(meta2$label)

# ---------------- faceted 2x2 block data ----------------
phy <- phy_frac %>% filter(identity >= 70 & identity <= 100) %>%
  left_join(meta2, by = "phylum")

phy_long <- bind_rows(
  phy %>% transmute(panel = "ASV count", identity, frac = 1 - frac_n_lt, label),
  phy %>% transmute(panel = "Reads", identity, frac = 1 - frac_r_lt, label),
  phy %>% transmute(panel = "Prevalence", identity, frac = 1 - frac_s_lt,
                    label)
)

tot_long <- bind_rows(
  tot_frac %>% filter(identity >= 70 & identity <= 100) %>%
    transmute(panel = "All phyla", metric = "ASV count",
              identity, frac = 1 - frac_n_lt),
  tot_frac %>% filter(identity >= 70 & identity <= 100) %>%
    transmute(panel = "All phyla", metric = "Reads",
              identity, frac = 1 - frac_r_lt),
  tot_frac %>% filter(identity >= 70 & identity <= 100) %>%
    transmute(panel = "All phyla", metric = "Prevalence",
              identity, frac = 1 - frac_s_lt)
)
tot_long$metric <- factor(tot_long$metric,
                          levels = c("ASV count", "Reads", "Prevalence"))
tot97 <- tot_frac %>% filter(abs(identity - 97) < 1e-9)
metric_ge97 <- c("ASV count" = 1 - tot97$frac_n_lt[1],
                 "Reads" = 1 - tot97$frac_r_lt[1],
                 "Prevalence" = 1 - tot97$frac_s_lt[1])
metric_labs <- sprintf("%s (%.1f%%)", names(metric_ge97),
                       100 * metric_ge97)
panel_levels <- c("All phyla", "ASV count", "Reads", "Prevalence")

# invisible geom_col in phylum curves carries the fill guide that is shared
# with the stacked bar (same geom/colour/linewidth so aplot merges them)
pfacet <- ggplot() +
  geom_hline(yintercept = 97, linetype = "dashed", colour = "grey40",
             linewidth = 0.35) +
  geom_line(data = phy_long, aes(x = frac, y = identity, colour = label),
            linewidth = 0.7) +
  scale_colour_manual(values = PAL, guide = "none") +
  geom_col(data = phy_long, aes(x = frac, y = 70, fill = label),
           width = 0, alpha = 0) +
  scale_fill_manual(values = PAL, breaks = legend_levels, name = NULL,
                    drop = FALSE,
                    guide = guide_legend(order = 2, direction = "vertical",
                                         override.aes = list(alpha = 1))) +
  new_scale_colour() +
  geom_line(data = tot_long, aes(x = frac, y = identity, colour = metric),
            linewidth = 0.9) +
  scale_colour_manual(values = MET_COL,
                      labels = metric_labs,
                      name = "Metric",
                      guide = guide_legend(order = 1, direction = "vertical")) +
  # each metric panel is headed by a strip tinted with that metric's colour, so
  # the panel and the combined curve it contains are read as the same thing
  ggh4x::facet_wrap2(~ factor(panel, levels = panel_levels), ncol = 2,
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = c("grey95", unname(MET_STRIP[panel_levels[-1]])),
        colour = NA),
      text_x = ggh4x::elem_list_text(
        size = 16, face = "bold",
        colour = c("grey20", unname(MET_TEXT[panel_levels[-1]]))))) +
  scale_x_continuous(breaks = seq(0, 1, 0.5),
                     labels = seq(0, 100, 50), limits = c(0, 1),
                     name = "Proportion (%)") +
  scale_y_continuous(breaks = seq(70, 100, 10),
                     name = "Best identity to eHOMD reference (%)") +
  coord_cartesian(ylim = c(70, 100)) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

# ---------------- two 100% stacked bars ----------------
bar70 <- meta2 %>% arrange(desc(n_70_97)) %>%
  mutate(x = 1, label2 = factor(as.character(label), levels = as.character(label)),
         frac = n_70_97 / sum(n_70_97))
bargt <- meta2 %>% arrange(desc(n_gt97)) %>%
  mutate(x = 2, label2 = factor(as.character(label), levels = as.character(label)),
         frac = n_gt97 / sum(n_gt97))

annotate_top5 <- function(d) {
  d %>% arrange(desc(frac)) %>%
    # ggplot stacks from the LAST factor level upwards, i.e. smallest share
    # at the bottom; compute segment midpoints accordingly
    mutate(cum = cumsum(frac),
           ymid = 1 - cum + frac / 2) %>%
    filter(frac >= 0.08) %>%
    mutate(pct = sprintf("%.1f%%", 100 * frac))
}
top5_lab <- bind_rows(annotate_top5(bar70), annotate_top5(bargt))

p5 <- ggplot() +
  geom_col(data = bar70, aes(x = x, y = frac, fill = label2), width = 0.85) +
  geom_col(data = bargt, aes(x = x, y = frac, fill = label2), width = 0.85) +
  geom_text(data = top5_lab, aes(x = x, y = ymid, label = pct),
            colour = "white", size = 3.3, fontface = "bold") +
  scale_fill_manual(values = PAL, breaks = legend_levels, name = NULL,
                    drop = FALSE,
                    guide = "none") +
  scale_x_continuous(breaks = 1:2,
                     labels = c("≥70% and\n<97%", "≥97%"),
                     limits = c(0.5, 2.5),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1),
                     expand = expansion(mult = c(0, 0)),
                     name = "Proportion") +
  theme_bw(base_size = 17) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = 12),
        axis.title.x = element_blank())

# ---------------- composite ----------------
combo <- aplot::plot_list(pfacet, p5, ncol = 2,
                          widths = c(2.7, 0.8),
                          guides = "collect") &
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(0.48, "cm"),
        legend.text = element_text(size = 11.8))
combo <- combo +
  patchwork::plot_annotation(theme = theme(
    plot.margin = margin(2, 12, 2, 12)
  ))

fig2 <- file.path(FIGOUT, "fig_v3_homd_phylum_identity_2x2_stack")
ggsave(paste0(fig2, ".png"), combo, width = 9, height = 5.5, dpi = 300)
ggsave(paste0(fig2, ".pdf"), combo, width = 9, height = 5.5, device = cairo_pdf)
message("[13i] faceted 2x2 + 70-97% stacked bar written")
