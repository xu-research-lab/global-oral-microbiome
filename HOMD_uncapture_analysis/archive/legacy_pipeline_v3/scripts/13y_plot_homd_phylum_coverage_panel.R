#!/usr/bin/env Rscript
# pipeline_v3: panel B, second take -- eHOMD coverage of the oral ASV catalogue.
#
# The first take (13i, fig_v3_homd_phylum_identity_2x2_stack) drew every phylum as
# one of fifteen coloured curves on a swapped axis pair, with the three weightings
# as a second family of curves in a shared 2x2 grid.  Readers had to be told how to
# read it, and fifteen hues is past the point where any of them can be told apart,
# so the one comparison that matters -- the same phylum under the three weightings --
# was exactly the one the encoding could not support.
#
# This take splits the panel by the question each part answers, and lets position
# carry the phyla so that colour is free to carry something with few levels:
#
#   p1  how much of the catalogue reaches the boundary, by weighting  (3 bars)
#   p2  the same quantity per phylum, ranked                          (15 x 3 bars)
#   p3  which phyla make up the two identity bands                   (5 x 2 bars)
#
# Phyla are never encoded by colour here; they are rows.  Colour carries the
# weighting (p1, p2) and the identity band (p3), so every colour scale has at most
# three levels.
#
# Run 13r_homd_lt_fractions.py and 13u_homd_lt97_ge70_phylum.py first.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")

# ---- palette ---------------------------------------------------------------
# The weightings keep the grey / blue / orange already used for them elsewhere in
# Figure 1, so the existing figure legend stays true; the two identity bands keep
# the red / green pair the rest of the figure already uses for below / at the
# species boundary.  Both trios and the pair were checked for colour-vision
# separation (worst pair CVD dE 16.1 and 9.6, normal-vision 17.2 and 32.1).
MET_COL <- c("ASVs" = "#4D4D4D", "Reads" = "#0072B2", "Subjects" = "#D55E00")
BAND_COL <- c("70–97% identity" = "#B2182B", "≥97% identity" = "#1A9850")
COL_TRACK <- "#EDEDED"          # the unfilled 100% track behind each bar
COL_TXT   <- "#1A1A1A"
COL_MUTED <- "#6E6E6E"
COL_GRID  <- "#E4E4E4"

# ---- data ------------------------------------------------------------------
cov <- read.delim(file.path(PIPE_V3, "13d_homd_phylum_coverage.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)
phy_frac <- read.delim(file.path(PIPE_V3, "13r_homd_phylum_lt_fractions.tsv"),
                       sep = "\t", header = TRUE, check.names = FALSE,
                       stringsAsFactors = FALSE)
tot_frac <- read.delim(file.path(PIPE_V3, "13r_homd_total_lt_fractions.tsv"),
                       sep = "\t", header = TRUE, check.names = FALSE,
                       stringsAsFactors = FALSE)
band70 <- read.delim(file.path(PIPE_V3, "13u_homd_lt97_ge70_phylum.tsv"),
                     sep = "\t", header = TRUE, stringsAsFactors = FALSE)
bandgt <- read.delim(file.path(PIPE_V3, "13u_homd_gt97_phylum.tsv"),
                     sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# every fraction is read off the identity grid at the species boundary
at97 <- function(df) df[abs(df$identity - 97) < 1e-9, , drop = FALSE]

tot97 <- at97(tot_frac)
METRICS <- data.frame(
  metric = factor(names(MET_COL), levels = names(MET_COL)),
  ge97   = c(1 - tot97$frac_n_lt, 1 - tot97$frac_r_lt, 1 - tot97$frac_s_lt))
METRICS$lab <- sprintf("%.1f%%", 100 * METRICS$ge97)

p97 <- at97(phy_frac) %>%
  transmute(phylum,
            ASVs     = 1 - frac_n_lt,
            Reads    = 1 - frac_r_lt,
            Subjects = 1 - frac_s_lt) %>%
  left_join(cov[, c("phylum", "n_asvs_total")], by = "phylum")
stopifnot(nrow(p97) == 15, !anyNA(p97))

# phylum order: best-covered first, by the primary (ASV) weighting
p97$phylum <- factor(p97$phylum, levels = p97$phylum[order(p97$ASVs)])

PHY_LONG <- bind_rows(lapply(names(MET_COL), function(m) {
  data.frame(phylum = p97$phylum, metric = m, ge97 = p97[[m]],
             stringsAsFactors = FALSE)
}))
PHY_LONG$metric <- factor(PHY_LONG$metric, levels = names(MET_COL))

# the two identity bands, as a share of each band's ASVs; the four phyla that
# carry the story are named and the rest are pooled, as the bars cannot hold
# fifteen labels
NAMED <- c("Bacillota", "Pseudomonadota", "Bacteroidota", "Actinomycetota")
band_share <- function(cnt, col) {
  cnt <- cnt[order(-cnt[[col]]), ]
  top <- cnt[cnt$phylum %in% NAMED, ]
  top <- top[order(match(top$phylum, NAMED)), ]
  rest <- sum(cnt[[col]]) - sum(top[[col]])
  data.frame(phylum = c(top$phylum, "Other phyla"),
             n = c(top[[col]], rest))
}
b70 <- data.frame(band = names(BAND_COL)[1],
                  band_share(band70, "n_asvs_70_97"), stringsAsFactors = FALSE)
bgt <- data.frame(band = names(BAND_COL)[2],
                  band_share(bandgt, "n_asvs_gt97"), stringsAsFactors = FALSE)
BANDS <- rbind(b70, bgt)
BANDS$share  <- BANDS$n / ave(BANDS$n, BANDS$band, FUN = sum)
BANDS$lab    <- sprintf("%.1f%%", 100 * BANDS$share)
BANDS$phylum <- factor(BANDS$phylum, levels = c(NAMED, "Other phyla"))
BANDS$band   <- factor(BANDS$band, levels = names(BAND_COL))
BANDS$solid  <- BANDS$phylum != "Other phyla"   # pooled category drawn pale

# ---- p1 | the catalogue as a whole ------------------------------------------
# one bar per weighting: how much of that weighting reaches the species boundary.
p1 <- ggplot(METRICS) +
  geom_col(aes(x = 100, y = metric), fill = COL_TRACK, width = 0.62) +
  geom_col(aes(x = 100 * ge97, y = metric, fill = metric), width = 0.62) +
  # the value is set against the end of the fill rather than its middle, so the
  # number marks where the bar stops
  geom_text(aes(x = 100 * ge97 - 1.5, y = metric, label = lab), hjust = 1,
            colour = "white", fontface = "bold", size = 5.2) +
  scale_fill_manual(values = MET_COL, guide = "none") +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "The catalogue as a whole",
       subtitle = "share of each weighting at or above the 97% species boundary") +
  theme_minimal(base_size = 18.3) +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 18.3, face = "bold", colour = COL_TXT,
                                  margin = margin(b = 1)),
        plot.subtitle = element_text(size = 12.4, colour = COL_MUTED,
                                     margin = margin(b = 6)))

# ---- p2 | by phylum --------------------------------------------------------
# Three weightings per phylum, dodged.  Only the primary weighting is numbered:
# the row is far too shallow to carry three labels, and one number per row in a
# fixed column right of the axis reads as a ranking table next to the bars.
# p2 and p3 share this range so that a bar of a given percentage is the same
# length in both; the 18 units past 100 are the value gutter.
GUTTER <- 103                     # where the value column sits, in x units
XMAX   <- 118

PHY_LAB <- PHY_LONG %>% filter(metric == "ASVs") %>%
  mutate(lab = sprintf("%.1f", 100 * ge97))

p2 <- ggplot(PHY_LONG, aes(y = phylum, x = 100 * ge97, fill = metric)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.24) +
  geom_text(data = PHY_LAB, aes(x = GUTTER, y = phylum, label = lab),
            hjust = 0, size = 4.0, colour = COL_TXT) +
  scale_fill_manual(values = MET_COL, name = NULL) +
  scale_x_continuous(limits = c(0, XMAX), breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "The same, phylum by phylum",
       subtitle = "the number is the share of that phylum's ASVs",
       x = "share at or above the 97% species boundary (%)") +
  theme_minimal(base_size = 15) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = COL_GRID, linewidth = 0.3),
        axis.text.y = element_text(colour = COL_TXT),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11.7, colour = COL_MUTED,
                                    margin = margin(t = 4)),
        plot.title = element_text(size = 15.6, face = "bold", colour = COL_TXT,
                                  margin = margin(b = 1)),
        plot.subtitle = element_text(size = 12.4, colour = COL_MUTED,
                                     margin = margin(b = 6)),
        legend.position = "bottom",
        legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 12.4),
        legend.margin = margin(t = -4))

# ---- p3 | where the gap sits -----------------------------------------------
# the same four phyla plus the pooled remainder, as a share of each band
p3 <- ggplot(BANDS, aes(y = phylum, x = 100 * share, fill = band,
                        alpha = solid)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.55) +
  scale_fill_manual(values = BAND_COL, name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5), guide = "none") +
  scale_x_continuous(limits = c(0, XMAX), breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0))) +
  # the two bars of a row are dodged along y, so each label can sit at its own
  # bar's end without meeting the other
  geom_text(aes(x = 100 * share, y = phylum, label = lab, group = band),
            position = position_dodge(width = 0.72), hjust = -0.15,
            size = 3.1, colour = COL_TXT) +
  labs(title = "Which phyla fill each band",
       subtitle = "share of the ASVs in that band",
       x = "share of the band's ASVs (%)") +
  theme_minimal(base_size = 13.7) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = COL_GRID, linewidth = 0.3),
        axis.text.y = element_text(colour = COL_TXT),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9, colour = COL_MUTED,
                                    margin = margin(t = 4)),
        plot.title = element_text(size = 15, face = "bold", colour = COL_TXT,
                                  margin = margin(b = 1)),
        plot.subtitle = element_text(size = 12.4, colour = COL_MUTED,
                                     margin = margin(b = 6)),
        legend.position = "bottom",
        legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 11.5),
        legend.margin = margin(t = -4))

# ---- composite -------------------------------------------------------------
panelB <- p1 / (p2 | p3) +
  plot_layout(heights = c(1, 3.5)) +
  plot_annotation(theme = theme(plot.margin = margin(2, 6, 2, 2)))

out <- file.path(FIGOUT, "fig_v3_panelB_coverage_v2")
ggsave(paste0(out, ".png"), panelB, width = 9, height = 5.5, dpi = 300,
       bg = "white")
ggsave(paste0(out, ".pdf"), panelB, width = 9, height = 5.5,
       device = cairo_pdf, bg = "white")
message("[13y] panel B v2 written: ", out, ".png | .pdf")
