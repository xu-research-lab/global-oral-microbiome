#!/usr/bin/env Rscript
# pipeline_v3: per-species eHOMD coverage as a ranked dot-line plot.
#
# One row per eHOMD species (n = 801), ranked by the ASV-weighted share of that
# species' own ASVs whose best hit reaches 97% identity.  Three weightings are
# drawn on each row -- ASVs (grey), reads (blue) and subject-level prevalence
# (orange), the same three the phylum panel uses -- so a row read across shows
# whether a species that is well covered by ASV count is also well covered by
# the reads and subjects that carry it.
#
# The species names cannot be printed: 801 rows do not leave room for them.  The
# row identity is carried by three heatmap tracks instead, assembled with aplot
# along the y axis -- the phylum the species belongs to (in the colours of the
# composition panel), its number of eHOMD reference sequences, and the number of
# ASVs mapped to it.  A handful of species still carry their names, as text
# beside their point: the extremes of the four richness/coverage quadrants of the
# previous version of this panel, so the same organisms stay identifiable.
#
# Run 13m_homd_species_metric_coverage.py first.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(aplot)
  library(ggrepel)
  library(scales)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

# the phylum colours are the composition panel's own, so a colour in the track
# and a segment in that panel are the same taxon
source(file.path(PIPE_V3, "scripts", "13z2_panel_homd_band_composition.R"))
NAMED <- c("Bacillota", "Bacteroidota", "Pseudomonadota",
           "Actinomycetota", "Fusobacteriota")

MET_COL <- c("ASVs" = "#4D4D4D", "Reads" = "#0072B2", "Prevalence" = "#D55E00")
COL_TXT <- "#1A1A1A"
LOW     <- "#F5F5F5"
HI_REF  <- "#1B7F79"          # eHOMD reference sequences
HI_ASV  <- "#3B6EA5"          # ASVs mapped

# ---- data ------------------------------------------------------------------
s <- read.delim(file.path(PIPE_V3, "13m_homd_species_metric_coverage.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE,
                stringsAsFactors = FALSE)
stopifnot(nrow(s) == 801, !anyNA(s$phylum))

s$phylum_grp <- factor(ifelse(s$phylum %in% NAMED, s$phylum, "Other phyla"),
                       levels = names(PHY_COL))
# rank by the ASV-weighted share, best covered at the top.  ggplot draws factor
# level 1 at the bottom, so the levels run from worst to best.
ORD <- s$species[order(s$frac_asv_ge97)]
s$sp_f <- factor(s$species, levels = ORD)
s$row  <- as.integer(s$sp_f)

# every row must carry a point for each weighting, or geom_path would join a
# metric across a gap that is not there
LONG <- bind_rows(lapply(names(MET_COL), function(m) {
  data.frame(sp_f  = s$sp_f,
             row   = s$row,
             metric = m,
             prop  = 100 * s[[c(ASVs = "frac_asv_ge97", Reads = "frac_reads_ge97",
                                Prevalence = "frac_prev_ge97")[[m]]]],
             stringsAsFactors = FALSE)
}))
LONG$metric <- factor(LONG$metric, levels = names(MET_COL))
# geom_line re-sorts each group by its x value before drawing; the ranking is the
# row order, so the path is drawn with geom_path from data already in row order
LONG <- LONG[order(LONG$metric, LONG$row), ]

# ---- raw values and the trend along the ranking -----------------------------
# Reads and prevalence jump from row to row -- neighbouring species can sit at
# 100% and at 40% in turn -- so joining the rows directly drew a hatch of
# near-horizontal strokes rather than anything that could be followed.  The two
# are therefore shown twice over: every species keeps its own value as a faint
# point, and a running median over RUNMED neighbouring ranks gives the trend
# ("species of comparable ASV coverage"), which smooths the noise without
# inventing a value for any individual species.  The ASV curve is the ranking
# itself, monotone by construction, and is drawn exactly.
RUNMED <- 41
# a running median left visible steps at every window boundary; a 20% trimmed
# running mean keeps the robustness and gives a line that can actually be read
roll_trim <- function(x, k = RUNMED, trim = 0.2) {
  h <- k %/% 2
  n <- length(x)
  vapply(seq_len(n), function(i) {
    mean(x[max(1, i - h):min(n, i + h)], trim = trim)
  }, numeric(1))
}
LONG$prop_s <- NA_real_
for (m in c("Reads", "Prevalence")) {
  i <- LONG$metric == m
  LONG$prop_s[i] <- roll_trim(LONG$prop[i])
}
# the ASV path uses its exact value; the other two use the trend
LINE <- LONG
LINE$xval <- ifelse(is.na(LINE$prop_s), LINE$prop, LINE$prop_s)
ASV_D <- LINE[LINE$metric == "ASVs", ]
TREND <- LINE[LINE$metric != "ASVs", ]
OTHER <- LONG[LONG$metric != "ASVs", ]

# ---- the named species ------------------------------------------------------
# The labelled set is the one the previous quadrant plot singled out: the three
# most ASV-rich species in each high-richness quadrant, the three least ASV-rich
# in each low-richness one, and the two best-referenced low-richness/high-gap
# species.
q <- read.delim(file.path(PIPE_V3, "13j_homd_species_coverage_fullname.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE,
                stringsAsFactors = FALSE)
q$frac_lt97 <- 1 - pmin(pmax(q$frac_ge97, 0), 1)
q$quad <- with(q, ifelse(frac_lt97 >= 0.5 & n_asvs_total >= 5000, "HH",
                  ifelse(frac_lt97 <  0.5 & n_asvs_total >= 5000, "HL",
                  ifelse(frac_lt97 >= 0.5, "LH", "LL"))))
lab <- bind_rows(
  q %>% filter(quad == "HH") %>% arrange(desc(n_asvs_total)) %>% slice_head(n = 3),
  q %>% filter(quad == "HL") %>% arrange(desc(n_asvs_total)) %>% slice_head(n = 3),
  q %>% filter(quad == "LH") %>% arrange(n_asvs_total)        %>% slice_head(n = 3),
  q %>% filter(quad == "LL") %>% arrange(n_asvs_total)        %>% slice_head(n = 3),
  q %>% filter(quad == "LH") %>% arrange(desc(n_homd_refseq)) %>% slice_head(n = 2)
) %>% distinct(species, .keep_all = TRUE)
lab <- lab[lab$species %in% s$species, ]
lab <- merge(lab[, c("species", "frac_ge97")], s[, c("species", "sp_f", "row")],
             by = "species")
lab$x <- 100 * lab$frac_ge97

# ---- the ranked plot --------------------------------------------------------
p_main <- ggplot() +
  geom_vline(xintercept = c(0, 25, 50, 75, 100), colour = "#E4E4E4",
             linewidth = 0.3) +
  # every species' own reads and prevalence, faint
  geom_point(data = OTHER, aes(x = prop, y = row, colour = metric),
             size = 0.20, alpha = 0.30) +
  # the trend of each weighting along the ranking
  geom_path(data = TREND, aes(x = xval, y = row, colour = metric),
            linewidth = 0.55, alpha = 0.90) +
  # the ASV curve is the ranking itself and is drawn exactly and heaviest
  geom_path(data = ASV_D, aes(x = prop, y = row, colour = metric),
            linewidth = 0.60, alpha = 0.95) +
  geom_point(data = ASV_D, aes(x = prop, y = row, colour = metric),
             size = 0.45, alpha = 0.95) +
  geom_text_repel(data = lab, aes(x = x, y = row, label = species),
                  inherit.aes = FALSE, size = 2.45, colour = COL_TXT,
                  fontface = "italic", max.overlaps = Inf, seed = 42,
                  min.segment.length = 0, box.padding = 0.55, force = 4,
                  max.iter = 1e5, segment.size = 0.25, segment.colour = "grey35",
                  bg.colour = "white", bg.r = 0.16) +
  scale_colour_manual(values = MET_COL, name = NULL, limits = names(MET_COL),
                      guide = guide_legend(
                        override.aes = list(size = 1.5, linewidth = 0.6,
                                            alpha = 1, shape = 16))) +
  scale_x_continuous(breaks = seq(0, 100, 25), limits = c(0, 100),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.004, 0.004))) +
  labs(x = "Proportion of each weighting at ≥97% identity (%)", y = NULL,
       title = "eHOMD species ranked by ASV-weighted coverage",
       subtitle = "801 species, best covered at the top") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        plot.title    = element_text(size = 9.5, face = "bold", colour = COL_TXT),
        plot.subtitle = element_text(size = 8, colour = "grey35",
                                     margin = margin(b = 3)),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(colour = COL_TXT),
        axis.title.x = element_text(colour = COL_TXT),
        # 801 species: there is no room for names, so neither names nor ticks
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(2, 4, 2, 2))

# ---- the tracks -------------------------------------------------------------
# Each track is a one-column tile strip on the same y scale; aplot aligns their
# panels with the ranked plot's, so a row in a track and a row in the plot are
# the same species.
track_theme <- function() {
  theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          axis.text  = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 5.6, hjust = 0.5, colour = COL_TXT,
                                    margin = margin(b = 1)),
          plot.margin = margin(2, 1, 2, 1),
          legend.position = "none")
}
p_phy <- ggplot(s, aes(x = 1, y = row, fill = phylum_grp)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = PHY_COL, drop = FALSE, name = "phylum") +
  scale_y_continuous(expand = expansion(mult = c(0.004, 0.004))) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "phylum") +
  # the phylum legend is the only key this panel carries, and it is collected
  # with the two colourbars at the right
  theme(legend.position = "right") +
  track_theme() + theme(legend.position = "right")

# a sequential track: log10 of the count, with the colourbar labelling the count
# itself rather than its logarithm
seq_track <- function(v, name, high, breaks) {
  ggplot(data.frame(row = s$row, value = pmax(v, 1)), aes(x = 1, y = row, fill = value)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_gradient(trans = "log10", low = LOW, high = high, name = name,
                        breaks = breaks,
                        labels = label_number(scale_cut = cut_short_scale()),
                        guide = guide_colourbar(barwidth = unit(0.40, "cm"),
                                                barheight = unit(2.0, "cm"),
                                                title.position = "top")) +
    scale_y_continuous(expand = expansion(mult = c(0.004, 0.004))) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = name) +
    track_theme()
}
p_refs <- seq_track(s$n_homd_refseq, "refs", HI_REF, c(1, 10, 100))
p_asvs <- seq_track(s$n_asvs_total,  "ASVs", HI_ASV, c(1, 10, 100, 1000, 10000))

# ---- assemble ---------------------------------------------------------------
# Tracks on the left, ranked plot on the right; widths in inches.  The legends
# are collected to the right of the plot.
W_FIG <- 7.4
H_FIG <- 11.0
p <- aplot::plot_list(p_phy, p_refs, p_asvs, p_main, ncol = 4,
                      widths = c(0.32, 0.44, 0.44, 4.6),
                      guides = "collect") &
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm"))

out <- file.path(FIGOUT, "fig_v3_homd_species_metric_dotline")
ggsave(paste0(out, ".png"), p, width = W_FIG, height = H_FIG, dpi = 300,
       bg = "white", limitsize = FALSE)
ggsave(paste0(out, ".pdf"), p, width = W_FIG, height = H_FIG,
       device = cairo_pdf, bg = "white", limitsize = FALSE)
message("[13m] species dot-line panel written: ", out, ".png | .pdf")
