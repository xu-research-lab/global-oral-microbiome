#!/usr/bin/env Rscript
# pipeline_v3: panel B, third take -- eHOMD coverage of the oral ASV catalogue.
#
# Three parts, each drawn with the mark its question wants:
#
#   p1  the three weightings as cumulative curves, identity on x, with the value
#       each curve takes at the 97% species boundary called out        (curves)
#   p2  the same quantity per phylum, the three weightings joined across the 15
#       phyla by a stepped line each                                     (steps)
#   p3  the composition of each identity band -- who is in it, not how big it is
#       -- under each of the three weightings                            (stacks)
#
# p1 and p2 run side by side along the top; p3 is three facets wide and takes the
# full width below them, which is the only way its bars are not squeezed.
#
# p3 is built by 13z2_panel_homd_band_composition.R, which the standalone preview
# 13z2_plot_homd_phylum_band_composition.R also sources.
#
# Run 13r_homd_lt_fractions.py and 13z2_homd_phylum_band_composition.py first.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
source("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/scripts/13z2_panel_homd_band_composition.R")

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

# ---- palette ---------------------------------------------------------------
# grey / blue / orange for the three weightings, as elsewhere in Figure 2; the
# composition's facet strips are filled with the same three, so a facet and its
# curve read as the same quantity.  The two identity bands no longer carry a
# colour of their own -- the composition labels them on the axis instead.
MET_COL <- c("ASVs" = "#4D4D4D", "Reads" = "#0072B2", "Prevalence" = "#D55E00")
COL_TXT    <- "#1A1A1A"
COL_GRID   <- "#E4E4E4"
COL_REF    <- "#8C8C8C"        # the 97% guide
SZ_AXT <- 11.7                 # axis-title size, shared by the three panels

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
band_comp <- band_comp_read(file.path(PIPE_V3,
                                      "13z2_homd_phylum_band_composition.tsv"))

# every fraction is read off the identity grid at the species boundary
at97 <- function(df) df[abs(df$identity - 97) < 1e-9, , drop = FALSE]

tot97 <- at97(tot_frac)
# the boundary figures are drawn on p1 rather than repeated in the legend
MET_LAB <- c(ASVs = "ASVs matched", Reads = "Reads represented",
             Prevalence = "Samples containing matched ASVs")
metric_scale <- function() {
  scale_colour_manual(values = MET_COL, labels = MET_LAB, name = NULL)
}
# shared chrome: black labels, a heavy black frame, recessive grid
house <- function(base = 15) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = COL_GRID, linewidth = 0.3),
          panel.border = element_rect(colour = "black", fill = NA,
                                      linewidth = 0.8),
          axis.text = element_text(colour = COL_TXT),
          axis.title = element_text(size = SZ_AXT, colour = COL_TXT),
          legend.text = element_text(colour = COL_TXT))
}

# ---- p1 | the catalogue as a whole ------------------------------------------
# cumulative curves: for each identity on x, the share of the weighting whose
# best hit is at least that identity.  Identity on x is the way round a reader
# expects a distribution to run.
CUM <- bind_rows(lapply(names(MET_COL), function(m) {
  data.frame(metric = m,
             identity = tot_frac$identity,
             prop = 1 - tot_frac[[c(ASVs = "frac_n_lt", Reads = "frac_r_lt",
                                    Prevalence = "frac_s_lt")[[m]]]])
}))
CUM$metric <- factor(CUM$metric, levels = names(MET_COL))
CUM <- CUM[CUM$identity >= 70 & CUM$identity <= 100, ]

# where each curve sits on the species boundary.  The point keeps the weighting's
# colour and the number stays in ink, so identity is never colour alone.  The
# values are printed without the per-cent sign: with it they overran the margin.
CROSS <- data.frame(metric = factor(names(MET_COL), levels = names(MET_COL)),
                    identity = 97,
                    val = 100 * c(1 - tot97$frac_n_lt, 1 - tot97$frac_r_lt,
                                  1 - tot97$frac_s_lt))
CROSS$lab <- sprintf("%.1f", CROSS$val)

# the axis is opened past 100 so the three callouts sit in the margin beside the
# boundary line instead of on top of the curves.  At the panel's final canvas
# width 103 was not enough and the numbers were cut off at the frame.
X1MAX <- 106

p1 <- ggplot(CUM, aes(x = identity, y = 100 * prop, colour = metric)) +
  geom_vline(xintercept = 97, linetype = "dashed", colour = COL_REF,
             linewidth = 0.35) +
  # linewidth is in mm, not points: 0.75 mm is a 2.1 pt stroke on the canvas,
  # which on the page (the panel is shown at a third of its canvas) still read
  # as a heavy black band against 4-pt type
  geom_line(linewidth = 0.45) +
  # shape 16: the default glyph draws the theme's 0.5 mm border round every
  # point, which at these sizes is as thick as the point is wide
  geom_point(data = CROSS, aes(x = identity, y = val, colour = metric),
             size = 2.3, shape = 16, inherit.aes = FALSE) +
  geom_text(data = CROSS, aes(x = identity + 1.0, y = val, label = lab),
            hjust = 0, size = 3.8, colour = COL_TXT, inherit.aes = FALSE) +
  metric_scale() +
  scale_x_continuous(breaks = seq(70, 100, 10), limits = c(70, X1MAX),
                     expand = expansion(mult = c(0.01, 0))) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  # the dashed line at 97% is left unlabelled: the caption names it
  labs(x = "Best-hit identity to eHOMD (%)",
       y = "Proportion\nat ≥ x (%)") +
  house() +
  # the key now spans the top row rather than sitting above p1 alone (see the
  # assembly below), so the three weightings fit on one line
  theme(axis.title.y = element_text(margin = margin(r = 4)),
        legend.position = "top") +
  # the key has to show the stroke as well as the point.  It was square, 0.38 cm
  # on a side, and the 2-mm point filled it, so the legend read as three dots;
  # it is widened and the key glyphs are drawn smaller than the plotted ones.
  guides(colour = guide_legend(keywidth = unit(0.52, "cm"),
                               keyheight = unit(0.26, "cm"),
                               nrow = 1, byrow = TRUE,
                               override.aes = list(size = 1.8,
                                                   linewidth = 0.45)))

# ---- p2 | by phylum, the three weightings ------------------------------------
# One lollipop per phylum per weighting: a horizontal segment from the axis out to
# the value, with the point on its end.  Rows are ranked by the primary weighting,
# so the row order is that ranking and the other two show how far it depends on
# which weighting is used.
b97 <- at97(phy_frac) %>%
  transmute(phylum,
            ASVs       = 1 - frac_n_lt,
            Reads      = 1 - frac_r_lt,
            Prevalence = 1 - frac_s_lt)
ORD <- b97$phylum[order(b97$ASVs)]          # last level sits at the top
PHY <- bind_rows(lapply(names(MET_COL), function(m) {
  data.frame(phylum = factor(b97$phylum, levels = ORD), metric = m, prop = b97[[m]],
             stringsAsFactors = FALSE)
}))
PHY$metric <- factor(PHY$metric, levels = names(MET_COL))
stopifnot(nrow(b97) == 15)
# All three segments of a row start at zero, so at the same y they would lie on
# top of one another and only the longest colour would show.  They are dodged
# inside the row instead; the axis ticks stay at the row's centre.  The offset is
# written into y rather than left to position_dodge, which moves y but not yend
# and would leave each segment pointing at a different row.
DODGE <- 0.24
PHY$yi <- as.integer(PHY$phylum)            # 15 at the top
PHY$yd <- PHY$yi + (as.integer(PHY$metric) - 2) * DODGE
PHY <- PHY[order(PHY$yd), ]

p2 <- ggplot(PHY, aes(x = 100 * prop, y = yd, colour = metric)) +
  geom_segment(aes(x = 0, xend = 100 * prop, yend = yd), linewidth = 0.45) +
  geom_point(size = 2.2, shape = 16) +
  metric_scale() +
  scale_x_continuous(breaks = seq(0, 100, 25), limits = c(0, 100),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = seq_along(ORD), labels = ORD,
                     expand = expansion(mult = c(0.035, 0.035))) +
  # the title has to fit the half-width panel, not the whole canvas
  labs(x = "Proportion at ≥97%\nidentity (%)") +
  house() +
  theme(panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  guides(colour = "none")          # the key over the top row already says this

# ---- p3 | who is in each identity band ---------------------------------------
# Stacked rather than dodged, because the phyla of a band sum to 100 and one bar
# per band says that in a way two bars side by side cannot.  The five largest
# phyla by ASV count are named and the other ten pooled, so the same six segments
# stand in every bar and the bars can be read straight across; the same six are
# used under all three weightings, so the facets stay comparable.
p3 <- band_comp_panel(band_comp)

# ---- assembly ---------------------------------------------------------------
# p1 and p2 along the top, p3 across the full width beneath them: p3 is three
# facets wide and needs the whole width to keep its bars off the floor.  p3's six
# phylum colours stay under its bars; only p1's three-curve key is collected, and
# only across the top row.
#
# No B1/B2/B3 tags: their row of margin above the top two plots is space the
# plots can have instead, and the sub-panels read in order on their own.  The
# figure's own letter B is placed by the compositor, outside this canvas.
# The three-curve key belongs to p1 alone, but at the top of a half-width panel
# its three entries did not fit on one line.  Collected over the top row it is
# centred above both plots and has the full width, so it reads on one line.
top_row <- (p1 | p2) + plot_layout(guides = "collect") &
  theme(legend.position = "top")

panelB <- top_row / p3 +
  plot_layout(heights = c(1.3, 1), guides = "keep") &
  # the six phylum keys under the bars have to fit the width of the three
  # facets, which they only do at a smaller key and tighter spacing; this also
  # sets the size of p1's three-curve key, so the two legends stay consistent
  theme(plot.margin = margin(2, 4, 2, 2),
        legend.key.size = unit(0.38, "cm"),
        legend.text = element_text(size = 9.0, colour = COL_TXT),
        legend.spacing.x = unit(0.12, "cm"),
        legend.margin = margin(b = -2))

out <- file.path(FIGOUT, "fig_v3_panelB_coverage_v4")
# The canvas is the size B is finally shown at, not the size it was designed at:
# drawn at 9 in and pasted into a 2.4-in column, its 15 pt type came out at 4 pt.
# At 6.43 in the type keeps its point size on the page, and the plot areas take
# the space the shrinking frees.
ggsave(paste0(out, ".png"), panelB, width = 6.43, height = 7.71, dpi = 300,
       bg = "white")
ggsave(paste0(out, ".pdf"), panelB, width = 6.43, height = 7.71,
       device = cairo_pdf, bg = "white")
message("[13z] panel B v4 written: ", out, ".png | .pdf")
