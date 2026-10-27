#!/usr/bin/env Rscript
# Shared builder for the identity-band composition panel.
#
# One facet per weighting (ASVs / Reads / Prevalence), and inside each facet a
# pair of 100% stacked bars: the <97% band against the >=97% band.  The question
# the panel answers is not "how big is each band" -- the curves above already say
# that -- but "who is in it": whether the poorly matched tail is the same phyla
# as the well matched bulk, or a different set altogether.
#
# The five phyla with the most ASVs are named and the other ten pooled, on ASV
# counts for all three facets, so the three facets stay comparable.  Stack order
# is held fixed across every bar -- top to bottom by ASV rank, pooled remainder
# last -- so a phylum sits at the same place in all six bars and the bars can be
# read straight across.
#
# Sourced by 13z2_plot_homd_phylum_band_composition.R, which renders it on its
# own, and by 13z_plot_homd_phylum_coverage_panel_v3.R, which places it inside
# panel B, so the two cannot drift apart.
#
# Run 13z2_homd_phylum_band_composition.py first.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(ggh4x)   # strip_themed(): one strip fill per facet
})

# ---- palette ---------------------------------------------------------------
# Six series is past the point where a categorical set can clear the floors on
# every pair, so the gate applied here is the one for stacks and bars: every
# ADJACENT pair in the stack.  This order was picked by enumerating the
# documented slots against the pooled grey and keeping only the orderings that
# clear both hard gates -- worst adjacent dE 11.0 under simulated colour-vision
# deficiency (floor 8), worst pair of any two dE 16.3 for normal vision (floor
# 15).  Stack direction does not change which pairs are adjacent, so the order
# is also the legend order.
#
# Three of the six sit below 3:1 against white: magenta, and yellow worst of
# all.  The relief rule therefore applies -- every segment worth reading carries
# its own number.  The segments sit directly on one another with no separating
# rule: a white hairline broke each bar into blocks and read as a gap rather
# than as one stack.
PHY_COL <- c(
  "Bacillota"       = "#2a78d6",
  "Bacteroidota"    = "#4a3aa7",
  "Pseudomonadota"  = "#e87ba4",
  "Actinomycetota"  = "#008300",
  "Fusobacteriota"  = "#eda100",
  "Other phyla"     = "#787878")
PHY_ORD <- names(PHY_COL)          # level 1 draws at the TOP of the stack

COL_TXT  <- "#1A1A1A"
COL_GRID <- "#E4E4E4"
INK_LT   <- "#FFFFFF"
INK_DK   <- "#1A1A1A"
LAB_MIN  <- 0.050                  # below this a share gets no number
LAB_SM   <- 0.090                  # below this the number drops a size

# the two identity ranges, as the rest of the panel writes them
BAND_LAB <- c("<97% identity" = "70–97% id",
              ">=97% identity" = "≥97% id")

# the three facet strips wear the same three colours the weightings wear on the
# curves above, so a facet and its curve are read as the same quantity
MET_COL <- c("ASVs" = "#4D4D4D", "Reads" = "#0072B2", "Prevalence" = "#D55E00")
# what the facet strips say; the colours keep the short names as their keys
MET_LAB <- c(ASVs = "ASVs matched", Reads = "Reads represented",
             Prevalence = "Samples matched")

# ---- data ------------------------------------------------------------------
band_comp_read <- function(path) {
  comp <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  comp$metric     <- factor(comp$metric, levels = names(MET_COL),
                            labels = unname(MET_LAB[names(MET_COL)]))
  comp$band       <- factor(comp$band, levels = names(BAND_LAB),
                            labels = unname(BAND_LAB))
  comp$phylum_grp <- factor(comp$phylum_grp, levels = PHY_ORD)
  stopifnot(!anyNA(comp$phylum_grp), nrow(comp) == 3 * 2 * length(PHY_ORD))
  comp
}

band_comp_panel <- function(comp, base_size = 15) {
  # ggplot draws factor level 1 at the TOP, so each segment's box runs from the
  # running total of the levels above it down to the running total including it
  place <- function(d) {
    d <- d[order(d$phylum_grp), ]
    above <- c(0, head(cumsum(d$share), -1))
    d$ymax <- 1 - above
    d$ymin <- d$ymax - d$share
    d$ymid <- (d$ymax + d$ymin) / 2
    d
  }
  BARS <- comp %>%
    group_by(metric, band) %>%
    group_modify(~ place(.x)) %>%
    ungroup()

  # label ink is whichever of white / near-black reads better on the segment, so
  # the pale yellow never carries white text
  lum <- function(hex) {
    v <- strtoi(substring(hex, c(2, 4, 6), c(3, 5, 7)), 16L) / 255
    v <- ifelse(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055)^2.4)
    sum(c(0.2126, 0.7152, 0.0722) * v)
  }
  contrast <- function(a, b) {
    la <- lum(a); lb <- lum(b)
    (max(la, lb) + 0.05) / (min(la, lb) + 0.05)
  }
  BARS$ink <- vapply(BARS$phylum_grp, function(p) {
    fill <- PHY_COL[[as.character(p)]]
    if (contrast(fill, INK_LT) >= contrast(fill, INK_DK)) INK_LT else INK_DK
  }, character(1))

  # segment numbers are absolute point sizes on the device, so they follow the
  # rest of the panel down when the panel is set smaller
  k <- base_size / 15
  LAB <- BARS %>% filter(share >= LAB_MIN) %>%
    mutate(txt = sprintf("%.1f", 100 * share),
           sz  = ifelse(share < LAB_SM, 2.9 * k, 3.5 * k))

  ggplot(BARS, aes(x = band, y = 100 * share, fill = phylum_grp)) +
    geom_col(width = 0.62) +
    geom_text(data = LAB, aes(y = 100 * ymid, label = txt, colour = ink,
                              size = sz),
              fontface = "bold", show.legend = FALSE) +
    scale_size_identity() +
    # two rows under the bars, left to right in stack order.  At the panel's
    # final canvas width one row of six keys no longer fits: the last key was
    # cut off at the canvas edge, so the keys wrap to two rows of three.
    scale_fill_manual(values = PHY_COL, breaks = PHY_ORD, name = NULL,
                      drop = FALSE,
                      guide = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_colour_identity() +
    scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100),
                       expand = expansion(mult = c(0, 0)),
                       labels = function(v) paste0(v, "%")) +
    # each strip is filled with its own weighting's colour
    ggh4x::facet_wrap2(~ metric, nrow = 1,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = unname(MET_COL), colour = NA),
        # one point smaller than before: "Samples containing" ran edge to edge
        # in its strip, so the longest label now has a margin
        text_x = ggh4x::elem_list_text(colour = "white", face = "bold",
                                       size = base_size - 3))) +
    # every bar is normalised to its own range, so the axis names the division the
    # bars make: relative to everything in that range, this much is this phylum
    labs(y = "Relative proportion of\neach phylum (%)", x = NULL) +
    theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = COL_GRID, linewidth = 0.3),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.text = element_text(colour = COL_TXT),
          axis.text.x = element_text(colour = COL_TXT, size = base_size - 4,
                                     angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = base_size - 2.5, colour = COL_TXT),
          # the legend sits under the bars, in the order the segments are stacked
          legend.position = "bottom",
          # one line of six keys is wider than three facets unless the keys and
          # the gaps between them are kept small; at the old sizes the last key
          # was cut off at the canvas edge
          legend.key.size = unit(0.38, "cm"),
          legend.text = element_text(size = base_size - 5.5),
          legend.margin = margin(b = 2),
          legend.box.spacing = unit(4, "pt"),
          legend.spacing.x = unit(0.12, "cm"),
          panel.spacing = unit(0.35, "lines"))
}
