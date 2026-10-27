#!/usr/bin/env Rscript
# Per-species best-hit identity curves, the species-level counterpart of the
# genus panel.
#
#   x  best-hit identity to eHOMD, 97-100 %
#   y  share of the species' ASVs with a best hit >= 97 % that reach that
#      identity
#
# 97 % is the species boundary, so the panel starts there: every curve starts at
# 100 % at x = 97 (the species' own >=97 % matches are the denominator) and asks
# how those matches are spread across identity -- how much of what identifies to
# the species is close enough to be the same strain, and how much only reaches
# the species boundary.
#
# The dashed line at 99 % is the strain-level threshold: 97-99 % separates
# species, 99-100 % is conventionally the same strain.
#
# One curve per species.  The five species with the highest mean share over the
# interval are named on the right, the five with the lowest on the left; the
# rest stay grey.  A species needs at least MIN_N ASVs at >=97 % identity to be
# named -- a handful of ASVs gives a step function, not a curve.
#
# The panel on the right annotates the y axis: each species is read once, at
# the dashed line, and the bars count how many species sit in each 5-point band
# of the y axis at 99 % identity.  It is attached with aplot so the two panels
# share one y scale.
#
# Run scripts/13z4_homd_species_identity_curves.py first.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(aplot)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
IN      <- file.path(PIPE_V3, "13z4_homd_species_identity_curves.tsv")

COL_TXT  <- "#1A1A1A"
COL_GRID <- "#E4E4E4"
COL_REF  <- "#8C8C8C"
COL_STRAIN <- "#4D4D4D"   # the 99 % line, dark enough to read over the curves
COL_OTHR <- "#9A9A9A"
COL_HI   <- "#B2182B"      # the flattest curves, named on the right
COL_LO   <- "#2166AC"      # the earliest to fall away, named on the left

d <- read.delim(IN, sep = "\t", stringsAsFactors = FALSE)
d$y <- 100 * d$frac

# ---- which species to name --------------------------------------------------
# Ranked by the mean share over 97-100: the curve's position over the whole
# interval.  Every curve starts at 100 %, so the value at a single x cannot
# order them.
N_LAB <- 5
MIN_N <- 50
stats <- d %>% group_by(species) %>%
  summarise(auc = mean(y), n = first(n_asvs), .groups = "drop") %>%
  filter(n >= MIN_N) %>%
  arrange(auc)
hi <- tail(stats$species, N_LAB)        # flattest, labelled on the right
lo <- head(stats$species, N_LAB)        # fall away first, labelled on the left
d$side <- ifelse(d$species %in% hi, "highest",
                 ifelse(d$species %in% lo, "lowest", "other"))
d$side <- factor(d$side, levels = c("highest", "lowest", "other"))
message(sprintf("[13z4] naming %d of %d species with >= %d ASVs at >=97%%",
                N_LAB * 2, nrow(stats), MIN_N))

# Several names run to 40 characters ("Clostridiales [F1 G2] bacterium
# HMT-402"); on one line that would take a margin wider than the curves.  Long
# names break at the space nearest the middle, so both lines stay short.
wrap <- function(s) {
  if (nchar(s) <= 20) return(s)
  # never split inside a bracketed placeholder ("[F1 G2]"): break after the
  # closing bracket if there is one, otherwise at the space nearest the middle
  br <- gregexpr("\\]", s)[[1]]
  br <- if (br[1] > 0) tail(br, 1) else -1
  k <- if (br > 0) {
    m <- regexpr(" ", substring(s, br + 1))
    if (m > 0) br + m else -1
  } else {
    pos <- gregexpr(" ", s)[[1]]
    if (pos[1] < 0) -1 else pos[which.min(abs(pos - nchar(s) / 2))]
  }
  if (k < 0) return(s)
  paste0(substr(s, 1, k - 1), "\n", substr(s, k + 1, nchar(s)))
}
# The highest curves separate as they bend down through 99-100 %; the lowest
# have already opened up by 97.5 %.  Anchoring each name in its own group's
# separating region keeps every leader short.
lab_hi <- d %>% filter(side == "highest", abs(identity - 100) < 1e-9) %>%
  transmute(species, side, x = identity, y, lab = vapply(species, wrap, ""))
lab_lo <- d %>% filter(side == "lowest", abs(identity - 97.5) < 1e-9) %>%
  transmute(species, side, x = identity, y, lab = vapply(species, wrap, ""))
stopifnot(nrow(lab_hi) == N_LAB, nrow(lab_lo) == N_LAB)

XL <- 95.60                             # label margin, left  (names end at 97.1)
XR <- 101.35                            # label margin, right (names start at 100.1)

p <- ggplot(d, aes(x = identity, y = y)) +
  # the unnamed species go first so a named curve is never hidden under them
  geom_line(data = filter(d, side == "other"), aes(group = species),
            colour = COL_OTHR, linewidth = 0.18, alpha = 0.35) +
  geom_line(data = filter(d, side != "other"), aes(group = species, colour = side),
            linewidth = 0.6) +
  # drawn after the curves, in a darker grey than they are: the 99 % column is
  # where the curves stack up, and a light line under them is invisible
  geom_vline(xintercept = 99, linetype = "dashed", colour = COL_STRAIN,
             linewidth = 0.45) +
  scale_colour_manual(values = c(highest = COL_HI, lowest = COL_LO),
                      guide = "none") +
  # the two blank bands are the naming margins; the curves still run 97-100
  scale_x_continuous(breaks = 97:100, limits = c(XL, XR),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  geom_text_repel(data = lab_hi, aes(x = x, y = y, label = lab, colour = side),
                  hjust = 0, nudge_x = 0.08, direction = "y",
                  size = 2.9, fontface = "italic", seed = 1, lineheight = 0.9,
                  min.segment.length = 0, segment.size = 0.25,
                  segment.colour = "grey65", box.padding = 0.18,
                  point.padding = 0.08, max.overlaps = Inf) +
  geom_text_repel(data = lab_lo, aes(x = x, y = y, label = lab, colour = side),
                  hjust = 1, nudge_x = -0.4, direction = "y",
                  size = 2.9, fontface = "italic", seed = 1, lineheight = 0.9,
                  min.segment.length = 0, segment.size = 0.25,
                  segment.colour = "grey65", box.padding = 0.18,
                  point.padding = 0.08, max.overlaps = Inf) +
  labs(x = "Sequencing identity threshold (%)",
       y = "Proportion of total (%)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = COL_GRID, linewidth = 0.3),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.text = element_text(colour = COL_TXT),
        axis.title = element_text(size = 12, colour = COL_TXT),
        legend.position = "none",
        plot.margin = margin(4, 6, 4, 4))

# ---- the y-axis histogram ---------------------------------------------------
# Every species contributes one point: where its curve sits at the dashed line.
# The bands are 5 points of the y axis wide, and the top band takes 100 %.
BIN <- 5
NB  <- 100 / BIN
y99 <- d$y[abs(d$identity - 99) < 1e-9]
stopifnot(length(y99) == length(unique(d$species)))
hv <- data.frame(lo = seq(0, 100 - BIN, BIN))
hv$hi <- hv$lo + BIN
hv$n  <- as.integer(table(factor(pmin(floor(y99 / BIN), NB - 1),
                                 levels = seq_len(NB) - 1)))
xmax <- max(hv$n)
BRK  <- pretty(c(0, xmax), 4)
BRK  <- BRK[BRK <= xmax * 1.04]          # drop the tick the limit cuts off

ph <- ggplot(hv) +
  geom_rect(aes(xmin = 0, xmax = n, ymin = lo, ymax = hi),
            fill = "#B4B4B4", colour = NA) +
  scale_x_continuous(breaks = BRK,
                     expand = expansion(mult = c(0, 0.04))) +
  # the same limits and expansion as the curves, so aplot can align the panels
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Number of species", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = COL_GRID, linewidth = 0.3),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.text = element_text(colour = COL_TXT),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 12, colour = COL_TXT),
        legend.position = "none",
        plot.margin = margin(4, 4, 4, 3))

# aplot stacks them on the shared y axis rather than patching two independent
# plots, which keeps the 0-100 % scale identical in both panels
p <- insert_right(p, ph, width = 0.26)

out <- file.path(FIGOUT, "fig_v3_homd_species_identity_curves")
# the curves need the naming margins on both sides; the histogram sits outside
ggsave(paste0(out, ".png"), p, width = 9.6, height = 4.8, dpi = 300, bg = "white")
ggsave(paste0(out, ".pdf"), p, width = 9.6, height = 4.8, device = cairo_pdf,
       bg = "white")
message("[13z4] species identity curves written: ", out, ".png | .pdf")
message("[13z4] highest: ", paste(hi, collapse = ", "))
message("[13z4] lowest:  ", paste(lo, collapse = ", "))
