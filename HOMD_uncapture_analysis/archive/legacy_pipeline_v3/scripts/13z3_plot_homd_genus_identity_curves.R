#!/usr/bin/env Rscript
# Per-genus best-hit identity curves, the genus-level counterpart of the
# top-left panel of panel B.
#
#   x  best-hit identity to eHOMD, 95-100 %
#   y  share of the genus's ASVs with a best hit >= 95 % that reach that identity
#
# Because the denominator is the genus's own >=95 % matches, every curve starts
# at 100 % at x = 95 and the panel asks how those matches are spread across
# identity, not how many of them the genus has (that is panel B's question).
#
# One curve per genus.  The five genera with the highest mean share over the
# interval are named on the right, the five with the lowest on the left; the
# rest stay grey.  A genus needs at least MIN_N ASVs at >=95 % identity to be
# named -- a genus represented by one or two ASVs produces a step function that
# says nothing about the genus.
#
# The panel on the right annotates the y axis: each genus is read once, at the
# dashed line (97 % identity), and the bars count how many genera sit in each
# 5-point band of the y axis at that threshold.  It is attached with aplot so
# the two panels share one y scale.
#
# Run scripts/13z3_homd_genus_identity_curves.py first.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(aplot)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
IN      <- file.path(PIPE_V3, "13z3_homd_genus_identity_curves.tsv")

COL_TXT  <- "#1A1A1A"
COL_GRID <- "#E4E4E4"
COL_REF  <- "#8C8C8C"
COL_OTHR <- "#9A9A9A"
COL_HI   <- "#B2182B"      # the flattest curves, named on the right
COL_LO   <- "#2166AC"      # the earliest to fall away, named on the left

d <- read.delim(IN, sep = "\t", stringsAsFactors = FALSE)
d$y <- 100 * d$frac

# ---- which genera to name ---------------------------------------------------
# Ranked by the mean share over 95-100: the curve's position over the whole
# interval.  Every curve starts at 100 %, so the value at a single x cannot
# order them.
N_LAB <- 5
MIN_N <- 20
# n_asvs is the genus's ASV count in the >=95 % universe (constant down its
# 51 rows), not the number of rows in this table
stats <- d %>% group_by(genus) %>%
  summarise(auc = mean(y), n = first(n_asvs), .groups = "drop") %>%
  filter(n >= MIN_N) %>%
  arrange(auc)
hi <- tail(stats$genus, N_LAB)          # flattest, labelled on the right
lo <- head(stats$genus, N_LAB)          # fall away first, labelled on the left
d$side <- ifelse(d$genus %in% hi, "highest",
                 ifelse(d$genus %in% lo, "lowest", "other"))
d$side <- factor(d$side, levels = c("highest", "lowest", "other"))

# The highest curves lie on top of one another until they plunge at 99-100 %,
# so they are read at the right-hand end of the axis; the lowest have already
# separated half a point into the interval.  Anchoring each name in its own
# group's separating region keeps every leader short.
wrap <- function(g) ifelse(nchar(g) > 18, sub(" \\[", "\n[", g), g)
lab_hi <- d %>% filter(side == "highest", abs(identity - 100) < 1e-9) %>%
  transmute(genus, side, x = identity, y, lab = wrap(genus))
lab_lo <- d %>% filter(side == "lowest", abs(identity - 95.5) < 1e-9) %>%
  transmute(genus, side, x = identity, y, lab = wrap(genus))
stopifnot(nrow(lab_hi) == N_LAB, nrow(lab_lo) == N_LAB)

XL <- 92.95                             # label margin, left  (names end at 94.85)
XR <- 102.60                            # label margin, right (names start at 100.14)

p <- ggplot(d, aes(x = identity, y = y)) +
  geom_vline(xintercept = 97, linetype = "dashed", colour = COL_REF,
             linewidth = 0.35) +
  # the unnamed genera go first so a named curve is never hidden under them
  geom_line(data = filter(d, side == "other"), aes(group = genus),
            colour = COL_OTHR, linewidth = 0.22, alpha = 0.40) +
  geom_line(data = filter(d, side != "other"), aes(group = genus, colour = side),
            linewidth = 0.65) +
  scale_colour_manual(values = c(highest = COL_HI, lowest = COL_LO),
                      guide = "none") +
  # the two blank bands are the naming margins; the curves still run 95-100
  scale_x_continuous(breaks = 95:100, limits = c(XL, XR),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  geom_text_repel(data = lab_hi, aes(x = x, y = y, label = lab, colour = side),
                  hjust = 0, nudge_x = 0.14, direction = "y",
                  size = 2.9, fontface = "italic", seed = 1,
                  min.segment.length = 0, segment.size = 0.25,
                  segment.colour = "grey65", box.padding = 0.18,
                  point.padding = 0.08, max.overlaps = Inf) +
  # pushed all the way into the left margin: at -0.12 the names were written
  # across the curves they were pointing at
  geom_text_repel(data = lab_lo, aes(x = x, y = y, label = lab, colour = side),
                  hjust = 1, nudge_x = -0.65, direction = "y",
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
# Every genus contributes one point: where its curve sits at the dashed line.
# The bands are 5 points of the y axis wide, and the top band takes 100 %.
BIN <- 5
NB  <- 100 / BIN
y97 <- d$y[abs(d$identity - 97) < 1e-9]
stopifnot(length(y97) == length(unique(d$genus)))
hv <- data.frame(lo = seq(0, 100 - BIN, BIN))
hv$hi <- hv$lo + BIN
hv$n  <- as.integer(table(factor(pmin(floor(y97 / BIN), NB - 1),
                                 levels = seq_len(NB) - 1)))

ph <- ggplot(hv) +
  geom_rect(aes(xmin = 0, xmax = n, ymin = lo, ymax = hi),
            fill = "#B4B4B4", colour = NA) +
  scale_x_continuous(breaks = pretty(c(0, max(hv$n)), 3),
                     expand = expansion(mult = c(0, 0.04))) +
  # the same limits and expansion as the curves, so aplot can align the panels
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Number of genera", y = NULL) +
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

out <- file.path(FIGOUT, "fig_v3_homd_genus_identity_curves")
# 6.9 in for the curves (each naming margin has to hold an 18-character name at
# this type size) plus the histogram on the right
ggsave(paste0(out, ".png"), p, width = 8.3, height = 4.6, dpi = 300, bg = "white")
ggsave(paste0(out, ".pdf"), p, width = 8.3, height = 4.6, device = cairo_pdf,
       bg = "white")
message("[13z3] genus identity curves written: ", out, ".png | .pdf")
message("[13z3] highest: ", paste(hi, collapse = ", "))
message("[13z3] lowest:  ", paste(lo, collapse = ", "))
