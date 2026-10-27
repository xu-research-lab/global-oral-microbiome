#!/usr/bin/env Rscript
# The per-genus identity curves as a square panel, for the composite figure.
#
# The same curves as 13z3_plot_homd_genus_identity_curves.R, with the ten named
# genera marked by the numbers 1-10 and their names written in a key strip above
# the panel.  Spelling the names at the ends of the curves cost a margin on each
# side wider than the curves themselves; a number is two characters, so the
# curves can use the width, and the panel comes out square.  The right-hand
# histogram of the standalone version is left out.
#
#   1-5  the five flattest curves        (red, named at the right end)
#   6-10 the five that fall away first   (blue, named at the left)
#
# Run scripts/13z3_homd_genus_identity_curves.py first.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(patchwork)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
IN      <- file.path(PIPE_V3, "13z3_homd_genus_identity_curves.tsv")

COL_TXT  <- "#1A1A1A"
COL_GRID <- "#E4E4E4"
COL_REF  <- "#8C8C8C"
COL_OTHR <- "#9A9A9A"
COL_HI   <- "#B2182B"      # the flattest curves
COL_LO   <- "#2166AC"      # the earliest to fall away

d <- read.delim(IN, sep = "\t", stringsAsFactors = FALSE)
d$y <- 100 * d$frac

# ---- which genera to number -------------------------------------------------
N_LAB <- 5
MIN_N <- 20
stats <- d %>% group_by(genus) %>%
  summarise(auc = mean(y), n = first(n_asvs), .groups = "drop") %>%
  filter(n >= MIN_N) %>%
  arrange(auc)
# 1 is the flattest curve, 6 the one that falls away first
hi <- rev(tail(stats$genus, N_LAB))
lo <- head(stats$genus, N_LAB)
NUM <- c(setNames(seq_along(hi), hi), setNames(N_LAB + seq_along(lo), lo))
d$side <- ifelse(d$genus %in% hi, "highest",
                 ifelse(d$genus %in% lo, "lowest", "other"))
d$side <- factor(d$side, levels = c("highest", "lowest", "other"))

# numbers sit in the margins, with a leader back to the curve where the end of a
# curve is crowded by its neighbours
lab_hi <- d %>% filter(side == "highest", abs(identity - 100) < 1e-9) %>%
  transmute(side, x = identity, y, lab = unname(NUM[genus]))
lab_lo <- d %>% filter(side == "lowest", abs(identity - 95.5) < 1e-9) %>%
  transmute(side, x = identity, y, lab = unname(NUM[genus]))
stopifnot(nrow(lab_hi) == N_LAB, nrow(lab_lo) == N_LAB)

XL <- 94.30                             # margin wide enough for a two-digit label
XR <- 100.45

p <- ggplot(d, aes(x = identity, y = y)) +
  # the unnamed genera go first so a numbered curve is never hidden under them
  geom_line(data = filter(d, side == "other"), aes(group = genus),
            colour = COL_OTHR, linewidth = 0.22, alpha = 0.40) +
  geom_line(data = filter(d, side != "other"), aes(group = genus, colour = side),
            linewidth = 0.65) +
  scale_colour_manual(values = c(highest = COL_HI, lowest = COL_LO),
                      guide = "none") +
  scale_x_continuous(breaks = 95:100, limits = c(XL, XR),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  geom_text_repel(data = lab_hi, aes(x = x, y = y, label = lab, colour = side),
                  hjust = 0, nudge_x = 0.06, direction = "y",
                  size = 3.3, fontface = "bold", seed = 1,
                  min.segment.length = 0, segment.size = 0.25,
                  segment.colour = "grey60", box.padding = 0.10,
                  point.padding = 0.05, max.overlaps = Inf) +
  geom_text_repel(data = lab_lo, aes(x = x, y = y, label = lab, colour = side),
                  hjust = 1, nudge_x = -0.62, direction = "y",
                  size = 3.3, fontface = "bold", seed = 1,
                  min.segment.length = 0, segment.size = 0.25,
                  segment.colour = "grey60", box.padding = 0.10,
                  point.padding = 0.05, max.overlaps = Inf) +
  geom_vline(xintercept = 97, linetype = "dashed", colour = COL_REF,
             linewidth = 0.4) +
  labs(x = "Sequencing identity threshold (%)",
       y = "Proportion of total (%)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = COL_GRID, linewidth = 0.3),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.text = element_text(colour = COL_TXT),
        axis.title = element_text(size = 13, colour = COL_TXT),
        legend.position = "none",
        plot.margin = margin(2, 4, 2, 2))

# ---- the histogram is not drawn ---------------------------------------------
# The curves carry the panel alone.  The bar chart of how many genera sit in each
# band of the y axis sat to the right of them; it is not part of this figure.
main <- p

# ---- the key strip ----------------------------------------------------------
# Ten names on two lines would need a line wider than the panel to stay
# readable, so they are set as two columns of five: 1-5 in the colour of the
# flattest curves, 6-10 in the colour of the ones that fall away first.
MAP <- data.frame(
  grp = factor(rep(c("highest", "lowest"), each = N_LAB),
               levels = c("highest", "lowest")),
  x   = rep(c(0.02, 0.52), each = N_LAB),
  y   = rep(rev(seq_len(N_LAB)), 2),
  txt = c(sprintf("%d %s", NUM[hi], hi), sprintf("%d %s", NUM[lo], lo)))

strip <- ggplot(MAP, aes(x = x, y = y, label = txt, colour = grp)) +
  geom_text(size = 4.3, fontface = "italic", hjust = 0, vjust = 0.5) +
  scale_colour_manual(values = c(highest = COL_HI, lowest = COL_LO),
                      guide = "none") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(0)) +
  scale_y_continuous(limits = c(0.4, 5.6), expand = expansion(0)) +
  theme_void() +
  theme(plot.margin = margin(2, 0, 2, 0))

# patchwork rather than aplot here: aplot's insert_top re-uses the main panel's
# x scale for the strip it inserts, which would place the key off the axis
whole <- strip / main + plot_layout(heights = c(0.24, 1))

out <- file.path(FIGOUT, "fig_v3_homd_genus_identity_curves_numbered")
ggsave(paste0(out, ".png"), whole, width = 5.8, height = 5.8, dpi = 300,
       bg = "white")
ggsave(paste0(out, ".pdf"), whole, width = 5.8, height = 5.8,
       device = cairo_pdf, bg = "white")
message("[13z3n] numbered genus panel written: ", out, ".png | .pdf")
message("[13z3n] key highest: ", paste(sprintf("%d %s", NUM[hi], hi), collapse = "; "))
message("[13z3n] key lowest:  ", paste(sprintf("%d %s", NUM[lo], lo), collapse = "; "))
