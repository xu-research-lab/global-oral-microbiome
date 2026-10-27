#!/usr/bin/env Rscript
# pipeline_v3: species-level eHOMD coverage split by coverage gap.
# x = fraction of ASVs <97% identity (eHOMD coverage gap);
# y = No. derep ASVs mapped (log10, ASV richness);
# the species are divided at x = 0.5 into a low-gap and a high-gap half; point
# size = eHOMD refseq count; the three richest species of each half are
# labelled.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(ggrepel)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

s <- read.delim(file.path(PIPE_V3, "13j_homd_species_coverage_fullname.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
s$frac_ge97 <- pmin(pmax(s$frac_ge97, 0), 1)
s$frac_lt97 <- 1 - s$frac_ge97

X0 <- 0.5
# One binary split only: at or below the 97% boundary on most ASVs, or not.
# The richness split that used to cut the panel in two is not drawn -- it made
# four colours for what is one question (does eHOMD cover this species?), and
# the species below it were named in a corner the reader cannot see.
s$side <- factor(ifelse(s$frac_lt97 >= X0, "High coverage gap", "Low coverage gap"),
                 levels = c("Low coverage gap", "High coverage gap"))
# red is the gap side, as everywhere else in the figure
side_cols <- c("Low coverage gap" = "#2166AC", "High coverage gap" = "#B2182B")

# ---- labels: the five richest species on each side ----
# Five per side is what the top half of the old two-way split carried, so the
# panel names the same species it used to; they are the extremes of ASV
# richness, which is what y measures.  The lower half now shares their side
# instead of having a set of its own, which leaves room for the two extra.
N_LAB <- 5
lab <- bind_rows(
  s %>% filter(side == "High coverage gap") %>% arrange(desc(n_asvs_total)) %>% slice_head(n = N_LAB),
  s %>% filter(side == "Low coverage gap")  %>% arrange(desc(n_asvs_total)) %>% slice_head(n = N_LAB)
) %>% ungroup()

# ---- half-panel summary annotation (safe positions, white label boxes) ----
qsum <- s %>% group_by(side) %>%
  summarise(n_sp = n(), n_asv = sum(n_asvs_total), .groups = "drop") %>%
  mutate(txt = sprintf("n=%d species\n%.1f%% ASVs", n_sp, 100 * n_asv / sum(n_asv)))
# One box per half, in the empty floor below the sparsest species: with ten
# labels the top of the panel is where they are repelled to, and a box there
# sat on top of them.
qpos <- data.frame(
  side  = levels(s$side),
  x     = c(0.02, 0.98),
  y     = 10^1.05,
  hjust = c(0, 1),
  vjust = c(0, 0),
  stringsAsFactors = FALSE
)
qsum <- qsum %>% left_join(qpos, by = "side")

# ---- main scatter ----
p_main <- ggplot(s, aes(x = frac_lt97, y = n_asvs_total)) +
  geom_vline(xintercept = X0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  # shape 16, not the default glyph: ggplot2 4 draws default points with the
  # theme's borderwidth (0.5 mm), so every bubble carried a rim as thick as a
  # quarter of its diameter
  geom_point(aes(colour = side, size = n_homd_refseq), alpha = 0.7, shape = 16) +
  # the labels are white boxes with no frame at all: at this scale any rule
  # round the text read as a heavy outline, and the fill alone separates the
  # label from the points.  The leader line is left as a hairline.
  geom_label_repel(data = lab, aes(label = species, colour = side),
                   size = 4.2, max.overlaps = Inf, segment.size = 0.10,
                   segment.colour = "grey60", segment.alpha = 0.9,
                   linewidth = 0,
                   min.segment.length = 0, box.padding = 0.55,
                   # no xlim/ylim: pinned to a bound the bottom three labels
                   # were forced onto one line and overlapped; free, ggrepel
                   # separates them (the lower y floor gives them the room)
                   force = 3, max.time = 5, max.iter = 1e5,
                   fill = "white", alpha = 0.85, fontface = "italic",
                   seed = 7) +
  geom_label(data = qsum, aes(x = x, y = y, label = txt, hjust = hjust, vjust = vjust),
             size = 4.6, colour = "grey25", fill = "white", alpha = 0.85,
             linewidth = 0, lineheight = 0.9) +
  # the title is two lines: on one line it is wider than the canvas and the
  # tail ("...identity") is cut off at the device edge
  scale_x_continuous(labels = percent_format(), limits = c(0, 1),
                     expand = expansion(mult = c(0.03, 0.03)),
                     name = "Proportion of ASVs with <97% best-hit identity to the species") +
  # the floor sits below the sparsest species so the labels for those rows have
  # somewhere to go: pinned to the old 50 floor they were forced into one line
  # and overlapped
  # the floor is below the sparsest species and the ceiling above the key, so
  # neither the summary boxes nor the labels crowd the tick range
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale()),
                limits = c(10^0.95, 10^6.4),
                breaks = c(100, 1e3, 1e4, 1e5, 1e6),
                name = "Number of best-hit ASVs (log10)") +
  scale_colour_manual(values = side_cols, name = NULL, guide = "none") +
  scale_size_continuous(range = c(0.8, 7), guide = "none") +
  # The canvas is drawn at 6.79 in and shown at 714 px in the composite, so a
  # point here is worth 0.35 of one on the page.  At 16 pt the axis titles came
  # out 19 px tall on the page against 13-15 px for the titles of the other
  # panels, i.e. the panel's type read a third larger than theirs.  12 pt puts
  # the titles at the same size as theirs and leaves the x title the same size
  # as the y title, which is what it is set to.
  theme_bw(base_size = 12) +
  # labels are allowed to sit outside the frame, and the margin is widened so
  # the ones at the edges land on white rather than on the canvas border
  theme(axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        plot.margin = margin(4, 16, 6, 4))

# The size key is drawn inside the panel rather than left to a legend: aplot puts
# any legend it finds at the right of the assembled block, which strands it in the
# gap between the density strip and the schematic column.  The key sits in the
# empty band above the highest species (the tallest point reaches ~3e4) and above
# the row the labels are repelled into, so the reference counts do not run into
# the species names.
p_key <- data.frame(n = c(1, 10, 50, 178), x = c(0.30, 0.44, 0.58, 0.72))
KEY_Y <- 10^6.07
p_main <- p_main +
  annotate("text", x = 0.29, y = 10^6.27,
           label = "eHOMD reference sequences per species",
           hjust = 0, size = 3.8, colour = "grey30") +
  geom_point(data = p_key, aes(x = x, y = KEY_Y, size = n), colour = "grey30",
             shape = 16) +
  geom_text(data = p_key, aes(x = x, y = KEY_Y / 2.4, label = n),
            size = 3.6, colour = "grey30")

# ---- the panel: the scatter alone -----------------------------------------
# The canvas is the size C is finally shown at (the 9.5-in design canvas divided
# by 1.4), so the type keeps its point size on the page instead of being shrunk
# with the image.
W_FIG  <- 6.79
H_MAIN <- 6.14
p <- p_main

fig <- file.path(FIGOUT, "fig_v3_homd_species_quadrant")
ggsave(paste0(fig, ".png"), p, width = W_FIG, height = H_MAIN, dpi = 300)
ggsave(paste0(fig, ".pdf"), p, width = W_FIG, height = H_MAIN, device = cairo_pdf)
message("[13l] species scatter panel written")

# ============================================================
# Alternative: species share scatter (fraction-of-total view)
#   x = n ASVs of species with best hit <97% / total ASVs <97%
#   y = n ASVs of species with best hit >=97% / total ASVs >=97%
#   log10 axes (species with zero >=97% ASVs are excluded)
#   + linear trend line in log-log space, Pearson R and P value
# ============================================================
s_full <- s %>% mutate(n_lt97 = pmax(n_asvs_total - n_asvs_ge97, 0))
tot_lt97 <- sum(s_full$n_lt97)
tot_ge97 <- sum(s_full$n_asvs_ge97)
sp <- s_full %>% filter(n_lt97 > 0, n_asvs_ge97 > 0)
sp$x_share <- sp$n_lt97 / tot_lt97
sp$y_share <- sp$n_asvs_ge97 / tot_ge97

# correlation and trend line are computed on log10-transformed shares
ct <- cor.test(log10(sp$x_share), log10(sp$y_share), method = "pearson")
rval  <- unname(ct$estimate)
pval  <- ct$p.value

# two dashed lines parallel to the diagonal that enclose ~95% of the points:
# offset = 95th percentile of the signed perpendicular distance to x=y in
# log-log space (central/symmetric envelope)
logx <- log10(sp$x_share)
logy <- log10(sp$y_share)
perp <- (logy - logx) / sqrt(2)          # signed distance from diagonal
off95 <- as.numeric(quantile(abs(perp), 0.95))
sp$outside <- abs(perp) > off95
lab10 <- sp %>% filter(outside) %>%
  arrange(desc(abs((log10(y_share) - log10(x_share)) / sqrt(2)))) %>%
  slice_head(n = 10)
xx <- 10^seq(log10(1e-6), log10(0.1), length.out = 300)
band <- data.frame(x = xx,
                   y_lo = xx * 10^(-off95 * sqrt(2)),
                   y_hi = xx * 10^(off95 * sqrt(2)))

p_alt <- ggplot(sp, aes(x = x_share, y = y_share)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40", linewidth = 0.5) +
  geom_line(data = band, aes(x = x, y = y_lo),
            linetype = "dashed", colour = "grey20", linewidth = 0.4) +
  geom_line(data = band, aes(x = x, y = y_hi),
            linetype = "dashed", colour = "grey20", linewidth = 0.4) +
  geom_point(data = sp %>% filter(!outside),
             size = 1.4, alpha = 0.5, colour = "#1f77b4") +
  geom_point(data = sp %>% filter(outside),
             size = 1.4, alpha = 0.7, colour = "#d62728") +
  geom_smooth(method = "lm", se = TRUE, colour = "#d62728",
              fill = "#d6272833", linewidth = 0.8) +
  geom_text_repel(data = lab10, aes(label = species), size = 2.1,
                  max.overlaps = Inf, seed = 42, colour = "grey15",
                  segment.colour = "grey30", segment.size = 0.2,
                  min.segment.length = 0, box.padding = 0.25) +
  annotate("text", x = 2e-6, y = 0.05,
           label = sprintf("R = %.3f\nP = %.2e", rval, pval),
           hjust = 0, vjust = 1, size = 3.5, colour = "grey15") +
  scale_x_log10(labels = percent_format(), limits = c(1e-6, 0.1),
                oob = scales::squish,
                name = "Fraction of <97%-identity ASVs (per species)") +
  scale_y_log10(labels = percent_format(), limits = c(1e-6, 0.1),
                oob = scales::squish,
                name = "Fraction of >=97%-identity ASVs (per species)") +
  coord_fixed() +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

fig_alt <- file.path(FIGOUT, "fig_v3_homd_species_share_ge_lt")
ggsave(paste0(fig_alt, ".png"), p_alt, width = 6, height = 6, dpi = 300)
ggsave(paste0(fig_alt, ".pdf"), p_alt, width = 6, height = 6)
message("[13l] species share scatter with trend line written")
