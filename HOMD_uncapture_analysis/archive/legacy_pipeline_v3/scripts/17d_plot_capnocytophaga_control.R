#!/usr/bin/env Rscript
# Control for cluster 21907: the part of Capnocytophaga that eHOMD does cover.
#
# One row of four panels, all on the OSCC cohort (913 cases / 431 HC):
#   a  how many >=97%-identity Capnocytophaga ASVs each sample carries
#   b  their combined relative abundance, cases vs controls, with the test
#   c  whether any of them is detected (stacked bars, Fisher P)
#   d  the same bars for cluster 21907, the uncaptured Capnocytophaga
#
# Panels c and d are the same statistic on the two halves of the genus, drawn
# the same way, so the contrast is read off one row.  Each bar is split into
# detected and not detected; the group sizes sit under the bars and the raw
# between-group Fisher test above them.
#
# Numbers come from scripts/17d_capnocytophaga_control.py; the cohort, the
# detection rule and the colours are 17b_plot_cluster_prevalence.R's.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIG  <- file.path(PIPE, "analysis_report", "figs")
SAMP <- file.path(PIPE, "17d_capnocytophaga_control_samples.tsv")
PS   <- file.path(PIPE, "17d_capnocytophaga_control_persample.tsv")
REL  <- file.path(PIPE, "14_cluster_relabund_dedup.tsv")
META <- file.path(PIPE, "14_cluster_sample_metadata.tsv")

COL      <- c(Disease = "#B2182B", HC = "#4393C3")
COL_ABS  <- "#D9D9D9"
COL_TXT  <- "#1A1A1A"
XL       <- c(Disease = "OSCC", HC = "HC")

ps <- read.delim(PS, stringsAsFactors = FALSE)
ps$group <- factor(ps$group, levels = c("Disease", "HC"))
ns <- tapply(ps$sample_key, ps$group, length)
n_det <- tapply(ps$n_asv > 0, ps$group, sum)
xlab2 <- c(Disease = sprintf("%s\n(n = %s)", XL["Disease"], comma(ns["Disease"])),
           HC      = sprintf("%s\n(n = %s)", XL["HC"], comma(ns["HC"])))
xlab3 <- c(Disease = sprintf("%s\n(%s detected)", XL["Disease"],
                             comma(n_det["Disease"])),
           HC      = sprintf("%s\n(%s detected)", XL["HC"], comma(n_det["HC"])))

base_th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 8, colour = COL_TXT, lineheight = 0.95),
        axis.text.y = element_text(size = 8, colour = COL_TXT),
        axis.title.y = element_text(size = 8.6, colour = COL_TXT),
        plot.title = element_text(size = 9.2, face = "bold", colour = COL_TXT),
        plot.subtitle = element_text(size = 7.2, colour = "grey35",
                                     lineheight = 1.02, margin = margin(b = 4)),
        plot.margin = margin(4, 7, 4, 5))

# ---- a: ASVs detected per sample ------------------------------------------
p_a <- ggplot(ps, aes(group, n_asv, fill = group)) +
  geom_boxplot(width = 0.5, linewidth = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0.2, size = 0.4, alpha = 0.25,
              show.legend = FALSE) +
  scale_fill_manual(values = COL, guide = "none") +
  scale_x_discrete(labels = xlab2) +
  scale_y_continuous(limits = c(-0.5, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "ASVs detected per sample",
       title = "ASVs per sample",
       subtitle = "eHOMD-covered Capnocytophaga") +
  base_th

# ---- b: relative abundance -------------------------------------------------
ab <- ps %>% filter(rel_abund > 0) %>% mutate(y = 100 * rel_abund)
wt <- wilcox.test(y ~ group, data = ab)
p_lab <- sprintf("Mann-Whitney P = %.1e", wt$p.value)
ybar <- max(ab$y) * 2.2

p_b <- ggplot(ab, aes(group, y, fill = group)) +
  geom_boxplot(width = 0.5, linewidth = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.25, show.legend = FALSE) +
  # the bracket is drawn in data units on the log axis; the drop is a factor,
  # so the two feet are equal drops on a log scale
  annotate("segment", x = 1, xend = 2, y = ybar, yend = ybar, linewidth = 0.3) +
  annotate("segment", x = 1, xend = 1, y = ybar / 1.8, yend = ybar, linewidth = 0.3) +
  annotate("segment", x = 2, xend = 2, y = ybar / 1.8, yend = ybar, linewidth = 0.3) +
  annotate("text", x = 1.5, y = ybar * 1.25, label = p_lab, size = 2.4,
           colour = "grey20") +
  scale_fill_manual(values = COL, guide = "none") +
  scale_x_discrete(labels = xlab3) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10),
                labels = c("0.01", "0.1", "1", "10"),
                limits = c(NA, ybar * 1.8),
                expand = expansion(mult = c(0.04, 0))) +
  labs(x = NULL, y = "Relative abundance (% of reads, log10)",
       title = "Relative abundance",
       subtitle = "Covered ASVs, detected samples") +
  base_th

# ---- c/d: stacked bars, group sizes on the axis only ----------------------
fmt_p <- function(p) if (p < 1e-3) sprintf("%.1e", p) else sprintf("%.2f", p)

bar_panel <- function(d, title, subtitle) {
  d$absent <- d$n - d$detected
  d$p_det <- d$detected / d$n
  # the two arms are independent samples, so the raw contingency test is Fisher's;
  # no covariates enter here, this figure is descriptive
  ft <- fisher.test(matrix(c(d$detected[1], d$absent[1],
                             d$detected[2], d$absent[2]), nrow = 2))
  lab <- paste0("Fisher P = ", fmt_p(ft$p.value))

  ggplot(d) +
    geom_col(aes(group, 1), width = 0.62, fill = COL_ABS) +
    geom_col(aes(group, p_det, fill = group), width = 0.62) +
    annotate("text", x = 1.5, y = 1.045, label = lab, size = 2.6,
             colour = "grey20", vjust = 0) +
    scale_fill_manual(values = COL, guide = "none") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       breaks = c(0, 0.5, 1), limits = c(0, 1.10),
                       expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(labels = xlab2) +
    labs(x = NULL, y = "Samples (%)", title = title, subtitle = subtitle) +
    base_th
}

s <- read.delim(SAMP, stringsAsFactors = FALSE)
s$group <- factor(s$group, levels = c("Disease", "HC"))
s$detected <- s$n_detected
s$n <- s$n_samples
p_c <- bar_panel(s, "Detection of covered ASVs",
                 "Any ASV, ≥97% HOMD match")

# cluster 21907, computed exactly as in 17b_plot_cluster_prevalence.R
rel <- read.delim(REL, sep = "\t", row.names = 1, check.names = FALSE)
meta <- read.delim(META, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(sample_key = "character"))
st <- meta$status_cleaned
di <- which(st == "OSCC")
hi <- which(st == "HC" & meta$bioproject_id %in% meta$bioproject_id[di])
v <- rel["c21907", ]
d4 <- data.frame(
  group = factor(c("Disease", "HC"), levels = c("Disease", "HC")),
  n = c(length(di), length(hi)),
  detected = c(sum(v[di] > 0), sum(v[hi] > 0)),
  stringsAsFactors = FALSE)
p_d <- bar_panel(d4, "Detection of phylotype 21907",
                 "Best HOMD match <97% (uncaptured)")

# ---- assemble --------------------------------------------------------------
p <- p_a + p_b + p_c + p_d + plot_layout(nrow = 1) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

out <- file.path(FIG, "fig_oscc_capnocytophaga_covered_prevalence")
ggsave(paste0(out, ".png"), p, width = 10.6, height = 3.9, dpi = 300, bg = "white")
ggsave(paste0(out, ".pdf"), p, width = 10.6, height = 3.9, device = cairo_pdf,
       bg = "white")

message(sprintf("[17d] covered: disease %d/%d, HC %d/%d; 21907: disease %d/%d, HC %d/%d",
                s$n_detected[1], s$n_samples[1], s$n_detected[2], s$n_samples[2],
                d4$detected[1], d4$n[1], d4$detected[2], d4$n[2]))
message(sprintf("[17d] relative abundance, detected samples: %s -> %s.png | .pdf",
                p_lab, out))
