#!/usr/bin/env Rscript
# Prevalence of one supported cluster, disease arm versus HC, as two stacked
# bars.  Each bar is the cohort's samples split into those in which the cluster
# was detected (any reads at all) and those in which it was not; the two
# segments of a bar sum to 100%.
#
# The cohort is the one used by 16_cluster_differential_wilcox.R and
# 17_maaslin3_disease.R: every sample of the disease status, plus the HC
# samples of the same bioprojects.  It is the full cohort -- prevalence needs no
# covariates -- so the bars describe the samples as they were collected.
#
# Usage: 17b_plot_cluster_prevalence.R [cluster_id] [disease_status]
#        defaults to 21907 (Capnocytophaga) in OSCC
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIG  <- file.path(PIPE, "analysis_report", "figs")
REL  <- file.path(PIPE, "14_cluster_relabund_dedup.tsv")
META <- file.path(PIPE, "14_cluster_sample_metadata.tsv")
ANN  <- file.path(PIPE, "14_supported_cluster_annotation.tsv")

a <- commandArgs(trailingOnly = TRUE)
CLUSTER <- if (length(a) >= 1) a[1] else "c21907"
STATUS  <- if (length(a) >= 2) a[2] else "OSCC"

COL      <- c(Disease = "#B2182B", HC = "#4393C3")
COL_ABS  <- "#D9D9D9"
COL_TXT  <- "#1A1A1A"

rel <- as.matrix(read.delim(REL, sep = "\t", row.names = 1, check.names = FALSE))
storage.mode(rel) <- "double"
meta <- read.delim(META, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(sample_key = "character"))
stopifnot(identical(colnames(rel), meta$sample_key))
stopifnot(CLUSTER %in% rownames(rel))

st <- meta$status_cleaned
di <- which(st == STATUS)
hi <- which(st == "HC" & meta$bioproject_id %in% meta$bioproject_id[di])
stopifnot(length(di) > 0, length(hi) > 0)

v <- rel[CLUSTER, ]
d <- data.frame(
  group    = factor(c("Disease", "HC"), levels = c("Disease", "HC")),
  n        = c(length(di), length(hi)),
  detected = c(sum(v[di] > 0), sum(v[hi] > 0)),
  stringsAsFactors = FALSE)
d$p_det <- d$detected / d$n
d$absent <- d$n - d$detected
d$lab_det <- sprintf("%s/%s (%.1f%%)", comma(d$detected), comma(d$n),
                     100 * d$p_det)
d$lab_abs <- sprintf("%s/%s (%.1f%%)", comma(d$absent), comma(d$n),
                     100 * (1 - d$p_det))

ann <- read.delim(ANN, sep = "\t", stringsAsFactors = FALSE,
                  colClasses = c(cluster_id = "character"))
nm <- ann$annotation[ann$cluster_id == CLUSTER]
if (!length(nm) || is.na(nm) || !nzchar(nm)) nm <- "unclassified"

# The HC segment is often a couple of percent tall, far too little to hold a
# label, so the detected count is written above the bar in the group colour and
# the absent count inside the grey segment, which is always large.
lab_det <- data.frame(group = d$group,
                      y = pmax(d$p_det, 0.02) + 0.055,
                      lab = d$lab_det, col = unname(COL[as.character(d$group)]))
lab_abs <- data.frame(group = d$group,
                      y = d$p_det + (1 - d$p_det) / 2,
                      lab = d$lab_abs)

p <- ggplot(d) +
  geom_col(aes(group, 1), width = 0.62, fill = COL_ABS) +
  geom_col(aes(group, p_det, fill = group), width = 0.62) +
  geom_text(data = lab_det, aes(group, y, label = lab, colour = col),
            size = 3.4, fontface = "bold", vjust = 0) +
  geom_text(data = lab_abs, aes(group, y, label = lab),
            size = 3.2, colour = "grey35") +
  scale_fill_manual(values = COL, guide = "none") +
  scale_colour_identity() +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1), expand = expansion(mult = c(0, 0.06))) +
  scale_x_discrete(labels = c(Disease = sprintf("%s\n(n = %s)", STATUS, comma(d$n[1])),
                               HC = sprintf("HC\n(n = %s)", comma(d$n[2])))) +
  labs(x = NULL, y = "Samples (%)",
       title = sprintf("%s (cluster %s)", nm, CLUSTER)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 9.5, colour = COL_TXT,
                                   lineheight = 0.95),
        axis.text.y = element_text(size = 9, colour = COL_TXT),
        axis.title.y = element_text(size = 10, colour = COL_TXT),
        plot.title = element_text(size = 11, face = "bold", colour = COL_TXT),
        plot.subtitle = element_text(size = 8, colour = "grey35",
                                     margin = margin(b = 6)),
        plot.margin = margin(6, 8, 6, 6))

W <- 3.5; H <- 4.0
out <- file.path(FIG, sprintf("fig_%s_%s_prevalence", tolower(STATUS), CLUSTER))
ggsave(paste0(out, ".png"), p, width = W, height = H, dpi = 300, bg = "white")
ggsave(paste0(out, ".pdf"), p, width = W, height = H, device = cairo_pdf, bg = "white")
message(sprintf("[17b] %s (cluster %s) in %s: disease %s, HC %s -> %s.png | .pdf",
                nm, CLUSTER, STATUS, d$lab_det[1], d$lab_det[2], out))
