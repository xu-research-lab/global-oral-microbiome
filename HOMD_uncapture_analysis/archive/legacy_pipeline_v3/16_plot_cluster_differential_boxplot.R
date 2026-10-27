#!/usr/bin/env Rscript
# Grouped boxplots of relative abundance for the clusters that differ most
# consistently between disease and HC.
#
# The 15 cohorts are those of 16_cluster_differential_wilcox.R (the largest by
# sample count).  Only significant (FDR < 0.05) cohort x cluster pairs in which
# the cluster is higher in the disease arm are drawn, and only for clusters in
# which that happens in at least two cohorts -- nothing else appears in the
# figure.  The corresponding rows are written to
# 16_cluster_differential_wilcox_sig.tsv; the complete 486 x 15 result is in
# 16_cluster_differential_wilcox.tsv.
#
# The clusters are rare: in most samples their relative abundance is exactly
# zero, so a linear axis would collapse every box onto zero.  Abundance is
# therefore drawn as percent of reads on a log10 axis, with zeros placed at a
# floor of 0.0001% (the tick labels are the real percentages).  Significance is
# taken from the rank-sum test and Benjamini-Hochberg correction; a star is
# drawn above a pair when FDR < 0.05, red when the cluster is higher in the
# disease arm and blue when it is higher in HC.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIG  <- file.path(PIPE, "analysis_report", "figs")
TESTS <- file.path(PIPE, "16_cluster_differential_wilcox.tsv")
REL   <- file.path(PIPE, "14_cluster_relabund_dedup.tsv")
META  <- file.path(PIPE, "14_cluster_sample_metadata.tsv")
ANN   <- file.path(PIPE, "14_supported_cluster_annotation.tsv")

COL   <- c("Disease" = "#B2182B", "HC" = "#4393C3")
COL_TXT <- "#1A1A1A"
FLOOR <- 1e-4          # percent; zeros are drawn here

tests <- read.delim(TESTS, sep = "\t", stringsAsFactors = FALSE)
tests$cluster_id <- as.character(tests$cluster_id)

# Direction is read from the change in detection rate: with so many zeros the
# mean-based fold change follows a few high values, while the rank test and the
# boxes follow how often the cluster is present at all.
tests$dir <- ifelse(tests$prev_disease > tests$prev_hc, "Disease", "HC")

# only the disease-higher pairs, and only clusters that repeat across cohorts
keep <- tests %>%
  filter(p_adj < 0.05, dir == "Disease") %>%
  count(cluster_id, name = "n_sig") %>%
  filter(n_sig >= 2) %>%
  arrange(desc(n_sig), cluster_id)
stopifnot(nrow(keep) > 0)

ann <- read.delim(ANN, sep = "\t", stringsAsFactors = FALSE)
ann$cluster_id <- as.character(ann$cluster_id)
ann_map <- setNames(ann$annotation, ann$cluster_id)
nm <- ann_map[keep$cluster_id]
nm[is.na(nm) | !nzchar(nm)] <- "unclassified"
# genus (or unclassified) on one line and the cluster id on the next, so a bar
# can be traced back to the table
lab <- setNames(sprintf("%s\n(%s)", nm, keep$cluster_id), keep$cluster_id)

meta <- read.delim(META, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(sample_key = "character"))
rel <- as.matrix(read.delim(REL, sep = "\t", row.names = 1, check.names = FALSE))
storage.mode(rel) <- "double"
stopifnot(identical(colnames(rel), meta$sample_key))

st <- meta$status_cleaned
cohorts <- tests %>%
  distinct(cohort, n_disease, n_hc) %>%
  arrange(desc(n_disease + n_hc))
cohorts$lab <- sprintf("%s\n%d vs %d", cohorts$cohort,
                       cohorts$n_disease, cohorts$n_hc)

# one row per sample per cluster, only for the cohorts and clusters in the panel
long <- vector("list", nrow(cohorts))
for (k in seq_len(nrow(cohorts))) {
  d <- cohorts$cohort[k]
  di <- which(st == d)
  hi <- which(st == "HC" & meta$bioproject_id %in% meta$bioproject_id[di])
  idx <- c(di, hi)
  grp <- c(rep("Disease", length(di)), rep("HC", length(hi)))
  vals <- rel[keep$cluster_id, idx, drop = FALSE]
  long[[k]] <- data.frame(
    cohort = d,
    cluster_id = rep(keep$cluster_id, times = length(idx)),
    cluster_pos = rep(seq_along(keep$cluster_id), times = length(idx)),
    group = rep(grp, each = length(keep$cluster_id)),
    y = FLOOR + 100 * as.vector(vals),
    stringsAsFactors = FALSE
  )
}
long <- bind_rows(long)

# keep only the pairs that passed correction
long$key <- paste(long$cohort, long$cluster_id, sep = "|")
tests$key <- paste(tests$cohort, tests$cluster_id, sep = "|")
sig_keys <- tests$key[tests$p_adj < 0.05 & tests$dir == "Disease"]
long <- long[long$key %in% sig_keys, ]

long$cluster_id <- factor(long$cluster_id, levels = keep$cluster_id)
long$group <- factor(long$group, levels = c("Disease", "HC"))
long$cohort <- factor(long$cohort, levels = cohorts$cohort)
long$cohort_lab <- factor(cohorts$lab[match(long$cohort, cohorts$cohort)],
                          levels = cohorts$lab)

# significance marks
sig <- tests %>%
  filter(cluster_id %in% keep$cluster_id, p_adj < 0.05, dir == "Disease") %>%
  mutate(mark = cut(p_adj, c(-Inf, 0.001, 0.01, 0.05), labels = c("***", "**", "*")),
         cohort_lab = cohorts$lab[match(cohort, cohorts$cohort)],
         cohort_lab = factor(cohort_lab, levels = cohorts$lab))

p <- ggplot(long, aes(x = cluster_id, y = y, fill = group, group = group)) +
  geom_boxplot(width = 0.72, orientation = "x", linewidth = 0.28,
               outlier.size = 0.22,
               outlier.alpha = 0.18, outlier.stroke = 0,
               position = position_dodge(width = 0.78)) +
  geom_text(data = sig, aes(x = cluster_id, y = Inf, label = mark, colour = dir),
            inherit.aes = FALSE, vjust = 1.25, size = 2.5, fontface = "bold") +
  facet_wrap(~cohort_lab, ncol = 3, scales = "free_x") +
  scale_fill_manual(values = COL, name = NULL, drop = FALSE) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_discrete(labels = lab, expand = expansion(add = 0.7)) +
  scale_y_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10),
                labels = c("0", "0.001", "0.01", "0.1", "1", "10"),
                expand = expansion(mult = c(0.02, 0.22))) +
  labs(x = NULL, y = "Relative abundance (% of reads, log scale)",
       title = "Clusters higher in disease than in HC",
       subtitle = sprintf(paste("Only the %d cohort x cluster pairs that are significant",
                                "(FDR < 0.05) and higher in the disease arm; clusters",
                                "repeated in at least two cohorts. Zeros at 0.0001%% floor",
                                sep = "\n"), nrow(sig))) +
  theme_bw(base_size = 7.4) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#EDEDED", linewidth = 0.25),
    panel.border = element_rect(colour = "grey35", linewidth = 0.35),
    strip.background = element_rect(fill = "#F2F2F2", colour = "grey35",
                                    linewidth = 0.3),
    strip.text = element_text(size = 6.1, colour = COL_TXT, lineheight = 0.95),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.2,
                               lineheight = 0.9, colour = COL_TXT),
    axis.text.y = element_text(size = 5.8, colour = COL_TXT),
    axis.title = element_text(size = 7, colour = COL_TXT),
    plot.title = element_text(size = 8.6, face = "bold", colour = COL_TXT),
    plot.subtitle = element_text(size = 6.6, colour = "grey35",
                                 margin = margin(b = 4)),
    legend.position = "top",
    legend.key.size = unit(0.26, "cm"),
    legend.text = element_text(size = 6.8),
    legend.margin = margin(b = -4),
    panel.spacing.x = unit(0.10, "in"),
    panel.spacing.y = unit(0.16, "in"),
    plot.margin = margin(4, 6, 4, 4)
  )

W <- 7.4
# facets with no significant pair are dropped, so the height follows what is left
N_FACET <- n_distinct(long$cohort)
H <- 1.95 * ceiling(N_FACET / 3) + 2.0
# the significant rows on their own
sigout <- file.path(PIPE, "16_cluster_differential_wilcox_sig.tsv")
write.table(tests[tests$p_adj < 0.05, setdiff(names(tests), "key")], sigout,
            sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("[16] significant rows written: %s (%d rows)",
                sigout, sum(tests$p_adj < 0.05)))

out <- file.path(FIG, "fig_cluster_differential_boxplot")
ggsave(paste0(out, ".png"), p, width = W, height = H, dpi = 300, bg = "white")
ggsave(paste0(out, ".pdf"), p, width = W, height = H, device = cairo_pdf, bg = "white")
message(sprintf("[16] wrote %s.png | .pdf (%.1f x %.1f in); %d clusters, %d sig marks, %d facets",
                out, W, H, nrow(keep), nrow(sig), N_FACET))
