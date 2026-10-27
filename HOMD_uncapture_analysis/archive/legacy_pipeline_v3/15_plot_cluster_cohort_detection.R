#!/usr/bin/env Rscript
# Detection of the 486 supported clusters across disease cohorts.
#
# One stacked bar per cluster, one facet per disease cohort with >100 subjects.
# The three segments are the samples of the cohort in which the cluster was
# detected in the disease arm, detected in the HC arm, or not detected at all,
# as fractions of the cohort's samples.  Clusters are held in one global order,
# by their pooled detection rate in the disease arms, so the same column is the
# same cluster in every facet.
#
# ICU has no HC samples and is kept (its HC segment is simply empty); the
# catch-all status 'Other' is not a disease and is dropped.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIG  <- file.path(PIPE, "analysis_report", "figs")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

d <- read.delim(file.path(PIPE, "15_cluster_cohort_detection.tsv"),
                sep = "\t", stringsAsFactors = FALSE)
d <- d[d$cohort != "Other", ]

COL <- c("Disease" = "#B2182B", "HC" = "#4393C3", "Not detected" = "#C9C9C9")
COL_TXT <- "#1A1A1A"

# ---- segments ---------------------------------------------------------------
# The three counts partition the cohort's samples (disease + HC), so each
# segment is drawn as its share of that common total and the bar sums to 100%.
# The per-arm rates are in 15_cluster_cohort_detection.tsv (rate_disease,
# rate_hc); here the share also reflects how many samples each arm contributes.
stopifnot(abs(d$share_disease + d$share_hc + d$share_absent - 1) < 1e-6)
d$s_disease <- d$share_disease
d$s_hc      <- d$share_hc
d$s_absent  <- d$share_absent

# ---- ordering ---------------------------------------------------------------
# pooled share of disease-positive samples, so the bars descend in the red
# segment.  The ranking is the same in every facet.
pool <- d %>%
  group_by(cluster_id) %>%
  summarise(rate = sum(n_pos_disease) / sum(n_total), .groups = "drop") %>%
  arrange(desc(rate), cluster_id)
# highest pooled disease share at the left edge
d$cluster_id <- factor(d$cluster_id, levels = pool$cluster_id)

# facets from the largest cohort down
coh <- d %>% distinct(cohort, n_subjects_disease) %>%
  arrange(desc(n_subjects_disease))
d$cohort <- factor(d$cohort, levels = coh$cohort)

# strip label: cohort name, then the samples the bar percentages are taken over,
# split into the two arms (the segments are shares of exactly these samples)
coh_n <- d %>% distinct(cohort, n_samples_disease, n_samples_hc)
COH_LAB <- setNames(
  sprintf("%s\nn = %s (%s disease / %s HC)", coh_n$cohort,
          comma(coh_n$n_samples_disease + coh_n$n_samples_hc),
          comma(coh_n$n_samples_disease), comma(coh_n$n_samples_hc)),
  as.character(coh_n$cohort))

# ---- long segments ----------------------------------------------------------
long <- bind_rows(
  data.frame(d[, c("cohort", "cluster_id")],
             segment = "Disease", prop = d$s_disease),
  data.frame(d[, c("cohort", "cluster_id")],
             segment = "HC", prop = d$s_hc),
  data.frame(d[, c("cohort", "cluster_id")],
             segment = "Not detected", prop = d$s_absent)
)
# the first level of a stacked bar sits at the top, so the order is written
# bottom-up as disease, HC, absent
long$segment <- factor(long$segment,
                       levels = c("Not detected", "HC", "Disease"))

p <- ggplot(long, aes(x = cluster_id, y = prop, fill = segment)) +
  geom_col(width = 1, position = position_stack()) +
  facet_grid(cohort ~ ., switch = "y", labeller = as_labeller(COH_LAB)) +
  scale_fill_manual(values = COL, breaks = c("Disease", "HC", "Not detected"),
                    name = NULL, drop = FALSE) +
  # the shares are written at 8 significant digits, so a bar's three segments
  # can sum to a hair over 1; squishing rather than dropping those rows keeps
  # the top of the grey segment in every bar
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, 0.5, 1),
                     limits = c(0, 1.0001), oob = squish,
                     expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(add = 0)) +
  labs(x = "Supported cluster (n = 486, ordered by detection in the disease arms)",
       y = "Samples (%)") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.35),
    panel.spacing.y = unit(0.11, "in"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6.2, colour = COL_TXT),
    axis.title = element_text(size = 7.4, colour = COL_TXT),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 0, size = 6.2,
                                     lineheight = 0.95,
                                     colour = COL_TXT, face = "bold"),
    strip.background = element_blank(),
    legend.position = "top",
    legend.key.size = unit(0.28, "cm"),
    legend.text = element_text(size = 6.8),
    legend.margin = margin(b = -2),
    plot.margin = margin(3, 6, 3, 3)
  )

W <- 7.4
H <- 0.40 * nlevels(d$cohort) + 1.10
out <- file.path(FIG, "fig_cluster_disease_cohort_prevalence")
ggsave(paste0(out, ".png"), p, width = W, height = H, dpi = 300, bg = "white")
ggsave(paste0(out, ".pdf"), p, width = W, height = H, device = cairo_pdf, bg = "white")
message(sprintf("[15] %d cohorts, %d clusters; written %s.png | .pdf (%.1f x %.1f in)",
                nlevels(d$cohort), nlevels(d$cluster_id), out, W, H))

# ---- companion: the leading clusters, where the names fit -------------------
# 486 bars leave no room for names, so the same segments are drawn again for
# the 30 clusters with the largest pooled disease detection.  The order and the
# three segments are unchanged, so this is the left edge of the panel above,
# enlarged and labelled.
N_TOP <- 30
TOP <- as.character(head(pool$cluster_id, N_TOP))
ann <- read.delim(file.path(PIPE, "14_supported_cluster_annotation.tsv"),
                  sep = "\t", stringsAsFactors = FALSE)
ann_map <- setNames(ann$annotation, as.character(ann$cluster_id))
nm <- ann_map[TOP]
nm[is.na(nm) | !nzchar(nm)] <- "unclassified"
lvl <- sprintf("%s (%s)", nm, TOP)

lc <- long[as.character(long$cluster_id) %in% TOP, ]
lc$cluster_id <- factor(as.character(lc$cluster_id), levels = TOP)
lc$label <- factor(as.character(lc$cluster_id), levels = TOP, labels = lvl)

p_top <- ggplot(lc, aes(x = label, y = prop, fill = segment)) +
  geom_col(width = 0.88, position = position_stack()) +
  facet_grid(cohort ~ ., switch = "y", labeller = as_labeller(COH_LAB)) +
  scale_fill_manual(values = COL, breaks = c("Disease", "HC", "Not detected"),
                    name = NULL, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, 0.5, 1),
                     limits = c(0, 1.0001), oob = squish,
                     expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  labs(x = sprintf("Supported cluster (top %d, ordered by detection in the disease arms)",
                   N_TOP),
       y = "Samples (%)") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.35),
    panel.spacing.y = unit(0.11, "in"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.4,
                               colour = COL_TXT),
    axis.ticks.x = element_line(colour = "grey50", linewidth = 0.25),
    axis.text.y = element_text(size = 6.2, colour = COL_TXT),
    axis.title = element_text(size = 7.4, colour = COL_TXT),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 0, size = 6.2,
                                     lineheight = 0.95,
                                     colour = COL_TXT, face = "bold"),
    strip.background = element_blank(),
    legend.position = "top",
    legend.key.size = unit(0.28, "cm"),
    legend.text = element_text(size = 6.8),
    legend.margin = margin(b = -2),
    plot.margin = margin(3, 6, 3, 3)
  )

out2 <- file.path(FIG, "fig_cluster_disease_cohort_prevalence_top30")
ggsave(paste0(out2, ".png"), p_top, width = W, height = H, dpi = 300, bg = "white")
ggsave(paste0(out2, ".pdf"), p_top, width = W, height = H, device = cairo_pdf, bg = "white")
message(sprintf("[15] companion written %s.png | .pdf", out2))
