#!/usr/bin/env Rscript
# Figure for the cluster-21907 cohort sensitivity analysis.
#
# Two rows, one per cohort definition, two columns, one per analysis:
#
#   left    abundance among samples in which the cluster is detected (boxplot)
#   right   prevalence over every sample of the cohort (stacked bars)
#
# The MaAsLin3 coefficient for the disease term is written into each panel, so
# the raw bars and the adjusted test can be read against each other -- which is
# the point of the figure: in the first row the bars differ sharply and the
# study-adjusted test is null.
#
# Run scripts/17e_21907_cohort_sensitivity_fit.R first.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIG  <- file.path(PIPE, "analysis_report", "figs")
SMP  <- file.path(PIPE, "17e_21907_cohort_sensitivity_samples.tsv")
RES  <- file.path(PIPE, "17e_21907_cohort_sensitivity_results.tsv")

COL     <- c(Disease = "#B2182B", HC = "#4393C3")
COL_ABS <- "#D9D9D9"
COL_TXT <- "#1A1A1A"
XL      <- c(Disease = "OSCC", HC = "HC")

smp <- read.delim(SMP, stringsAsFactors = FALSE)
res <- read.delim(RES, stringsAsFactors = FALSE)
smp$group <- factor(smp$group, levels = c("Disease", "HC"))
smp$scenario <- factor(smp$scenario, levels = c("S1", "S2"))

SCEN_LAB <- c(
  S1 = "All cohorts carrying the cluster",
  S2 = "Cohorts with both OSCC and HC")

# MaAsLin3 line for one scenario x model
stat_lab <- function(k, model) {
  r <- res[res$scenario == k & res$model == model, ]
  if (!nrow(r)) return("MaAsLin3: not fitted")
  # the confidence interval carries the precision of the estimate: in row 1 the
  # abundance coefficient is large but its interval spans zero
  sprintf("MaAsLin3: β = %.2f (95%% CI %.2f to %.2f), p = %s",
          r$coef, r$coef - 1.96 * r$stderr, r$coef + 1.96 * r$stderr,
          format.pval(r$pval_individual, digits = 2))
}

n_by <- function(d) {
  t <- table(d$group)
  c(Disease = sum(d$group == "Disease"), HC = sum(d$group == "HC"))
}

# ---- abundance: samples where the cluster is detected ----------------------
p_abund <- function(k) {
  d <- smp[smp$scenario == k & smp$detected & smp$rel_abund > 0, ]
  n <- n_by(d)
  xl <- c(Disease = sprintf("%s\n(n = %s)", XL["Disease"], comma(n["Disease"])),
          HC = sprintf("%s\n(n = %s)", XL["HC"], comma(n["HC"])))
  ggplot(d, aes(group, rel_abund, fill = group)) +
    geom_boxplot(width = 0.45, outlier.shape = NA, linewidth = 0.3,
                 colour = "grey25", alpha = 0.55) +
    geom_jitter(width = 0.12, height = 0, size = 0.7, alpha = 0.35,
                colour = "#2B2B2B", shape = 16) +
    # one range for both rows: the two abundance panels are read against each
    # other, and a shared scale is what makes that possible
    scale_y_log10(labels = label_number(scale_cut = cut_short_scale()),
                  limits = c(1.5e-5, 3e-2), oob = scales::squish,
                  breaks = c(1e-5, 1e-4, 1e-3, 1e-2),
                  name = "Relative abundance (log10)") +
    scale_x_discrete(labels = xl) +
    scale_fill_manual(values = COL, guide = "none") +
    labs(title = "Abundance (detected samples only)",
         subtitle = stat_lab(k, "abundance"), x = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(colour = COL_TXT),
          axis.text.x = element_text(size = 9),
          axis.title = element_text(size = 9.5, colour = COL_TXT),
          plot.title = element_text(size = 10, face = "bold", colour = COL_TXT),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          plot.margin = margin(4, 6, 4, 4))
}

# ---- prevalence: every sample in the cohort --------------------------------
p_prev <- function(k) {
  d <- smp[smp$scenario == k, ]
  s <- d %>% group_by(group) %>%
    summarise(n = n(), det = sum(detected), .groups = "drop") %>%
    mutate(p = det / n, absent = n - det,
           lab = sprintf("%.1f%%\n%s/%s", 100 * p, comma(det), comma(n)))
  s$group <- factor(s$group, levels = c("Disease", "HC"))
  ggplot(s) +
    geom_col(aes(group, 1), width = 0.55, fill = COL_ABS) +
    geom_col(aes(group, p, fill = group), width = 0.55) +
    geom_text(aes(group, pmin(p + 0.05, 0.92), label = lab,
                  colour = group), size = 3, fontface = "bold",
              vjust = 0, lineheight = 0.9, show.legend = FALSE) +
    scale_colour_manual(values = COL, guide = "none") +
    scale_fill_manual(values = COL, guide = "none") +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.02)),
                       name = "Samples with the cluster (%)") +
    scale_x_discrete(labels = c(Disease = "OSCC", HC = "HC")) +
    labs(title = "Prevalence (all samples)",
         subtitle = stat_lab(k, "prevalence"), x = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(colour = COL_TXT),
          axis.text.x = element_text(size = 9),
          axis.title = element_text(size = 9.5, colour = COL_TXT),
          plot.title = element_text(size = 10, face = "bold", colour = COL_TXT),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          plot.margin = margin(4, 6, 4, 4))
}

# ---- assembly ---------------------------------------------------------------
row <- function(k) {
  n <- length(unique(smp$bioproject_id[smp$scenario == k]))
  n_both <- length(unique(smp$bioproject_id[smp$scenario == k &
                                              smp$group == "HC"]))
  ttl <- if (k == "S1")
    sprintf("%s  (%d projects, %d with controls)", SCEN_LAB[[k]], n, n_both) else
    sprintf("%s  (%d projects)", SCEN_LAB[[k]], n)
  (p_abund(k) | p_prev(k)) +
    plot_annotation(title = ttl,
                    theme = theme(plot.title = element_text(
                      size = 10.5, face = "bold", colour = COL_TXT,
                      margin = margin(2, 0, 2, 2))))
}

# wrap_elements: a row's own title is dropped when a patchwork carrying it is
# nested into another one, so each row is wrapped as a single element first
fig <- wrap_elements(row("S1")) / wrap_elements(row("S2")) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

out <- file.path(FIG, "fig_oscc_21907_cohort_sensitivity")
ggsave(paste0(out, ".png"), fig, width = 7.4, height = 6.6, dpi = 300,
       bg = "white")
ggsave(paste0(out, ".pdf"), fig, width = 7.4, height = 6.6,
       device = cairo_pdf, bg = "white")
message("[17e] figure written: ", out, ".png | .pdf")
