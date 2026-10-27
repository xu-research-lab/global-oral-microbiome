#!/usr/bin/env Rscript
# Figures for the restricted-cohort re-fit of the head-and-neck and periodontal
# significant clusters.
#
# One figure per disease, two panels side by side on a shared row order:
#   left   abundance among the samples in which the cluster was detected
#   right  prevalence over every sample of the restricted cohort
# Rows are ordered by the original q, so the clusters that were strongest in the
# full cohort sit at the top and can be followed across both panels.  The
# Benjamini-Hochberg q of the restricted fit is written at the right of each
# panel; it is recomputed here across every cluster carried into the same model
# of the same disease, because 17f_ctrl_only_fit.R splits the clusters into one
# MaAsLin3 call per set of detecting bioprojects and its own correction is then
# made within a fragment of the family.
#
# Run scripts/17f_ctrl_only_fit.R first.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIG  <- file.path(PIPE, "analysis_report", "figs")
REL  <- file.path(PIPE, "14_cluster_relabund_dedup.tsv")
META <- file.path(PIPE, "14_cluster_sample_metadata.tsv")
RES0 <- file.path(PIPE, "17_maaslin3_all_results.tsv")
FIT  <- file.path(PIPE, "17f_ctrl_only_results.tsv")
ANN  <- file.path(PIPE, "14_supported_cluster_annotation.tsv")

COL     <- c(Disease = "#B2182B", HC = "#4393C3")
COL_TXT <- "#1A1A1A"

# Left out of the figure although they pass the correction.  Both sit just
# under q = 0.05 (c21907 0.032, s6635 0.011) while their pooled prevalence bars
# differ by less than two percentage points, because the model's within-project
# contrast is carried by one study and diluted by the other: s6635 is 46.2 % vs
# 30.9 % in PRJNA901226 but 7/499 vs 7/495 in PRJEB37445.  They stay in
# 17f_ctrl_only_results.tsv and in the text; they are not drawn, and they are
# kept out of the correction family so the q values belong to the rows shown.
EXCLUDE <- list(`Head and Neck Cancer` = c("c21907", "s6635"))

rel <- as.matrix(read.delim(REL, sep = "\t", row.names = 1, check.names = FALSE))
storage.mode(rel) <- "double"
meta <- read.delim(META, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(sample_key = "character"))
res0 <- read.delim(RES0, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(feature = "character"))
fit  <- read.delim(FIT, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(feature = "character"))
ann  <- read.delim(ANN, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(cluster_id = "character"))
ann_map <- setNames(ann$annotation, ann$cluster_id)

fmt_p <- function(p) if (is.na(p)) "n.a." else
  if (p < 0.001) sprintf("%.1e", p) else sprintf("%.3f", p)

lab_of <- function(id) {
  nm <- ann_map[[id]]
  if (is.null(nm) || is.na(nm) || !nzchar(nm)) nm <- "unclassified"
  sprintf("%s (%s)", nm, id)
}

plot_disease <- function(d) {
  # rows: clusters significant in the original cohort, strongest first, minus
  # (a) the clusters that occur in a single bioproject -- the question is which
  # clusters separate cases from controls across data sets, and a cluster found
  # in one study has no across-study contrast to contribute -- and (b) the ones
  # that do not survive the correction in either model.  A row is kept when the
  # abundance model or the prevalence model still reaches q < 0.05 on the
  # cohorts that have controls; a cluster that does neither has nothing left to
  # show and would only fill the panel with grey.
  ord <- res0 %>% filter(cohort == d, !is.na(qval_individual),
                         qval_individual < 0.05) %>%
    group_by(feature) %>% summarise(q = min(qval_individual), .groups = "drop") %>%
    arrange(q)
  ids <- ord$feature
  pa0 <- fit[fit$disease == d & fit$model == "abundance", ]
  pp0 <- fit[fit$disease == d & fit$model == "prevalence", ]
  n_det <- pp0$n_projects_detected[match(ids, pp0$feature)]
  single <- ids[is.na(n_det) | n_det < 2]
  ids <- ids[!(ids %in% single)]
  if (length(single)) {
    message(sprintf("[17f] %s: dropping %d cluster(s) detected in one bioproject only",
                    d, length(single)))
  }
  n_multi <- length(ids)                 # found in two or more bioprojects
  excl <- EXCLUDE[[d]]
  if (length(excl)) {
    gone <- intersect(ids, excl)
    ids  <- setdiff(ids, excl)
    message(sprintf("[17f] %s: leaving out %d cluster(s) by hand: %s",
                    d, length(gone), paste(gone, collapse = ", ")))
  }
  # BH over every cluster carried into the same model of this disease, not over
  # the per-detecting-set MaAsLin3 calls: the sets fragment the family (24 of
  # the HNC clusters were corrected in a call of their own) while the panels are
  # read as one set of tests.
  q_of <- function(p) {
    p <- as.numeric(p); q <- rep(NA_real_, length(p))
    k <- !is.na(p); q[k] <- p.adjust(p[k], "BH"); q
  }
  q_ab <- setNames(q_of(pa0$pval_individual[match(ids, pa0$feature)]), ids)
  q_pr <- setNames(q_of(pp0$pval_individual[match(ids, pp0$feature)]), ids)
  sig <- (!is.na(q_ab) & q_ab < 0.05) | (!is.na(q_pr) & q_pr < 0.05)
  dropped <- ids[!sig]
  ids <- ids[sig]
  if (length(dropped)) {
    message(sprintf("[17f] %s: dropping %d cluster(s) with q >= 0.05 in both models: %s",
                    d, length(dropped), paste(dropped, collapse = ", ")))
  }
  rows <- rev(ids)                       # ggplot puts the first level at the bottom
  labs <- setNames(vapply(rows, lab_of, character(1)), rows)

  # ---- abundance: detected samples only ----
  smp <- read.delim(file.path(PIPE, "17f_ctrl_only_samples.tsv"),
                    stringsAsFactors = FALSE)
  smp <- smp[smp$disease == d, ]
  md <- meta[match(smp$sample_key, meta$sample_key), ]
  grp <- smp$condition

  long <- bind_rows(lapply(ids, function(id) {
    v <- 100 * rel[id, smp$sample_key]
    k <- v > 0
    data.frame(cluster = id, group = grp[k], y = v[k])
  }))
  long$cluster <- factor(long$cluster, levels = rows)
  long$group <- factor(long$group, levels = c("Disease", "HC"))

  pa_lab <- data.frame(cluster = factor(rows, levels = rows),
                       lab = paste0("q = ", vapply(q_ab[rows], fmt_p, "")))

  # the axis is set so that the data occupy its first 80 %: the last empty fifth
  # is the column the p values go in, and ending it there keeps the highest
  # tick inside the panel instead of on its border, where it ran into the left
  # edge of the prevalence panel
  lo <- min(long$y, na.rm = TRUE) * 0.75
  hi <- lo * 10^(log10(max(long$y, na.rm = TRUE) / lo) / 0.8)
  brk <- 10^(-6:6); brk <- brk[brk >= lo & brk <= hi]

  n_d <- sum(grp == "Disease"); n_h <- sum(grp == "HC")
  p_ab <- ggplot(long, aes(y = cluster, x = y, fill = group)) +
    geom_boxplot(width = 0.62, linewidth = 0.25, outlier.shape = NA,
                 position = position_dodge(width = 0.72),
                 orientation = "y") +
    geom_point(position = position_jitterdodge(jitter.height = 0.16,
                                               dodge.width = 0.72),
               size = 0.3, alpha = 0.45, stroke = 0, show.legend = FALSE) +
    # x = Inf puts the p at the panel's right edge whatever the data range is,
    # and hjust > 1 holds it just inside the border
    geom_text(data = pa_lab, aes(y = cluster, x = Inf, label = lab),
              inherit.aes = FALSE, hjust = 1.1, size = 2.1, colour = "grey30") +
    scale_fill_manual(values = COL, name = NULL) +
    # the right-hand fifth of the axis is left empty for the p column; on a log
    # axis a multiplicative expansion would add decades, so the limits are set
    # on the data instead
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale()),
                  limits = c(lo, hi), breaks = brk,
                  expand = expansion(mult = 0)) +
    scale_y_discrete(labels = labs) +
    labs(x = "Relative abundance (% of reads)", y = NULL,
         title = "Abundance",
         subtitle = "Samples in which the cluster was detected\nMaAsLin3 q (BH, restricted cohort)") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 6.4, colour = COL_TXT),
          axis.text.x = element_text(size = 7, colour = COL_TXT),
          axis.title.x = element_text(size = 7.6, colour = COL_TXT),
          plot.title = element_text(size = 9.5, face = "bold", colour = COL_TXT),
          plot.subtitle = element_text(size = 7, colour = "grey35",
                                       lineheight = 1.05, margin = margin(b = 4)),
          legend.position = "top",
          legend.key.size = unit(0.26, "cm"),
          legend.text = element_text(size = 7),
          legend.margin = margin(b = -4))

  # ---- prevalence: the samples the model was fitted on ----
  # Each cluster is drawn on the same samples as its p value: the samples of the
  # bioprojects in which it was detected.  Averaged over the whole cohort, a
  # cluster found in one of five studies had its bars diluted by the four
  # studies where it never occurred, which is not where the contrast is made.
  pr_of <- setNames(pp0$projects, pp0$feature)
  prev <- bind_rows(lapply(ids, function(id) {
    prj <- strsplit(pr_of[[id]], ",", fixed = TRUE)[[1]]
    k <- smp$bioproject_id %in% prj
    v <- rel[id, smp$sample_key] > 0
    data.frame(cluster = id,
               group = c("Disease", "HC"),
               det = c(sum(v[k & grp == "Disease"]), sum(v[k & grp == "HC"])),
               n = c(sum(k & grp == "Disease"), sum(k & grp == "HC")))
  }))
  prev$p <- prev$det / prev$n
  prev$cluster <- factor(prev$cluster, levels = rows)
  prev$group <- factor(prev$group, levels = c("Disease", "HC"))

  # The axis stops short of 100 % and keeps the last quarter for the p values.
  # The limit follows the data so the longest bar cannot run into that strip:
  # the longest bar is put at 72 % of the axis, rounded up to the next 5 %, with
  # 40 % as the floor.
  XMAX <- max(40, 5 * ceiling(max(100 * prev$p) / 0.72 / 5))

  pp_lab <- data.frame(cluster = factor(rows, levels = rows),
                       lab = paste0("q = ", vapply(q_pr[rows], fmt_p, "")))

  p_pr <- ggplot(prev, aes(y = cluster, x = 100 * p, fill = group)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62,
             orientation = "y") +
    geom_text(data = pp_lab, aes(y = cluster, x = Inf, label = lab),
              inherit.aes = FALSE, hjust = 1.1, size = 2.1, colour = "grey30") +
    scale_fill_manual(values = COL, name = NULL) +
    scale_x_continuous(limits = c(0, XMAX), breaks = seq(0, XMAX, 10),
                       labels = function(v) paste0(v, "%"),
                       expand = expansion(mult = 0)) +
    scale_y_discrete(labels = labs) +
    labs(x = "Samples with the cluster (%)", y = NULL,
         title = "Prevalence",
         subtitle = "Samples of the detecting bioprojects\nMaAsLin3 q (BH, restricted cohort)") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 7, colour = COL_TXT),
          axis.title.x = element_text(size = 7.6, colour = COL_TXT),
          plot.title = element_text(size = 9.5, face = "bold", colour = COL_TXT),
          plot.subtitle = element_text(size = 7, colour = "grey35",
                                       lineheight = 1.05, margin = margin(b = 4)),
          legend.position = "top",
          legend.key.size = unit(0.26, "cm"),
          legend.text = element_text(size = 7),
          legend.margin = margin(b = -4))

  n_proj <- unique(fit$n_projects[fit$disease == d])
  fig <- p_ab | p_pr
  fig <- fig + plot_annotation(
    title = sprintf("%s: clusters re-fitted on the cohorts that have controls", d),
    subtitle = sprintf(paste0("%d bioprojects with both cases and controls; %s cases / %s HC.",
                              "\nOf the %d clusters significant in the full cohort, %d occur in two or",
                              " more bioprojects and %d are shown; each shown row reaches q < 0.05",
                              "\nin at least one model (BH across the %d shown, per model). Each row",
                              " is tested, and its bars drawn, on the bioprojects where it was",
                              "\ndetected."),
                       n_proj[1], comma(n_d), comma(n_h),
                       length(ord$feature), n_multi, length(ids), length(ids)),
    theme = theme(plot.title = element_text(size = 10, face = "bold",
                                            colour = COL_TXT),
                  plot.subtitle = element_text(size = 7.6, colour = "grey35")))

  # The pair was drawn 10 in wide; with the prevalence axis no longer running
  # to 100 % the two panels had more width than the data, so the canvas is
  # narrowed to the width the row labels and two half-width panels need.
  H <- max(3.5, 0.16 * length(ids) + 1.9)
  W <- 7.2
  out <- file.path(FIG, sprintf("fig_%s_ctrl_only",
                                gsub("[^A-Za-z0-9]+", "_", tolower(d))))
  ggsave(paste0(out, ".png"), fig, width = W, height = H, dpi = 300, bg = "white")
  ggsave(paste0(out, ".pdf"), fig, width = W, height = H, device = cairo_pdf,
         bg = "white")
  message(sprintf("[17f] %-22s %d clusters -> %s.png | .pdf", d, length(ids), out))
}

for (d in c("Head and Neck Cancer", "Periodontal Disease")) plot_disease(d)
