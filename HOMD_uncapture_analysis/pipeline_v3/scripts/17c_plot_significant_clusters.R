#!/usr/bin/env Rscript
# Figures for the clusters MaAsLin3 called significant, one disease cohort at a
# time.
#
# A cluster significant in the abundance model is drawn as grouped boxplots of
# relative abundance, disease versus HC, one box pair per cluster.  The abundance
# model is fitted to log relative abundance of the samples in which the cluster
# was detected, so the boxes are restricted to those samples; the detection count
# of each arm is written under the cluster name, and prevalence is what the other
# figure is for.
#
# A cluster significant in the prevalence model is drawn as two stacked bars,
# disease and HC, each summing to 100% and split into detected and not detected.
# One facet per cluster, ordered by prevalence in the disease arm.
#
# The samples are the ones the models were fitted to: the disease arm plus the
# HC samples of the same bioprojects, reduced to those complete on the covariates
# (age, sex, site, sequencing region, continent), exactly as 17_maaslin3_disease.R
# does.  The counts under the bars and in the boxplot axes are therefore the ones
# behind the reported p-values, not the raw cohort.
#
# Usage: 17c_plot_significant_clusters.R [cohort ...]
#        defaults to the cohorts that have any significant cluster
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(tidyr)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIG  <- file.path(PIPE, "analysis_report", "figs")
REL  <- file.path(PIPE, "14_cluster_relabund_dedup.tsv")
META <- file.path(PIPE, "14_cluster_sample_metadata.tsv")
RES  <- file.path(PIPE, "17_maaslin3_all_results.tsv")
ANN  <- file.path(PIPE, "14_supported_cluster_annotation.tsv")

COV <- c("age_cleaned", "sex_cleaned", "site_cleaned",
         "sequencing_region_cleaned", "continent_group")
COL     <- c(Disease = "#B2182B", HC = "#4393C3")
COL_ABS <- "#D9D9D9"
COL_TXT <- "#1A1A1A"
FLOOR   <- 1e-4        # percent; only reached if a detected sample rounds to 0

args <- commandArgs(trailingOnly = TRUE)
cohorts <- if (length(args)) args else
  c("Colorectal Cancer", "Head and Neck Cancer", "Lung Cancer",
    "Periodontal Disease")

meta <- read.delim(META, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(sample_key = "character"))
rel <- as.matrix(read.delim(REL, sep = "\t", row.names = 1, check.names = FALSE))
storage.mode(rel) <- "double"
stopifnot(identical(colnames(rel), meta$sample_key))

res <- read.delim(RES, sep = "\t", stringsAsFactors = FALSE,
                  colClasses = c(feature = "character"))
ann <- read.delim(ANN, sep = "\t", stringsAsFactors = FALSE,
                  colClasses = c(cluster_id = "character"))
ann_map <- setNames(ann$annotation, ann$cluster_id)

# ---- the samples each model was fitted to ----------------------------------
model_samples <- function(d) {
  di <- which(meta$status_cleaned == d)
  hi <- which(meta$status_cleaned == "HC" &
                meta$bioproject_id %in% meta$bioproject_id[di])
  md <- meta[c(di, hi), ]
  for (v in COV) md[[v]][md[[v]] == ""] <- NA
  md[complete.cases(md[, COV]), ]
}

clab <- function(id) {
  vapply(id, function(i) {
    nm <- ann_map[i]
    if (length(nm) == 0 || is.na(nm) || !nzchar(nm)) nm <- "unclassified"
    # a facet one sixth of the canvas wide cannot hold "Prochlorococcus
    # marinus_32507" on one line; wrapped, it takes three and stays readable
    sprintf("%s\n(%s)", paste(strwrap(nm, 19), collapse = "\n"), i)
  }, character(1), USE.NAMES = FALSE)
}

# ===========================================================================
# prevalence: two bars per cluster, detected / not detected
# ===========================================================================
plot_prevalence <- function(d, ids, md) {
  grp <- ifelse(md$status_cleaned == "HC", "HC", "Disease")
  s <- lapply(ids, function(id) {
    v <- rel[id, md$sample_key]
    data.frame(cluster = id,
               group = factor(c("Disease", "HC"), levels = c("Disease", "HC")),
               n = c(sum(grp == "Disease"), sum(grp == "HC")),
               det = c(sum(v[grp == "Disease"] > 0), sum(v[grp == "HC"] > 0)),
               stringsAsFactors = FALSE)
  })
  s <- bind_rows(s)
  s$p <- s$det / s$n
  s$lab <- sprintf("%s/%s (%.1f%%)", comma(s$det), comma(s$n), 100 * s$p)
  ord <- s %>% filter(group == "Disease") %>% arrange(desc(p)) %>% pull(cluster)
  s$cluster <- factor(s$cluster, levels = ord)
  s$facet <- factor(clab(as.character(s$cluster)), levels = clab(ord))

  # Labels sit above their bar so thin segments still carry their count.  In the
  # large cohorts the counts run to "1,306/3,398 (14.1%)", wider than the gap
  # between the two bars of a facet, so the HC label is lifted clear of the
  # disease one -- above it by a fixed amount, not to a fixed height, or a tall
  # disease bar would put the two back on the same line.
  s <- s %>% group_by(cluster) %>%
    mutate(p_d = p[group == "Disease"]) %>% ungroup()
  s$y <- ifelse(s$group == "Disease",
                s$p + 0.035,
                pmax(s$p + 0.035, s$p_d + 0.115))
  s$col <- unname(COL[as.character(s$group)])

  ncol <- min(6, length(ids))
  nrow <- ceiling(length(ids) / ncol)
  n_d <- sum(grp == "Disease"); n_h <- sum(grp == "HC")
  p <- ggplot(s) +
    geom_col(aes(group, 1), width = 0.68, fill = COL_ABS) +
    geom_col(aes(group, p, fill = group), width = 0.68) +
    geom_text(aes(group, y, label = lab, colour = col), size = 1.95,
              fontface = "bold", vjust = 0) +
    facet_wrap(~facet, ncol = ncol) +
    scale_fill_manual(values = COL, guide = "none") +
    scale_colour_identity() +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1), breaks = c(0, 0.5, 1),
                       expand = expansion(mult = c(0, 0.10))) +
    labs(x = NULL, y = "Samples (%)",
         title = sprintf("%s: prevalence", d),
         subtitle = paste0("Clusters with a prevalence difference\n",
                           sprintf("Complete-case samples: %s disease, %s HC\n",
                                   comma(n_d), comma(n_h)),
                           "Colour = cluster detected, grey = not detected")) +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(colour = "grey45", linewidth = 0.3),
          strip.background = element_rect(fill = "#F2F2F2", colour = "grey45",
                                          linewidth = 0.3),
          strip.text = element_text(size = 5.9, lineheight = 0.92,
                                    colour = COL_TXT),
          axis.text.x = element_text(size = 5.4, colour = COL_TXT),
          axis.text.y = element_text(size = 5.2, colour = "grey35"),
          axis.title.y = element_text(size = 6.4, colour = COL_TXT),
          plot.title = element_text(size = 8.4, face = "bold", colour = COL_TXT),
          plot.subtitle = element_text(size = 6.2, colour = "grey35",
                                       lineheight = 1.05, margin = margin(b = 5)),
          panel.spacing.x = unit(0.06, "in"),
          panel.spacing.y = unit(0.20, "in"),
          plot.margin = margin(5, 7, 5, 5))

  # the row has to be tall enough for the strip (three lines of a wrapped
  # species name) plus the pair of tick labels the next row's strip would
  # otherwise sit on
  W <- max(3.4, ncol * 1.18 + 0.45)
  H <- nrow * 1.32 + 1.15
  out <- file.path(FIG, sprintf("fig_%s_prevalence_sig",
                                gsub("[^A-Za-z0-9]+", "_", tolower(d))))
  ggsave(paste0(out, ".png"), p, width = W, height = H, dpi = 300, bg = "white")
  ggsave(paste0(out, ".pdf"), p, width = W, height = H, device = cairo_pdf,
         bg = "white")
  message(sprintf("[17c] %-24s prevalence: %d clusters -> %s.png | .pdf",
                  d, length(ids), out))
}

# ===========================================================================
# abundance: grouped boxplots among the samples where the cluster was detected
# ===========================================================================
plot_abundance <- function(d, ids, md) {
  grp <- ifelse(md$status_cleaned == "HC", "HC", "Disease")
  long <- lapply(ids, function(id) {
    v <- 100 * rel[id, md$sample_key]
    keep <- v > 0
    data.frame(cluster = id, group = grp[keep], y = v[keep],
               stringsAsFactors = FALSE)
  })
  long <- bind_rows(long)
  long$cluster <- factor(long$cluster, levels = ids)
  long$group <- factor(long$group, levels = c("Disease", "HC"))

  # detection counts per arm, which is the set of samples actually modelled
  cnt <- long %>% count(cluster, group) %>%
    tidyr::complete(cluster, group, fill = list(n = 0))
  lab <- vapply(ids, function(id) {
    nd <- cnt$n[cnt$cluster == id & cnt$group == "Disease"]
    nh <- cnt$n[cnt$cluster == id & cnt$group == "HC"]
    sprintf("%s\n%d vs %d detected", clab(id), nd, nh)
  }, character(1))
  names(lab) <- ids

  n_d <- sum(grp == "Disease"); n_h <- sum(grp == "HC")
  p <- ggplot(long, aes(cluster, y, fill = group)) +
    geom_boxplot(width = 0.62, linewidth = 0.3, outlier.shape = NA,
                 position = position_dodge(width = 0.7)) +
    # every sample is shown: 8 to 36 per arm is few enough that the box alone
    # would hide the spread it summarises
    geom_point(position = position_jitterdodge(jitter.width = 0.20, dodge.width = 0.7),
               size = 0.35, alpha = 0.45, stroke = 0, show.legend = FALSE) +
    scale_fill_manual(values = COL, name = NULL) +
    scale_x_discrete(labels = lab) +
    scale_y_log10(labels = label_number(accuracy = 1e-4, trim = TRUE),
                  expand = expansion(mult = c(0.04, 0.10))) +
    labs(x = NULL, y = "Relative abundance (% of reads, log scale)",
         title = sprintf("%s: abundance", d),
         subtitle = paste0("Clusters with an abundance difference\n",
                           sprintf("Detected samples only; complete-case set: %s disease, %s HC",
                                   comma(n_d), comma(n_h)))) +
    theme_bw(base_size = 8.5) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 6.2, colour = COL_TXT,
                                    lineheight = 0.95),
          axis.text.y = element_text(size = 6.4, colour = COL_TXT),
          axis.title.y = element_text(size = 7.2, colour = COL_TXT),
          plot.title = element_text(size = 8.6, face = "bold", colour = COL_TXT),
          plot.subtitle = element_text(size = 6.4, colour = "grey35",
                                       margin = margin(b = 5)),
          legend.position = "top",
          legend.key.size = unit(0.28, "cm"),
          legend.text = element_text(size = 6.6),
          legend.margin = margin(b = -3),
          plot.margin = margin(5, 8, 4, 4))

  W <- max(2.8, 0.95 * length(ids) + 1.1)
  H <- 3.5
  out <- file.path(FIG, sprintf("fig_%s_abundance_sig",
                                gsub("[^A-Za-z0-9]+", "_", tolower(d))))
  ggsave(paste0(out, ".png"), p, width = W, height = H, dpi = 300, bg = "white")
  ggsave(paste0(out, ".pdf"), p, width = W, height = H, device = cairo_pdf,
         bg = "white")
  message(sprintf("[17c] %-24s abundance: %d clusters -> %s.png | .pdf",
                  d, length(ids), out))
}

# ---------------------------------------------------------------------------
for (d in cohorts) {
  sig <- res %>% filter(cohort == d, !is.na(qval_individual),
                        qval_individual < 0.05)
  if (!nrow(sig)) {
    message(sprintf("[17c] %-24s no significant cluster", d)); next
  }
  # a cluster significant in both models is drawn once, as abundance: that is
  # the model that says something about how much of it there is
  ab <- sig %>% filter(model == "abundance") %>% arrange(qval_individual)
  pr <- sig %>% filter(model == "prevalence", !feature %in% ab$feature) %>%
    arrange(qval_individual)

  md <- model_samples(d)
  message(sprintf("[17c] %-24s %d disease / %d HC samples in the model",
                  d, sum(md$status_cleaned != "HC"),
                  sum(md$status_cleaned == "HC")))
  if (nrow(ab)) plot_abundance(d, ab$feature, md)
  if (nrow(pr)) plot_prevalence(d, pr$feature, md)
}
