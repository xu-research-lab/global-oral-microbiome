#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# figure3_panels_abc_0910.R
#
# Renders Figure 3 panels (replacing the phylum stacked-bar supplementary
# figures with Nature/Cell-style designs):
#
#   A : site x phylum bubble matrix
#       x = sampling site, y = top-N phyla, dot size = mean relative
#       abundance, dot color = prevalence (% samples with abundance > 0).
#       Answers "which phyla define which oral niche".
#
#   B : age trajectory
#       x = age group (ordered), y = mean relative abundance, one line per
#       top phylum, faceted by sampling site. Shows the developmental
#       trajectory of the oral microbiome from infancy to old age.
#
#   C : continent x phylum bubble matrix
#       x = continent, y = top-N phyla, dot size = mean relative abundance,
#       dot color = prevalence. (Statistical annotation intentionally omitted
#       for now; add MaAsLin2 / ANCOM-BC later if needed.)
#
# plus a combined figure of A/B/C.
#
# Input paths follow figure1_render_panels_0910.R (metadata + phylum table).
# Outputs go to <script_dir>/results/  (i.e. age-continent-site/results).
#
# Run once on the server:
#     Rscript figure3_panels_abc_0910.R
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---- Locate this script's directory (works with Rscript and source()) -----
find_script_dir <- function() {
  if (interactive()) return(normalizePath(getwd()))
  args <- commandArgs(trailingOnly = FALSE)
  f <- args[grep("^--file=", args)]
  if (length(f) == 1) return(dirname(normalizePath(sub("^--file=", "", f))))
  normalizePath(getwd())
}
script_dir <- find_script_dir()

# ============================================================================
# CONFIG
# ============================================================================
# Production tree. Its input/ folder holds only biom / metadata / taxonomy,
# so the phylum table keeps its original server path. Absent this tree (local
# checkout) everything falls back to script-relative folders. Env overrides:
# FIG1_METADATA, FIG1_PHYLUM_TABLE, FIG3_PANEL_DIR.
SERVER_ROOT <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-map-2/04-patch-pipeline-map/03-patch-pipeline-map"
SERVER_OUT_ROOT <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-pipeline-map/output"
LOCAL_BASE <- normalizePath(file.path(script_dir, ".."))
on_server <- dir.exists(file.path(SERVER_ROOT, "input")) &&
             dir.exists(SERVER_OUT_ROOT)
BASE_DIR <- if (on_server) SERVER_ROOT else LOCAL_BASE

metadata_path <- Sys.getenv("FIG1_METADATA",
                            unset = file.path(BASE_DIR, "input",
                                              "metadata_oral_qc.tsv"))
abund_file <- Sys.getenv("FIG1_PHYLUM_TABLE",
                         unset = if (on_server)
  file.path(SERVER_OUT_ROOT, "phylum-relative-table",
            "phylum-relative-table.tsv")
else file.path(LOCAL_BASE, "output", "phylum-relative-table",
               "phylum-relative-table.tsv"))

# Outputs go to age-continent-site/results (this script's own results folder).
PANEL_DIR <- Sys.getenv("FIG3_PANEL_DIR", unset = file.path(script_dir, "results"))

# ---- Tunable parameters -----------------------------------------------------
TOP_N   <- 10L          # number of phyla shown in every panel
FIG_DPI <- 300L
BASE_SZ <- 13           # base font size (pt)

# ============================================================================
# Small helpers
# ============================================================================
# Metadata column resolution: prefer the *_cleaned column, fall back to the
# raw column only when no clean variant exists.
pick_meta_col <- function(md, preferred, fallback = NULL) {
  if (preferred %in% names(md)) return(preferred)
  if (!is.null(fallback) && fallback %in% names(md)) return(fallback)
  preferred
}

# Parent phylum ("p__..." segment) from a qiime2-collapsed label. Vectorised
# (runs over ~8M melted rows, so avoid per-element vapply). NOTE: regmatches()
# drops no-match positions, so use substr() + ifelse() to preserve length.
extract_parent_phylum <- function(label) {
  m  <- regexpr("p__[^;]+", label)
  ml <- attr(m, "match.length")
  out <- ifelse(m > 0, substr(label, m, m + ml - 1L), NA_character_)
  sub("^p__", "", out)
}

# Strip the GTDB numeric strain suffix: "Bacillota_A_368345" -> "Bacillota_A".
clean_phylum_name <- function(x) {
  sub("_[0-9]+$", "", x)
}

# Global top-N phyla by summed relative abundance across all samples.
global_top_phyla <- function(long_dt, n = 10) {
  tot <- long_dt[, .(Total = sum(Abundance, na.rm = TRUE)), by = Phylum]
  tot <- tot[!is.na(Phylum) & Phylum != ""][order(-Total)]
  as.character(head(tot$Phylum, n))
}

# Simple Title Case ("oral mucosa" -> "Oral Mucosa") without stringr.
to_title <- function(x) {
  vapply(strsplit(x, " "), function(parts) {
    paste0(toupper(substr(parts, 1, 1)), substr(parts, 2, nchar(parts)),
           collapse = " ")
  }, character(1))
}

# Categorical palette for the phylum lines (NPG / Nature-style, muted).
phylum_color_map <- function(top_phyla) {
  pal <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
           "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
           "#5F559B", "#A20056", "#808180", "#1B1919", "#E69F00")
  setNames(rep(pal, length.out = length(top_phyla)), top_phyla)
}

# Shared minimal theme (Nature/Cell aesthetic: no gridlines, no box fills).
theme_fig3 <- function() {
  theme_classic(base_size = BASE_SZ, base_family = "sans") +
    theme(
      axis.text         = element_text(color = "#303030"),
      axis.title        = element_text(color = "#303030", face = "bold"),
      axis.text.x       = element_text(angle = 30, hjust = 1, vjust = 1),
      legend.title      = element_text(color = "#303030", face = "bold",
                                      size = BASE_SZ - 3),
      legend.text       = element_text(color = "#303030", size = BASE_SZ - 4),
      legend.key        = element_blank(),
      strip.text        = element_text(color = "#303030", face = "bold"),
      strip.background  = element_blank(),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      panel.grid        = element_blank()
    )
}

# ============================================================================
# Data preparation
# ============================================================================
prepare_data <- function() {
  if (!file.exists(metadata_path)) stop("metadata not found: ", metadata_path)
  if (!file.exists(abund_file))    stop("phylum table not found: ", abund_file)

  md <- data.table::fread(metadata_path, header = TRUE, sep = "\t")
  col_site <- pick_meta_col(md, "site_cleaned",    "site")
  col_ageg <- pick_meta_col(md, "age_cleaned_group", "age_group")
  col_cont <- pick_meta_col(md, "continent_group", "continent")

  rel <- data.table::fread(abund_file, header = TRUE, sep = "\t", skip = 1)
  data.table::setnames(rel, 1, "Taxon")
  sample_cols <- setdiff(names(rel), "Taxon")
  rel_long <- data.table::melt(rel, id.vars = "Taxon", measure.vars = sample_cols,
                               variable.name = "sample_id", value.name = "Abundance")
  rel_long[, sample_id := as.character(sample_id)]
  rel_long[, Phylum := clean_phylum_name(extract_parent_phylum(Taxon))]
  rel_long <- rel_long[!is.na(Phylum) & Phylum != ""]

  meta_sub <- md[, .(sample_id, site = get(col_site), age = get(col_ageg),
                     continent = get(col_cont))]
  meta_sub[, sample_id := as.character(sample_id)]
  merged <- merge(rel_long, meta_sub, by = "sample_id", all.x = FALSE)

  # Clean site / age / continent (normalise the full-width "＞" character).
  merged[, site      := tolower(trimws(as.character(site)))]
  merged[, age       := gsub("＞", ">", as.character(age))]
  merged[, continent := trimws(as.character(continent))]

  merged <- merged[!is.na(site) & site != "" &
                     !(site %in% c("na", "other", "others",
                                   "other (specify)", "others (specify)"))]
  merged <- merged[!is.na(age) & age != "" & age != "NA"]
  merged <- merged[!is.na(continent) & continent != "" & continent != "NA"]

  top_phyla <- global_top_phyla(merged, TOP_N)
  merged <- merged[Phylum %in% top_phyla]
  merged[, Phylum := factor(Phylum, levels = top_phyla)]

  list(merged = merged, top_phyla = top_phyla)
}

# Mean relative abundance + prevalence per (group x phylum).
group_stats <- function(merged, by_cols) {
  merged[, .(
    Mean_Abundance = mean(Abundance, na.rm = TRUE),
    Prevalence     = mean(Abundance > 0, na.rm = TRUE),
    N              = .N
  ), by = c(by_cols, "Phylum")]
}

# ============================================================================
# Panel A : site x phylum bubble matrix
# ============================================================================
render_panel_A_site <- function(dat, out_dir) {
  merged <- dat$merged
  site_order  <- merged[, .N, by = site][order(-N)]$site
  site_labels <- to_title(site_order)

  stats <- group_stats(merged, "site")
  stats[, site   := factor(site,   levels = site_order)]
  stats[, Phylum := factor(Phylum, levels = rev(levels(merged$Phylum)))]

  p <- ggplot(stats, aes(x = site, y = Phylum,
                         size = Mean_Abundance, color = Prevalence)) +
    geom_point(alpha = 0.85) +
    scale_x_discrete(labels = site_labels) +
    scale_size_continuous(range = c(1, 10), trans = "sqrt",
                          breaks = c(0.01, 0.05, 0.1, 0.2, 0.3),
                          labels = scales::percent_format(accuracy = 1)) +
    scale_color_viridis_c(option = "C", limits = c(0, 1),
                          breaks = c(0, 0.25, 0.5, 0.75, 1),
                          labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = NULL,
         size = "Mean relative\nabundance", color = "Prevalence") +
    theme_fig3()

  ggplot2::ggsave(file.path(out_dir, "fig3_panelA_site_bubble.png"), p,
                  width = 7, height = 6, dpi = FIG_DPI, bg = "white")
  data.table::fwrite(stats, file.path(out_dir, "fig3_panelA_site_bubble_stats.csv"))
  p
}

# ============================================================================
# Panel B : age trajectory (faceted by site)
# ============================================================================
render_panel_B_age <- function(dat, out_dir) {
  merged <- dat$merged
  site_order  <- merged[, .N, by = site][order(-N)]$site
  site_labels <- to_title(site_order)
  age_order   <- c("0-6 months", "6-36 months", "3-12 years",
                   "12-18 years", "18-60 years", ">60 years")
  age_short   <- c("0-6 mo", "6-36 mo", "3-12 y", "12-18 y", "18-60 y", ">60 y")

  stats <- group_stats(merged, c("site", "age"))
  stats[, site   := factor(site,   levels = site_order)]
  stats[, age    := factor(age,    levels = age_order)]
  stats[, Phylum := factor(Phylum, levels = levels(merged$Phylum))]

  p <- ggplot(stats, aes(x = age, y = Mean_Abundance,
                         color = Phylum, group = Phylum)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) +
    facet_wrap(~ site, ncol = 3, scales = "free_y") +
    scale_x_discrete(labels = age_short) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values = phylum_color_map(levels(merged$Phylum)),
                       name = "Phylum") +
    labs(x = "Age group", y = "Mean relative abundance") +
    theme_fig3() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.text  = element_text(size = BASE_SZ - 1, face = "bold"))

  ggplot2::ggsave(file.path(out_dir, "fig3_panelB_age_trajectory.png"), p,
                  width = 14, height = 7, dpi = FIG_DPI, bg = "white")
  data.table::fwrite(stats, file.path(out_dir, "fig3_panelB_age_trajectory_stats.csv"))
  p
}

# ============================================================================
# Panel C : continent x phylum bubble matrix
# ============================================================================
render_panel_C_continent <- function(dat, out_dir) {
  merged <- dat$merged
  cont_order <- merged[, .N, by = continent][order(-N)]$continent

  stats <- group_stats(merged, "continent")
  stats[, continent := factor(continent, levels = cont_order)]
  stats[, Phylum    := factor(Phylum,    levels = rev(levels(merged$Phylum)))]

  p <- ggplot(stats, aes(x = continent, y = Phylum,
                         size = Mean_Abundance, color = Prevalence)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(1, 10), trans = "sqrt",
                          breaks = c(0.01, 0.05, 0.1, 0.2, 0.3),
                          labels = scales::percent_format(accuracy = 1)) +
    scale_color_viridis_c(option = "C", limits = c(0, 1),
                          breaks = c(0, 0.25, 0.5, 0.75, 1),
                          labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = NULL,
         size = "Mean relative\nabundance", color = "Prevalence") +
    theme_fig3()

  ggplot2::ggsave(file.path(out_dir, "fig3_panelC_continent_bubble.png"), p,
                  width = 7, height = 6, dpi = FIG_DPI, bg = "white")
  data.table::fwrite(stats, file.path(out_dir, "fig3_panelC_continent_bubble_stats.csv"))
  p
}

# ============================================================================
# Combined figure
# ============================================================================
render_combined <- function(p_a, p_b, p_c, out_dir) {
  # Layout: site (A) and continent (C) bubble matrices side by side on top,
  # age trajectory (B) full width below. Tags follow visual order:
  # A = site, B = continent, C = age.
  combined <- (p_a | p_c) / p_b +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 18, face = "bold"))

  ggplot2::ggsave(file.path(out_dir, "fig3_panels_ABC_combined.png"), combined,
                  width = 14, height = 13, dpi = FIG_DPI, bg = "white")
  invisible(NULL)
}

# ============================================================================
# Main
# ============================================================================
main <- function() {
  if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)
  cat("Script dir :", script_dir, "\n")
  cat("PANEL_DIR  :", PANEL_DIR, "\n")
  cat("metadata   :", metadata_path, "\n")
  cat("phylum tab :", abund_file, "\n")

  dat <- prepare_data()
  cat("samples    :", data.table::uniqueN(dat$merged$sample_id), "\n")
  cat("top phyla  :", paste(dat$top_phyla, collapse = ", "), "\n")

  p_a <- render_panel_A_site(dat, PANEL_DIR)
  cat("panel A (site bubble) done\n")
  p_b <- render_panel_B_age(dat, PANEL_DIR)
  cat("panel B (age trajectory) done\n")
  p_c <- render_panel_C_continent(dat, PANEL_DIR)
  cat("panel C (continent bubble) done\n")

  render_combined(p_a, p_b, p_c, PANEL_DIR)
  cat("combined figure done\n")
  cat("All outputs written to:", PANEL_DIR, "\n")
}

main()
