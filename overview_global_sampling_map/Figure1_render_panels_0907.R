#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# figure1_render_panels_0907.R
#
# Renders the Figure 1 panels (and the two supplementary stacked-bar figures)
# into PANEL_DIR as PNG files. This is the HEAVY analysis step: it reads the
# 80,788-sample metadata, the phylum/class/order relative-abundance tables and
# the rarefaction table, and can take ~15 minutes. Run it ONCE.
#
# After the panels exist, assemble the composite figure quickly and repeatedly
# with the companion script:
#     Rscript figure1_composite_all_in_one_0907.R
# (that script only reads the panel PNGs and lays them out - it is fast).
#
# Panels / data inputs:
#   A    : pipeline funnel                      (ggplot2, ggforce)      [unchanged]
#   B    : sequencing-depth histogram           (tidyverse, scales)    [NEW]
#   C    : global sampling world map            (tidyverse, sf, ...)    [unchanged]
#   D    : core-taxa prevalence phylum/class/order (tidyverse,data.table) [NEW]
#   E    : top-5 phylum relative-abundance density (tidyverse,data.table) [NEW]
#   F    : alpha diversity (Shannon) histogram   (tidyverse, scales)   [NEW]
#   G    : rarefaction / species-accumulation curves  (tidyverse, mgcv, scales) [unchanged]
#   supp : phylum x continent, phylum x age stacked bars  (tidyverse, data.table)
#
# All phylum-colored panels (D, E, stacked bars) share one global top-N phylum
# list and one color map (see helper section). New panels use all 80,788
# samples in input/metadata_oral_qc.tsv (the dedup'd 51,186 count appears only
# inside the funnel narrative).
#
# Input/output paths and env overrides are documented in the CONFIG block
# below; on machines without the server tree all defaults fall back to
# script-relative folders automatically.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(png)   # readPNG
  library(grid)  # grid.raster / grid.text / gpar
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
# Production tree. Its  input/  folder holds only biom / metadata / taxonomy,
# so the rarefaction table and the phylum table keep their original server
# paths below. Outputs go to .../overview_figure_1/results/. Absent this tree
# (local checkout) everything falls back to script-relative folders. Env
# overrides: FIG1_BASE_DIR, FIG1_PANEL_DIR, FIG1_METADATA, FIG1_RARE_TABLE,
# FIG1_PHYLUM_TABLE, FIG1_AUTO_INSTALL.
SERVER_ROOT <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-map-2/04-patch-pipeline-map/03-patch-pipeline-map"
# Server-side output root for the collapsed relative tables (phylum/class/order).
# These are produced by phylum_relative_collapse-260802.sh and
# taxa_relative_collapse_multi_level.sh into .../02-patch-pipeline-map/output/.
SERVER_OUT_ROOT <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-pipeline-map/output"
LOCAL_BASE  <- normalizePath(file.path(script_dir, ".."))
# True only on the production server, where BOTH the input tree and the
# collapsed-table output root exist. A local checkout mirrors the input/ folder
# but not the server output root, so it correctly falls back to local paths.
on_server   <- dir.exists(file.path(SERVER_ROOT, "input")) &&
               dir.exists(SERVER_OUT_ROOT)
if (on_server) {
  BASE_DIR  <- SERVER_ROOT
  PANEL_DIR <- file.path(SERVER_ROOT, "overview_figure_1", "results")
} else {
  BASE_DIR  <- LOCAL_BASE
  PANEL_DIR <- file.path(script_dir, "results")
}
BASE_DIR    <- Sys.getenv("FIG1_BASE_DIR",  unset = BASE_DIR)
PANEL_DIR   <- Sys.getenv("FIG1_PANEL_DIR", unset = PANEL_DIR)
AUTO_INSTALL <- as.logical(Sys.getenv("FIG1_AUTO_INSTALL", unset = "FALSE"))

metadata_path <- Sys.getenv("FIG1_METADATA",
                            unset = file.path(BASE_DIR, "input",
                                              "metadata_oral_qc.tsv"))
rare_file     <- Sys.getenv("FIG1_RARE_TABLE",
                            unset = if (on_server)
  "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/01-left-pipeline-right-map/output/rarefaction_pcofg.tsv"
else file.path(LOCAL_BASE, "output", "rarefaction_pcofg.tsv"))
abund_file    <- Sys.getenv("FIG1_PHYLUM_TABLE",
                            unset = if (on_server)
  file.path(SERVER_OUT_ROOT, "phylum-relative-table",
            "phylum-relative-table.tsv")
else file.path(LOCAL_BASE, "output", "phylum-relative-table",
               "phylum-relative-table.tsv"))
# Class/order collapsed relative tables: generated on the server by
# taxa_relative_collapse_multi_level.sh (LEVELS="3 4") using the same
# phylum_relative_collapse-260802.sh convention (writes into the parent
# pipeline-map project's output/<level>-relative-table/). Defaults therefore
# sit next to the phylum table; override via FIG1_CLASS_TABLE / FIG1_ORDER_TABLE.
class_file    <- Sys.getenv("FIG1_CLASS_TABLE",
                            unset = if (on_server)
  file.path(SERVER_OUT_ROOT, "class-relative-table",
            "class-relative-table.tsv")
else file.path(LOCAL_BASE, "output", "class-relative-table",
               "class-relative-table.tsv"))
order_file    <- Sys.getenv("FIG1_ORDER_TABLE",
                            unset = if (on_server)
  file.path(SERVER_OUT_ROOT, "order-relative-table",
            "order-relative-table.tsv")
else file.path(LOCAL_BASE, "output", "order-relative-table",
               "order-relative-table.tsv"))

# Canonical panel filenames. Keys A/B/C are referenced by render functions kept
# verbatim (funnel / rarefaction / map); new panels and the two supplementary
# stacked-bar panels use their own keys. Composite letters are mapped
# separately via COMPOSITE_ROWS.
panel_png <- c(
  A     = "pipeline_two_stage_funnel_HD.png",
  B     = "rarefaction_curve.png",
  C     = "global_sampling_world_map_final_cropped_260904.png",
  DEPTH = "seq_depth_histogram.png",
  PREV  = "prevalence_phyla_class_order.png",
  DENS  = "top5_phylum_density.png",
  SHAN  = "shannon_diversity.png",
  CONT  = "supp_phylum_continent_site_stacked.png",
  AGE   = "supp_phylum_age_site_stacked.png"
)

# Composite figure settings
compo_base    <- "Figure1_combined_260907"
compo_w       <- 16L
compo_dpi     <- 300L
compo_mg      <- 0.14       # outer margin (inches)
compo_gap     <- 0.18       # gap between panels within a row (inches)
compo_row_gap <- 0.18       # gap between rows (inches)
compo_tag_inset <- c(x = 0.12, y = 0.10)
compo_tag_size  <- 30

# Row layout for the composite. Each entry is a named vector mapping the
# composite letter (e.g. "A") to a panel_png key (e.g. "A"/"B"/"DEPTH"/...).
composite_rows <- list(
  c(A = "A",     B = "DEPTH", G = "B"),
  c(C = "C",     E = "DENS",  F = "SHAN"),
  c(D = "PREV")
)

# ============================================================================
# Small helpers
# ============================================================================
# Returns NULL when every package is available, else a message.
load_pkgs <- function(pkgs, auto_install = FALSE) {
  pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(pkgs) == 0) return(invisible(NULL))
  if (auto_install) {
    install.packages(pkgs, repos = "http://cran.us.r-project.org")
    pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(pkgs) == 0) return(invisible(NULL))
  }
  paste("missing package(s):", paste(pkgs, collapse = ", "))
}
load_libs <- function(pkgs, auto_install = FALSE) {
  ok <- load_pkgs(pkgs, auto_install)
  if (!is.null(ok)) return(ok)
  suppressPackageStartupMessages(lapply(pkgs, library, character.only = TRUE))
  invisible(NULL)
}

# ============================================================================
# Shared helpers for all phylum-colored panels (D, E, supplementary stacked bars)
# ============================================================================
# Fixed qualitative palette shared by every phylum color assignment; the first
# n colors are tied (by rank) to the global top-n phyla, so panels always agree.
core_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                  "#FFFF33", "#A65628", "#F781BF", "#999999")
other_color  <- "#CCCCCC"

# Terminal taxon name from a qiime2-collapsed label, e.g.
#   "d__Bacteria; p__Firmicutes_X; c__Clostridia_Y"  ->  "Clostridia_Y"
#   "d__Bacteria; p__Bacillota_A_1"                  ->  "Bacillota_A_1"
last_taxon_name <- function(label) {
  vapply(strsplit(as.character(label), ";\\s*"), function(parts) {
    sub("^[dpcogfs]__\\s*", "", tail(parts, 1))
  }, character(1))
}

# Parent phylum ("p__..." segment) from a collapsed label, or NA if missing.
extract_parent_phylum <- function(label) {
  vapply(as.character(label), function(x) {
    m <- regmatches(x, regexpr("p__[^;]+", x))
    if (length(m) == 0) NA_character_ else sub("^p__", "", m)
  }, character(1))
}

# Metadata column resolution: always prefer the *_cleaned / *_cleaned_group
# column; fall back to the raw column only when no clean variant exists in the
# given metadata table. Returns the name to use (the preferred name even when
# both are absent, so callers keep failing loudly in that case).
pick_meta_col <- function(md, preferred, fallback = NULL) {
  if (preferred %in% names(md)) return(preferred)
  if (!is.null(fallback) && fallback %in% names(md)) return(fallback)
  preferred
}

# Read a collapsed relative-abundance table (top line starts with '#') and
# return a long data.table: Taxon, sample_id, Abundance, parent_phylum,
# Taxon_final.
read_collapsed_rel_table <- function(file) {
  dt <- data.table::fread(file, header = TRUE, sep = "\t", skip = 1)
  sample_cols <- setdiff(names(dt), names(dt)[1])
  data.table::setnames(dt, names(dt)[1], "Taxon")
  long <- data.table::melt(dt, id.vars = "Taxon", measure.vars = sample_cols,
                           variable.name = "sample_id", value.name = "Abundance")
  long[, sample_id := as.character(sample_id)]
  long[, parent_phylum := extract_parent_phylum(Taxon)]
  long[, Taxon_final := last_taxon_name(Taxon)]
  long <- long[!is.na(Taxon_final) & Taxon_final != "" & !is.na(Abundance)]
  long
}

# Fast, melt-free prevalence for a collapsed relative-abundance table.
# Prevalence of a taxon is a pure row statistic (fraction of samples with
# abundance > 0), so this avoids melting over the ~80k sample columns, which is
# the dominant cost of the long-format approach for the class/order tables.
# Averaging over taxa that collapse to the same (parent, final) name matches
# the melt-based aggregation exactly.
prevalence_collapsed_table <- function(file) {
  dt <- data.table::fread(file, header = TRUE, sep = "\t", skip = 1)
  sample_cols <- setdiff(names(dt), names(dt)[1])
  data.table::setnames(dt, names(dt)[1], "Taxon")
  out <- data.table::data.table(
    parent_phylum = extract_parent_phylum(dt$Taxon),
    Taxon_final   = last_taxon_name(dt$Taxon),
    Count         = rowSums(dt[, ..sample_cols] > 0, na.rm = TRUE)
  )
  out <- out[!is.na(Taxon_final) & Taxon_final != ""]
  out[, .(Count = mean(Count)), by = .(parent_phylum, Taxon_final)]
}

# Global top-n phyla by summed relative abundance across all samples.
global_top_phyla <- function(long_dt, n = 9) {
  tot <- long_dt[, .(Total = sum(Abundance, na.rm = TRUE)), by = parent_phylum]
  tot <- tot[!is.na(parent_phylum) & parent_phylum != ""][order(-Total)]
  as.character(head(tot$parent_phylum, n))
}

# Named color vector: one color per top phylum (by rank), plus "Other" grey.
phylum_color_map <- function(top_phyla) {
  pal <- rep(core_palette, length.out = length(top_phyla))
  c(setNames(pal, top_phyla), Other = other_color)
}

# ============================================================================
# Panel A : pipeline funnel (merged from pipeline_two_stage_funnel.R)
# ============================================================================
render_panel_A_funnel <- function(out_dir) {
  if (!is.null(err <- load_libs(c("ggplot2", "ggforce"), AUTO_INSTALL))) return(err)

  palette <- c(blue = "#E9F1F8", orange = "#F9EDE0", yellow = "#FBF4E2", green = "#EAF3E9")
  zone_title_pal <- c(blue = "#2F5D7E", orange = "#A5672A", yellow = "#8A6D1F", green = "#3E6B4F")

  node_fill   <- "white"
  node_border <- "#9AA3AB"
  arrow_col   <- "#7E8A94"
  bullet_col  <- "#8A94A0"
  txt_main    <- "#262626"
  txt_sub     <- "#3F3F3F"
  txt_detail  <- "#6B6B6B"
  txt_excl    <- "#B75D69"

  # Master font scale for panel A text (A_TXT inside step nodes, A_TXT_SUB for
  # zone titles / side bullets).
  A_TXT     <- as.numeric(Sys.getenv("FIG1_A_TXT_SCALE",     unset = "1.40"))
  A_TXT_SUB <- as.numeric(Sys.getenv("FIG1_A_TXT_SCALE_SUB", unset = "1.25"))

  rect_xy <- function(cx, cy, w, h, grp = NULL) {
    d <- data.frame(
      x = c(cx - w/2, cx + w/2, cx + w/2, cx - w/2),
      y = c(cy - h/2, cy - h/2, cy + h/2, cy + h/2)
    )
    if (!is.null(grp)) d$grp <- grp
    d
  }
  arrowhead_down <- function(cx, cy, a = 3.2, b = 1.9)
    data.frame(x = c(cx - b, cx + b, cx), y = c(cy, cy, cy - a))

  # Canvas & zones (x domain widened so the enlarged text fits inside its boxes)
  zones <- data.frame(
    name = c("1 . Sample retrieval from NCBI SRA", "2 . Metadata harmonization",
             "3 . Filtering and sequencing-data linkage", "4 . Deduplicated resource"),
    xmin = c(-40, 63, -40, -40), xmax = c(63, 160, 160, 160),
    ymin = c(82, 82, 17, -20), ymax = c(127.5, 127.5, 80, 15),
    key  = c("blue", "orange", "yellow", "green")
  )
  zones$fill <- palette[zones$key]
  zones$tcol <- zone_title_pal[zones$key]
  zones$ty <- zones$ymax - 2.2

  zones_poly <- do.call(rbind, lapply(seq_len(nrow(zones)), function(i)
    rect_xy((zones$xmin[i] + zones$xmax[i])/2, (zones$ymin[i] + zones$ymax[i])/2,
            zones$xmax[i] - zones$xmin[i], zones$ymax[i] - zones$ymin[i], paste0("z", i))))
  zones_poly$fill <- rep(zones$fill, each = 4)

  # Top modules (parallel tasks)
  bullet <- data.frame(
    x = c(rep(6.5, 4), rep(94.0, 4)),
    y = c(121, 118, 115, 112)
  )
  row_txt <- data.frame(
    x = c(rep(9.5, 4), rep(97.0, 4)),
    y = c(121, 118, 115, 112),
    label = c("116,643 samples from NCBI", "334 BioProjects from NCBI",
              "16S rRNA amplicons · DADA2", "All body sites (oral + extra-oral)",
              "Country (47 countries)", "Body site (7 oral categories)",
              "Health status (healthy + diseased)", "Country · site · platform · age · sex")
  )

  top_pill <- data.frame(id = c("tblue", "torange"), cx = c(11.5, 111), cy = c(95.5, 92),
                         w = c(74, 84), h = c(29, 21))
  top_poly <- do.call(rbind, lapply(seq_len(nrow(top_pill)), function(i)
    rect_xy(top_pill$cx[i], top_pill$cy[i], top_pill$w[i], top_pill$h[i], top_pill$id[i])))

  top_txt <- do.call(rbind, list(
    data.frame(x = 11.5,
               y = 95.5 + c(11.2, 8.0, 4.8, 1.6, -1.6, -4.8, -8.0, -11.2),
               label = c("116,643 samples", "after Step 1 · Sample retrieval",
                         "116,624 samples", "after Step 2 · SRA download validation",
                         "114,959 samples", "after Step 3 · Quality filtering & DADA2",
                         "Excluded: low-quality or denoising failures", "−1,665 samples"),
               size = c(5.2, 3.8, 5.2, 3.8, 5.2, 3.8, 3.5, 3.5),
               fontface = c("bold", "plain", "bold", "plain", "bold", "plain", "italic", "italic"),
               color = c(txt_main, txt_sub, txt_main, txt_sub, txt_main, txt_sub, txt_detail, txt_excl)),
    data.frame(x = 111.5,
               y = 92 + c(3.6, 0, -3.6),
               label = c("Harmonized attribute table", "country · site · status · age · sex",
                         "standardized metadata for linkage"),
               size = c(5.2, 3.8, 3.5), fontface = c("bold", "plain", "italic"),
               color = c(txt_main, txt_sub, txt_detail))
  ))

  row_arrow <- data.frame(x = 111, y = 108, xend = 111, yend = 101.0)

  merge_seg <- data.frame(x = c(11.5, 111, 11.5, 60), y = c(82, 83.25, 80, 80),
                          xend = c(11.5, 111, 111, 60), yend = c(80, 80, 80, 73.0))
  merge_head <- arrowhead_down(60, 73.0, a = 2.2, b = 1.3)
  merge_dot  <- data.frame(x = c(11.5, 111), y = c(80, 80))

  # Funnel nodes
  step_pill <- data.frame(
    id = paste0("s", 1:4),
    cx = rep(60, 4), cy = c(64, 44, 24, -3),
    w = rep(120, 4), h = c(17.5, 17.5, 17.5, 34)
  )
  step_poly <- do.call(rbind, lapply(seq_len(nrow(step_pill)), function(i)
    rect_xy(step_pill$cx[i], step_pill$cy[i], step_pill$w[i], step_pill$h[i], step_pill$id[i])))

  step_line <- do.call(rbind, list(
    data.frame(x = 60, y = 64 + c(6.0, 2.0, -2.0, -6.0),
               label = c("113,233 samples", "after Step 4 · Chimera Filtering",
                         "Excluded: >30% of oral samples with a chimera rate >25%",
                         "−1,726 samples · −22 studies"),
               size = c(5.2, 4.2, 3.8, 3.8), fontface = c("bold", "plain", "italic", "italic"),
               color = c(txt_main, txt_sub, txt_detail, txt_excl)),
    data.frame(x = 60, y = 44 + c(6.0, 2.0, -2.0, -6.0),
               label = c("108,691 samples", "after Step 5 · Closed-Reference Mapping",
                         "Excluded: oral aggregate GG2 closed-reference mapping rate <70%",
                         "−4,542 samples · −24 studies"),
               size = c(5.2, 4.2, 3.8, 3.8), fontface = c("bold", "plain", "italic", "italic"),
               color = c(txt_main, txt_sub, txt_detail, txt_excl)),
    data.frame(x = 60, y = 24 + c(6.0, 2.0, -2.0, -6.0),
               label = c("80,788 samples", "after Step 6 · Site Filtering",
                         "Excluded: non-oral samples (site_group ≠ oral)",
                         "−27,903 samples"),
               size = c(5.2, 4.2, 3.8, 3.8), fontface = c("bold", "plain", "italic", "italic"),
               color = c(txt_main, txt_sub, txt_detail, txt_excl)),
    data.frame(x = 60, y = -3 + c(14.375, 9.75, 5.125, 0.5, -5.125, -9.75, -14.375),
               label = c("51,186 samples", "after Step 7 · Deduplication",
                         "Excluded: >1 sample per subject × site",
                         "−29,602 samples", "deduplicated resource · 288 studies",
                         "42,224 subjects · 61,235 features",
                         "47 countries · 7 oral site categories"),
               size = c(5.2, 4.2, 3.8, 3.8, 4.2, 3.8, 3.8),
               fontface = c("bold", "plain", "italic", "italic", "plain", "italic", "italic"),
               color = c(txt_main, txt_sub, txt_detail, txt_excl, txt_sub, txt_detail, txt_detail))
  ))

  final_div <- data.frame(x = 36, y = -5.25, xend = 84, yend = -5.25)

  # Apply the text-scale factor to the in-panel text (headlines/captions/notes)
  top_txt$size   <- top_txt$size  * A_TXT
  step_line$size <- step_line$size * A_TXT
  step_arrow <- data.frame(x = rep(60, 3), y = c(55.25, 35.25, 15.25),
                           xend = rep(60, 3), yend = c(52.75, 32.75, 14.20))

  # Left-align zone bullets with the "." of each "N . Title" step label,
  # measured at render time so the offset adapts to device/fonts.
  zone_font_pt <- 4.0 * A_TXT_SUB * ggplot2::.pt
  grid::pushViewport(grid::viewport(gp = grid::gpar(fontsize = zone_font_pt, fontface = "bold")))
  pref_blue_mm   <- as.numeric(grid::convertWidth(grid::stringWidth("1 "), "mm", valueOnly = TRUE))
  pref_orange_mm <- as.numeric(grid::convertWidth(grid::stringWidth("2 "), "mm", valueOnly = TRUE))
  grid::popViewport()
  mg_mm       <- 6 / ggplot2::.pt
  dev_w_mm    <- 13 * 25.4 - 2 * mg_mm
  dev_h_mm    <- 10.6 * 25.4 - 2 * mg_mm
  xrange      <- (160 - (-40)) * 1.1
  yrange      <- (134 - (-22)) * 1.1
  panel_w_mm  <- min(dev_w_mm, dev_h_mm * xrange / yrange)
  units_per_mm <- xrange / panel_w_mm
  dot_blue   <- (-40 + 2.5) + pref_blue_mm   * units_per_mm
  dot_orange <- (63  + 2.5) + pref_orange_mm * units_per_mm
  row_txt$x  <- c(rep(dot_blue, 4), rep(dot_orange, 4))
  bullet$x   <- row_txt$x - 3.0

  p2 <- ggplot2::ggplot() +
    ggforce::geom_shape(data = zones_poly, aes(x = x, y = y, group = grp),
                        fill = zones_poly$fill, radius = unit(0.35, "cm"), color = NA) +
    ggplot2::geom_text(data = zones, aes(x = xmin + 2.5, y = ty, label = name, color = tcol),
                       hjust = 0, vjust = 1, fontface = "bold", size = 4.0 * A_TXT_SUB) +
    ggplot2::scale_color_identity() +
    ggplot2::geom_point(data = bullet, aes(x = x, y = y), shape = 16, size = 1.6, color = bullet_col) +
    ggplot2::geom_text(data = row_txt, aes(x = x, y = y, label = label), hjust = 0, size = 4.0 * A_TXT_SUB, color = txt_sub) +
    ggplot2::geom_segment(data = row_arrow, aes(x = x, y = y, xend = xend, yend = yend),
                          color = arrow_col, linewidth = 1.0,
                          arrow = ggplot2::arrow(length = unit(0.18, "cm"), type = "closed")) +
    ggplot2::geom_segment(data = merge_seg, aes(x = x, y = y, xend = xend, yend = yend),
                          color = arrow_col, linewidth = 1.0) +
    ggplot2::geom_polygon(data = merge_head, aes(x = x, y = y), fill = arrow_col) +
    ggplot2::geom_point(data = merge_dot, aes(x = x, y = y), shape = 21, size = 1.8,
                        fill = arrow_col, color = arrow_col) +
    ggforce::geom_shape(data = top_poly, aes(x = x, y = y, group = grp),
                        fill = node_fill, color = node_border, radius = unit(0.35, "cm"), linewidth = 0.8) +
    ggplot2::geom_text(data = top_txt, aes(x = x, y = y, label = label),
                       size = top_txt$size, fontface = top_txt$fontface, color = top_txt$color) +
    ggforce::geom_shape(data = step_poly, aes(x = x, y = y, group = grp),
                        fill = node_fill, color = node_border, radius = unit(0.35, "cm"), linewidth = 0.8) +
    ggplot2::geom_text(data = step_line, aes(x = x, y = y, label = label),
                       size = step_line$size, fontface = step_line$fontface, color = step_line$color) +
    ggplot2::geom_segment(data = final_div, aes(x = x, y = y, xend = xend, yend = yend),
                          color = node_border, linewidth = 0.5) +
    ggplot2::geom_segment(data = step_arrow, aes(x = x, y = y, xend = xend, yend = yend),
                          color = arrow_col, linewidth = 1.0,
                          arrow = ggplot2::arrow(length = unit(0.18, "cm"), type = "closed")) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(-40, 160), ylim = c(-22, 134)) +
    ggplot2::theme_void() + ggplot2::theme(plot.margin = margin(6, 6, 6, 6, "pt"))

  ggplot2::ggsave(file.path(out_dir, panel_png["A"]), plot = p2, width = 13, height = 10.6,
                  device = "png", dpi = 600, bg = "white")
  invisible(NULL)
}

# ============================================================================
# Panel B : rarefaction curves (merged from rarefaction_curve_plot.R)
# ============================================================================
render_panel_B_rarefaction <- function(out_dir) {
  if (!is.null(err <- load_libs(c("tidyverse", "mgcv", "scales"), AUTO_INSTALL))) return(err)
  if (!file.exists(rare_file))
    return(sprintf("rarefaction table not found: %s", rare_file))

  base_sz <- 22
  text_sz <- base_sz / .pt

  ggplot2::update_geom_defaults("text",  list(family = "sans"))
  ggplot2::update_geom_defaults("label", list(family = "sans"))
  ggplot2::theme_set(ggplot2::theme_minimal(base_size = base_sz, base_family = "sans"))

  df_rare <- read.table(rare_file, header = TRUE, sep = "\t") %>%
    pivot_longer(cols = c("Phylum", "Class", "Order", "Family", "Genus"),
                 names_to = "Taxonomy", values_to = "Unique_taxa") %>%
    mutate(Taxonomy = factor(Taxonomy,
                             levels = c("Genus", "Family", "Order", "Class", "Phylum")))

  p_e <- ggplot(df_rare, aes(x = Samples, y = Unique_taxa, color = Taxonomy)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = 0.7) +
    scale_color_manual(values = c(
      "Genus"  = "#5C6BC0", "Family" = "#FFB300", "Order" = "#EF5350",
      "Class"  = "#26C6DA", "Phylum" = "#66BB6A")) +
    scale_x_continuous(breaks = seq(0, 80000, by = 20000), limits = c(0, 98000),
                       labels = c("0", "20K", "40K", "60K", "80K"), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000),
                       labels = c("0", "1K", "2K", "3K", "4K", "5K"),
                       expand = expansion(mult = c(0.1, 0.06))) +
    annotate("text", x = 81000, y = 4900, label = "Genus", color = "#5C6BC0",
             fontface = "bold", size = base_sz / 3.8, hjust = 0) +
    annotate("text", x = 81000, y = 1400, label = "Family", color = "#FFB300",
             fontface = "bold", size = base_sz / 3.8, hjust = 0) +
    annotate("text", x = 81000, y = 800, label = "Order", color = "#EF5350",
             fontface = "bold", size = base_sz / 3.8, hjust = 0) +
    annotate("text", x = 81000, y = 350, label = "Class", color = "#26C6DA",
             fontface = "bold", size = base_sz / 3.8, hjust = 0) +
    annotate("text", x = 81000, y = 5, label = "Phylum", color = "#66BB6A",
             fontface = "bold", size = base_sz / 3.8, hjust = 0) +
    labs(x = "Samples", y = "Observed taxa") +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = base_sz - 4, base_family = "sans") +
    theme(
      axis.text.x     = element_text(size = base_sz * 0.72, face = "bold", color = "#303030"),
      axis.text.y     = element_text(size = base_sz * 0.72, face = "bold", color = "#303030"),
      axis.title.y    = element_text(size = base_sz * 0.85, face = "bold", color = "#303030",
                                     margin = margin(r = 6)),
      axis.title.x    = element_text(size = base_sz * 0.85, face = "bold", color = "#303030",
                                     margin = margin(t = 10)),
      legend.position = "none",
      plot.margin     = margin(5, 5, 5, 5, unit = "pt"),
      plot.background = element_blank(),
      panel.background = element_blank()
    )

  ggplot2::ggsave(file.path(out_dir, panel_png["B"]), plot = p_e, width = 7, height = 5.5,
                  device = "png", dpi = 600, bg = "white")
  invisible(NULL)
}

# ============================================================================
# Panel B : sequencing-depth histogram (new; 80,788-sample metadata)
#   x = merged reads per sample on a log10 axis, dashed line = median
# ============================================================================
render_panel_depth <- function(out_dir) {
  if (!is.null(err <- load_libs(c("tidyverse", "scales"), AUTO_INSTALL))) return(err)
  if (!file.exists(metadata_path))
    return(sprintf("metadata table not found: %s", metadata_path))
  base_sz <- as.numeric(Sys.getenv("FIG1_DEPTH_BASE_SZ", unset = "20"))

  depth <- readr::read_tsv(metadata_path, show_col_types = FALSE) %>%
    filter(!is.na(total_reads), total_reads > 0)
  if (nrow(depth) == 0) return("no usable total_reads in metadata")
  med_reads <- median(depth$total_reads)

  p <- ggplot(depth, aes(x = total_reads)) +
    geom_histogram(bins = 60, fill = "#5A4A42", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = med_reads, linetype = "dashed", color = "#B75D69", linewidth = 0.7) +
    annotate("text", x = med_reads, y = Inf, vjust = 1.6, hjust = 0.5,
             label = paste0("median = ", scales::comma(round(med_reads))),
             color = "#B75D69", fontface = "bold", size = base_sz / 3.5) +
    scale_x_log10(limits = c(1e4, NA),
                  breaks = 10^(4:6),
                  labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "Merged reads", y = "Number of samples") +
    theme_classic(base_size = base_sz - 4, base_family = "sans") +
    theme(
      axis.text         = element_text(size = base_sz * 0.72, face = "bold", color = "#303030"),
      axis.title        = element_text(size = base_sz * 0.82, face = "bold", color = "#303030"),
      plot.background   = element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA)
    )

  ggplot2::ggsave(file.path(out_dir, panel_png["DEPTH"]), plot = p,
                  width = 7, height = 5.2, device = "png", dpi = 600, bg = "white")
  invisible(NULL)
}

# ============================================================================
# Panel F : alpha diversity (Shannon index) histogram
# ============================================================================
render_panel_shannon <- function(out_dir) {
  if (!is.null(err <- load_libs(c("tidyverse", "scales"), AUTO_INSTALL))) return(err)
  if (!file.exists(metadata_path))
    return(sprintf("metadata table not found: %s", metadata_path))
  base_sz <- as.numeric(Sys.getenv("FIG1_SHAN_BASE_SZ", unset = "20"))

  shan <- readr::read_tsv(metadata_path, show_col_types = FALSE) %>%
    filter(!is.na(Shannon))
  if (nrow(shan) == 0) return("no usable Shannon values in metadata")
  med_shan <- median(shan$Shannon)

  p <- ggplot(shan, aes(x = Shannon)) +
    geom_histogram(bins = 60, fill = "#3E6B4F", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = med_shan, linetype = "dashed", color = "#B75D69", linewidth = 0.7) +
    annotate("text", x = med_shan, y = Inf, vjust = 1.6, hjust = 0.5,
             label = sprintf("median = %.2f", med_shan),
             color = "#B75D69", fontface = "bold", size = base_sz / 3.5) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "Shannon diversity", y = "Number of samples") +
    theme_classic(base_size = base_sz - 4, base_family = "sans") +
    theme(
      axis.text         = element_text(size = base_sz * 0.72, face = "bold", color = "#303030"),
      axis.title        = element_text(size = base_sz * 0.82, face = "bold", color = "#303030"),
      plot.background   = element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA)
    )

  ggplot2::ggsave(file.path(out_dir, panel_png["SHAN"]), plot = p,
                  width = 7, height = 5.2, device = "png", dpi = 600, bg = "white")
  invisible(NULL)
}

# ============================================================================
# Panel E : top-5 phylum relative-abundance density (pseudo-log y axis)
#   five coloured curves = global top-5 phyla, grey line = all other phyla
# ============================================================================
render_panel_density <- function(out_dir) {
  if (!is.null(err <- load_libs(c("tidyverse", "data.table"), AUTO_INSTALL))) return(err)
  if (!file.exists(abund_file))
    return(sprintf("phylum table not found: %s", abund_file))
  base_sz <- as.numeric(Sys.getenv("FIG1_DENS_BASE_SZ", unset = "20"))

  rel_long <- read_collapsed_rel_table(abund_file)
  if (nrow(rel_long) == 0) return("phylum table is empty after cleaning")

  top9 <- global_top_phyla(rel_long, 9)
  top5 <- head(top9, 5)
  colmap <- phylum_color_map(top9)

  top_lines   <- rel_long[parent_phylum %in% top5,
                          .(sample_id, Abundance, Group = parent_phylum)]
  other_lines <- rel_long[!(parent_phylum %in% top5),
                          .(Abundance = sum(Abundance, na.rm = TRUE),
                            Group     = "Others"),
                          by = sample_id]
  dens_data <- rbindlist(list(top_lines, other_lines), use.names = TRUE)
  dens_data[, Group := factor(Group, levels = c(top5, "Others"))]
  dens_colors <- c(setNames(colmap[top5], top5), "Others" = "#9E9E9E")

  p <- ggplot(dens_data, aes(x = Abundance, color = Group)) +
    # key_glyph = "path" overrides ggplot2 4.x's default hollow-rectangle key
    # (GeomDensity inherits GeomArea/GeomRibbon) so each phylum is drawn as a
    # coloured line segment in the legend.
    geom_density(linewidth = 0.9, key_glyph = "path") +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                       limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(trans = "pseudo_log", labels = scales::label_number()) +
    scale_color_manual(values = dens_colors, name = "Phylum") +
    labs(x = "Relative abundance",
         y = "Samples") +
    theme_classic(base_size = base_sz - 4, base_family = "sans") +
    theme(
      legend.position  = "right",
      legend.title     = element_text(size = base_sz * 0.78, face = "bold"),
      legend.text      = element_text(size = base_sz * 0.7),
      axis.text        = element_text(size = base_sz * 0.72, face = "bold", color = "#303030"),
      axis.title       = element_text(size = base_sz * 0.82, face = "bold", color = "#303030"),
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )

  ggplot2::ggsave(file.path(out_dir, panel_png["DENS"]), plot = p,
                  width = 8, height = 5.5, device = "png", dpi = 600, bg = "white")
  invisible(NULL)
}

# ============================================================================
# Panel D : core-taxa prevalence, phylum / class / order
#   Three independent horizontal-bar plots (one per level) stitched side by
#   side with patchwork. Bars are sorted by descending sample count (longest at
#   the top). Only the top-5 phyla get bright colours; every other phylum is
#   grey. Class and order bars inherit the colour of their parent phylum, so
#   the reader can see at a glance which phylum each class/order belongs to.
#   Taxon names are drawn inside the bars; bars carry a thin outline, the panel
#   background is warm white with vertical-only gridlines and a solid frame.
# ============================================================================
render_panel_prevalence <- function(out_dir) {
  if (!is.null(err <- load_libs(c("tidyverse", "data.table", "patchwork"),
                                AUTO_INSTALL))) return(err)
  if (!file.exists(abund_file))
    return(sprintf("phylum table not found: %s", abund_file))
  base_sz <- as.numeric(Sys.getenv("FIG1_PREV_BASE_SZ", unset = "15"))

  tab_files <- c(Phylum = abund_file)
  if (file.exists(class_file)) {
    tab_files <- c(tab_files, Class = class_file)
  } else {
    warning("class relative table not found; Class panel skipped: ", class_file)
  }
  if (file.exists(order_file)) {
    tab_files <- c(tab_files, Order = order_file)
  } else {
    warning("order relative table not found; Order panel skipped: ", order_file)
  }

  # Top-9 phyla by total relative abundance, computed row-wise (no melt).
  # Only the top-5 get bright colours; the 6th-10th phyla are grey.
  phy_tab <- data.table::fread(abund_file, header = TRUE, sep = "\t", skip = 1)
  sample_cols_phy <- setdiff(names(phy_tab), names(phy_tab)[1])
  data.table::setnames(phy_tab, names(phy_tab)[1], "Taxon")
  phy_tab[, parent_phylum := extract_parent_phylum(Taxon)]
  phy_tab[, Total := rowSums(.SD, na.rm = TRUE), .SDcols = sample_cols_phy]
  top9   <- as.character(phy_tab[!is.na(parent_phylum) & parent_phylum != ""]
                         [order(-Total), head(parent_phylum, 9)])
  top5   <- head(top9, 5)
  fill_values <- c(setNames(core_palette[seq_along(top5)], top5),
                   Other = other_color)

  prev_tabs <- lapply(seq_along(tab_files), function(i) {
    lvl  <- names(tab_files)[i]
    sub  <- prevalence_collapsed_table(tab_files[[i]])
    if (nrow(sub) == 0) return(NULL)
    sub <- sub[order(-Count)]
    sub <- sub[seq_len(min(10, nrow(sub)))]
    sub[, Level := lvl]
    sub
  })
  prev_tabs <- rbindlist(prev_tabs, use.names = TRUE)

  if (nrow(prev_tabs) == 0) return("no prevalence data after cleaning")

  # Colour by parent phylum: top-5 phyla get their panel colour, everything
  # else (6th-10th phyla and all non-top-5 classes/orders) is grey.
  prev_tabs[, ColorKey := ifelse(parent_phylum %in% top5, parent_phylum, "Other")]
  prev_tabs[, ColorKey := factor(ColorKey, levels = c(top5, "Other"))]
  # Taxon names sit inside the bars: white on the bright phylum colours, dark
  # grey on the light "Other" bars.
  prev_tabs[, TextColor := ifelse(ColorKey == "Other", "#303030", "white")]

  level_levels <- c("Phylum", "Class", "Order")[c("Phylum", "Class", "Order") %in%
                                                  unique(prev_tabs$Level)]

  make_plot <- function(lvl) {
    sub <- prev_tabs[Level == lvl]
    # Longest bar at the top: on a discrete y-axis ggplot2 draws the FIRST level
    # at the bottom and the LAST at the top, so levels run shortest -> longest.
    sub[, Taxon_final := factor(Taxon_final,
                                levels = sub$Taxon_final[order(sub$Count)])]
    ggplot(sub, aes(x = Count, y = Taxon_final, fill = ColorKey)) +
      geom_col(width = 0.72, color = "#4A4A4A", linewidth = 0.3) +
      geom_text(aes(x = Count / 2, label = Taxon_final, color = TextColor),
                size = base_sz * 0.38, fontface = "bold") +
      scale_color_identity() +
      scale_fill_manual(values = fill_values, name = "Phylum",
                        breaks = c(top5, "Other"),
                        labels = c(top5, "Others")) +
      scale_x_continuous(breaks = seq(0, 80000, by = 20000),
                         labels = function(x) ifelse(x == 0, "0",
                                                     paste0(x / 1000, "K")),
                         expand = expansion(mult = c(0, 0.04))) +
      labs(x = "Samples", y = lvl) +
      theme_classic(base_size = base_sz - 4, base_family = "sans") +
      theme(
        axis.title.y       = element_text(size = base_sz, face = "bold",
                                          color = "#303030", angle = 90,
                                          margin = margin(r = 6)),
        panel.background   = element_rect(fill = "#F8F5F0", color = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "#E0DCD5", linewidth = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.border       = element_rect(color = "#4A4A4A", fill = NA,
                                          linewidth = 0.5),
        axis.line          = element_blank(),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        axis.text.x        = element_text(size = base_sz * 0.75, face = "bold",
                                          color = "#303030"),
        axis.title.x       = element_text(size = base_sz, face = "bold",
                                          color = "#303030"),
        legend.position    = "top",
        legend.title       = element_text(size = base_sz, face = "bold"),
        legend.text        = element_text(size = base_sz * 0.85),
        legend.key.size    = unit(3.2, "mm"),
        plot.background    = element_rect(fill = "transparent", color = NA)
      )
  }

  plots <- lapply(level_levels, make_plot)
  # guides="collect" without an outer override lets the assembled patchwork
  # inherit theme_classic()'s default legend.position="right", which pushes the
  # collected legend into its own right-hand column and squeezes the three
  # panels into ~2/3 of the canvas. The `&` operator forces the assembly-level
  # theme so the single shared legend sits on top and the panels fill the width.
  combined <- (patchwork::wrap_plots(plots, ncol = length(level_levels),
                                     guides = "collect")) &
    theme(legend.position = "top")

  n_facets <- length(level_levels)
  w_out <- if (n_facets == 1) 8 else if (n_facets == 2) 12 else 16
  h_out <- if (n_facets == 1) 3.2 else if (n_facets == 2) 3.6 else 4.5
  ggplot2::ggsave(file.path(out_dir, panel_png["PREV"]), plot = combined,
                  width = w_out, height = h_out, device = "png", dpi = 600, bg = "white")
  invisible(NULL)
}

# ============================================================================
# Panel C : global sampling world map (merged from global_sampling_world_map_enhanced.R)
# ============================================================================
render_panel_C_world_map <- function(out_dir) {
  if (!is.null(err <- load_libs(c("tidyverse", "sf", "countrycode", "rnaturalearth",
                                  "rnaturalearthdata", "ggrepel", "scales",
                                  "patchwork", "ggtext"), AUTO_INSTALL))) return(err)
  if (!file.exists(metadata_path))
    return(sprintf("metadata table not found: %s", metadata_path))

  sf_use_s2(FALSE)

  # base_sz = master font-size knob for panel C.
  base_sz  <- as.numeric(Sys.getenv("FIG1_C_BASE_SZ", unset = "22"))
  text_sz  <- base_sz / .pt

  theme_map <- theme_void(base_size = base_sz, base_family = "sans") +
    theme(
      legend.title   = element_text(size = base_sz * 0.8, face = "bold", color = "#303030"),
      legend.text    = element_text(size = base_sz * 0.7, color = "#303030"),
      legend.position = "top",
      legend.box      = "horizontal",
      legend.spacing.x = unit(3, "mm"),
      legend.key      = element_blank(),
      plot.background  = element_rect(fill = "#F8F5F0", color = NA),
      panel.background = element_rect(fill = "#F8F5F0", color = NA),
      plot.margin      = margin(5, 5, 5, 5, "pt")
    )

  continent_colors <- c(
    "NAM" = "#FDD835", "SA" = "#AB47BC", "EU" = "#4CAF50",
    "AS"  = "#FF9800", "AF" = "#795548", "OA" = "#F48FB1"
  )
  continent_order <- c("NAM", "SA", "EU", "AS", "AF", "OA")

  site_cols_ordered <- c("Saliva", "Tooth", "Oral Mucosa", "Oral Swab", "Tongue", "Gingiva", "Others")
  site_colors <- c(
    "Saliva"      = "#08306B",
    "Tooth"       = "#E7298A",
    "Oral Mucosa" = "#666666",
    "Oral Swab"   = "#E6AB02",
    "Tongue"      = "#A6761D",
    "Gingiva"     = "#CAB2D6",
    "Others"      = "#117733"
  )

  # Read and clean metadata. Always prefer the *_cleaned columns; fall back to
  # the raw columns only when no clean variant exists in this metadata file.
  md_raw <- read_tsv(metadata_path, show_col_types = FALSE)
  col_site    <- pick_meta_col(md_raw, "site_cleaned",    "site")
  col_country <- pick_meta_col(md_raw, "country_cleaned", "country")
  col_age     <- pick_meta_col(md_raw, "age_cleaned",     "age")
  col_bmi     <- pick_meta_col(md_raw, "bmi_cleaned",     "bmi")
  col_sex     <- pick_meta_col(md_raw, "sex_cleaned",     "sex")

  metadata <- tibble::tibble(
    BioProject      = md_raw[["bioproject_id"]],
    subject_id      = md_raw[[pick_meta_col(md_raw, "subject_id_cleaned", "subject_id")]],
    Country_Cleaned = md_raw[[col_country]],
    Site_Cleaned    = md_raw[[col_site]]
  ) %>%
    mutate(
      Site_Cleaned = str_to_title(Site_Cleaned),
      Site_Cleaned = case_when(Site_Cleaned == "Other" ~ "Others", TRUE ~ Site_Cleaned)
    ) %>%
    filter(!(Site_Cleaned %in% c("Nasal", "Tonsils")), !is.na(Country_Cleaned),
           Country_Cleaned != "") %>%
    mutate(iso3c = countrycode(Country_Cleaned, origin = "country.name", destination = "iso3c"))

  south_america_iso3c <- c("ARG", "BOL", "BRA", "CHL", "COL", "ECU", "GUY", "PRY",
                           "PER", "SUR", "URY", "VEN", "FLK", "GUF")

  meta_clean <- metadata %>%
    mutate(
      Temp_Continent = countrycode(iso3c, "iso3c", "continent"),
      Continent = case_when(
        Temp_Continent == "Americas" & (iso3c %in% south_america_iso3c) ~ "SA",
        Temp_Continent == "Americas" & !(iso3c %in% south_america_iso3c) ~ "NAM",
        Temp_Continent == "Europe"   ~ "EU",
        Temp_Continent == "Asia"     ~ "AS",
        Temp_Continent == "Africa"   ~ "AF",
        Temp_Continent == "Oceania"  ~ "OA",
        iso3c == "ATA"               ~ "Antarctica",
        TRUE                         ~ Temp_Continent
      ),
      age_num = suppressWarnings(as.numeric(md_raw[[col_age]])),
      bmi_num = suppressWarnings(as.numeric(md_raw[[col_bmi]])),
      sex_label = case_when(
        tolower(md_raw[[col_sex]]) %in% c("f", "female") ~ "F",
        tolower(md_raw[[col_sex]]) %in% c("m", "male")   ~ "M",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!Continent %in% c("Other", "Antarctica")) %>%
    mutate(Continent = factor(Continent, levels = continent_order))

  # Continent-level summary
  continent_summary <- meta_clean %>%
    group_by(Continent) %>%
    summarise(
      n_countries = n_distinct(Country_Cleaned),
      n_bp   = n_distinct(BioProject),
      n_subj = n_distinct(subject_id),
      n_samp = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Full_Name = case_when(
        Continent == "NAM" ~ "North America (NA)",
        Continent == "SA"  ~ "South America (SA)",
        Continent == "EU"  ~ "Europe (EU)",
        Continent == "AS"  ~ "Asia (AS)",
        Continent == "AF"  ~ "Africa (AF)",
        Continent == "OA"  ~ "Oceania (OA)"
      ),
      card_label = paste0(
        "<span style='font-size: ", round(base_sz * 0.70), "pt;'><b>", Full_Name, "</b></span><br>",
        "<span style='font-size: ", round(base_sz * 0.52), "pt;'>",
        "Countries = ", n_countries, "<br>",
        "Bioprojects = ", n_bp, "<br>",
        "Subjects = ", scales::comma(n_subj), "<br>",
        "Samples = ", scales::comma(n_samp),
        "</span>"
      )
    )

  # Absolute safe-zone coordinates for country info boxes
  continent_coords <- data.frame(
    Continent = c("NAM", "SA", "EU", "AS", "AF", "OA"),
    lon_pie   = c(-175, -135,  -40,  145,   -40,   75),
    lat_pie   = c(  30,  -30,   55,   10,  -40,  -35)
  )
  continent_summary <- continent_summary %>%
    left_join(continent_coords, by = "Continent") %>%
    mutate(
      lon_lab = if_else(
        Continent == "EU",
        lon_pie[Continent == "AF"] + 9.5,
        lon_pie + 9.5
      )
    )

  # Country-level summary (BioProject count per country, for gray points)
  country_summary <- meta_clean %>%
    group_by(Country_Cleaned, Continent) %>%
    summarise(
      n_samples = n(),
      n_studies = n_distinct(BioProject),
      .groups = "drop"
    ) %>%
    mutate(iso3c = countrycode(Country_Cleaned, "country.name", "iso3c"))

  world_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") %>%
    filter(name != "Antarctica") %>%
    mutate(iso3c = iso_a3) %>%
    mutate(
      Temp_Continent = countrycode(iso3c, "iso3c", "continent"),
      Continent = case_when(
        Temp_Continent == "Americas" & (iso3c %in% south_america_iso3c) ~ "SA",
        Temp_Continent == "Americas" & !(iso3c %in% south_america_iso3c) ~ "NAM",
        Temp_Continent == "Europe"   ~ "EU",
        Temp_Continent == "Asia"     ~ "AS",
        Temp_Continent == "Africa"   ~ "AF",
        Temp_Continent == "Oceania"  ~ "OA",
        name == "Somaliland"         ~ "AF",
        name == "N. Cyprus"          ~ "AS",
        name == "Kosovo"             ~ "EU",
        TRUE                         ~ NA_character_
      )
    ) %>%
    select(-Temp_Continent)

  # Randomly distributed gray points: one per BioProject within each country
  set.seed(42)
  bp_points <- lapply(seq_len(nrow(country_summary)), function(i) {
    cs <- country_summary[i, ]
    idx <- which(world_sf$iso3c == cs$iso3c)
    if (length(idx) == 0) return(NULL)
    pts <- tryCatch(st_sample(world_sf$geometry[idx], size = cs$n_studies), error = function(e) NULL)
    if (is.null(pts) || length(pts) == 0) {
      if (cs$iso3c == "BRB") {
        return(data.frame(
          iso3c = "BRB", lon = rep(-59.5, cs$n_studies),
          lat = rep(13.2, cs$n_studies), n_studies = cs$n_studies
        ))
      }
      return(NULL)
    }
    coords <- st_coordinates(pts)
    data.frame(
      iso3c = cs$iso3c,
      lon   = coords[, 1],
      lat   = coords[, 2],
      n_studies = cs$n_studies
    )
  }) %>%
    bind_rows()

  # Pie charts per continent as grobs
  site_continent_data <- meta_clean %>%
    filter(!is.na(Site_Cleaned)) %>%
    count(Continent, Site_Cleaned) %>%
    group_by(Continent) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(Site_Cleaned = factor(Site_Cleaned, levels = site_cols_ordered))

  pie_grobs <- lapply(continent_order, function(cont) {
    df <- site_continent_data %>% filter(Continent == cont)
    if (nrow(df) == 0) return(grid::nullGrob())
    p <- ggplot(df, aes(x = "", y = prop, fill = Site_Cleaned)) +
      geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = site_colors, limits = site_cols_ordered) +
      theme_void() +
      theme(legend.position = "none", plot.background = element_blank(),
            panel.background = element_blank())
    ggplotGrob(p)
  })

  r_pie <- 7.5
  pie_layers <- lapply(seq_along(continent_order), function(i) {
    cont <- continent_order[i]
    coords <- continent_coords %>% filter(Continent == cont)
    annotation_custom(
      grob = pie_grobs[[i]],
      xmin = coords$lon_pie - r_pie, xmax = coords$lon_pie + r_pie,
      ymin = coords$lat_pie - r_pie, ymax = coords$lat_pie + r_pie
    )
  })

  # Core map rendering
  p_map_base <- ggplot() +
    geom_sf(data = world_sf, aes(fill = Continent), color = "#D5D0C5", linewidth = 0.3) +
    scale_fill_manual(values = scales::alpha(continent_colors, 0.6),
                      na.value = "#EFEBDD", guide = "none") +
    geom_point(data = bp_points, aes(x = lon, y = lat, size = n_studies),
               shape = 16, color = "#0D47A1", alpha = 1) +
    scale_size_continuous(
      name   = "Number of samples",
      range  = c(2, 6),
      limits = c(0, 100),
      breaks = c(25, 50, 100),
      labels = c("25", "50", "100"),
      guide  = guide_legend(override.aes = list(color = "#0D47A1", alpha = 1,
                                                fill = "transparent", shape = 16))
    ) +
    pie_layers +
    geom_richtext(
      data = continent_summary,
      aes(x = lon_lab, y = lat_pie, label = card_label),
      hjust = 0, vjust = 0.5,
      fill = "#F8F5F0", color = "#303030",
      label.color = "grey75",
      label.padding = unit(0.4, "lines"),
      label.r = unit(0.15, "lines"),
      lineheight = 1.3,
      alpha = 0.95
    ) +
    coord_sf(xlim = c(-180, 180), ylim = c(-55, 75), expand = FALSE, clip = "off") +
    theme_map

  # Transparent legend
  p_site_legend <- ggplot(data.frame(Site = factor(site_cols_ordered, levels = site_cols_ordered)),
                          aes(x = 1, y = 1, fill = Site)) +
    geom_point(shape = 21, size = 0, color = NA, stroke = 0) +
    scale_fill_manual(
      name = "Sampling Site",
      values = site_colors,
      limits = site_cols_ordered,
      guide = guide_legend(override.aes = list(size = 5, color = "white",
                                               stroke = 0.5, alpha = 1))
    ) +
    theme_void() +
    theme(
      legend.position = "left",
      legend.title = element_text(size = base_sz * 0.85, face = "bold", color = "#303030"),
      legend.text = element_text(size = base_sz * 0.75, color = "#303030"),
      legend.key.size = unit(5, "mm"),
      legend.background = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )

  # Composite map
  p_final_map <- p_map_base +
    inset_element(
      p_site_legend,
      left = 0.01, bottom = 0.05, right = 0.20, top = 0.35,
      align_to = "full", clip = FALSE
    ) +
    plot_annotation(
      theme = theme(
        plot.background  = element_rect(fill = "#F8F5F0", color = NA),
        panel.background = element_rect(fill = NA, color = NA)
      )
    )

  ggplot2::ggsave(file.path(out_dir, panel_png["C"]), plot = p_final_map,
                  width = 18, height = 7.2, device = "png", dpi = 600, bg = "#F8F5F0")
  invisible(NULL)
}

# ============================================================================
# Panel D : phylum-level relative abundance, stacked bars by continent x site
#   (merged from backup/Figure1_composite_v2.R, "Panel c")
#   X = continent (NA/SA/EU/AS/AF/OA), faceted by site, fill = top-9 phyla + Others
# ============================================================================
render_panel_D_phylum_continent <- function(out_dir) {
  if (!is.null(err <- load_libs(c("tidyverse", "data.table"), AUTO_INSTALL))) return(err)
  if (!file.exists(abund_file) || !file.exists(metadata_path))
    return(sprintf("inputs not found (need %s and %s)", abund_file, metadata_path))

  # base_sz = master font-size knob for panel D.
  base_sz <- as.numeric(Sys.getenv("FIG1_D_BASE_SZ", unset = "26"))

  metadata_raw <- readr::read_tsv(metadata_path, show_col_types = FALSE)
  rel_dt <- data.table::fread(abund_file, header = TRUE, sep = "\t", skip = 1)
  data.table::setnames(rel_dt, 1, "Phylum")
  phylum_cols <- setdiff(names(rel_dt), "Phylum")
  rel_long <- data.table::melt(rel_dt, id.vars = "Phylum", measure.vars = phylum_cols,
                               variable.name = "sample_id", value.name = "Abundance")
  rel_long[, sample_id := as.character(sample_id)]
  rel_long[, Phylum_clean := gsub("^d__Bacteria;\\s*p__", "", Phylum)]
  rel_long[, Phylum_clean := gsub("^p__", "", Phylum_clean)]
  rel_long <- rel_long[!is.na(Phylum_clean) & Phylum_clean != ""]

  col_site <- pick_meta_col(metadata_raw, "site_cleaned", "site")
  meta_sub <- as.data.table(metadata_raw)[, .(sample_id, x_site = get(col_site),
                                              continent_group)]
  data.table::setnames(meta_sub, "x_site", "site_cleaned")
  meta_sub[, sample_id := as.character(sample_id)]
  merged_dt <- merge(rel_long, meta_sub, by = "sample_id", all.x = FALSE)

  target_continents    <- c("Africa", "Asia", "Europe", "North America",
                            "Oceania", "South America")

  merged_dt <- merged_dt[!is.na(site_cleaned) & site_cleaned != ""]
  merged_dt <- merged_dt[!is.na(continent_group) & continent_group != "" &
                         continent_group != "NA"]
  merged_dt[, site_cleaned    := tolower(trimws(as.character(site_cleaned)))]
  merged_dt[, continent_group := trimws(as.character(continent_group))]
  # Show every oral site category present (ordered by sample count); the
  # "Other"/unclassified category is left out of the figure.
  merged_dt <- merged_dt[!(site_cleaned %in% c("other", "others", "other",
                                               "others (specify)", "other (specify)"))]
  merged_dt <- merged_dt[continent_group %in% target_continents]
  site_targets   <- merged_dt[, .N, by = site_cleaned][order(-N)]$site_cleaned
  site_labels    <- stringr::str_to_title(site_targets)
  site_label_map <- setNames(site_labels, site_targets)
  merged_dt[, site_display := site_label_map[site_cleaned]]
  merged_dt[, site_display := factor(site_display, levels = site_labels)]

  # One shared global-top-9 / color map for every phylum-colored panel.
  rel_long[, parent_phylum := Phylum_clean]
  global_top9    <- global_top_phyla(rel_long, 9)
  unified_colors <- phylum_color_map(global_top9)
  merged_dt[, Phylum_group := ifelse(Phylum_clean %in% global_top9,
                                     Phylum_clean, "Other")]
  plot_data <- merged_dt[, .(Mean_Abundance = mean(Abundance, na.rm = TRUE)),
                         by = .(site_display, continent_group, Phylum_group)]
  plot_data[, Mean_Abundance := Mean_Abundance / sum(Mean_Abundance),
            by = .(site_display, continent_group)]

  # Within each site facet, order continents so the dominant phylum decreases from
  # left to right; the dominant phylum is drawn at the bar bottom.
  dominant_phylum <- global_top9[1]
  combined_levels <- character(0)
  for (s in site_labels) {
    sub_c <- plot_data[site_display == s & Phylum_group == dominant_phylum]
    if (nrow(sub_c) == 0) {
      ordered_conts <- intersect(target_continents,
                                 unique(plot_data[site_display == s]$continent_group))
    } else {
      ordered_conts <- sub_c[order(-Mean_Abundance)]$continent_group
    }
    combined_levels <- c(combined_levels, paste(s, ordered_conts, sep = "|"))
  }
  plot_data[, facet_var := factor(paste(site_display, continent_group, sep = "|"),
                                  levels = combined_levels)]
  strip_facet_prefix <- function(x) {
    raw_names <- sapply(x, function(lab) gsub("^[^|]*\\|", "", as.character(lab)))
    abbr_map <- c("North America" = "NA", "South America" = "SA",
                  "Europe" = "EU", "Asia" = "AS",
                  "Africa" = "AF", "Oceania" = "OA")
    unname(abbr_map[raw_names])
  }

  other_phyla     <- setdiff(global_top9, dominant_phylum)
  # Stack order: "Other" at the very top, remaining phyla in descending
  # abundance, dominant phylum at the bar bottom. ggplot2 draws the FIRST factor
  # level at the TOP of a stacked bar, so levels are listed top-to-bottom.
  stack_order     <- c("Other", rev(other_phyla), dominant_phylum)
  plot_data[, Phylum_group := factor(Phylum_group, levels = stack_order)]

  p_d <- ggplot(plot_data, aes(x = facet_var, y = Mean_Abundance, fill = Phylum_group)) +
    geom_bar(stat = "identity", position = "stack", width = 0.90) +
    facet_grid(. ~ site_display, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = unified_colors, name = "Phylum",
                      breaks = c(dominant_phylum, other_phyla, "Other"),
                      labels = c(dominant_phylum, other_phyla, "Other")) +
    scale_x_discrete(labels = strip_facet_prefix) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Relative Abundance") +
    guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
    theme_classic(base_size = base_sz - 4, base_family = "sans") +
    theme(
      axis.text.x      = element_text(size = base_sz - 6, face = "bold", color = "#303030",
                                      angle = 30, hjust = 1, vjust = 1),
      axis.text.y      = element_text(size = base_sz - 5, face = "bold", color = "#303030"),
      axis.title.y     = element_text(size = base_sz - 4, face = "bold", color = "#303030",
                                      margin = margin(r = 4)),
      legend.position  = "top",
      legend.title     = element_text(size = base_sz - 6, face = "bold", color = "#303030"),
      legend.text      = element_text(size = base_sz - 7, color = "#303030"),
      legend.key.size  = unit(3.0, "mm"),
      legend.margin    = margin(b = -5),
      strip.text       = element_text(size = base_sz - 8, face = "bold", color = "#303030",
                                      margin = margin(b = 4, t = 4)),
      strip.background = element_rect(fill = "#ECEFF1", color = NA),
      panel.spacing    = unit(0.5, "lines"),
      axis.line        = element_line(color = "#303030", linewidth = 0.3),
      axis.ticks       = element_line(color = "#303030", linewidth = 0.2),
      plot.margin      = margin(2, 2, 2, 2, unit = "pt"),
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )

  ggplot2::ggsave(file.path(out_dir, panel_png["CONT"]), plot = p_d,
                  width = 20, height = 5.7, device = "png", dpi = 600, bg = "white")
  invisible(NULL)
}

# ============================================================================
# Panel E (supplementary) : phylum-level relative abundance, stacked bars by
#   age group x site. Shares the global top-9 and color map with the other
#   phylum-colored panels (D, E and the phylum x continent supp panel).
#   and the other stacked-bar panel. X = age group, faceted by site.
render_panel_E_phylum_age <- function(out_dir) {
  if (!is.null(err <- load_libs(c("tidyverse", "data.table"), AUTO_INSTALL))) return(err)
  if (!file.exists(abund_file) || !file.exists(metadata_path))
    return(sprintf("inputs not found (need %s and %s)", abund_file, metadata_path))

  # base_sz = master font-size knob for panel E (aligned with panel D).
  base_sz <- as.numeric(Sys.getenv("FIG1_E_BASE_SZ", unset = "26"))

  rel <- data.table::fread(abund_file, header = TRUE, sep = "\t", skip = 1)
  data.table::setnames(rel, 1, "Phylum")
  phylum_cols <- setdiff(names(rel), "Phylum")
  rel_long <- data.table::melt(rel, id.vars = "Phylum", measure.vars = phylum_cols,
                               variable.name = "sample_id", value.name = "Abundance")
  rel_long[, sample_id := as.character(sample_id)]
  rel_long[, Phylum_clean := gsub("^d__Bacteria;\\s*p__", "", Phylum)]
  rel_long[, Phylum_clean := gsub("^p__", "", Phylum_clean)]
  rel_long <- rel_long[!is.na(Phylum_clean) & Phylum_clean != ""]

  meta <- data.table::fread(metadata_path, header = TRUE, sep = "\t")
  col_site <- pick_meta_col(meta, "site_cleaned",       "site")
  col_ageg <- pick_meta_col(meta, "age_cleaned_group",  "age_group")
  required_cols <- c("sample_id", col_site, col_ageg)
  if (!all(required_cols %in% names(meta)))
    return(sprintf("metadata missing required columns: %s",
                   paste(setdiff(required_cols, names(meta)), collapse = ", ")))
  meta_sub <- meta[, .(sample_id, x_site = get(col_site), x_ageg = get(col_ageg))]
  data.table::setnames(meta_sub, c("x_site", "x_ageg"),
                       c("site_cleaned", "age_cleaned_group"))
  meta_sub[, sample_id := as.character(sample_id)]
  merged <- merge(rel_long, meta_sub, by = "sample_id", all.x = FALSE)

  age_order            <- c("0-6 months", "6-36 months", "3-12 years",
                            "12-18 years", "18-60 years", ">60 years")

  merged <- merged[!is.na(site_cleaned) & site_cleaned != ""]
  merged <- merged[!is.na(age_cleaned_group) & age_cleaned_group != ""]
  merged[, site_cleaned := tolower(trimws(as.character(site_cleaned)))]
  merged[, age_cleaned_group := gsub("＞", ">", age_cleaned_group)]
  # Show every oral site category present (ordered by sample count); the
  # "Other"/unclassified category is left out of the figure.
  merged <- merged[!(site_cleaned %in% c("other", "others", "other",
                                         "others (specify)", "other (specify)"))]
  merged <- merged[age_cleaned_group %in% age_order]
  if (nrow(merged) == 0)
    return("no data left after site/age filtering for panel E")

  site_targets   <- merged[, .N, by = site_cleaned][order(-N)]$site_cleaned
  site_labels    <- stringr::str_to_title(site_targets)
  site_label_map <- setNames(site_labels, site_targets)
  merged[, site_display := site_label_map[site_cleaned]]
  merged[, site_display := factor(site_display, levels = site_labels)]
  merged[, age_cleaned_group := factor(age_cleaned_group, levels = age_order)]

  rel_long[, parent_phylum := Phylum_clean]
  global_top9    <- global_top_phyla(rel_long, 9)
  unified_colors <- phylum_color_map(global_top9)
  merged[, Phylum_group := ifelse(Phylum_clean %in% global_top9, Phylum_clean, "Other")]
  plot_data <- merged[, .(Mean_Abundance = mean(Abundance, na.rm = TRUE)),
                      by = .(site_display, age_cleaned_group, Phylum_group)]
  plot_data[, Mean_Abundance := Mean_Abundance / sum(Mean_Abundance),
            by = .(site_display, age_cleaned_group)]

  dominant_phylum <- global_top9[1]
  other_phyla     <- setdiff(global_top9, dominant_phylum)
  # Stack order: "Other" at the very top, remaining phyla in descending
  # abundance, dominant phylum at the bar bottom (matching the continent panel).
  # ggplot2 draws the FIRST factor level at the TOP of a stacked bar.
  stack_order     <- c("Other", rev(other_phyla), dominant_phylum)
  plot_data[, Phylum_group := factor(Phylum_group, levels = stack_order)]

  p_e <- ggplot(plot_data, aes(x = age_cleaned_group, y = Mean_Abundance,
                               fill = Phylum_group)) +
    geom_bar(stat = "identity", position = "stack", width = 0.95) +
    facet_grid(. ~ site_display, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = unified_colors, name = "Phylum",
                      breaks = c(dominant_phylum, other_phyla, "Other"),
                      labels = c(dominant_phylum, other_phyla, "Other")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Age Group", y = "Relative Abundance") +
    theme_classic(base_size = base_sz - 4, base_family = "sans") +
    theme(
      axis.text.x      = element_text(size = base_sz - 6, face = "bold", color = "#303030",
                                      angle = 45, hjust = 1, vjust = 1),
      axis.text.y      = element_text(size = base_sz - 5, face = "bold", color = "#303030"),
      axis.title.x     = element_text(size = base_sz - 4, face = "bold", color = "#303030",
                                      margin = margin(t = 6)),
      axis.title.y     = element_text(size = base_sz - 4, face = "bold", color = "#303030",
                                      margin = margin(r = 4)),
      legend.position  = "right",
      legend.title     = element_text(size = base_sz - 6, face = "bold", color = "#303030"),
      legend.text      = element_text(size = base_sz - 7, color = "#303030"),
      legend.key.size  = unit(3.0, "mm"),
      strip.text       = element_text(size = base_sz - 8, face = "bold", color = "#303030",
                                      margin = margin(b = 4, t = 4)),
      strip.background = element_rect(fill = "#ECEFF1", color = NA),
      panel.spacing    = unit(0.5, "lines"),
      axis.line        = element_line(color = "#303030", linewidth = 0.3),
      axis.ticks       = element_line(color = "#303030", linewidth = 0.2),
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )

  ggplot2::ggsave(file.path(out_dir, panel_png["AGE"]), plot = p_e,
                  width = 16, height = 5.7, device = "png", dpi = 600, bg = "white")
  invisible(NULL)
}

# ============================================================================
# Main
# ============================================================================
render_all_panels <- function() {
  if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)
  jobs <- list(
    A     = list(fn = render_panel_A_funnel,            name = "pipeline funnel"),
    B     = list(fn = render_panel_B_rarefaction,       name = "rarefaction curve"),
    C     = list(fn = render_panel_C_world_map,         name = "global sampling map",
                 quiet = TRUE),
    DEPTH = list(fn = render_panel_depth,               name = "sequencing-depth histogram"),
    PREV  = list(fn = render_panel_prevalence,          name = "core-taxa prevalence"),
    DENS  = list(fn = render_panel_density,             name = "top-5 phylum density"),
    SHAN  = list(fn = render_panel_shannon,             name = "alpha diversity (Shannon)"),
    CONT  = list(fn = render_panel_D_phylum_continent,  name = "phylum x continent stacked (supp)"),
    AGE   = list(fn = render_panel_E_phylum_age,        name = "phylum x age stacked (supp)")
  )
  for (k in names(jobs)) {
    j <- jobs[[k]]
    cat("== rendering panel", k, "-", j$name, "... ")
    # quiet suppresses the harmless sf "planar"/countrycode/geom_point
    # warnings of the map panel (errors still surface via tryCatch).
    run <- function() {
      if (isTRUE(j$quiet)) suppressMessages(suppressWarnings(j$fn(PANEL_DIR)))
      else j$fn(PANEL_DIR)
    }
    ans <- tryCatch(run(), error = function(e) paste("error:", conditionMessage(e)))
    cat(if (is.null(ans)) "OK\n" else paste("SKIPPED:", ans, "\n", sep = " "))
  }
  invisible(NULL)
}

cat("Script dir :", script_dir, "\n")
cat("BASE_DIR   :", BASE_DIR, "\n")
cat("PANEL_DIR  :", PANEL_DIR, "\n")
cat("\n== Regenerating panels ==\n")
render_all_panels()
cat("\nDone. Panels written to:", PANEL_DIR, "\n")
cat("Now assemble the composite (fast, repeatable) with:\n")
cat("  Rscript figure1_composite_all_in_one_0907.R\n")
