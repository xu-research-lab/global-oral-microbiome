#!/usr/bin/env Rscript
# pipeline_v3: split old supported-catalog profile (12i) into TWO figures
#   Fig A: multi-member clusters (one row per CLUSTER)
#          y labels = species annotation, phylum panel removed,
#          continuous annotation-distance heatmap added,
#          rightmost = boxplot + jitter of per-ASV subject prevalence,
#          sorted by max subject prevalence (desc)
#   Fig B: singleton top N by subject prevalence (one row per ASV,
#          rightmost = subject prevalence bar chart, sorted desc)
#   In addition, Fig A and Fig B are combined into ONE row figure
#   (aplot::plot_list) with a shared legend area at the bottom.
# Layout: left label column (right-aligned species names, wide enough for the
#         longest annotation) + tile panels; colour scales are forced to the
#         same limits in both figures so aplot can collect their legends.
suppressPackageStartupMessages({
  library(ggplot2)
  library(aplot)
  library(cowplot)    # get_legend(): the two-line legend strip of Fig A
  library(patchwork)  # wrap_elements(): stacking that strip above the row
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

MIN_MEMBERS <- 2   # >=2 = multi-member (interpreted); set to 3 for strictly >2
TOP_N       <- 40  # singleton top N by subject prevalence

# physical column widths (inches) used for both the two standalone figures
# and the final side-by-side figure
# The canvas is drawn at the size D is finally shown at, so that the type keeps
# its point size instead of being shrunk with the image: these widths are the
# 7.64-in canvas divided by 1.4.  The legends, whose type also does not shrink,
# sit on two lines instead of three.
# A's label column has to hold the longest annotation on the canvas, e.g.
# "Pseudoalteromonas luteoviolacea_G" (1.9 in at the row-label size), so it is
# widened and the tracks narrowed by the same amount: the canvas width is held
# fixed, otherwise D (displayed at a fixed page width) would shrink its type.
# The label column also carries the cluster ID after the annotation, which adds
# up to ten characters (" (c15805)"): the longest row of A is 45 characters
# ("Fructilactobacillus sanfranciscensis (c21350)") and of B 41, so both columns
# are widened enough to hold them.  The column is wider than that needs at the
# old row-label size because the labels are set larger now: the legend went from
# three lines to two, the space it gave up went into the rows, and the rows carry
# bigger type -- but the type widens the label column too, and at a canvas this
# narrow that is what actually limits the size.
LAB_W   <- c(A = 3.50, B = 2.60)   # annotation + cluster ID label columns
# the prevalence panel carries three tick labels on a log axis and needs most of
# an inch for them; the tile columns give up what it gains, so the canvas width
# (and with it the size the panel's type ends up at on the page) is unchanged
TILE_W  <- 0.177                    # width of ONE tile column / ONE site column
PREV_W  <- 0.95                     # rightmost prevalence panel
GAP_W   <- 0.21                     # blank gap between Fig A and Fig B

# geom_text size is mm, axis.text size is pt.  The row labels set the size of the
# label column, so they are the one that is tied to it; the labels under the tile
# tracks do not get the same treatment -- each stands in a 0.18-in column, and at
# the row-label size they ran into one another along the diagonal.
ROW_LAB_SZ_MM  <- 3.45
AXIS_LAB_SZ_PT <- ROW_LAB_SZ_MM * 72.27 / 25.4   # ~9.8 pt, prevalence panels
TILE_XLAB_SZ_PT <- 8.3                           # track names under the tiles
LEGEND_W       <- 5.30             # min canvas width that fits the top legend strip
                                   # (= sum(widths) exactly; the legend strip is
                                   # tightened below so that it fits inside it)

# shared colour-scale limits for panels that appear in BOTH figures:
# identical limits/breaks are required for aplot/patchwork to merge legends
DIST_RANGE <- c(1e-6, 2)       # annotation distance (log10 scale)
BP_RANGE   <- c(0.2, 2.5)      # log10 n_bioprojects
CONT_RANGE <- c(2, 6)          # n_continents
HOMD_RANGE <- c(40, 100)       # HOMD identity (%)
ASV_RANGE  <- c(0.3, 2.1)      # Fig A only: log10 n_derep_asvs

split_values <- function(x) {
  if (is.na(x) || x == "" || x == "NA" || x == "nan") return(character())
  vals <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
  vals[vals != "" & vals != "NA" & vals != "nan"]
}
nvals <- function(x) length(split_values(x))
# Row label: the annotation, then the cluster ID in brackets.  A cluster the
# annotation step left blank still has to be nameable, so it is called
# "unannotated" rather than dropped: the ID is what the tables are keyed on.
row_label <- function(ann, id) {
  a <- as.character(ann)
  a[is.na(a) | a == ""] <- "unannotated"
  sprintf("%s (%s)", a, as.character(id))
}
clean <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA" | x == "nan" | x == "Unassigned"] <- NA_character_
  x
}

message("[12j_v3] loading anchors + master ASV annotation + headline profile")
anchors <- read.delim(file.path(PIPE_V3, "12d_anchors_supported_v3.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
master  <- read.delim(file.path(PIPE_V3, "07_master_asv_table.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
prof    <- read.delim(file.path(PIPE_V3, "12d_headline_profile_v3.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

mcols <- c("derep_id", "annotation", "annotation_distance", "homd_identity", "annotation_in_homd")
master <- master[, intersect(mcols, colnames(master))]
master$derep_id <- as.character(master$derep_id)
prof$derep_id <- as.character(prof$derep_id)
anchors$anchor_tip <- as.character(anchors$anchor_tip)
prof <- merge(prof, master[, c("derep_id", "annotation_in_homd")], by = "derep_id", all.x = TRUE)
prof$annotation_in_homd[is.na(prof$annotation_in_homd)] <- FALSE

# ---- shared track styling --------------------------------------------------
# One neutral low end for every sequential track, so the heatmap reads as one
# family instead of seven different off-whites, and one high-end hue per track.
LOW <- "#F5F5F5"
TRACK <- list(
  dist   = list(track = "distance to known tips", name = "Distance to known tips",
                high = "#B2182B"),
  asv    = list(track = "ASV number", name = "ASV number", high = "#3B6EA5"),
  bp     = list(track = "Bioprojects", name = "Bioprojects",
                high = "#1B7F79"),
  cont   = list(track = "Continents", name = "Continents", high = "#E08214"),
  homd   = list(track = "eHOMD identity", name = "eHOMD identity (%)",
                high = "#1A9850"),
  tax    = list(track = "taxonomy in eHOMD", name = "Taxonomy in eHOMD",
                yes = "#8ECAE6"),
  site   = list(track = "sites", name = "Sites", yes = "#4A7C59"))
FLAG_OFF <- "#F2F2F2"

cb_guide <- function(width = 1.35) {
  # no frame at all: the bar is a gradient and the ticks mark its ends, so the
  # box added nothing.  A rule this size can only be one rendered pixel wide,
  # which at this scale still read heavier than the bar it enclosed.
  guide_colourbar(barwidth = unit(width, "cm"), barheight = unit(0.44, "cm"),
                  title.position = "top", title.hjust = 0.5,
                  ticks.colour = "grey55", ticks.linewidth = 0.1,
                  frame.colour = NA, frame.linewidth = 0)
}
lg_guide <- function(ncol = 1, nrow = NULL, ksize = 0.38) {
  guide_legend(ncol = ncol, keywidth = unit(ksize, "cm"), keyheight = unit(ksize, "cm"),
               nrow = nrow, override.aes = list(size = 0.1))
}

mk_theme_h <- function(base_sz, xlab_sz) {
  theme_minimal(base_size = base_sz) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = xlab_sz),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(2, 4, 2, 2))
}
mk_tile <- function(d, fill_var, scale_fun, base_sz = 8,
                    xlab_sz = TILE_XLAB_SZ_PT) {
  th <- mk_theme_h(base_sz, xlab_sz)
  p <- ggplot(d, aes(y = rank_f, x = track, fill = !!sym(fill_var))) +
    geom_tile(width = 1, height = 1) +
    scale_fun +
    th +
    theme(axis.text.y = element_blank(),
          legend.position = "top",
          legend.box = "horizontal",
          legend.title = element_text(size = 7.5),
          legend.text = element_text(size = 7, margin = margin(l = 2, r = 2)),
          legend.key.size = unit(0.36, "cm"),
          legend.spacing.x = unit(0.30, "cm"),
          legend.margin = margin(0, 0, 0, 0))
  p
}
# `title` / `subtitle` may be NULL, for figures whose heading is carried by the
# panel letter in the composite instead
mk_label_col <- function(d, title, subtitle, ylab_sz) {
  ggplot(d, aes(y = rank_f, x = 1, label = label)) +
    geom_text(hjust = 1, vjust = 0.5, size = ylab_sz, colour = "grey15") +
    scale_x_continuous(limits = c(-0.12, 1.05), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 6.2, lineheight = 1.05),
          plot.margin = margin(2, 2, 2, 2))
}

save_fig <- function(p, base, width, height) {
  ggsave(paste0(base, ".pdf"), p, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(paste0(base, ".png"), p, width = width, height = height, dpi = 300,
         limitsize = FALSE)
}

# Legends have to be *collected* into a strip above the row: a tile panel is only
# ~0.3 in wide, so a legend left inside its own panel has nowhere to go and
# spills out to the right of the figure.
assemble_row <- function(lab, panels, lab_w, widths) {
  p <- do.call(aplot::plot_list,
               c(list(lab), panels,
                 list(ncol = length(panels) + 1,
                      widths = c(lab_w, widths), guides = "collect")))
  p & theme(legend.position = "top",
            legend.box = "horizontal",
            legend.title = element_text(size = 7),
            legend.text = element_text(size = 6.5),
            legend.key.size = unit(0.32, "cm"),
            legend.spacing.x = unit(0.15, "cm"))
}

# ============================================================
# Fig A : multi-member clusters, one row per cluster
# ============================================================
build_cluster_figure <- function(anc, prof_all) {
  anc <- anc[anc$n_members >= MIN_MEMBERS & !is.na(anc$n_members), ]
  rep <- merge(anc, master, by.x = "anchor_tip", by.y = "derep_id", all.x = TRUE)
  stopifnot(nrow(rep) == nrow(anc))
  for (cc in c("n_members", "max_subject_prevalence", "n_bioprojects_total",
               "n_continents", "homd_identity", "annotation_distance")) {
    rep[[cc]] <- suppressWarnings(as.numeric(rep[[cc]]))
  }
  rep$annotation <- clean(rep$annotation)
  rep$annotation_in_homd[is.na(rep$annotation_in_homd)] <- FALSE
  rep$label <- row_label(rep$annotation, rep$cluster_id)

  rep <- rep[order(-rep$max_subject_prevalence, rep$cluster_id), ]
  rep$rank <- seq_len(nrow(rep))
  rep$rank_f <- factor(rep$rank, levels = rev(rep$rank))
  n <- nrow(rep)
  # Row spacing, not a height: the canvas is as tall as the ranks need.  The row
  # pitch is what the row labels need at ROW_LAB_SZ_MM (a label is ~0.9 of its
  # size in ink, so the pitch is one size plus a little lead), and 1.83 in covers
  # the two-line legend and the 45-degree site labels underneath.
  fig_h <- max(3.2, 1.83 + n * 0.158)

  # ---- panels ----
  dist_df <- data.frame(rank_f = rep$rank_f, track = TRACK$dist$track,
                        value = pmax(rep$annotation_distance, 1e-6))
  asv_df <- data.frame(rank_f = rep$rank_f, track = TRACK$asv$track,
                       value = log10(pmax(rep$n_members, 1)))
  bp_df  <- data.frame(rank_f = rep$rank_f, track = TRACK$bp$track,
                       value = log10(pmax(rep$n_bioprojects_total, 1)))
  cont_df <- data.frame(rank_f = rep$rank_f, track = TRACK$cont$track, value = rep$n_continents)
  homd_df <- data.frame(rank_f = rep$rank_f, track = TRACK$homd$track, value = rep$homd_identity)
  inhomd_df <- data.frame(rank_f = rep$rank_f, track = TRACK$tax$track,
                          in_homd = ifelse(rep$annotation_in_homd, "yes", "no"))
  prof_all$sites <- as.character(prof_all$sites)
  sites_by_cluster <- aggregate(sites ~ cluster_id, data = prof_all, FUN = function(x) {
    vals <- unique(unlist(lapply(x, split_values)))
    paste(vals, collapse = ";")
  })
  sites_map <- setNames(sites_by_cluster$sites, as.character(sites_by_cluster$cluster_id))
  rep$sites <- sites_map[as.character(rep$cluster_id)]
  site_order <- c("saliva", "tooth", "tongue", "oral mucosa", "oral swab", "gingiva", "Other")
  all_sites <- intersect(site_order, unique(unlist(lapply(rep$sites, split_values))))
  sites_mat <- do.call(rbind, lapply(seq_len(n), function(i) {
    vals <- split_values(rep$sites[i])
    data.frame(rank_f = rep$rank_f[i], track = all_sites, detected = all_sites %in% vals)
  }))
  sites_mat$track <- factor(sites_mat$track, levels = all_sites)
  sites_w <- nlevels(sites_mat$track) * TILE_W   # one site sub-column per site

  p_lab <- mk_label_col(rep, title = NULL, subtitle = NULL,
                        ylab_sz = ROW_LAB_SZ_MM)

  p_dist <- mk_tile(dist_df, "value",
                    scale_fill_gradient(trans = "log10", low = LOW,
                                        high = TRACK$dist$high,
                                        limits = DIST_RANGE,
                                        # only the two ends: three labels
                                        # 2 decades apart collide at this bar
                                        # width
                                        breaks = c(1e-4, 1),
                                        labels = c("1e-04", "1"),
                                        name = TRACK$dist$name,
                                        guide = cb_guide()))
  # ASV counts are shown on a log10 fill but the legend is labelled with the
  # counts themselves, so the track name can be plain "ASV number"
  p_asv  <- mk_tile(asv_df, "value",
                    scale_fill_gradient(low = LOW, high = TRACK$asv$high,
                                        limits = ASV_RANGE,
                                        breaks = log10(c(3, 10, 100)),
                                        labels = c("3", "10", "100"),
                                        name = TRACK$asv$name,
                                        guide = cb_guide()))
  p_bp   <- mk_tile(bp_df, "value",
                    scale_fill_gradient(low = LOW, high = TRACK$bp$high,
                                        limits = BP_RANGE,
                                        breaks = log10(c(3, 10, 100)),
                                        labels = c("3", "10", "100"),
                                        name = TRACK$bp$name,
                                        guide = cb_guide()))
  p_cont <- mk_tile(cont_df, "value",
                    scale_fill_gradient(low = LOW, high = TRACK$cont$high,
                                        limits = CONT_RANGE, breaks = 2:6,
                                        name = TRACK$cont$name,
                                        guide = lg_guide(nrow = 2, ncol = 3)))
  p_homd <- mk_tile(homd_df, "value",
                    scale_fill_gradient(low = LOW, high = TRACK$homd$high,
                                        limits = HOMD_RANGE,
                                        breaks = c(50, 70, 90),
                                        name = TRACK$homd$name,
                                        guide = cb_guide()))
  p_inhomd <- mk_tile(inhomd_df, "in_homd",
                      scale_fill_manual(values = c("yes" = TRACK$tax$yes,
                                                   "no" = FLAG_OFF),
                                        name = TRACK$tax$name, drop = FALSE,
                                        guide = lg_guide()))
  p_sites <- mk_tile(sites_mat, "detected",
                     scale_fill_manual(values = c(`FALSE` = FLAG_OFF,
                                                  `TRUE` = TRACK$site$yes),
                                       # the keys are logical, so without labels
                                       # the legend read "FALSE" / "TRUE"
                                       labels = c(`FALSE` = "no", `TRUE` = "yes"),
                                       name = TRACK$site$name, drop = FALSE,
                                       guide = lg_guide()))

  prof_all$subject_prevalence <- suppressWarnings(as.numeric(prof_all$subject_prevalence))
  prev_vals <- prof_all[prof_all$cluster_id %in% rep$cluster_id,
                        c("cluster_id", "derep_id", "subject_prevalence")]
  # subject_prevalence = 0 cannot be shown on a log10 axis; drop those rows
  prev_vals <- prev_vals[!is.na(prev_vals$subject_prevalence) &
                           prev_vals$subject_prevalence > 0, ]
  rank_map <- setNames(rep$rank_f, as.character(rep$cluster_id))
  prev_vals$rank_f <- rank_map[as.character(prev_vals$cluster_id)]
  p_box <- ggplot(prev_vals, aes(y = rank_f, x = subject_prevalence)) +
    # the box outline is a hairline as well: at the default 0.5 mm, scaled
    # down to the page, the box read as a dark frame round every row
    geom_boxplot(width = 0.7, outlier.shape = NA, colour = "#3F6F75",
                 fill = "#2F5D6255", linewidth = 0.25) +
    # the points are there to be counted.  `size` is in mm and they are drawn
    # as shape 16: the default glyph carries the theme's 0.5 mm border, which
    # is wider than these dots and turned each one into a blob
    geom_jitter(width = 0, height = 0.20, size = 0.45, alpha = 0.75,
                colour = "#2F5D62", shape = 16) +
    # the panel is less than an inch wide: four decade labels were wider than
    # the axis and ran into each other ("10" on "100", "100" on "1K"), so only
    # the outer decades are labelled and the axis keeps extra room at both ends
    # for the ones at the limits
    scale_x_log10(breaks = c(1, 100, 1000), labels = c("1", "100", "1K"),
                  expand = expansion(mult = c(0.04, 0.06))) +
    labs(y = NULL, x = "Prevalence") +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = AXIS_LAB_SZ_PT),
          axis.title.x = element_text(size = AXIS_LAB_SZ_PT),
          plot.margin = margin(2, 4, 2, 2))

  leaves <- list(p_lab, p_dist, p_asv, p_bp, p_cont, p_homd,
                 p_inhomd, p_sites, p_box)
  widths <- c(LAB_W[["A"]], rep(TILE_W, 6), sites_w, PREV_W)
  combined <- assemble_row(p_lab, leaves[-1], LAB_W[["A"]],
                           c(rep(TILE_W, 6), sites_w, PREV_W))

  # Seven legends on one line crowd each other: at this canvas width the strip
  # had to be squeezed to fit inside it.  patchwork collects them into a single
  # box and has no way to wrap that box, so the legends are lifted out of their
  # panels and set on two lines here instead.  The two lines take every inch the
  # three did plus the space they free: the first line is the three continuous
  # tracks and the discrete continents key, the second the three that are left.
  legs <- lapply(list(p_dist, p_asv, p_bp, p_cont, p_homd, p_inhomd, p_sites),
                 function(q) cowplot::get_legend(q))
  strip <- cowplot::plot_grid(
    cowplot::plot_grid(plotlist = legs[1:4], nrow = 1),
    cowplot::plot_grid(plotlist = legs[5:7], nrow = 1),
    ncol = 1)

  STRIP_H <- 1.30                       # the two legend lines, inches
  combined <- patchwork::wrap_elements(full = strip) /
    (combined & theme(legend.position = "none")) +
    patchwork::plot_layout(heights = c(STRIP_H, fig_h - STRIP_H))

  # the canvas no longer has to be wide enough for one long legend strip
  out_base <- file.path(FIGOUT, "fig_v3_supported_multi_asv_profile")
  save_fig(combined, out_base, width = sum(widths), height = fig_h)
  message(sprintf("[12j_v3] cluster figure written (%d clusters, %.1f in)", n, fig_h))
  list(combined = combined, leaves = leaves, widths = widths,
       n = n, fig_h = fig_h)
}

# ============================================================
# Fig B : singleton top N, one row per ASV
# ============================================================
build_singleton_figure <- function(prof_all) {
  prof_all$n_derep_asvs <- suppressWarnings(as.numeric(prof_all$n_derep_asvs))
  prof_all$subject_prevalence <- suppressWarnings(as.numeric(prof_all$subject_prevalence))
  prof_all$homd_identity <- suppressWarnings(as.numeric(prof_all$homd_identity))
  prof_all$annotation_distance <- suppressWarnings(as.numeric(prof_all$annotation_distance))
  prof_all$annotation <- clean(prof_all$annotation)
  prof_all$annotation_in_homd[is.na(prof_all$annotation_in_homd)] <- FALSE
  prof_all$n_cont <- vapply(prof_all$continents, nvals, integer(1))
  prof_all$label <- row_label(prof_all$annotation, prof_all$cluster_id)

  d <- prof_all[prof_all$n_derep_asvs == 1, ]
  d <- d[order(-d$subject_prevalence, d$derep_id), ][1:min(TOP_N, nrow(d)), ]
  d$rank <- seq_len(nrow(d))
  d$rank_f <- factor(d$rank, levels = rev(d$rank))
  n <- nrow(d)
  # 1.9 in above the ranks: 1.3 for the two-line legend strip and the rest for
  # the label column's title and subtitle
  fig_h <- max(3.4, 1.9 + n * 0.16)

  dist_df <- data.frame(rank_f = d$rank_f, track = TRACK$dist$track,
                        value = pmax(d$annotation_distance, 1e-6))
  bp_df  <- data.frame(rank_f = d$rank_f, track = TRACK$bp$track,
                       value = log10(pmax(d$n_bioprojects, 1)))
  cont_df <- data.frame(rank_f = d$rank_f, track = TRACK$cont$track, value = d$n_cont)
  homd_df <- data.frame(rank_f = d$rank_f, track = TRACK$homd$track, value = d$homd_identity)
  inhomd_df <- data.frame(rank_f = d$rank_f, track = TRACK$tax$track,
                          in_homd = ifelse(d$annotation_in_homd, "yes", "no"))
  site_order <- c("saliva", "tooth", "tongue", "oral mucosa", "oral swab", "gingiva", "Other")
  all_sites <- intersect(site_order, unique(unlist(lapply(d$sites, split_values))))
  sites_mat <- do.call(rbind, lapply(seq_len(n), function(i) {
    vals <- split_values(d$sites[i])
    data.frame(rank_f = d$rank_f[i], track = all_sites, detected = all_sites %in% vals)
  }))
  sites_mat$track <- factor(sites_mat$track, levels = all_sites)
  sites_w <- nlevels(sites_mat$track) * TILE_W   # one site sub-column per site

  p_lab <- mk_label_col(d,
                        title = sprintf("Singleton ASVs (top %d)", n),
                        subtitle = "Sorted by subject prevalence\nright: prevalence (bar)",
                        ylab_sz = ROW_LAB_SZ_MM)

  p_dist <- mk_tile(dist_df, "value",
                    scale_fill_gradient(trans = "log10", low = LOW,
                                        high = TRACK$dist$high,
                                        limits = DIST_RANGE,
                                        # only the two ends: three labels
                                        # 2 decades apart collide at this bar
                                        # width
                                        breaks = c(1e-4, 1),
                                        labels = c("1e-04", "1"),
                                        name = TRACK$dist$name,
                                        guide = cb_guide()))
  p_bp   <- mk_tile(bp_df, "value",
                    scale_fill_gradient(low = LOW, high = TRACK$bp$high,
                                        limits = BP_RANGE,
                                        breaks = log10(c(3, 10, 100)),
                                        labels = c("3", "10", "100"),
                                        name = TRACK$bp$name,
                                        guide = cb_guide()))
  p_cont <- mk_tile(cont_df, "value",
                    scale_fill_gradient(low = LOW, high = TRACK$cont$high,
                                        limits = CONT_RANGE, breaks = 2:6,
                                        name = TRACK$cont$name,
                                        guide = lg_guide(nrow = 2, ncol = 3)))
  p_homd <- mk_tile(homd_df, "value",
                    scale_fill_gradient(low = LOW, high = TRACK$homd$high,
                                        limits = HOMD_RANGE,
                                        breaks = c(50, 70, 90),
                                        name = TRACK$homd$name,
                                        guide = cb_guide()))
  p_inhomd <- mk_tile(inhomd_df, "in_homd",
                      scale_fill_manual(values = c("yes" = TRACK$tax$yes,
                                                   "no" = FLAG_OFF),
                                        name = TRACK$tax$name, drop = FALSE,
                                        guide = lg_guide()))
  p_sites <- mk_tile(sites_mat, "detected",
                     scale_fill_manual(values = c(`FALSE` = FLAG_OFF,
                                                  `TRUE` = TRACK$site$yes),
                                       # the keys are logical, so without labels
                                       # the legend read "FALSE" / "TRUE"
                                       labels = c(`FALSE` = "no", `TRUE` = "yes"),
                                       name = TRACK$site$name, drop = FALSE,
                                       guide = lg_guide()))
  p_bar <- ggplot(data.frame(rank_f = d$rank_f, subject_prevalence = d$subject_prevalence),
                  aes(y = rank_f, x = subject_prevalence)) +
    geom_col(fill = "#2F5D62", width = 0.8) +
    # the panel is under an inch wide: "1000" and "1500" were drawn on top of
    # one another, so the thousands are abbreviated
    scale_x_continuous(breaks = c(0, 500, 1000, 1500),
                       labels = c("0", "500", "1K", "1.5K"),
                       expand = expansion(mult = c(0.02, 0.06))) +
    labs(y = NULL, x = "Prevalence") +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = AXIS_LAB_SZ_PT),
          axis.title.x = element_text(size = AXIS_LAB_SZ_PT),
          plot.margin = margin(2, 4, 2, 2))

  leaves <- list(p_lab, p_dist, p_bp, p_cont, p_homd,
                 p_inhomd, p_sites, p_bar)
  widths <- c(LAB_W[["B"]], rep(TILE_W, 5), sites_w, PREV_W)
  combined <- assemble_row(p_lab, leaves[-1], LAB_W[["B"]],
                           c(rep(TILE_W, 5), sites_w, PREV_W))

  # Six keys do not fit on one line at this canvas width (5.7 in): the collected
  # box was wider than the canvas and was cut off at both ends, so the first key
  # lost its title and the last lost half its swatch.  As in figure A the legends
  # are lifted out of the panels and set on two lines of three.
  legs <- lapply(list(p_dist, p_bp, p_cont, p_homd, p_inhomd, p_sites),
                 function(q) cowplot::get_legend(q))
  strip <- cowplot::plot_grid(
    cowplot::plot_grid(plotlist = legs[1:3], nrow = 1),
    cowplot::plot_grid(plotlist = legs[4:6], nrow = 1),
    ncol = 1)
  STRIP_H <- 1.30                       # the two legend lines, inches
  combined <- patchwork::wrap_elements(full = strip) /
    (combined & theme(legend.position = "none")) +
    patchwork::plot_layout(heights = c(STRIP_H, fig_h - STRIP_H))

  out_base <- file.path(FIGOUT, "fig_v3_supported_singleton_top40_profile")
  save_fig(combined, out_base, width = max(sum(widths), LEGEND_W), height = fig_h)
  message(sprintf("[12j_v3] singleton figure written (%d ASVs, %.1f in)", n, fig_h))
  list(combined = combined, leaves = leaves, widths = widths,
       n = n, fig_h = fig_h)
}

# ---- run ----
figA <- build_cluster_figure(anchors, prof)
figB <- build_singleton_figure(prof)

# ---- final: Fig A | Fig B in one row with a shared legend ----
blank <- ggplot() + theme_void()
flat_leaves <- c(figA$leaves, list(blank), figB$leaves)
# enforce Fig A : Fig B = 1.1 : 1 so long species labels on the left have room
ratio_AB <- 1.1
sumA <- sum(figA$widths)
sumB <- sum(figB$widths)
scaleA <- ratio_AB * sumB / sumA
flat_widths <- c(figA$widths * scaleA, GAP_W, figB$widths)

combined_row <- do.call(
  aplot::plot_list,
  c(flat_leaves,
    list(ncol = length(flat_leaves), widths = flat_widths, guides = "collect"))
)
combined_row <- combined_row &
  theme(legend.position = "top",
        legend.box = "horizontal",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.32, "cm"),
        legend.spacing.x = unit(0.15, "cm"))

out_comb <- file.path(FIGOUT, "fig_v3_supported_asv_profiles_row")
save_fig(combined_row, out_comb, width = 12, height = 8)
message(sprintf("[12j_v3] combined row figure written (%.1f x %.1f in)",
                12, 8))
cat("Done!\n")
