#!/usr/bin/env Rscript
# ============================================================================
# Panel A: cohort provenance Sankey (BioProject-weighted).
# Ribbon width = BioProject count; shade = sample-count rank within group.
# Tuned for Cell Press full width (174 mm).
#
# Three outputs from one builder:
#   panel_A_sankey                 Disease -> Disease category -> Continent -> 16S region
#   panel_A_sankey_no_category     Disease -> Continent -> 16S region
#                                  (ordered and coloured by anatomical site)
#   panel_A_sankey_AnatomicalGroup Disease -> Anatomical site -> Continent -> 16S region
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggforce)
})

source("/beegfs/xuhuan/global_oral_tmp/panel_scripts_edit/pub_theme.R")

DATA_PATH <- file.path(WORK_DIR, "cohort_overview", "cohort_overview_data.csv")
if (!file.exists(DATA_PATH)) stop("Missing: ", DATA_PATH)

CATEGORY_MEMBERS <- list(
  "Oral and dental diseases" = c(
    "Periodontal Disease", "Dental Caries", "Aphthous", "Oral Lichen Planus"
  ),
  "Malignancies" = c(
    "OSCC", "Head and Neck Cancer", "Colorectal Cancer", "Lung Cancer",
    "Breast Cancer", "Esophageal Cancer", "Thyroid Cancer"
  ),
  "Infectious and immune-mediated diseases" = c("COVID-19", "HIV", "RA"),
  "Cardiometabolic and sleep-related disorders" = c(
    "Obesity", "Diabetes", "Hypertension", "Obstructive Sleep Apnea"
  ),
  "Gastrointestinal diseases" = c("IBD", "Gastroesophageal Reflux Disease")
)

CATEGORY_SHORT <- c(
  "Oral and dental diseases" = "Oral / dental",
  "Malignancies" = "Malignancies",
  "Infectious and immune-mediated diseases" = "Infectious / immune",
  "Cardiometabolic and sleep-related disorders" = "Cardiometabolic / sleep",
  "Gastrointestinal diseases" = "GI"
)

CATEGORY_COLOURS <- setNames(
  c("#2C7FB8", "#B2182B", "#7B3294", "#E08214", "#1B7837"),
  names(CATEGORY_MEMBERS)
)

# Anatomical site of disease. List order is the fixed top-to-bottom order.
ANATOMICAL_MEMBERS <- list(
  "Oral cavity" = c(
    "Aphthous", "Dental Caries", "Oral Lichen Planus",
    "Periodontal Disease", "OSCC"
  ),
  "Head/neck & upper airway" = c(
    "Head and Neck Cancer", "Obstructive Sleep Apnea", "Thyroid Cancer"
  ),
  "Gastrointestinal tract" = c(
    "Esophageal Cancer", "Gastroesophageal Reflux Disease",
    "Colorectal Cancer", "IBD"
  ),
  "Respiratory / thoracic" = c("Lung Cancer", "COVID-19", "Breast Cancer"),
  "Systemic" = c("Diabetes", "Hypertension", "Obesity", "HIV", "RA")
)

ANATOMICAL_SHORT <- c(
  "Oral cavity" = "Oral cavity",
  "Head/neck & upper airway" = "Head / neck &\nupper airway",
  "Gastrointestinal tract" = "GI tract",
  "Respiratory / thoracic" = "Respiratory /\nthoracic",
  "Systemic" = "Systemic"
)

# Same five colours as the category version; oral stays blue.
ANATOMICAL_COLOURS <- c(
  "Oral cavity" = "#2C7FB8",
  "Head/neck & upper airway" = "#7B3294",
  "Gastrointestinal tract" = "#1B7837",
  "Respiratory / thoracic" = "#B2182B",
  "Systemic" = "#E08214"
)

DISPLAY_NAME <- c(
  "Aphthous" = "RAS",
  "Colorectal Cancer" = "CRC",
  "Head and Neck Cancer" = "H&N Cancer",
  "Obstructive Sleep Apnea" = "OSA",
  "Gastroesophageal Reflux Disease" = "GERD",
  "Periodontal Disease" = "Periodontal",
  "Esophageal Cancer" = "Esophageal Ca",
  "Breast Cancer" = "Breast Ca",
  "Lung Cancer" = "Lung Ca",
  "Thyroid Cancer" = "Thyroid Ca",
  "Dental Caries" = "Dental Caries",
  "Oral Lichen Planus" = "OLP"
)

display_disease <- function(x) {
  out <- as.character(x)
  hit <- out %in% names(DISPLAY_NAME)
  out[hit] <- unname(DISPLAY_NAME[out[hit]])
  out
}

# Compact geometry for 174 mm width.
SANKEY_SEP <- 0.035
SANKEY_AXIS_W <- 0.10
LABEL_PT_OTHER <- 5.5
LABEL_PT_DISEASE <- 5.5

cohort <- read.csv(DATA_PATH, stringsAsFactors = FALSE, check.names = FALSE)

continent_levels <- cohort %>% count(Continent, sort = TRUE) %>% pull(Continent)
region_levels <- cohort %>% count(Region, sort = TRUE) %>% pull(Region)

# Used only by the compact by-disease figure.
CONTINENT_ABBREV <- c(
  "Asia" = "AS",
  "Europe" = "EU",
  "Africa" = "AF",
  "Oceania" = "OC",
  "North America" = "NA",
  "South America" = "SA"
)

META_PATH <- file.path(WORK_DIR, "metadata_oral_qc_dedup.tsv")
SITE_DISPLAY <- c(
  "saliva" = "Saliva",
  "Other" = "Other",
  "tooth" = "Tooth",
  "tongue" = "Tongue",
  "oral mucosa" = "Oral mucosa",
  "gingiva" = "Gingiva",
  "oral swab" = "Oral swab"
)

# Sampling-site share of each (disease, BioProject) pair: cases and healthy
# controls of that pair only. Fractions sum to 1 per pair, so a BioProject
# that mixes sites still contributes one BioProject of ribbon width.
load_site_fractions <- function() {
  meta <- read.delim(META_PATH, stringsAsFactors = FALSE, check.names = FALSE,
                     quote = "")
  meta <- meta[, c("sample_id", "bioproject_id", "status_cleaned", "site_cleaned")]
  for (col in names(meta)) meta[[col]] <- trimws(as.character(meta[[col]]))
  meta$status_cleaned[meta$status_cleaned %in% c("HC", "Healthy")] <- "Healthy Control"
  meta <- meta[!duplicated(meta$sample_id), , drop = FALSE]

  pairs <- distinct(cohort, Disease, bioproject_id)
  samples <- bind_rows(lapply(split(pairs, pairs$Disease), function(p) {
    d <- p$Disease[[1]]
    meta %>%
      filter(bioproject_id %in% p$bioproject_id,
             status_cleaned %in% c(d, "Healthy Control")) %>%
      mutate(Disease = d)
  }))
  unmapped <- setdiff(unique(samples$site_cleaned), names(SITE_DISPLAY))
  if (length(unmapped)) {
    stop("Sampling sites without a display name: ", paste(unmapped, collapse = ", "))
  }
  samples$Site <- unname(SITE_DISPLAY[samples$site_cleaned])

  frac <- samples %>%
    count(Disease, bioproject_id, Site, name = "n") %>%
    group_by(Disease, bioproject_id) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()
  check <- frac %>%
    group_by(Disease, bioproject_id) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    inner_join(cohort, by = c("Disease", "bioproject_id"))
  if (nrow(check) != nrow(cohort) || any(check$n != check$n_samples)) {
    stop("Site sample counts do not match cohort_overview_data.csv")
  }
  site_levels <- samples %>% count(Site, sort = TRUE) %>% pull(Site)
  list(frac = frac, levels = site_levels)
}

shade_ramp <- function(base, n) {
  if (n <= 1) return(base)
  dark <- grDevices::colorRampPalette(c(base, "#1A1A1A"))(10)[3]
  light <- grDevices::colorRampPalette(c(base, "white"))(10)[6]
  grDevices::colorRampPalette(c(dark, light))(n)
}

# Twenty medium-dark colours, ordered so a stride through the list jumps hue.
# Light yellows and greys are avoided: at alpha 0.72 they disappear on white.
DISEASE_PALETTE <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4",
  "#1199B8", "#F032E6", "#7CB518", "#0E8A7D", "#9A6324",
  "#800000", "#000075", "#E56AAA", "#6C4AB6", "#C4A000",
  "#008B8B", "#E09A4A", "#2E7D4F", "#5C4B8A", "#B85C38"
)

# Same hues and same order, with chroma clamped at 52 (HCL conversion via the
# colorspace package) so the ribbons read as a muted scientific palette rather
# than a high-saturation rainbow. Hue assignments -- and therefore the
# disease-to-colour mapping -- are unchanged.
DISEASE_PALETTE_MUTED <- c(
  "#AC5F68", "#68AD6E", "#5E68A5", "#D19279", "#78468D",
  "#2698B6", "#B976B4", "#8BB064", "#0E8A7D", "#966433",
  "#6D2323", "#000075", "#CD7BA3", "#675498", "#BBA15D",
  "#008B8B", "#D19F72", "#2E7D4F", "#5C4B8A", "#A66552"
)

# Every cohort disease must sit in exactly one group.
group_lookup <- function(group_members) {
  long <- stack(group_members)
  dup <- unique(as.character(long$values[duplicated(long$values)]))
  if (length(dup)) stop("Disease in more than one group: ", paste(dup, collapse = ", "))
  map <- setNames(as.character(long$ind), as.character(long$values))
  cohort_dz <- unique(cohort$Disease)
  missing <- setdiff(cohort_dz, names(map))
  extra <- setdiff(names(map), cohort_dz)
  if (length(missing)) stop("Diseases without a group: ", paste(missing, collapse = ", "))
  if (length(extra)) stop("Grouped diseases absent from cohort: ", paste(extra, collapse = ", "))
  map
}

# ============================================================================
# Stretched-column sankey variant.
#
# Standard ggforce gives every axis a total height proportional to the value
# flowing through it, so the Disease column towers over the others. In this
# variant every column is stretched to the total height of the first axis
# (Disease): inside a column the boxes keep their value proportions and share
# the height left after the gaps, the gaps are the same absolute size on all
# columns (the Disease column's gap), and a ribbon's end thickness is its
# value times the endpoint column's stretch factor. Ribbons therefore still
# tile every box exactly, while the two ends of a ribbon may differ when the
# two columns are stretched by different factors.
# ============================================================================

eq_df0 <- getFromNamespace("data_frame0", "ggplot2")
eq_unique0 <- getFromNamespace("unique0", "ggplot2")
eq_flip_data <- getFromNamespace("flip_data", "ggplot2")
eq_remove_group <- getFromNamespace("remove_group", "ggforce")
eq_complete_data <- getFromNamespace("complete_data", "ggforce")
eq_stat_diag <- getFromNamespace("StatDiagonalWide", "ggforce")
eq_geom_shape <- getFromNamespace("GeomShape", "ggforce")
eq_geom_axes <- getFromNamespace("GeomParallelSetsAxes", "ggforce")

# Per-axis box geometry stretched to a common column height. Returns the box
# frames (stacked from the bottom up in ggforce's own level order) plus, per
# axis, the multiplier that maps a flow value to drawn y units at that axis.
# `target` is the common column total (set once from the Disease column in
# build_sankey, so partial layers such as the label subplots agree).
eq_axes_and_scale <- function(data, sep, target) {
  axes_raw <- lapply(split(data, data$x), function(d) {
    splits <- split(d$value, as.character(d$split))
    splits <- splits[rev(order(match(names(splits), levels(d$split))))]
    splits <- splits[names(splits) != ".ggforce_missing"]
    d <- eq_df0(
      split = factor(names(splits)),
      value = sapply(splits, sum),
      x = d$x[1]
    )
    gap <- sum(d$value) * sep
    d$ymax <- (seq_len(nrow(d)) - 1) * gap + cumsum(d$value)
    d$ymin <- d$ymax - d$value
    d
  })
  # The value total is identical on every axis, so the gap is too.
  gap1 <- sum(axes_raw[[1]]$value) * sep
  scale <- vapply(axes_raw, function(ax) {
    (target - (nrow(ax) - 1) * gap1) / sum(ax$value)
  }, numeric(1))
  axes <- mapply(function(ax, sc) {
    ax$value <- ax$value * sc
    ax$ymax <- (seq_len(nrow(ax)) - 1) * gap1 + cumsum(ax$value)
    ax$ymin <- ax$ymax - ax$value
    ax
  }, axes_raw, scale, SIMPLIFY = FALSE)
  list(axes = do.call(vctrs::vec_rbind, axes), scale = scale, gap = gap1)
}

# Per-pair ribbon geometry: thickness at each end is the flow's value scaled
# by the endpoint column's factor, and flows stack inside each box from its
# top down, in ggforce's own flow-group order.
sankey_diag_data_eq <- function(data, axes_data, scale, groups, axis.width) {
  axes <- sort(eq_unique0(data$x))
  pair_list <- unname(lapply(seq_len(length(axes) - 1), function(i) {
      from <- data[data$x == axes[i], , drop = FALSE]
      to <- data[data$x == axes[i + 1], , drop = FALSE]
      flows <- split(seq_len(nrow(from)), list(from$group, from$split, to$split))
      flows <- flows[lengths(flows) != 0]
      rep_row <- vapply(flows, `[`, integer(1), 1)
      n <- length(flows)
      g <- from$group[rep_row]
      s_from <- as.character(from$split[rep_row])
      s_to <- as.character(to$split[rep_row])
      vals <- vapply(flows, function(ii) sum(from$value[ii]), numeric(1))
      # Flows into or out of a phantom (.ggforce_missing) slot have no drawn
      # box to tile; drop them.
      keep <- s_from != ".ggforce_missing" & s_to != ".ggforce_missing"
      if (!all(keep)) {
        g <- g[keep]; s_from <- s_from[keep]; s_to <- s_to[keep]
        vals <- vals[keep]
        n <- length(g)
      }
      th_from <- scale[[as.character(axes[i])]] * vals
      th_to <- scale[[as.character(axes[i + 1])]] * vals
      stack <- function(s, th, ax) {
        y <- numeric(length(s))
        for (sp in unique(s)) {
          ii <- which(s == sp)
          ii <- ii[order(g[ii])]
          y[ii] <- as.numeric(ax$ymax[as.character(ax$split) == sp]) -
            cumsum(th[ii])
        }
        y
      }
      y_from <- stack(s_from, th_from, axes_data[axes_data$x == axes[i], , drop = FALSE])
      y_to <- stack(s_to, th_to, axes_data[axes_data$x == axes[i + 1], , drop = FALSE])
      xf <- from$x[1] + axis.width / 2
      xt <- to$x[1] - axis.width / 2
      out <- eq_df0(
        x = c(rep(xf, 2 * n), rep(xt, 2 * n)),
        y = c(y_from, y_from + th_from, y_to, y_to + th_to),
        group = rep(seq_len(n), 4)
      )
      aes_cols <- groups[match(g, groups$group),
                         setdiff(names(groups), "group"), drop = FALSE]
      if (ncol(aes_cols) > 0) {
        out <- cbind(out, aes_cols[rep(seq_len(n), 4), , drop = FALSE])
      }
      out
    }
  ))
  # Ribbon group ids must be unique across axis pairs.
  n_pairs <- vapply(pair_list, function(d) nrow(d) / 4, numeric(1))
  group_offset <- c(0, cumsum(n_pairs)[-length(n_pairs)])
  pair_list <- Map(function(d, i) { d$group <- d$group + i; d },
                   pair_list, group_offset)
  do.call(vctrs::vec_rbind, pair_list)
}

StatParallelSetsEq <- ggproto(
  "StatParallelSetsEq", getFromNamespace("StatParallelSets", "ggforce"),
  compute_panel = function(data, scales, sep = 0.05, strength = 0.5, n = 100,
                           axis.width = 0, eq.target = numeric(0),
                           flipped_aes = FALSE) {
    data <- eq_flip_data(data, flipped_aes)
    data <- eq_remove_group(data)
    data <- eq_complete_data(data)
    cols <- c("group", "colour", "color", "fill", "size", "linewidth",
              "alpha", "linetype")
    data_groups <- do.call(vctrs::vec_rbind, unname(lapply(
      split(data[, names(data) %in% cols, drop = FALSE], data$group),
      function(d) do.call(eq_df0, lapply(d, function(x) na.omit(x)[1]))
    )))
    geom <- eq_axes_and_scale(data, sep, eq.target)
    diagonals <- sankey_diag_data_eq(data, geom$axes, geom$scale, data_groups, axis.width)
    diagonals <- eq_flip_data(diagonals, flipped_aes)
    eq_stat_diag$compute_panel(diagonals, scales, strength, n, flipped_aes)
  }
)

StatParallelSetsAxesEq <- ggproto(
  "StatParallelSetsAxesEq", getFromNamespace("StatParallelSetsAxes", "ggforce"),
  compute_panel = function(data, scales, sep = 0.05, axis.width = 0,
                           eq.target = numeric(0), flipped_aes = FALSE) {
    data <- eq_flip_data(data, flipped_aes)
    split_levels <- levels(data$split)
    data <- eq_remove_group(data)
    data <- eq_complete_data(data, FALSE)
    data_axes <- eq_axes_and_scale(data, sep, eq.target)$axes
    data_axes <- data_axes[data_axes$split != ".ggforce_missing", ]
    cols <- c("x", "split", "colour", "color", "fill", "size",
              "linewidth", "alpha", "linetype")
    aes <- data[, names(data) %in% cols]
    aes <- eq_unique0(aes)
    if (nrow(aes) != nrow(data_axes)) {
      cli::cli_abort("Axis aesthetics must be constant in each split")
    }
    data_axes$split <- factor(as.character(data_axes$split), levels = split_levels)
    aes$split <- factor(as.character(aes$split), levels = split_levels)
    data <- merge(data_axes, aes, by = c("x", "split"), all.x = TRUE, sort = FALSE)
    names(data)[names(data) == "split"] <- "label"
    data$y <- data$ymin + data$value / 2
    data$xmin <- data$x - axis.width / 2
    data$xmax <- data$x + axis.width / 2
    eq_flip_data(data, flipped_aes)
  }
)

geom_parallel_sets_eq <- function(mapping = NULL, data = NULL, stat = NULL,
    position = "identity", n = 100, na.rm = FALSE, orientation = NA,
    sep = 0.05, strength = 0.5, axis.width = 0, show.legend = NA,
    inherit.aes = TRUE, ...) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = StatParallelSetsEq,
    geom = eq_geom_shape, position = position, show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, orientation = orientation, n = n,
                  strength = strength, sep = sep, axis.width = axis.width, ...)
  )
}

geom_parallel_sets_axes_eq <- function(mapping = NULL, data = NULL, stat = NULL,
    position = "identity", na.rm = FALSE, orientation = NA, show.legend = NA,
    inherit.aes = TRUE, ...) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = StatParallelSetsAxesEq,
    geom = eq_geom_axes, position = position, show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, orientation = orientation, ...)
  )
}

# ---------------------------------------------------------------------------
# Crossing-reducing node ordering (weighted barycentric / Sugiyama heuristic).
#
# Aggregated rows connect one node on every axis at once, so the flows between
# two adjacent axes are just the rows read on those two columns. A crossing
# between two flows is therefore a pair of rows whose relative order on one
# axis is inverted on the next.
#
# sankey_crossings counts them per axis pair (unweighted). sankey_bary_order
# runs barycentric sweeps: each non-Disease node is re-sorted by the
# value-weighted mean position of its flow endpoints on the neighbouring axes,
# iterated until the ordering is stable. Ties keep the existing order and the
# Disease column (axis 1) always keeps its order, so the disease labels and
# their BioProject counts read top-down exactly as before.
# ---------------------------------------------------------------------------

sankey_crossings <- function(agg, axis_cols) {
  pos <- lapply(axis_cols, function(col) {
    stats::setNames(seq_along(levels(agg[[col]])), levels(agg[[col]]))
  })
  vapply(seq_len(length(axis_cols) - 1), function(i) {
    a <- unname(pos[[i]][as.character(agg[[axis_cols[[i]]]])])
    b <- unname(pos[[i + 1]][as.character(agg[[axis_cols[[i + 1]]]])])
    bb <- b[order(a)]
    cnt <- 0L
    for (k in seq_along(bb)) {
      cnt <- cnt + sum(bb[seq_len(k - 1L)] > bb[k])
    }
    cnt
  }, integer(1))
}

sankey_bary_order <- function(agg, axis_cols, n_rounds = 4) {
  n_axis <- length(axis_cols)
  cols <- lapply(axis_cols, function(col) as.character(agg[[col]]))
  lvls <- lapply(axis_cols, function(col) levels(agg[[col]]))
  w <- agg$value

  cross_total <- function() {
    pos <- lapply(lvls, function(l) stats::setNames(seq_along(l), l))
    tot <- 0L
    for (i in seq_len(n_axis - 1)) {
      bb <- unname(pos[[i + 1]][cols[[i + 1]]])[order(unname(pos[[i]][cols[[i]]]))]
      for (k in seq_along(bb)) tot <- tot + sum(bb[seq_len(k - 1L)] > bb[k])
    }
    tot
  }

  # One round: barycentric sweeps (each non-Disease node follows the weighted
  # mean position of its flow endpoints on the neighbouring axes), then a
  # hill-climb that accepts an adjacent swap on any non-Disease column only
  # when it strictly lowers the total crossing count.
  for (round_i in seq_len(n_rounds)) {
    previous <- cross_total()
    repeat {
      moved <- FALSE
      for (i in seq_len(n_axis)) {
        if (i == 1L) next # the Disease column keeps its order
        pos <- lapply(lvls, function(l) stats::setNames(seq_along(l), l))
        bary <- vapply(lvls[[i]], function(u) {
          rows <- which(cols[[i]] == u)
          if (!length(rows)) return(NA_real_)
          sides <- if (i == n_axis) i - 1L else c(i - 1L, i + 1L)
          num <- 0; den <- 0
          for (j in sides) {
            num <- num + sum(w[rows] * pos[[j]][cols[[j]][rows]])
            den <- den + sum(w[rows])
          }
          if (den == 0) NA_real_ else num / den
        }, numeric(1))
        # Ties and unused levels keep their current rank.
        ord <- order(bary, pos[[i]][lvls[[i]]], na.last = TRUE)
        if (!identical(ord, seq_along(ord))) {
          lvls[[i]] <- lvls[[i]][ord]
          moved <- TRUE
        }
      }
      if (!moved) break
    }
    repeat {
      best <- cross_total()
      improved <- FALSE
      for (i in 2:n_axis) {
        for (k in seq_len(length(lvls[[i]]) - 1)) {
          lvls[[i]][k:(k + 1)] <- lvls[[i]][(k + 1):k]
          if (cross_total() < best) {
            best <- cross_total()
            improved <- TRUE
          } else {
            lvls[[i]][k:(k + 1)] <- lvls[[i]][(k + 1):k] # swap back
          }
        }
      }
      if (!improved) break
    }
    if (cross_total() >= previous) break
  }
  for (i in seq_len(n_axis)) {
    agg[[axis_cols[[i]]]] <- factor(
      as.character(agg[[axis_cols[[i]]]]), levels = lvls[[i]]
    )
  }
  agg
}

build_sankey <- function(group_members, group_short, group_axis_lab,
                         show_group_axis, group_colours, stem,
                         fixed_order = FALSE, show_group_legend = FALSE,
                         abbreviate_cancers = FALSE,
                         colour_by = c("group_shade", "disease"),
                         compact = FALSE,
                         extra_labels = NULL,
                         tighter = FALSE,
                         show_site = FALSE,
                         rainbow = FALSE,
                         compact_width = NULL,
                         font_scale = 1,
                         bar_scale = 1,
                         compact_height = NULL,
                         equal_column_height = FALSE,
                         flow_alpha = 0.72,
                         muted_palette = FALSE,
                         optimize_order = FALSE,
                         node_fill = "grey88",
                         node_border = "grey55",
                         node_linewidth = 0.2,
                         label_colour = "grey10",
                         continent_order = NULL) {
  colour_by <- match.arg(colour_by)
  label_pt_disease <- if (compact) 7 else LABEL_PT_DISEASE
  # 6 pt sits in the gap above each box; 6.5 pt is taller than that gap.
  label_pt_other <- if (compact) 6 else LABEL_PT_OTHER
  # font_scale upsizes all text for figures that are later enlarged again
  # (e.g. placed on a presentation canvas); the gaps between boxes scale
  # with it so the bigger labels still fit.
  label_pt_disease <- label_pt_disease * font_scale
  label_pt_other <- label_pt_other * font_scale
  # On the doubled-width canvas the axes are drawn at half the original data
  # width so their printed width stays the same as the 174 mm layout.
  # Axis boxes keep the physical width they were tuned at (120 mm for the
  # compact figures) when the compact canvas is overridden wider.
  axis_base <- if (compact && !is.null(compact_width)) 120 / compact_width else 1
  axis_w <- (if (compact) 0.08 else SANKEY_AXIS_W / 2) * axis_base
  # bar_scale widens the axis boxes (and the ribbons meeting them) without
  # changing anything else; the label nudges scale with it so the gap between
  # a label and the bar edge stays the same.
  axis_w <- axis_w * bar_scale
  # Wider gaps on the compact figure so a label fits above each box.
  flow_sep <- if (compact) 0.036 * font_scale else SANKEY_SEP
  # The equal-height variant only redefines box geometry in the compact layout.
  if (equal_column_height && !compact) {
    stop("equal_column_height is only implemented for compact = TRUE")
  }
  label_fn <- if (abbreviate_cancers) disease_abbrev else display_disease
  if (!is.null(extra_labels)) {
    base_fn <- label_fn
    label_fn <- function(x) {
      out <- base_fn(x)
      hit <- out %in% names(extra_labels)
      out[hit] <- unname(extra_labels[out[hit]])
      out
    }
  }
  group_levels <- names(group_members)
  group_map <- group_lookup(group_members)

  dat <- cohort
  cont_levels_use <- as.character(continent_levels)
  region_levels_use <- as.character(region_levels)
  if (compact || show_site) {
    # Top to bottom. V3 is not in the requested list; it sits before V3-V4.
    region_order <- c(
      "V1-V2", "V2-V3", "V1-V3", "V3", "V3-V4",
      "V4", "V4-V5", "V6-V8", "Unknown"
    )
    missing <- setdiff(region_levels_use, region_order)
    if (length(missing)) {
      stop("Region missing from display order: ", paste(missing, collapse = ", "))
    }
    region_levels_use <- intersect(region_order, region_levels_use)
    cont <- as.character(dat$Continent)
    mapped <- CONTINENT_ABBREV[cont]
    dat$Continent <- ifelse(is.na(mapped), cont, unname(mapped))
    mapped_lv <- CONTINENT_ABBREV[cont_levels_use]
    cont_levels_use <- ifelse(is.na(mapped_lv), cont_levels_use, unname(mapped_lv))
  }
  dat$Group <- factor(unname(group_map[dat$Disease]), levels = group_levels)

  disease_stats <- dat %>%
    group_by(Group, Disease) %>%
    summarise(n_bp = n(), n_samples = sum(n_samples), .groups = "drop")

  sankey_group_levels <- if (fixed_order) {
    group_levels
  } else {
    dat %>%
      count(Group, name = "n_bp_pairs") %>%
      mutate(
        Group = as.character(Group),
        n_diseases = vapply(
          Group, function(g) length(group_members[[g]]), integer(1)
        )
      ) %>%
      arrange(desc(n_diseases), desc(n_bp_pairs), Group) %>%
      pull(Group)
  }

  disease_levels <- unlist(lapply(sankey_group_levels, function(g) {
    disease_stats %>%
      filter(as.character(Group) == g) %>%
      mutate(disp = display_disease(Disease)) %>%
      arrange(desc(n_bp), desc(n_samples), disp, Disease) %>%
      pull(Disease)
  }), use.names = FALSE)

  if (colour_by == "group_shade") {
    disease_colour_map <- unlist(lapply(group_levels, function(g) {
      df <- disease_stats %>%
        filter(as.character(Group) == g) %>%
        mutate(disp = display_disease(Disease)) %>%
        arrange(desc(n_samples), disp, Disease)
      stats::setNames(
        shade_ramp(group_colours[[g]], nrow(df)),
        as.character(df$Disease)
      )
    }), use.names = TRUE)
  } else {
    n_dz_col <- length(disease_levels)
    if (rainbow) {
      # Bright hue-rainbow, one equal-luminance hue per disease in axis
      # order; matches the classic multi-colour reference look.
      hues <- (seq_len(n_dz_col) - 1L) * 360 / n_dz_col + 15
      disease_colour_map <- stats::setNames(
        grDevices::hcl(h = hues, c = 95, l = 62),
        disease_levels
      )
      if (anyDuplicated(unname(disease_colour_map))) {
        stop("Disease colours are not unique")
      }
    } else {
      # Stride coprime to the palette length, so diseases that sit next to
      # each other on the axis do not receive neighbouring hues.
      pal <- if (muted_palette) DISEASE_PALETTE_MUTED else DISEASE_PALETTE
      if (n_dz_col > length(pal)) {
        stop("Need ", n_dz_col, " disease colours, palette has ", length(pal))
      }
      stride <- 7L
      idx <- ((seq_len(n_dz_col) - 1L) * stride) %% length(pal) + 1L
      disease_colour_map <- stats::setNames(pal[idx], disease_levels)
      if (anyDuplicated(unname(disease_colour_map))) {
        stop("Disease colours are not unique")
      }
    }
  }

  # Compact labels use a narrower space before "(" so the closing
  # parenthesis does not sit against the disease bar.
  label_fmt <- if (compact) {
    "%s<span style=\"font-size:2pt\"> </span>(%s)"
  } else {
    "%s (%s)"
  }
  disease_label_map <- stats::setNames(
    sprintf(
      label_fmt,
      label_fn(as.character(disease_stats$Disease)),
      format(disease_stats$n_bp, big.mark = ",", trim = TRUE)
    ),
    as.character(disease_stats$Disease)
  )

  # Fixed Continent display order, top to bottom. Only the level order
  # changes; boxes, labels and ribbons all follow the new stacking.
  if (!is.null(continent_order)) {
    missing <- setdiff(cont_levels_use, continent_order)
    if (length(missing)) {
      stop("Continent missing from display order: ", paste(missing, collapse = ", "))
    }
    cont_levels_use <- intersect(continent_order, cont_levels_use)
  }
  if (show_site) {
    site <- load_site_fractions()
    agg <- dat %>%
      inner_join(site$frac, by = c("Disease", "bioproject_id")) %>%
      group_by(Disease, Group, Continent, Site, Region) %>%
      summarise(value = sum(frac), .groups = "drop") %>%
      mutate(Site = factor(Site, levels = site$levels))
  } else {
    agg <- dat %>%
      group_by(Disease, Group, Continent, Region) %>%
      summarise(value = n(), .groups = "drop")
  }
  agg <- agg %>%
    mutate(
      DiseaseLabel = factor(
        unname(disease_label_map[as.character(Disease)]),
        levels = unname(disease_label_map[as.character(disease_levels)])
      ),
      GroupLabel = factor(
        unname(group_short[as.character(Group)]),
        levels = unname(group_short[sankey_group_levels])
      ),
      Continent = factor(Continent, levels = cont_levels_use),
      Region = factor(Region, levels = region_levels_use),
      Disease = factor(as.character(Disease), levels = disease_levels)
    )

  if (show_group_axis) {
    axis_cols <- c("DiseaseLabel", "GroupLabel", "Continent", "Region")
    axis_labs <- c("Disease", group_axis_lab, "Continent", "16S region")
  } else if (show_site) {
    axis_cols <- c("DiseaseLabel", "Continent", "Site", "Region")
    axis_labs <- c("Disease", "Continent", "Site", "16S region")
  } else {
    axis_cols <- c("DiseaseLabel", "Continent", "Region")
    axis_labs <- c("Disease", "Continent", "16S region")
  }
  n_axis <- length(axis_cols)
  all_levels <- unlist(lapply(axis_cols, function(col) levels(agg[[col]])),
                       use.names = FALSE)
  if (anyDuplicated(all_levels)) {
    stop("Category shared by two axes: ",
         paste(unique(all_levels[duplicated(all_levels)]), collapse = ", "))
  }
  agg <- agg %>%
    select(all_of(axis_cols), Disease, value) %>%
    as.data.frame()

  # Optional: re-sort the non-Disease columns by weighted barycentric sweeps
  # to reduce flow crossings between adjacent axes. Node weights, connections
  # and heights are untouched; only the level order of axes 2..n changes.
  if (optimize_order) {
    cross_before <- sankey_crossings(agg, axis_cols)
    agg <- sankey_bary_order(agg, axis_cols)
    cross_after <- sankey_crossings(agg, axis_cols)
    cat(sprintf(
      "%s: flow crossings by axis pair -- before %s | after %s (total %d -> %d)\n",
      stem,
      paste(cross_before, collapse = "/"),
      paste(cross_after, collapse = "/"),
      sum(cross_before), sum(cross_after)
    ))
  }

  gathered <- gather_set_data(agg, seq_len(n_axis))
  if (is.numeric(gathered$x) ||
      all(as.character(gathered$x) %in% as.character(seq_len(n_axis)))) {
    gathered$x <- factor(as.integer(as.character(gathered$x)),
                         levels = seq_len(n_axis), labels = axis_labs)
  } else {
    gathered$x <- factor(as.character(gathered$x),
                         levels = axis_cols, labels = axis_labs)
  }
  gathered$y <- factor(
    as.character(gathered$y),
    levels = unlist(lapply(axis_cols, function(col) levels(agg[[col]])),
                    use.names = FALSE)
  )
  gathered$Disease <- factor(as.character(gathered$Disease), levels = disease_levels)

  gathered_disease <- gathered[gathered$x == "Disease", , drop = FALSE]
  gathered_other <- gathered[gathered$x != "Disease", , drop = FALSE]

  # Stretched-column variant: common column total, taken from the Disease
  # column (axis 1) so every layer -- including the label layers, which only
  # hold a subset of the axes -- stretches by the same recipe. The value total
  # is identical on every axis, so only the box count differs per column.
  eq_target <- NULL
  if (equal_column_height) {
    v_tot <- sum(gathered$value[gathered$x == "Disease"])
    n_dis_boxes <- length(unique(gathered$y[gathered$x == "Disease"]))
    eq_target <- v_tot * (1 + (n_dis_boxes - 1) * flow_sep)
  }

  # Padding tuned for 4 axes (0.2 right).
  # Half the bar is 0.04; -0.035 put the closing parenthesis on the bar.
  label_nudge <- (if (compact) -0.048 else -0.12) * bar_scale
  # Measure the widest disease label once so the left gutter fits it exactly.
  grDevices::pdf(NULL)
  lab_mm <- max(vapply(unname(disease_label_map), function(s) {
    grob <- gridtext::richtext_grob(
      s,
      gp = grid::gpar(
        fontsize = label_pt_disease, fontface = "bold", fontfamily = FONT
      ),
      padding = grid::unit(0, "pt")
    )
    as.numeric(grid::convertWidth(grid::grobWidth(grob), "mm"))
  }, numeric(1))) + 0.3
  grDevices::dev.off()
  if (compact) {
    # Left pad is only as wide as the longest disease label. The 18 mm
    # margin used by the full-width figures is not added on top.
    fig_w <- if (!is.null(compact_width)) {
      compact_width
    } else if (show_site) 120 else if (tighter) 94 else 120
    margin_l <- 1.0
    margin_t <- if (tighter) 0.3 else 2.0
    margin_b <- if (tighter) 1.0 else 1.5
    # Half of the "16S region" axis title hangs past the last column.
    margin_r <- (if (tighter) 12 else 8) * font_scale
    panel_w <- fig_w - margin_l - margin_r
    pad_r <- 0.12
    pad_l <- lab_mm * ((n_axis - 1) + pad_r) / (panel_w - lab_mm) - label_nudge
    x_pad <- c(pad_l, pad_r)
  } else {
    # 2026-10-10: with the doubled canvas the fixed 1.55-unit gutter grew to
    # ~100 mm of empty space, so it is now sized from the measured label width.
    margin_l <- 2.0
    margin_r <- 4.0
    panel_w <- 2 * W_FULL_MM - margin_l - margin_r
    pad_r <- 0.12
    pad_l <- lab_mm * ((n_axis - 1) + pad_r) / (panel_w - lab_mm) - label_nudge
    x_pad <- c(pad_l, pad_r)
  }

  p <- ggplot(gathered, aes(x = x, id = id, split = y, value = value)) +
    {
      if (equal_column_height) {
        geom_parallel_sets_eq(
          aes(fill = Disease),
          alpha = flow_alpha,
          axis.width = axis_w,
          sep = flow_sep,
          eq.target = eq_target,
          n = 50,
          strength = 0.55
        )
      } else {
        geom_parallel_sets(
          aes(fill = Disease),
          alpha = flow_alpha,
          axis.width = axis_w,
          sep = flow_sep,
          n = 50,
          strength = 0.55
        )
      }
    } +
    {
      if (equal_column_height) {
        geom_parallel_sets_axes_eq(
          axis.width = axis_w,
          sep = flow_sep,
          eq.target = eq_target,
          fill = node_fill,
          colour = node_border,
          linewidth = node_linewidth
        )
      } else {
        geom_parallel_sets_axes(
          axis.width = axis_w,
          sep = flow_sep,
          fill = node_fill,
          colour = node_border,
          linewidth = node_linewidth
        )
      }
    } +
    {
      if (compact && equal_column_height) {
        # Names sit in the white gap above each box. Stretched columns use the
        # variant stat so their label positions match the stretched boxes.
        geom_text(
          data = gathered_other,
          stat = StatParallelSetsAxesEq,
          aes(label = after_stat(label), y = after_stat(ymax)),
          axis.width = axis_w,
          sep = flow_sep,
          eq.target = eq_target,
          hjust = 0.5,
          vjust = 0,
          nudge_y = 0.15 * font_scale,
          colour = label_colour,
          size = label_pt_other / ggplot2::.pt,
          lineheight = 0.9,
          fontface = "bold",
          family = FONT
        )
      } else if (compact) {
        # Names sit in the white gap above each box. Inside the box they
        # land on the ribbons and are hard to read.
        geom_text(
          data = gathered_other,
          stat = "parallel_sets_axes",
          aes(label = after_stat(label), y = after_stat(ymax)),
          axis.width = axis_w,
          sep = flow_sep,
          hjust = 0.5,
          vjust = 0,
          nudge_y = 0.15 * font_scale,
          colour = label_colour,
          size = label_pt_other / ggplot2::.pt,
          lineheight = 0.9,
          fontface = "bold",
          family = FONT
        )
      } else {
        geom_parallel_sets_labels(
          data = gathered_other,
          axis.width = axis_w,
          sep = flow_sep,
          colour = label_colour,
          size = label_pt_other / ggplot2::.pt,
          angle = 0,
          hjust = 0.5,
          lineheight = 0.9,
          fontface = "bold",
          family = FONT
        )
      }
    } +
    {
      if (compact) {
        ggtext::geom_richtext(
          data = gathered_disease,
          stat = if (equal_column_height) StatParallelSetsAxesEq else "parallel_sets_axes",
          aes(label = after_stat(label)),
          axis.width = axis_w,
          sep = flow_sep,
          eq.target = if (equal_column_height) eq_target,
          hjust = 1,
          nudge_x = label_nudge,
          colour = label_colour,
          size = label_pt_disease / ggplot2::.pt,
          fontface = "bold",
          family = FONT,
          fill = NA,
          label.colour = NA,
          label.size = 0,
          label.padding = grid::unit(0, "pt"),
          label.r = grid::unit(0, "pt")
        )
      } else {
        geom_parallel_sets_labels(
          data = gathered_disease,
          axis.width = axis_w,
          sep = flow_sep,
          colour = label_colour,
          size = label_pt_disease / ggplot2::.pt,
          angle = 0,
          hjust = 1,
          nudge_x = label_nudge,
          fontface = "bold",
          family = FONT
        )
      }
    } +
    scale_x_discrete(expand = expansion(add = x_pad)) +
    {
      if (tighter || !compact) {
        # Drop the default headroom above the top bar and below the bottom
        # bar so the axis titles sit close to the plot. The stretched-column
        # variant needs a little more top room: its tallest boxes now reach
        # the column top and the label of the topmost box must fit above.
        expand_bottom <- if (equal_column_height) 0.018 else 0.018
        expand_top <- if (equal_column_height) 0.032 else 0.002
        scale_y_continuous(expand = expansion(mult = c(expand_bottom, expand_top)))
      }
    } +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = disease_colour_map, guide = "none", drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_cell() +
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = if (tighter || !compact) {
        element_text(
          face = "bold", size = 7 * font_scale, colour = "black",
          margin = margin(t = 0.4 * font_scale, r = 0, b = 0, l = 0)
        )
      } else {
        element_text(
          face = "bold",
          size = (if (compact) 7 else BASE_PT) * font_scale,
          colour = "black"
        )
      },
      plot.margin = if (compact) {
        margin(margin_t, margin_r, margin_b, margin_l)
      } else {
        margin(2, margin_r, 2, margin_l)
      }
    )

  # Ribbons are filled per disease; this legend shows the group base colour.
  if (show_group_legend) {
    legend_df <- data.frame(
      x = factor("Disease", levels = axis_labs),
      y = 0,
      Group = factor(sankey_group_levels, levels = sankey_group_levels)
    )
    p <- p +
      geom_point(
        data = legend_df, aes(x = x, y = y, colour = Group),
        inherit.aes = FALSE, alpha = 0, size = 0
      ) +
      scale_colour_manual(
        values = group_colours[sankey_group_levels],
        name = group_axis_lab,
        guide = guide_legend(
          nrow = 1, title.position = "left",
          override.aes = list(alpha = 1, size = 3.2, shape = 15)
        )
      ) +
      theme(
        legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_text(size = 6.5, face = "bold"),
        legend.text = element_text(size = 6.2),
        legend.key.size = unit(3.2, "mm"),
        legend.margin = margin(4, 0, 0, 0)
      )
  }

  n_dz <- n_distinct(dat$Disease)
  if (compact) {
    # Same aspect as the 174 x 190 mm figure, with the empty left gutter removed.
    # The abbrev figure drops the top and bottom headroom and narrows the
    # gaps between Disease, Continent and 16S region.
    width_mm <- fig_w
    height_mm <- if (!is.null(compact_height)) compact_height else
      if (tighter) 116 else 131
  } else {
    # 2026-10-10: doubled width, halved height relative to the original
    # 174 x max(160, 20 + n_dz * 8.5) mm layout.
    width_mm <- 2 * W_FULL_MM
    height_mm <- min(MAX_H_MM, max(160, 20 + n_dz * 8.5)) / 2
  }
  save_panel(p, stem, width_mm = width_mm, height_mm = height_mm)
  cat(sprintf(
    "%s: %d BioProjects | %d diseases | %d axes | height %.0f mm\n",
    stem, n_distinct(dat$bioproject_id), n_dz, n_axis, height_mm
  ))
  invisible(p)
}

build_sankey(
  ANATOMICAL_MEMBERS, ANATOMICAL_SHORT, "Anatomical site",
  show_group_axis = FALSE, group_colours = ANATOMICAL_COLOURS,
  stem = "panel_A_sankey", fixed_order = TRUE,
  abbreviate_cancers = TRUE, colour_by = "disease", show_site = TRUE,
  rainbow = TRUE,
  extra_labels = c(
    "Periodontal" = "PD",
    "Dental Caries" = "DC",
    "Obesity" = "OB",
    "Diabetes" = "DB",
    "Hypertension" = "HTN"
  )
)
build_sankey(
  ANATOMICAL_MEMBERS, ANATOMICAL_SHORT, "Anatomical site",
  show_group_axis = FALSE, group_colours = ANATOMICAL_COLOURS,
  stem = "panel_A_sankey_no_category", fixed_order = TRUE,
  show_group_legend = TRUE, abbreviate_cancers = TRUE
)
build_sankey(
  ANATOMICAL_MEMBERS, ANATOMICAL_SHORT, "Anatomical site",
  show_group_axis = FALSE, group_colours = ANATOMICAL_COLOURS,
  stem = "panel_A_sankey_by_disease", fixed_order = TRUE,
  abbreviate_cancers = TRUE, colour_by = "disease", compact = TRUE
)
build_sankey(
  ANATOMICAL_MEMBERS, ANATOMICAL_SHORT, "Anatomical site",
  show_group_axis = FALSE, group_colours = ANATOMICAL_COLOURS,
  stem = "panel_A_sankey_by_disease_abbrev", fixed_order = TRUE,
  abbreviate_cancers = TRUE, colour_by = "disease", compact = TRUE,
  tighter = TRUE, show_site = TRUE, compact_width = 240, font_scale = 2.5,
  bar_scale = 4.5,
  # Stretched columns: every column gets the Disease column's total height,
  # with the boxes inside keeping their BioProject-count proportions.
  equal_column_height = TRUE, compact_height = 160,
  # Publication polish: muted low-saturation ribbon palette at alpha 0.65,
  # uniform light-grey nodes with a hairline border and #222222 label text.
  # Node order stays exactly as before (no crossing-reducing re-sort).
  flow_alpha = 0.65, muted_palette = TRUE, optimize_order = FALSE,
  node_fill = "#E6E6E6", node_border = "#BDBDBD", node_linewidth = 0.5,
  label_colour = "#222222",
  # Fixed top-to-bottom node order for the Continent column.
  continent_order = c("AS", "NA", "EU", "OC", "SA", "AF"),
  extra_labels = c(
    "Periodontal" = "PD",
    "Dental Caries" = "DC",
    "Hypertension" = "HTN"
  )
)
build_sankey(
  ANATOMICAL_MEMBERS, ANATOMICAL_SHORT, "Anatomical site",
  show_group_axis = TRUE, group_colours = ANATOMICAL_COLOURS,
  stem = "panel_A_sankey_AnatomicalGroup", fixed_order = TRUE
)
