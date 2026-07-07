# =========================================================================
# 1. Load required packages
# =========================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrastr)   # rasterise dense scatter so reds and blues stop swallowing each other
  library(Cairo)     # raster backend for ggrastr + PDF backend on headless HPC nodes
})

# =========================================================================
# 2. Define paths and variables
# =========================================================================
# Input paths for both Weighted and Unweighted PCoA coordinates
input_unweight <- "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/beta-diversity/input/unifrac-unweight-pcoa-R-input-0602.csv"
input_weight <- "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/beta-diversity/input/unifrac-weight-pcoa-R-input-0602.csv"
input_permanova_weight <- "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/beta-diversity/input/03.02-unifrac-weight-distance-matrix-permanova-0604.tsv"
input_permanova_unweight <- "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/beta-diversity/input/03.02-unifrac-unweight-distance-matrix-permanova-0604.tsv"
input_meta <- "/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/metadata_oral_qc.tsv"
out_dir <- "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/beta-diversity/results/results-color/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

color_vars <- c("age_cleaned_group", "bmi_cleaned_group", "sex_cleaned",
                "platform_cleaned", "sequencing_region_cleaned",
                "continent_group", "site_cleaned")

pair_colors <- c("#CC3311", "#0072B2")

# Assign one fixed, high-contrast color per group across all panels.
fallback_palette <- c(
  "#0072B2", # blue
  "#E69F00", # orange
  "#009E73", # green
  "#CC3311", # red
  "#7F3C8D", # purple
  "#8C564B", # brown
  "#808000", # olive
  "#7F7F7F", # gray
  "#6A3D9A", # deep violet
  "#A6761D", # dark gold-brown
  "#1B1B1B", # charcoal
  "#B8860B"  # dark goldenrod
)

# =========================================================================
# 3. Get Slurm Array Task ID and set independent log
# =========================================================================
task_id_str <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id_str == "") {
  task_id <- 1
  warning("SLURM_ARRAY_TASK_ID not detected; defaulting to task 1.")
} else {
  task_id <- as.integer(task_id_str)
}

current_var <- color_vars[task_id]

if (is.na(current_var)) {
  stop(paste("Task ID", task_id, "is outside the variable list. Check Slurm --array."))
}

log_file <- file.path(out_dir, paste0("pcoa_matrix_drawing_log_", current_var, ".txt"))

write_log <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- paste0("[", timestamp, "] ", msg)
  cat(full_msg, "\n")
  cat(full_msg, "\n", file = log_file, append = TRUE)
}

# =========================================================================
# 4. Read data and clean metadata
# =========================================================================
write_log(paste("=== Start: Matrix PCoA drawing | target variable:", current_var, "==="))
write_log("Reading input files...")
dt_pcoa_un <- fread(input_unweight)
dt_pcoa_wt <- fread(input_weight)
dt_perm_wt <- fread(input_permanova_weight)
dt_perm_un <- fread(input_permanova_unweight)
dt_meta <- fread(input_meta)

normalize_group_str <- function(x) {
  x <- gsub("＞", ">", x, fixed = TRUE)
  x <- gsub("＜", "<", x, fixed = TRUE)
  x <- gsub("≥", ">=", x, fixed = TRUE)
  x <- gsub("≧", ">=", x, fixed = TRUE)
  x <- trimws(x)
  x <- ifelse(x %in% c("F", "Female"), "female", x)
  x <- ifelse(x %in% c("M", "Male"), "male", x)
  x
}
dt_perm_un[, Group1 := normalize_group_str(Group1)]
dt_perm_un[, Group2 := normalize_group_str(Group2)]
dt_perm_wt[, Group1 := normalize_group_str(Group1)]
dt_perm_wt[, Group2 := normalize_group_str(Group2)]

dup_idx <- duplicated(names(dt_meta))
if (any(dup_idx)) {
  dt_meta <- dt_meta[, !dup_idx, with = FALSE]
}

if ("bioproject_id" %in% names(dt_meta) && !"bioproject" %in% names(dt_meta)) {
  setnames(dt_meta, "bioproject_id", "bioproject")
}

# Fix full-width and special character replacement logic
dt_meta[, age_cleaned_group := gsub("＞", ">", age_cleaned_group)]
dt_meta[, age_cleaned_group := gsub("＜", "<", age_cleaned_group)]
dt_meta[, bmi_cleaned_group := gsub("＜", "<", bmi_cleaned_group)]
dt_meta[, bmi_cleaned_group := gsub("≥", ">=", bmi_cleaned_group)]
dt_meta[, bmi_cleaned_group := gsub("≧", ">=", bmi_cleaned_group)]

# =========================================================================
# 5. Merge data
# =========================================================================
write_log("Merging data...")
plot_data_un <- merge(dt_pcoa_un, dt_meta, by.x = "SampleID", by.y = "sample_id", all.x = TRUE)
plot_data_wt <- merge(dt_pcoa_wt, dt_meta, by.x = "SampleID", by.y = "sample_id", all.x = TRUE)

rm(dt_meta, dt_pcoa_un, dt_pcoa_wt)
gc()

custom_orders <- list(
  age_cleaned_group = c("0-6 months", "6-36 months", "3-12 years", "12-18 years", "18-60 years", ">60 years"),
  bmi_cleaned_group = c("<18.5", "18.5-24.0", "24.0-28.0", ">=28.0"),
  continent_group = c("North America", "Asia", "Europe", "Africa", "Oceania", "South America"),
  platform_cleaned = c("Illumina MiSeq", "Illumina HiSeq", "Illumina NextSeq", "Illumina NovaSeq", "Illumina GAIIx", "454 GS FLX Titanium", "Ion Torrent"),
  site_cleaned = c("saliva", "tooth", "oral mucosa", "oral swab", "tongue", "gingiva", "other")
)

title_labels <- list(
  age_cleaned_group = "PCoA: age group",
  bmi_cleaned_group = "PCoA: bmi group",
  sex_cleaned = "PCoA: sex",
  platform_cleaned = "PCoA: platform",
  sequencing_region_cleaned = "PCoA: sequencing region",
  continent_group = "PCoA: continent",
  site_cleaned = "PCoA: site"
)

if (current_var %in% names(custom_orders)) {
  existing_levels_un <- unique(as.character(plot_data_un[[current_var]]))
  existing_levels_wt <- unique(as.character(plot_data_wt[[current_var]]))
  existing_levels <- unique(c(existing_levels_un, existing_levels_wt))
  existing_levels <- existing_levels[!is.na(existing_levels) & existing_levels != ""]
  final_levels <- c(intersect(custom_orders[[current_var]], existing_levels),
                    setdiff(existing_levels, custom_orders[[current_var]]))
  plot_data_un[[current_var]] <- factor(plot_data_un[[current_var]], levels = final_levels)
  plot_data_wt[[current_var]] <- factor(plot_data_wt[[current_var]], levels = final_levels)
} else {
  plot_data_un[[current_var]] <- as.factor(as.character(plot_data_un[[current_var]]))
  plot_data_wt[[current_var]] <- as.factor(as.character(plot_data_wt[[current_var]]))
}

# -------------------------------------------------------------------------
# Per-variable group-size filter
# Sequencing-region groups with fewer than this many samples are dropped
# from the matrix plot — they otherwise produce panels too sparse to
# interpret. PERMANOVA values shown in retained panels remain those
# pre-computed on the FULL data by the upstream script (unchanged).
# -------------------------------------------------------------------------
min_n_by_var <- list(sequencing_region_cleaned = 500)
if (current_var %in% names(min_n_by_var)) {
  min_n <- min_n_by_var[[current_var]]
  counts_un <- table(as.character(plot_data_un[[current_var]]))
  counts_wt <- table(as.character(plot_data_wt[[current_var]]))
  keep_un <- names(counts_un)[counts_un >= min_n]
  keep_wt <- names(counts_wt)[counts_wt >= min_n]
  keep_levels <- intersect(keep_un, keep_wt)
  dropped <- setdiff(union(names(counts_un), names(counts_wt)), keep_levels)
  if (length(dropped) > 0) {
    write_log(paste0("Filter ", current_var, ": dropping groups with n<", min_n,
                     " -> [", paste(dropped, collapse = ", "), "]"))
  }
  plot_data_un <- plot_data_un[as.character(plot_data_un[[current_var]]) %in% keep_levels, ]
  plot_data_wt <- plot_data_wt[as.character(plot_data_wt[[current_var]]) %in% keep_levels, ]
  plot_data_un[[current_var]] <- factor(as.character(plot_data_un[[current_var]]),
                                        levels = intersect(levels(plot_data_un[[current_var]]), keep_levels))
  plot_data_wt[[current_var]] <- factor(as.character(plot_data_wt[[current_var]]),
                                        levels = intersect(levels(plot_data_wt[[current_var]]), keep_levels))
}

# -------------------------------------------------------------------------
# Per-variable group exclusion (by name)
# For the site matrix the "other" bucket is a heterogeneous catch-all rather
# than a single biologically meaningful oral site, so it is dropped entirely:
# it gets no panel row/column and none of its samples are plotted. PERMANOVA
# values shown in retained panels remain those pre-computed on the FULL data
# by the upstream script (unchanged).
# -------------------------------------------------------------------------
drop_levels_by_var <- list(site_cleaned = c("other"))
if (current_var %in% names(drop_levels_by_var)) {
  drop_lv <- drop_levels_by_var[[current_var]]
  plot_data_un <- plot_data_un[!(as.character(plot_data_un[[current_var]]) %in% drop_lv), ]
  plot_data_wt <- plot_data_wt[!(as.character(plot_data_wt[[current_var]]) %in% drop_lv), ]
  plot_data_un[[current_var]] <- factor(as.character(plot_data_un[[current_var]]),
                                        levels = setdiff(levels(plot_data_un[[current_var]]), drop_lv))
  plot_data_wt[[current_var]] <- factor(as.character(plot_data_wt[[current_var]]),
                                        levels = setdiff(levels(plot_data_wt[[current_var]]), drop_lv))
  write_log(paste0("Exclude ", current_var, ": removing group(s) [",
                   paste(drop_lv, collapse = ", "), "] from the matrix"))
}

# Base theme
custom_theme <- theme_bw() + theme(
  text = element_text(face = "bold", size = 16), # Increased size
  plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), # Increased size
  axis.title = element_blank(), # Matrix plots look cleaner without repetitive axis titles
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "none", # Legends are hidden in subplots to keep matrix clean
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
  plot.margin = margin(2, 2, 2, 2)
)

# =========================================================================
# 6. Matrix Plot Functions
# =========================================================================

# Look up the PERMANOVA row for a given variable and unordered group pair.
lookup_permanova <- function(perm_dt, var_name, g1, g2) {
  g1_norm <- normalize_group_str(as.character(g1))
  g2_norm <- normalize_group_str(as.character(g2))
  matched_row <- perm_dt[
    Variable == var_name &
      ((Group1 == g1_norm & Group2 == g2_norm) |
       (Group1 == g2_norm & Group2 == g1_norm))
  ]
  if (nrow(matched_row) == 0) return(NULL)
  matched_row[1]
}

# Format a PERMANOVA row into a compact label.
format_perm_text <- function(perm_row) {
  if (is.null(perm_row) || nrow(perm_row) == 0) return(NA_character_)
  p_value <- as.numeric(perm_row$p_value)
  # Report P as-computed. With 999 permutations the floor value is 1/(999+1)
  # = 0.001 \u2014 writing "P < 0.001" would over-claim the test resolution.
  p_label <- if (is.finite(p_value)) sprintf("P = %.3f", p_value) else "P = NA"
  sprintf("R\u00b2 = %.4f\n%s", as.numeric(perm_row$R2), p_label)
}

# Generate Pairwise off-diagonal plot
generate_pairwise_plot <- function(data, var_name, g1, g2, type_label, perm_dt) {
  sub_data <- data[data[[var_name]] %in% c(g1, g2), ]
  if (nrow(sub_data) == 0) return(plot_spacer())

  sub_data[[var_name]] <- factor(sub_data[[var_name]], levels = c(g1, g2))
  group_counts <- sort(table(sub_data[[var_name]]), decreasing = TRUE)
  bottom_group <- names(group_counts)[1]                       # majority
  top_group <- names(group_counts)[length(group_counts)]       # minority
  bottom_data <- sub_data[sub_data[[var_name]] == bottom_group, ]
  top_data <- sub_data[sub_data[[var_name]] == top_group, ]

  # Stable colour per group — looked up by group NAME from the global
  # group_colors map populated from named_palette_by_var. This guarantees
  # the same group always uses the same colour across every panel of the
  # matrix, regardless of which two groups are being compared or which is
  # the majority / minority.
  bottom_color <- unname(group_colors[as.character(bottom_group)])
  top_color    <- unname(group_colors[as.character(top_group)])

  if (nrow(bottom_data) > 1) bottom_data <- bottom_data[sample.int(nrow(bottom_data)), ]
  if (nrow(top_data) > 1) top_data <- top_data[sample.int(nrow(top_data)), ]

  pt_size  <- 1.8
  # Aggressively lowered alpha to defeat colour-swamping in dense panels:
  # even after stacking ~100 majority points on the same pixel, alpha
  # accumulation stays well below saturation so the minority layer on top
  # remains visible. Both groups are translucent enough that their overlap
  # regions read as a true blend rather than one colour winning outright.
  pt_alpha_bottom <- 0.15
  pt_alpha_top    <- 0.45
  raster_dpi <- 600

  # Look up PERMANOVA early — we need its presence to decide whether to carve
  # out a top-right exclusion zone for the label box so points never overlap it.
  perm_row <- lookup_permanova(perm_dt, var_name, g1, g2)
  perm_text <- format_perm_text(perm_row)

  x_range <- range(sub_data$PC1, na.rm = TRUE)
  y_range <- range(sub_data$PC2, na.rm = TRUE)
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  if (!is.finite(x_span) || x_span == 0) x_span <- 1
  if (!is.finite(y_span) || y_span == 0) y_span <- 1

  # Drop the few points that would land underneath the PERMANOVA label box so
  # the box never visually overlaps any sample dot. The box itself sits in the
  # outer top-right padding; we only need to clear a thin strip right at the
  # inner top-right corner of the data range.
  if (!is.na(perm_text)) {
    excl_x_min <- x_range[2] - x_span * 0.12
    excl_y_min <- y_range[2] - y_span * 0.12
    bottom_data <- bottom_data[!(PC1 >= excl_x_min & PC2 >= excl_y_min)]
    top_data    <- top_data[!(PC1 >= excl_x_min & PC2 >= excl_y_min)]
  }

  # -----------------------------------------------------------------------
  # Summary geometry layer to defeat the "one-group-only" illusion in dense
  # panels. Even when the two faint point clouds fully overlap and occlude
  # each other, a 95% confidence ellipse + centroid for EACH group is drawn
  # on the very top layer (vector, never rasterised, so it is never buried by
  # any point). The ellipse encloses ~95% of that group's samples and is
  # computed from all of the group's points; its colour is the group's stable
  # palette colour, so the two ellipses' overlap area reads directly as the
  # between-group similarity. Guarded by a minimum point count because an
  # ellipse cannot be fit to a near-empty / collinear group.
  ellipse_layers <- list()
  if (nrow(bottom_data) >= 5) {
    ellipse_layers <- c(ellipse_layers, list(
      stat_ellipse(
        data = bottom_data, aes(x = PC1, y = PC2),
        colour = bottom_color, type = "t", level = 0.95, linewidth = 1.0
      )))
  }
  if (nrow(top_data) >= 5) {
    ellipse_layers <- c(ellipse_layers, list(
      stat_ellipse(
        data = top_data, aes(x = PC1, y = PC2),
        colour = top_color, type = "t", level = 0.95, linewidth = 1.0
      )))
  }

  p <- ggplot() +
    custom_theme +
    # Majority group drawn first — semi-transparent so its pixels never
    # saturate to fully opaque and bury the minority on top
    rasterise(
      geom_point(
        data = bottom_data,
        aes(x = PC1, y = PC2),
        colour = bottom_color,
        shape = 16, size = pt_size, alpha = pt_alpha_bottom, stroke = 0
      ),
      dpi = raster_dpi, dev = "cairo"
    ) +
    # Minority group drawn LAST — slightly higher alpha than the bottom layer
    # so it reads clearly, but still translucent enough to let the bottom
    # colour show through in zones where the two clouds overlap
    rasterise(
      geom_point(
        data = top_data,
        aes(x = PC1, y = PC2),
        colour = top_color,
        shape = 16, size = pt_size, alpha = pt_alpha_top, stroke = 0
      ),
      dpi = raster_dpi, dev = "cairo"
    ) +
    # 95% confidence ellipses (one per group) drawn on top of both point
    # clouds — vector lines, so they are never occluded by the dense scatter
    ellipse_layers +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    coord_cartesian(clip = "off")

  if (!is.na(perm_text)) {
    # Box pinned tight to the panel's top-right OUTER corner: tiny right/top
    # padding (×0.22) gives the box just enough room without leaving a big
    # empty quadrant, and the box itself (×0.18 offset) sits right at the
    # very corner of the expanded panel. Left/bottom padding kept small
    # (×0.10) so the scatter cloud occupies the panel centre.
    box_x <- x_range[2] + x_span * 0.18
    box_y <- y_range[2] + y_span * 0.18

    p <- p +
      expand_limits(
        x = c(x_range[1] - x_span * 0.10, x_range[2] + x_span * 0.22),
        y = c(y_range[1] - y_span * 0.10, y_range[2] + y_span * 0.22)
      ) +
      annotate(
        "label",
        x = box_x, y = box_y,
        label = perm_text,
        hjust = 1, vjust = 1,
        size = 3.8,
        fontface = "bold",
        color = "black",
        fill = "white",
        alpha = 0.92,
        label.padding = grid::unit(0.20, "lines"),
        label.r = grid::unit(0.10, "lines"),
        label.size = 0.25
      )
  }

  return(p)
}

# Generate Diagonal Plot (single group with text label coloured per group)
generate_diagonal_plot <- function(data, var_name, g, text_color = "black") {
  g_spaced <- gsub(",\\s*", ", ", as.character(g))
  wrapped_g <- paste(strwrap(g_spaced, width = 8), collapse = "\n")

  p <- ggplot() +
    custom_theme +
    annotate("text", x = 0.5,
             y = 0.5,
             label = wrapped_g, size = 8, fontface = "bold", color = text_color) + # Increased size
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)

  return(p)
}

# =========================================================================
# 7. Construct the N x N Matrix
# =========================================================================
all_levels <- levels(plot_data_un[[current_var]])
all_levels <- all_levels[!is.na(all_levels) & all_levels != ""]
n_levels <- length(all_levels)

named_palette_by_var <- list(
  # age: 6 groups — hues evenly spaced around the colour wheel (~60° apart),
  # each from a different family (deep blue / orange / green / red / purple / brown)
  age_cleaned_group = c(
    "0-6 months"  = "#1F77B4",  # strong blue
    "6-36 months" = "#FF7F0E",  # vivid orange
    "3-12 years"  = "#2CA02C",  # green
    "12-18 years" = "#D62728",  # red
    "18-60 years" = "#9467BD",  # purple
    ">60 years"   = "#8C564B"   # brown
  ),
  # bmi: 4 groups — cool to warm gradient with maximum hue separation
  bmi_cleaned_group = c(
    "<18.5"      = "#00A0B0",   # teal
    "18.5-24.0"  = "#6A8E3E",   # olive green
    "24.0-28.0"  = "#F2A300",   # amber
    ">=28.0"     = "#D62728"    # red
  ),
  # sex: 2-3 groups — magenta vs sky-blue (~180° apart)
  sex_cleaned = c(
    "female"  = "#E63946",      # magenta-red
    "male"    = "#1D9BD8",      # sky blue
    "Female"  = "#E63946",
    "Male"    = "#1D9BD8",
    "unknown" = "#6E6E6E",
    "Unknown" = "#6E6E6E"
  ),
  # platform: 7 distinct families from Set1 + custom tweaks
  platform_cleaned = c(
    "Illumina MiSeq"      = "#E41A1C",  # red
    "Illumina HiSeq"      = "#377EB8",  # blue
    "Illumina NextSeq"    = "#4DAF4A",  # green
    "Illumina NovaSeq"    = "#984EA3",  # purple
    "Illumina GAIIx"      = "#FF7F00",  # orange
    "454 GS FLX Titanium" = "#A65628",  # brown
    "Ion Torrent"         = "#F781BF"   # pink
  ),
  # sequencing region: 8 high-contrast hues (no overlap with each other)
  sequencing_region_cleaned = c(
    "V1-V2"       = "#E41A1C",   # red
    "V1-V3"       = "#377EB8",   # blue
    "V3-V4"       = "#4DAF4A",   # green
    "V4"          = "#FF7F00",   # orange
    "V4-V5"       = "#984EA3",   # purple
    "V5-V7"       = "#A65628",   # brown
    "V6-V8"       = "#F781BF",   # pink (same as platform's Ion Torrent) — distinct from the V3-V4 green
    "full-length" = "#666666"    # dark grey
  ),
  # continent: 6 highly distinct hues, no cool-cluster collisions
  continent_group = c(
    "North America" = "#E41A1C",   # red
    "Asia"          = "#FF7F00",   # orange
    "Europe"        = "#377EB8",   # blue
    "Africa"        = "#4DAF4A",   # green
    "Oceania"       = "#984EA3",   # purple
    "South America" = "#A65628"    # brown
  ),
  # oral sites: 7 highly distinct hues, no warm/cool cluster collisions
  site_cleaned = c(
    "saliva"      = "#E41A1C",     # red
    "tooth"       = "#377EB8",     # blue
    "oral mucosa" = "#4DAF4A",     # green
    "oral swab"   = "#FF7F00",     # orange
    "tongue"      = "#984EA3",     # purple
    "gingiva"     = "#A65628",     # brown
    "other"       = "#E7298C"      # magenta-pink — high contrast vs the other site hues
  )
)

current_named_palette <- named_palette_by_var[[current_var]]
if (is.null(current_named_palette)) {
  current_named_palette <- c()
}

group_colors <- setNames(rep(NA_character_, n_levels), all_levels)
matched_levels <- intersect(all_levels, names(current_named_palette))
group_colors[matched_levels] <- unname(current_named_palette[matched_levels])

missing_levels <- names(group_colors)[is.na(group_colors)]
if (length(missing_levels) > 0) {
  extra_pool <- unique(c(
    fallback_palette,
    grDevices::hcl.colors(max(length(missing_levels), 12), palette = "Dark 3")
  ))
  extra_pool <- setdiff(extra_pool, unname(group_colors[!is.na(group_colors)]))
  group_colors[missing_levels] <- extra_pool[seq_along(missing_levels)]
}

if (n_levels > 1) {
  write_log(paste("Building", n_levels, "x", n_levels, "Matrix plot..."))

  plot_list <- list()

  for (i in 1:n_levels) {
    for (j in 1:n_levels) {
      g_row <- all_levels[i]
      g_col <- all_levels[j]

      if (i == j) {
        # Diagonal: Single Group — label text coloured with the group's
        # assigned palette colour so the reader can read off colour <-> group
        # mapping straight from the matrix diagonal.
        p <- generate_diagonal_plot(plot_data_un, current_var, g_row,
                                    text_color = unname(group_colors[as.character(g_row)]))
      } else if (i > j) {
        # Lower Triangle: Unweighted
        p <- generate_pairwise_plot(plot_data_un, current_var, g_row, g_col, "Unweighted", dt_perm_un)
      } else {
        # Upper Triangle: Weighted
        p <- generate_pairwise_plot(plot_data_wt, current_var, g_row, g_col, "Weighted", dt_perm_wt)
      }

      plot_list <- append(plot_list, list(p))
    }
  }

  # Assemble using Patchwork
  write_log("Patching plots together...")
  subtitle_text <- paste0(
    "Upper Triangle: Weighted UniFrac | Lower Triangle: Unweighted UniFrac",
    "\nEllipses: 95% confidence ellipse per group"
  )
  if (current_var %in% names(min_n_by_var)) {
    subtitle_text <- paste0(subtitle_text,
                            "\nGroups with n < ", min_n_by_var[[current_var]],
                            " excluded from the matrix")
  }
  matrix_plot <- wrap_plots(plot_list, ncol = n_levels, nrow = n_levels) +
    plot_annotation(
      title = title_labels[[current_var]],
      subtitle = subtitle_text,
      theme = theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5), # Increased size
                    plot.subtitle = element_text(size = 18, hjust = 0.5)) # Increased size
    )

  # Save the output
  plot_size <- n_levels * 3 # 3 inches per subplot dynamically
  out_pdf <- file.path(out_dir, paste0("PCoA-Matrix-", current_var, "-0610.pdf"))
  out_png <- file.path(out_dir, paste0("PCoA-Matrix-", current_var, "-0610.png"))

  write_log("Saving PDF...")
  ggsave(out_pdf, plot = matrix_plot, width = plot_size, height = plot_size, limitsize = FALSE, device = cairo_pdf)

  write_log("Saving PNG...")
  ggsave(out_png, plot = matrix_plot, width = plot_size, height = plot_size, dpi = 300, device = "png", type = "cairo", bg = "white", limitsize = FALSE)

  write_log(paste("Matrix PCoA saved successfully at:", plot_size, "x", plot_size, "inches"))

} else {
  write_log(paste("Warning: variable", current_var, "has less than 2 valid groups; no matrix generated."))
}

write_log("=== Task finished ===")
