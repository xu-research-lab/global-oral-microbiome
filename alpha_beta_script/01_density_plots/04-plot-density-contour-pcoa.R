# =========================================================================
# 1. Load required packages
# =========================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(Cairo)
})

# =========================================================================
# 2. Define paths and variables
# =========================================================================
input_unweight <- "/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/beta_pcoa_out_20260909/unifrac_unweight/unifrac-unweight-pcoa-R-input-20260909.csv"
input_weight <- "/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/beta_pcoa_out_20260909/unifrac_weight/unifrac-weight-pcoa-R-input-20260909.csv"
input_permanova_weight <- "/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/beta_pcoa_out_20260909/unifrac_weight/02-unifrac-weight-distance-matrix-permanova-20260909.tsv"
input_permanova_unweight <- "/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/beta_pcoa_out_20260909/unifrac_unweight/02-unifrac-unweight-distance-matrix-permanova-20260909.tsv"
input_meta <- "/home/hexincheng/projects/hxc/oral_project/data/metadata_oral_qc_dedup.tsv"
out_dir <- "/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/beta_pcoa_out_20260909/pcoa_plots/density/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

color_vars <- c("age_cleaned_group", "bmi_cleaned_group", "sex_cleaned",
                "platform_cleaned", "sequencing_region_cleaned",
                "continent_group", "site_cleaned")

# Only two colours per pairwise panel: row group = red, column group = blue
color_row <- "#CC3311"
color_col <- "#0072B2"

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

log_file <- file.path(out_dir, paste0("pcoa_density_log_", current_var, ".txt"))

write_log <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- paste0("[", timestamp, "] ", msg)
  cat(full_msg, "\n")
  cat(full_msg, "\n", file = log_file, append = TRUE)
}

# =========================================================================
# 4. Read data and clean metadata
# =========================================================================
write_log(paste("=== Start: Density contour PCoA | target variable:", current_var, "==="))
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

dt_meta[, age_cleaned_group := gsub("＞", ">", age_cleaned_group)]
dt_meta[, age_cleaned_group := gsub("＜", "<", age_cleaned_group)]
dt_meta[, bmi_cleaned_group := gsub("＜", "<", bmi_cleaned_group)]
dt_meta[, bmi_cleaned_group := gsub("≥", ">=", bmi_cleaned_group)]
dt_meta[, bmi_cleaned_group := gsub("≧", ">=", bmi_cleaned_group)]

# Normalise site labels to lowercase. The dedup metadata stores the
# heterogeneous catch-all as "Other" (capital O), while custom_orders and
# drop_levels_by_var below use lowercase "other". The drop is case-sensitive,
# so without this "Other" survives and an extra catch-all panel leaks into
# the site matrix (7x7 instead of 6x6). Lowercasing is safe — every other
# site label is already lowercase.
dt_meta[, site_cleaned := tolower(site_cleaned)]

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
  text = element_text(face = "bold", size = 16),
  plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "none",
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
  plot.margin = margin(2, 2, 2, 2)
)

# =========================================================================
# 6. Helper Functions
# =========================================================================

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

format_perm_text <- function(perm_row) {
  if (is.null(perm_row) || nrow(perm_row) == 0) return(NA_character_)
  p_value <- as.numeric(perm_row$p_value)
  p_label <- if (is.finite(p_value)) sprintf("P = %.3f", p_value) else "P = NA"
  sprintf("R² = %.4f\n%s", as.numeric(perm_row$R2), p_label)
}

# -------------------------------------------------------------------------
# Density contour pairwise plot
# Two semi-transparent filled-contour layers (one per group) show where
# each group's mass actually sits. Fixed alpha per band; inner bands
# accumulate opacity through overlap (6 bands x 0.15 ~ 0.9 at centre).
# Only two colours per panel: row group = red, column group = blue.
# -------------------------------------------------------------------------
generate_density_pairwise <- function(data, var_name, g_row, g_col,
                                      perm_dt, row_color = color_row, col_color = color_col) {
  sub_data <- data[data[[var_name]] %in% c(g_row, g_col), ]
  if (nrow(sub_data) < 10) return(plot_spacer())

  data_row <- sub_data[sub_data[[var_name]] == g_row, ]
  data_col <- sub_data[sub_data[[var_name]] == g_col, ]
  if (nrow(data_row) < 5 || nrow(data_col) < 5) return(plot_spacer())

  perm_row <- lookup_permanova(perm_dt, var_name, g_row, g_col)
  perm_text <- format_perm_text(perm_row)

  p <- ggplot() + custom_theme

  # Row group density contours (red) — drawn first, background layer
  if (nrow(data_row) >= 5) {
    p <- p + stat_density_2d(
      data = data_row, aes(x = PC1, y = PC2),
      geom = "polygon", fill = row_color, color = row_color,
      alpha = 0.15, linewidth = 0.3, bins = 6, contour = TRUE
    )
  }

  # Column group density contours (blue) — drawn second, foreground layer
  if (nrow(data_col) >= 5) {
    p <- p + stat_density_2d(
      data = data_col, aes(x = PC1, y = PC2),
      geom = "polygon", fill = col_color, color = col_color,
      alpha = 0.15, linewidth = 0.3, bins = 6, contour = TRUE
    )
  }

  p <- p +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    coord_cartesian(clip = "off")

  # PERMANOVA annotation
  if (!is.na(perm_text)) {
    x_range <- range(sub_data$PC1, na.rm = TRUE)
    y_range <- range(sub_data$PC2, na.rm = TRUE)
    x_span <- diff(x_range)
    y_span <- diff(y_range)
    if (!is.finite(x_span) || x_span == 0) x_span <- 1
    if (!is.finite(y_span) || y_span == 0) y_span <- 1

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

# Diagonal: group label in black text
generate_diagonal_plot <- function(g) {
  g_spaced <- gsub(",\\s*", ", ", as.character(g))
  wrapped_g <- paste(strwrap(g_spaced, width = 8), collapse = "\n")

  p <- ggplot() +
    custom_theme +
    annotate("text", x = 0.5, y = 0.5,
             label = wrapped_g, size = 8, fontface = "bold", color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)

  return(p)
}

# =========================================================================
# 7. Build density contour matrix
# =========================================================================
all_levels <- levels(plot_data_un[[current_var]])
all_levels <- all_levels[!is.na(all_levels) & all_levels != ""]
n_levels <- length(all_levels)

if (n_levels > 1) {
  write_log(paste("Building", n_levels, "x", n_levels, "Density contour matrix..."))
  density_list <- list()

  for (i in 1:n_levels) {
    for (j in 1:n_levels) {
      g_row <- all_levels[i]
      g_col <- all_levels[j]

      if (i == j) {
        p <- generate_diagonal_plot(g_row)
      } else if (i > j) {
        p <- generate_density_pairwise(plot_data_un, current_var, g_row, g_col, dt_perm_un)
      } else {
        p <- generate_density_pairwise(plot_data_wt, current_var, g_row, g_col, dt_perm_wt)
      }

      density_list <- append(density_list, list(p))
    }
  }

  write_log("Patching density contour plots together...")
  density_subtitle <- paste0(
    "Upper Triangle: Weighted UniFrac | Lower Triangle: Unweighted UniFrac",
    "\nDensity contours: two semi-transparent layers per group pair (red = row, blue = column)"
  )
  if (current_var %in% names(min_n_by_var)) {
    density_subtitle <- paste0(density_subtitle,
                              "\nGroups with n < ", min_n_by_var[[current_var]],
                              " excluded from the matrix")
  }
  density_plot <- wrap_plots(density_list, ncol = n_levels, nrow = n_levels) +
    plot_annotation(
      title = title_labels[[current_var]],
      subtitle = density_subtitle,
      theme = theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 18, hjust = 0.5))
    )

  plot_size <- n_levels * 3
  density_pdf <- file.path(out_dir, paste0("PCoA-Density-", current_var, "-20260909.pdf"))
  density_png <- file.path(out_dir, paste0("PCoA-Density-", current_var, "-20260909.png"))

  write_log("Saving density contour PDF...")
  ggsave(density_pdf, plot = density_plot, width = plot_size, height = plot_size,
         limitsize = FALSE, device = cairo_pdf)

  write_log("Saving density contour PNG...")
  ggsave(density_png, plot = density_plot, width = plot_size, height = plot_size,
         dpi = 300, device = "png", type = "cairo", bg = "white", limitsize = FALSE)

  write_log(paste("Density contour matrix saved at:", plot_size, "x", plot_size, "inches"))

} else {
  write_log(paste("Warning: variable", current_var, "has less than 2 valid groups; no matrix generated."))
}

write_log("=== Task finished ===")