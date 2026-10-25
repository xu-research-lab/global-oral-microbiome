suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

find_script_dir <- function() {
  if (interactive()) return(normalizePath(getwd()))
  args <- commandArgs(trailingOnly = FALSE)
  f <- args[grep("^--file=", args)]
  if (length(f) == 1) return(dirname(normalizePath(sub("^--file=", "", f))))
  normalizePath(getwd())
}
script_dir <- find_script_dir()


SERVER_ROOT <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-map-2/04-patch-pipeline-map/03-patch-pipeline-map"
SERVER_METADATA_DEDUP <- file.path(
  SERVER_ROOT, "input", "metadata_oral_qc_dedup.tsv"
)
SERVER_OUT_DIR <- file.path(
  SERVER_ROOT, "overview_figure_1", "Fig1-0926", "results", "Fig1"
)

LOCAL_BASE <- normalizePath(file.path(script_dir, "..", ".."))
on_server  <- dir.exists(file.path(SERVER_ROOT, "input"))
BASE_DIR   <- if (on_server) SERVER_ROOT else LOCAL_BASE

PREP_DIR <- file.path(script_dir, "input")
input_permanova_unweight <- Sys.getenv(
  "FIG1_PERMANOVA_TABLE",
  unset = file.path(PREP_DIR,
                    "02-unifrac-unweight-distance-matrix-permanova-20260909.tsv")
)
input_meta <- Sys.getenv(
  "FIG1_METADATA",
  unset = if (on_server) {
    SERVER_METADATA_DEDUP
  } else {
    file.path(LOCAL_BASE, "input", "metadata_oral_qc_dedup.tsv")
  }
)
out_dir <- Sys.getenv(
  "FIG1_PERMANOVA_OUT_DIR",
  unset = if (on_server) SERVER_OUT_DIR else file.path(script_dir, "results", "Fig1")
)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

heatmap_vars <- c("age_cleaned_group", "continent_group", "site_cleaned",
                  "sequencing_region_cleaned")

OUT_PREFIX <- "Fig1_permanovaR2"
OUT_DATE   <- "20260909"
FILE_KEY <- c(age_cleaned_group = "age", continent_group = "continent",
              site_cleaned = "site", sequencing_region_cleaned = "seqregion")

log_file <- file.path(out_dir, paste0(OUT_PREFIX, "_log_", OUT_DATE, ".txt"))
write_log <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- paste0("[", timestamp, "] ", msg)
  cat(full_msg, "\n")
  cat(full_msg, "\n", file = log_file, append = TRUE)
}

write_log("=== Start: PERMANOVA R2 lower-triangle heatmaps (Unweighted UniFrac) ===")
dt_perm <- fread(input_permanova_unweight)
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
dt_perm[, Group1 := normalize_group_str(Group1)]
dt_perm[, Group2 := normalize_group_str(Group2)]
dt_perm <- dt_perm[Variable %in% heatmap_vars]

dup_idx <- duplicated(names(dt_meta))
if (any(dup_idx)) dt_meta <- dt_meta[, !dup_idx, with = FALSE]
if ("bioproject_id" %in% names(dt_meta) && !"bioproject" %in% names(dt_meta)) {
  setnames(dt_meta, "bioproject_id", "bioproject")
}
dt_meta[, age_cleaned_group := gsub("＞", ">", age_cleaned_group)]
dt_meta[, age_cleaned_group := gsub("＜", "<", age_cleaned_group)]
dt_meta[, bmi_cleaned_group := gsub("＜", "<", bmi_cleaned_group)]
dt_meta[, bmi_cleaned_group := gsub("≥", ">=", bmi_cleaned_group)]
dt_meta[, bmi_cleaned_group := gsub("≧", ">=", bmi_cleaned_group)]

dt_meta[, site_cleaned := tolower(site_cleaned)]

custom_orders <- list(
  age_cleaned_group = c("0-6 months", "6-36 months", "3-12 years", "12-18 years", "18-60 years", ">60 years"),
  continent_group = c("North America", "Asia", "Europe", "Africa", "Oceania", "South America"),
  site_cleaned = c("saliva", "tooth", "oral mucosa", "oral swab", "tongue", "gingiva", "other")
)
title_labels <- list(
  age_cleaned_group = "Age group",
  continent_group = "Continent",
  site_cleaned = "Site",
  sequencing_region_cleaned = "16S region"
)

display_label <- function(var, g) {
  if (var == "age_cleaned_group") {
    m <- c("0-6 months" = "0-6 mo", "6-36 months" = "6-36 mo",
           "3-12 years" = "3-12 yr", "12-18 years" = "12-18 yr",
           "18-60 years" = "18-60 yr", ">60 years" = ">60 yr")
    if (g %in% names(m)) return(unname(m[g]))
  }
  if (var == "continent_group") {
    m <- c("North America" = "NA", "Asia" = "AS", "Europe" = "EU",
           "Africa" = "AF", "Oceania" = "OC", "South America" = "SA")
    if (g %in% names(m)) return(unname(m[g]))
  }
  if (var == "site_cleaned") {
    m <- c("saliva" = "Saliva", "tooth" = "Tooth", "oral mucosa" = "Oral Mucosa",
           "oral swab" = "Oral Swab", "tongue" = "Tongue", "gingiva" = "Gingiva")
    if (g %in% names(m)) return(unname(m[g]))
  }
  g
}

min_n_by_var <- list(sequencing_region_cleaned = 500)
drop_levels_by_var <- list(site_cleaned = c("other"))

get_levels <- function(var) {
  vals <- as.character(dt_meta[[var]])
  vals <- vals[!is.na(vals) & vals != ""]
  counts <- sort(table(vals), decreasing = TRUE)
  if (var %in% names(custom_orders)) {
    lv <- c(intersect(custom_orders[[var]], names(counts)),
            setdiff(names(counts), custom_orders[[var]]))
  } else {
    lv <- sort(names(counts))
  }
  if (var %in% names(min_n_by_var)) {
    lv <- lv[counts[lv] >= min_n_by_var[[var]]]
  }
  if (var %in% names(drop_levels_by_var)) {
    lv <- setdiff(lv, drop_levels_by_var[[var]])
  }
  lv
}

levels_by_var <- sapply(heatmap_vars, get_levels, simplify = FALSE)
for (v in heatmap_vars) {
  write_log(paste0(v, ": ", length(levels_by_var[[v]]), " groups [",
                   paste(levels_by_var[[v]], collapse = ", "), "]"))
}

make_r2_matrix <- function(var, lv) {
  n <- length(lv)
  m <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(lv, lv))
  perm_sub <- dt_perm[Variable == var]
  key <- function(a, b) paste(sort(c(a, b)), collapse = " || ")
  perm_sub[, pair_key := mapply(function(a, b) key(a, b), Group1, Group2)]
  perm_keyed <- perm_sub[!duplicated(pair_key), .(pair_key, R2)]
  for (i in 1:n) {
    for (j in 1:n) {
      if (j >= i) next
      hit <- perm_keyed[ pair_key == key(lv[i], lv[j]) ]
      if (nrow(hit) > 0) m[i, j] <- as.numeric(hit$R2[1])
    }
  }
  m
}

r2_by_var <- sapply(heatmap_vars, function(v) make_r2_matrix(v, levels_by_var[[v]]), simplify = FALSE)

all_r2 <- unlist(lapply(r2_by_var, function(m) m[lower.tri(m)]))
fill_lim <- c(0, max(all_r2, na.rm = TRUE))
write_log(paste0("Shared fill scale: 0 to ", sprintf("%.4f", fill_lim[2])))

seq_low  <- "#cde2fb"
seq_high <- "#0d366b"
diag_fill <- "#f0efec"
ink <- "#1b1b1b"

SZ_TITLE       <- 24
SZ_SUBTITLE    <- 13
SZ_AXIS        <- 13
SZ_LEG_TITLE   <- 15
SZ_LEG_TEXT    <- 13
SZ_PANEL_TITLE <- 16
SZ_CELL        <- 3.4
SZ_DIAG        <- 2.8
STANDALONE_W   <- 6.4;  STANDALONE_H <- 6.0
COMBINED_W     <- 13.0; COMBINED_H   <- 4.5

custom_theme <- theme_bw() + theme(
  text = element_text(face = "bold", colour = ink),
  plot.title = element_text(size = SZ_TITLE, hjust = 0.5, face = "bold"),
  plot.subtitle = element_text(size = SZ_SUBTITLE, hjust = 0.5, colour = "#555555"),
  axis.title = element_blank(),
  axis.text.x = element_text(size = SZ_AXIS, angle = 45, hjust = 1, colour = ink),
  axis.text.y = element_text(size = SZ_AXIS, colour = ink),
  axis.ticks = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
  panel.grid = element_blank(),
  legend.title = element_text(size = SZ_LEG_TITLE),
  legend.text = element_text(size = SZ_LEG_TEXT),
  plot.margin = margin(10, 10, 6, 10)
)

n_panels <- 0
panel_list <- list()
for (v in heatmap_vars) {
  lv <- levels_by_var[[v]]
  r2m <- r2_by_var[[v]]
  n <- length(lv)
  lv_disp <- vapply(lv, function(g) display_label(v, g), character(1))
  if (n < 2) {
    write_log(paste("Skip", v, ": fewer than 2 groups."))
    next
  }
  if (anyNA(r2m[lower.tri(r2m)])) {
    write_log(paste("Warning:", v, "has missing PERMANOVA pairs; those cells stay blank."))
  }

  cells <- list()
  diag_cells <- list()
  for (i in 1:n) {
    for (j in 1:n) {
      if (j < i) {
        cells[[length(cells) + 1]] <- data.frame(x = j, y = n - i + 1, R2 = r2m[i, j])
      } else if (j == i) {
        diag_cells[[length(diag_cells) + 1]] <- data.frame(x = j, y = n - i + 1, label = lv[i])
      }
    }
  }
  cells <- do.call(rbind, cells)
  diag_cells <- do.call(rbind, diag_cells)
  cells$text_col <- ifelse(cells$R2 >= 0.5 * fill_lim[2], "white", ink)

  diag_cells$label_wrapped <- vapply(diag_cells$label, function(g) {
    g <- display_label(v, g)
    if (v == "site_cleaned" && grepl(" ", g)) return(gsub(" ", "\n", g))
    paste(strwrap(gsub(",\\s*", ", ", g), width = 10), collapse = "\n")
  }, character(1))

  p <- ggplot() +
    geom_tile(data = cells, aes(x = x, y = y, fill = R2),
              colour = "white", linewidth = 0.9) +
    geom_tile(data = diag_cells, aes(x = x, y = y),
              fill = diag_fill, colour = "white", linewidth = 0.9) +
    geom_text(data = cells, aes(x = x, y = y, label = sprintf("%.3f", R2),
                                colour = text_col),
              size = SZ_CELL, fontface = "bold", show.legend = FALSE) +
    geom_text(data = diag_cells, aes(x = x, y = y, label = label_wrapped),
              size = SZ_DIAG, fontface = "bold", colour = ink, lineheight = 0.9) +
    scale_colour_identity() +
    scale_fill_gradient(low = seq_low, high = seq_high, limits = fill_lim,
                        name = "R²") +
    scale_x_continuous(breaks = seq_len(n), labels = lv_disp,
                       expand = expansion(add = 0.15)) +
    scale_y_continuous(breaks = seq_len(n), labels = rev(lv_disp),
                       expand = expansion(add = 0.15)) +
    coord_equal() +
    custom_theme

  p_ind <- p +
    labs(
      title = paste0("PERMANOVA R²: ", title_labels[[v]]),
      subtitle = paste0(
        "Unweighted UniFrac | Lower-triangle pairwise R²,\n",
        "shared colour scale across the four panels\n",
        "All pairwise P = 0.001 (999 permutations)"
      )
    )

  panel_list[[length(panel_list) + 1]] <- p +
    labs(title = title_labels[[v]]) +
    theme(plot.title = element_text(size = SZ_PANEL_TITLE, hjust = 0.5,
                                    face = "bold", colour = ink))

  out_stem <- paste0(OUT_PREFIX, "_", FILE_KEY[[v]], "_", OUT_DATE)
  out_pdf <- file.path(out_dir, paste0(out_stem, ".pdf"))
  out_png <- file.path(out_dir, paste0(out_stem, ".png"))
  ggsave(out_pdf, plot = p_ind, width = STANDALONE_W, height = STANDALONE_H,
         device = cairo_pdf)
  ggsave(out_png, plot = p_ind, width = STANDALONE_W, height = STANDALONE_H,
         dpi = 300, device = "png", type = "cairo", bg = "white")
  n_panels <- n_panels + 1
  write_log(paste0("Saved heatmap for ", v, " (", n, " x ", n, ") -> ", basename(out_png)))
}

if (length(panel_list) > 0) {
  write_log("Patching 1x4 combined heatmap...")
  combined <- wrap_plots(panel_list, ncol = 4, nrow = 1, guides = "collect") &
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  comb_stem <- paste0(OUT_PREFIX, "_combined_", OUT_DATE)
  comb_pdf <- file.path(out_dir, paste0(comb_stem, ".pdf"))
  comb_png <- file.path(out_dir, paste0(comb_stem, ".png"))
  ggsave(comb_pdf, plot = combined, width = COMBINED_W, height = COMBINED_H,
         device = cairo_pdf)
  ggsave(comb_png, plot = combined, width = COMBINED_W, height = COMBINED_H,
         dpi = 300, device = "png", type = "cairo", bg = "white")
  write_log(paste0("Saved combined 1x4 heatmap (", basename(comb_png), ")."))
}

tidy_out <- rbindlist(lapply(heatmap_vars, function(v) {
  lv <- levels_by_var[[v]]; m <- r2_by_var[[v]]; n <- length(lv)
  rows <- list()
  for (i in 1:n) for (j in 1:n) {
    if (j < i) rows[[length(rows) + 1]] <-
      data.frame(Variable = v, Group1 = lv[i], Group2 = lv[j], R2 = m[i, j])
  }
  rbindlist(rows)
}))
fwrite(tidy_out, file.path(out_dir, paste0(OUT_PREFIX, "_pairs_", OUT_DATE, ".tsv")), sep = "\t")
write_log(paste0("Saved tidy R2 table with ", nrow(tidy_out), " pairs."))
write_log(paste0("=== Finished: ", n_panels, " heatmaps generated ==="))