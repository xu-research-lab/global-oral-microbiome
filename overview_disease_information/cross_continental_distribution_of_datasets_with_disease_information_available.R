suppressPackageStartupMessages({
  library(png)
  library(grid)
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
LOCAL_BASE  <- normalizePath(file.path(script_dir, "..", ".."))
on_server   <- dir.exists(file.path(SERVER_ROOT, "input"))
if (on_server) {
  BASE_DIR  <- SERVER_ROOT
} else {
  BASE_DIR  <- LOCAL_BASE
}

PANEL_DIR <- file.path(script_dir, "results", "Fig1")
BASE_DIR    <- Sys.getenv("FIG1_BASE_DIR",  unset = BASE_DIR)
PANEL_DIR   <- Sys.getenv("FIG1_PANEL_DIR", unset = PANEL_DIR)
AUTO_INSTALL <- as.logical(Sys.getenv("FIG1_AUTO_INSTALL", unset = "FALSE"))

PREP_DIR <- Sys.getenv("FIG1_PREP_DIR", unset = file.path(script_dir, "input_0928"))
top5_metadata_path <- Sys.getenv(
  "FIG1_TOP5_METADATA",
  unset = file.path(PREP_DIR, "metadata_oral_qc_matched_0928.tsv")
)

panel_png <- c(
  TOP5D  = "fig1_top5_disease_continent_0928.png"
)

save_plot_png_svg <- function(plot, png_path, width, height, bg = "white") {
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height,
                  device = "png", dpi = 600, bg = bg)
  invisible(png_path)
}

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


render_panel_disease_top5_continent <- function(out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(err <- load_libs(c("tidyverse", "patchwork"), AUTO_INSTALL))) return(err)
  if (!file.exists(top5_metadata_path))
    return(sprintf("matched metadata table not found: %s", top5_metadata_path))

  base_sz <- as.numeric(Sys.getenv("FIG1_TOP5_BASE_SZ", unset = "15"))
  cont_order <- c("North America", "South America", "Europe", "Asia",
                  "Africa", "Oceania")
  cont_abbr <- c(
    "North America" = "NA",
    "South America" = "SA",
    "Europe" = "EU",
    "Asia" = "AS",
    "Africa" = "AF",
    "Oceania" = "OC"
  )
  status_abbr_map <- c(
    "dental caries" = "DC",
    "periodontal disease" = "PD",
    "oral mucositis" = "OM",
    "oral squamous cell carcinoma" = "OSCC",
    "oscc" = "OSCC",
    "lung cancer" = "LC",
    "head and neck cancer" = "HNC",
    "colorectal cancer" = "CRC",
    "breast cancer" = "BC",
    "acute myeloid leukemia" = "AML",
    "preterm" = "PTB",
    "sjogren s syndrome" = "SS",
    "inflammatory bowel disease" = "IBD",
    "ibd" = "IBD",
    "aphthous" = "AS",
    "chagas disease" = "CD",
    "nonmalignant breast disease" = "NBD",
    "rheumatoid arthritis" = "RA",
    "covid 19" = "COVID-19",
    "hiv" = "HIV",
    "icu" = "ICU"
  )
  normalize_status <- function(x) {
    x <- tolower(trimws(x))
    x <- gsub("[[:punct:]]+", " ", x)
    gsub("[[:space:]]+", " ", trimws(x))
  }
  abbreviate_status <- function(x) {
    key <- normalize_status(x)
    out <- unname(status_abbr_map[key])
    ifelse(is.na(out), trimws(x), out)
  }

  md_raw <- readr::read_tsv(top5_metadata_path, show_col_types = FALSE)
  required_cols <- c("bioproject_id", "continent_group", "status_cleaned")
  missing_cols <- setdiff(required_cols, names(md_raw))
  if (length(missing_cols) > 0) {
    return(sprintf("top-5 metadata is missing required column(s): %s",
                   paste(missing_cols, collapse = ", ")))
  }
  cat(">=2-BioProject disease x continent metadata:", top5_metadata_path, "\n")

  md <- md_raw %>%
    transmute(
      bioproject_id = trimws(as.character(bioproject_id)),
      continent_group = trimws(as.character(continent_group)),
      status_cleaned = trimws(as.character(status_cleaned))
    ) %>%
    filter(
      !is.na(continent_group), continent_group %in% cont_order,
      !is.na(bioproject_id),
      !toupper(bioproject_id) %in% c("", "NA", "N/A"),
      !is.na(status_cleaned), status_cleaned != ""
    ) %>%
    mutate(
      continent_group = factor(continent_group, levels = cont_order),
      status_abbr = abbreviate_status(status_cleaned)
    )

  status_annotated_rows <- md

  disease_rows <- md %>%
    filter(!tolower(trimws(status_cleaned)) %in% c("hc", "other", ""))

  continent_totals <- status_annotated_rows %>%
    distinct(continent_group, bioproject_id) %>%
    count(continent_group, name = "n_total_projects") %>%
    tidyr::complete(
      continent_group = factor(cont_order, levels = cont_order),
      fill = list(n_total_projects = 0L)
    )
  cat("Continent BioProject totals:\n")
  print(as.data.frame(continent_totals))

  disease_counts <- disease_rows %>%
    group_by(continent_group, status_abbr) %>%
    summarise(
      n_bioprojects = dplyr::n_distinct(bioproject_id),
      n_samples     = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(continent_group, desc(n_bioprojects), desc(n_samples), status_abbr) %>%
    group_by(continent_group) %>%
    slice_head(n = 3) %>%
    ungroup() %>%
    left_join(continent_totals, by = "continent_group") %>%
    filter(n_bioprojects >= 2) %>%
    arrange(continent_group, desc(n_bioprojects), desc(n_samples), status_abbr)
  if (nrow(disease_counts) == 0)
    return("no disease data with >=2 BioProjects after filtering")

  missing_disease_continents <- setdiff(
    cont_order, unique(as.character(disease_counts$continent_group))
  )
  if (length(missing_disease_continents) > 0) {
    warning("No disease rows for continent(s): ",
            paste(missing_disease_continents, collapse = ", "))
  }

  hm_cols <- disease_counts %>%
    group_by(status_abbr) %>%
    summarise(total = sum(n_bioprojects), .groups = "drop") %>%
    arrange(desc(total), status_abbr) %>%
    pull(status_abbr)

  plot_continents <- continent_totals %>%
    filter(continent_group %in% unique(as.character(disease_counts$continent_group))) %>%
    arrange(desc(n_total_projects)) %>%
    pull(continent_group) %>%
    as.character()
  row_order <- plot_continents
  row_labs <- unname(cont_abbr[row_order])
  row_levels <- rev(row_labs)
  scale_y_shared <- ggplot2::scale_y_discrete(limits = row_levels,
                                              expand = expansion(add = 0.6))

  grid <- disease_counts %>%
    select(continent_group, status_abbr, n_bioprojects) %>%
    tidyr::complete(
      continent_group = factor(row_order, levels = row_order),
      status_abbr     = hm_cols,
      fill = list(n_bioprojects = 0L)
    ) %>%
    mutate(
      status_abbr = factor(status_abbr, levels = hm_cols),
      cont_lab    = factor(unname(cont_abbr[as.character(continent_group)]),
                           levels = row_levels)
    )

  max_n  <- max(grid$n_bioprojects)
  lim_hi <- max(12, max_n)
  size_range <- c(5, 15)
  size_breaks <- c(2, 5, 8, 12)
  col_expand <- 0.22

  bubble_txt <- grid %>%
    dplyr::filter(n_bioprojects > 0) %>%
    dplyr::mutate(txt_col = ifelse(n_bioprojects >= 5, "white", "#16334A"))

  p_matrix <- ggplot(grid, aes(x = status_abbr, y = cont_lab)) +
    geom_vline(xintercept = seq_along(hm_cols), linetype = "dashed",
               colour = "#E4EBF2", linewidth = 0.3) +
    geom_hline(yintercept = seq_along(row_levels), linetype = "dashed",
               colour = "#E4EBF2", linewidth = 0.3) +
    geom_point(data = dplyr::filter(grid, n_bioprojects > 0),
               aes(size = n_bioprojects, fill = n_bioprojects),
               shape = 21, colour = "white", stroke = 0.6) +
    geom_text(data = bubble_txt,
              aes(label = n_bioprojects, colour = txt_col),
              size = (base_sz * 0.52) / ggplot2::.pt, fontface = "bold",
              show.legend = FALSE) +
    scale_size_continuous(range = size_range, limits = c(2, lim_hi),
                          breaks = size_breaks) +
    scale_fill_gradient(low = "#8FB8D6", high = "#0A3052",
                        limits = c(2, lim_hi), breaks = size_breaks,
                        name = "BioProjects", guide = "none") +
    scale_colour_identity() +
    scale_x_discrete(limits = hm_cols, expand = expansion(add = col_expand)) +
    scale_y_shared +
    guides(size = "none") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = base_sz, base_family = "sans") +
    theme(
      axis.text.x       = element_text(size = base_sz * 0.82, face = "bold",
                                       color = "#303030", angle = 45,
                                       hjust = 1, vjust = 1,
                                       margin = margin(t = 2)),
      axis.text.y       = element_text(size = base_sz * 0.9, face = "bold",
                                       color = "#303030",
                                       margin = margin(r = 3)),
      panel.grid        = element_blank(),
      legend.position   = "none",
      plot.margin       = margin(2, 2, 4, 4),
      plot.background   = element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA)
    )

  top_df <- disease_counts %>%
    group_by(status_abbr) %>%
    summarise(global_bp = sum(n_bioprojects), .groups = "drop") %>%
    tidyr::complete(status_abbr = hm_cols, fill = list(global_bp = 0)) %>%
    mutate(status_abbr = factor(status_abbr, levels = hm_cols))

  top_max <- max(top_df$global_bp, na.rm = TRUE)
  legend_df <- data.frame(
    x = c(3.30, 3.85, 4.40, 4.95),
    y = top_max * 1.06,
    bp = size_breaks,
    point_size = scales::rescale(size_breaks, to = size_range,
                                 from = c(2, lim_hi))
  )

  p_top <- ggplot(top_df, aes(x = status_abbr, y = global_bp)) +
    geom_col(width = 0.55, fill = "#1F5F8B") +
    geom_text(aes(label = global_bp), vjust = -0.35,
              size = (base_sz * 0.6) / ggplot2::.pt, colour = "#303030") +
    geom_point(data = legend_df, aes(x = x, y = y), inherit.aes = FALSE,
               shape = 21, fill = "#5B7B9A", colour = "white",
               stroke = 0.5, size = legend_df$point_size) +
    geom_text(data = legend_df,
              aes(x = x, y = y - top_max * 0.16, label = bp),
              inherit.aes = FALSE,
              size = (base_sz * 0.5) / ggplot2::.pt, colour = "#303030") +
    annotate("text", x = 3.20, y = top_max * 1.20,
             label = "BioProjects", hjust = 0,
             size = (base_sz * 0.55) / ggplot2::.pt, fontface = "bold",
             colour = "#303030") +
    scale_x_discrete(limits = hm_cols, expand = expansion(add = col_expand)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.55))) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = base_sz, base_family = "sans") +
    theme(
      plot.margin      = margin(4, 2, 0, 4),
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )

  right_df <- data.frame(
    cont_lab = factor(row_labs, levels = row_levels),
    total    = continent_totals$n_total_projects[
      match(row_order, as.character(continent_totals$continent_group))]
  )
  p_right <- ggplot(right_df, aes(x = total, y = cont_lab)) +
    geom_col(width = 0.62, fill = "#5B7B9A") +
    geom_text(aes(label = total), hjust = -0.18,
              size = (base_sz * 0.6) / ggplot2::.pt, colour = "#303030") +
    scale_y_shared +
    scale_x_continuous(expand = expansion(mult = c(0, 0.4))) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = base_sz, base_family = "sans") +
    theme(
      axis.text.y       = element_blank(),
      axis.text.x       = element_text(size = base_sz * 0.82, angle = 45,
                                       hjust = 1, vjust = 1, colour = NA,
                                       margin = margin(t = 2)),
      panel.grid        = element_blank(),
      plot.margin       = margin(2, 4, 4, 0),
      plot.background   = element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA)
    )

  p <- p_top + patchwork::plot_spacer() + p_matrix + p_right +
    patchwork::plot_layout(design = "AB\nCD",
                           widths = c(1, 0.22), heights = c(0.32, 1))

  save_plot_png_svg(p, file.path(out_dir, panel_png["TOP5D"]),
                    width = 6.0, height = 5.2)
  ggplot2::ggsave(file.path(out_dir, "fig1_top5_disease_continent_0928.pdf"),
                  plot = p, width = 6.0, height = 5.2, device = "pdf")

  hm_bp_label <- setNames(
    paste0(row_labs, " (", continent_totals$n_total_projects[
      match(row_order, as.character(continent_totals$continent_group))], "bp)"),
    row_labs
  )
  hm_grid <- grid %>%
    dplyr::mutate(fill_val = ifelse(n_bioprojects >= 2, n_bioprojects, NA_real_))
  p_hm <- ggplot(hm_grid, aes(x = status_abbr, y = cont_lab)) +
    geom_tile(aes(fill = fill_val), colour = "#D9E2EC", linewidth = 0.5) +
    geom_text(data = dplyr::filter(hm_grid, n_bioprojects >= 2),
              aes(label = n_bioprojects), colour = "white",
              size = (base_sz * 0.85) / ggplot2::.pt, fontface = "bold") +
    scale_fill_gradientn(colours = c("white", "#8FB8D6", "#0A3052"),
                         values = c(0, 2 / lim_hi, 1),
                         na.value = "white", name = "BioProjects",
                         limits = c(0, lim_hi), breaks = c(0, 4, 8, 12)) +
    scale_x_discrete(limits = hm_cols, expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_discrete(limits = row_levels, labels = hm_bp_label,
                     expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = base_sz, base_family = "sans") +
    theme(
      axis.text.x     = element_text(size = base_sz * 0.85, face = "bold",
                                     color = "#303030", angle = 45,
                                     hjust = 1, vjust = 1, margin = margin(t = 4)),
      axis.text.y     = element_text(size = base_sz * 0.9, face = "bold",
                                     color = "#303030"),
      panel.grid      = element_blank(),
      legend.position = "right",
      legend.title    = element_text(size = base_sz * 0.75, face = "bold",
                                     color = "#303030"),
      legend.text     = element_text(size = base_sz * 0.7, color = "#303030"),
      plot.margin     = margin(6, 6, 6, 6)
    )
  save_plot_png_svg(p_hm, file.path(out_dir, "fig1_bp2_disease_continent_0928.png"),
                    width = 5.5, height = 4.2)

  write.table(disease_counts,
              file.path(out_dir, "fig1_top5_disease_continent_stats_0928.csv"),
              sep = ",", quote = TRUE, row.names = FALSE)
  invisible(NULL)
}

cat("Script dir :", script_dir, "\n")
cat("BASE_DIR   :", BASE_DIR, "\n")
cat("PANEL_DIR  :", PANEL_DIR, "\n")
if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)

cat("\n== rendering panel g - top-5 disease x continent ... ")
ans <- tryCatch(render_panel_disease_top5_continent(PANEL_DIR),
                error = function(e) paste("error:", conditionMessage(e)))
cat(if (is.null(ans)) "OK\n" else paste("SKIPPED:", ans, "\n", sep = " "))
cat("\nDone. Panel written to:", PANEL_DIR, "\n")
