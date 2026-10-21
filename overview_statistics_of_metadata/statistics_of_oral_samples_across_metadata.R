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

FIG1_TEXT_SCALE  <- as.numeric(Sys.getenv("FIG1_TEXT_SCALE",  unset = "1.18"))
FIG1_PANEL_SCALE <- as.numeric(Sys.getenv("FIG1_PANEL_SCALE", unset = "0.88"))
fig1_base_sz <- function(env_var, default) {
  configured <- Sys.getenv(env_var, unset = "")
  base <- if (nzchar(configured)) as.numeric(configured) else default
  base * FIG1_TEXT_SCALE
}
fig1_panel_in <- function(x) x * FIG1_PANEL_SCALE

fig1_site_levels <- c(
  "Saliva", "Tooth", "Oral Mucosa", "Oral Swab",
  "Tongue", "Gingiva", "Others"
)
fig1_site_colors <- c(
  "Saliva"      = "#0072B2",
  "Tooth"       = "#D55E00",
  "Oral Mucosa" = "#009E73",
  "Oral Swab"   = "#CC79A7",
  "Tongue"      = "#E69F00",
  "Gingiva"     = "#6A3D9A",
  "Others"      = "#56B4E9"
)

PREP_DIR <- Sys.getenv("FIG1_PREP_DIR", unset = file.path(script_dir, "input_0928"))

panel_png <- c(
  COHORT = "cohort_composition_0928.png"
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


render_panel_cohort <- function(out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(err <- load_libs(c("tidyverse", "patchwork"), AUTO_INSTALL))) return(err)
  base_sz <- fig1_base_sz("FIG1_COHORT_BASE_SZ", 18)
  bar_col <- "#2F5D7E"
  txt_col <- "#303030"

  read_prep <- function(f) {
    p <- file.path(PREP_DIR, f)
    if (!file.exists(p)) return(NULL)
    readr::read_tsv(p, show_col_types = FALSE)
  }

  vbar <- function(d, title, bar_width = 0.76, x_angle = 36,
                   x_text_size = 0.66) {
    if (is.null(d) || nrow(d) == 0)
      return(ggplot() + labs(title = title) + theme_void(base_size = base_sz - 4))
    d <- d %>%
      arrange(desc(n)) %>%
      mutate(label = factor(label, levels = label))
    y_breaks <- pretty(c(0, max(d$n, na.rm = TRUE)), n = 4)
    y_breaks <- y_breaks[y_breaks >= 0]
    y_max <- max(y_breaks)
    y_label <- function(x) {
      ifelse(x == 0, "0",
             ifelse(x %% 1000 == 0, paste0(x / 1000, "K"),
                    as.character(x)))
    }
    ggplot(d, aes(x = label, y = n)) +
      geom_col(fill = bar_col, width = bar_width, color = NA) +
      scale_y_continuous(
        breaks = y_breaks,
        limits = c(0, y_max),
        labels = y_label,
        expand = expansion(mult = c(0, 0.02))) +
      labs(title = title, x = NULL, y = NULL) +
      theme_classic(base_size = base_sz - 4, base_family = "sans") +
      theme(
        plot.title   = element_text(size = base_sz * 1.05, face = "bold",
                                    color = txt_col, hjust = 0.5,
                                    margin = margin(b = 3)),
        axis.text.x  = element_text(size = base_sz * x_text_size, face = "bold",
                                    color = "#303030",
                                    angle = x_angle,
                                    hjust = 1,
                                    vjust = if (x_angle == 90) 0.5 else 1),
        axis.text.y  = element_text(size = base_sz * 0.72, face = "bold",
                                    color = "#303030"),
        axis.line.x  = element_line(color = "#B0B0B0", linewidth = 0.4),
        axis.line.y  = element_line(color = "#B0B0B0", linewidth = 0.4),
        axis.ticks.x = element_line(color = "#B0B0B0"),
        axis.ticks.y = element_line(color = "#B0B0B0"),
        plot.margin  = margin(2, 3, 2, 2)
      )
  }

  donut_pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
                 "#6A3D9A", "#56B4E9", "#8C564B", "#B79F00", "#4D4D4D")

  vdonut <- function(d, title, n_total = NULL, color_values = NULL) {
    if (is.null(d) || nrow(d) == 0)
      return(ggplot() + labs(title = title) + theme_void(base_size = base_sz - 4))
    r_in  <- 1.55
    r_out <- 3.15
    d <- d %>%
      filter(!is.na(label), n > 0) %>%
      arrange(desc(n)) %>%
      mutate(prop = n / sum(n),
             ymax = cumsum(prop),
             ymin = cumsum(prop) - prop,
             mid  = (ymax + ymin) / 2,
             pct  = 100 * prop,
             lab  = as.character(label))
    tot <- if (is.null(n_total)) sum(d$n) else n_total
    lv  <- d$lab
    col_map <- if (is.null(color_values)) {
      setNames(donut_pal[seq_len(nrow(d))], lv)
    } else {
      vals <- unname(color_values[lv])
      if (anyNA(vals)) {
        stop("Missing colour(s) for: ",
             paste(lv[is.na(vals)], collapse = ", "))
      }
      setNames(vals, lv)
    }
    leg_map <- setNames(d$lab, lv)
    d$lab_f <- factor(d$lab, levels = lv)

    ggplot(d, aes(ymax = ymax, ymin = ymin, xmax = r_out, xmin = r_in,
                  fill = lab_f)) +
      geom_rect(color = "white", linewidth = 0.6) +
      coord_polar(theta = "y", start = 0) +
      xlim(0, r_out) +
      scale_fill_manual(values = col_map, breaks = lv,
                        labels = unname(leg_map[lv]), name = NULL) +
      labs(title = title) +
      theme_void(base_size = base_sz - 4, base_family = "sans") +
      theme(
        plot.title       = element_text(size = base_sz * 0.92, face = "bold",
                                        color = txt_col, hjust = 0.5,
                                        margin = margin(b = 1)),
        legend.position  = "right",
        legend.text.position = "right",
        legend.text      = element_text(size = base_sz * 0.62, color = "#303030"),
        legend.key.width = unit(2.4, "mm"),
        legend.key.height = unit(2.4, "mm"),
        legend.key       = element_rect(fill = NA, colour = NA),
        legend.spacing.y = unit(0.15, "mm"),
        legend.margin    = margin(0, 0, 0, 0),
        plot.margin      = margin(0, 1, 0, 0)
      )
  }

  site_d <- read_prep("cohort_site.tsv")
  p_site    <- vdonut(site_d, "Site", color_values = fig1_site_colors)

  p_sex     <- vdonut(
    read_prep("cohort_sex.tsv") %>%
      mutate(label = ifelse(label == "F", "Female",
                     ifelse(label == "M", "Male", label))),
    "Sex")

  p_region  <- vdonut(read_prep("cohort_region.tsv"),   "16S region")
  p_plat    <- vdonut(read_prep("cohort_platform.tsv"), "Platform")
  country_abbr <- c(
    "United States" = "USA", "China" = "CHN", "Canada" = "CAN",
    "Netherlands" = "NLD", "Spain" = "ESP", "Sweden" = "SWE",
    "Australia" = "AUS", "Japan" = "JPN", "Ghana" = "GHA"
  )

  country_d <- read_prep("cohort_country.tsv")
  if (!is.null(country_d)) country_d <- country_d %>% slice_max(n, n = 7)
  p_country <- vdonut(
    country_d %>% mutate(label = unname(country_abbr[label])),
    "Country")

  pdens <- function(d, xcol, xbreaks, title) {
    if (is.null(d) || nrow(d) == 0)
      return(ggplot() + labs(title = title) + theme_void(base_size = base_sz - 4))
    ggplot(d, aes(x = .data[[xcol]])) +
      geom_density(aes(y = after_stat(count)), fill = bar_col, color = bar_col,
                   alpha = 0.55, linewidth = 0.6) +
      scale_x_continuous(breaks = xbreaks, expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number(scale = 1/1000, suffix = "K")(x)),
                         expand = expansion(mult = c(0, 0.16))) +
      labs(title = title, x = NULL, y = NULL) +
      theme_classic(base_size = base_sz - 4, base_family = "sans") +
      theme(
        plot.title   = element_text(size = base_sz * 1.05, face = "bold",
                                    color = txt_col, hjust = 0.5,
                                    margin = margin(b = 3)),
        axis.text.x  = element_text(size = base_sz * 0.8, face = "bold",
                                    color = "#303030", angle = 45,
                                    hjust = 1, vjust = 1),
        axis.text.y  = element_text(size = base_sz * 0.72, face = "bold",
                                    color = "#303030"),
        axis.line.x  = element_line(color = "#B0B0B0", linewidth = 0.4),
        axis.line.y  = element_line(color = "#B0B0B0", linewidth = 0.4),
        plot.margin  = margin(3, 4, 2, 2)
      )
  }

  agev <- read_prep("cohort_agevals.tsv")
  p_age <- pdens(agev, "age", seq(0, 100, 20), "Age")

  bmiv <- read_prep("cohort_bmivals.tsv")
  p_bmi <- pdens(bmiv, "bmi", seq(0, 50, 10), "BMI")

  pairs <- read_prep("cohort_pairs.tsv")
  if (!is.null(pairs) && nrow(pairs) > 0) {
    n_pairs <- pairs$n_pairs[1]
    n_total <- if ("n_total" %in% names(pairs)) pairs$n_total[1] else NA_real_
  } else {
    n_pairs <- NA
    n_total <- NA_real_
  }
  if (is.na(n_total)) n_total <- n_pairs
  pairs_d <- data.frame(label = c("Paired", "Not paired"),
                        n     = c(n_pairs, max(0, n_total - n_pairs)))
  p_pairs <- vdonut(pairs_d, "Oral-gut paired")

  combined <- patchwork::wrap_plots(
    list(p_age, p_bmi, p_site, p_sex,
         p_country, p_region, p_plat, p_pairs),
    ncol = 4, nrow = 2) +
    patchwork::plot_layout(ncol = 4, nrow = 2, widths = rep(1, 4),
                           heights = rep(1, 2)) &
    theme(plot.margin = margin(1, 1, 1, 1))

  save_plot_png_svg(combined, file.path(out_dir, panel_png["COHORT"]),
                    width = fig1_panel_in(14), height = fig1_panel_in(5.8))
  invisible(NULL)
}


cat("Script dir :", script_dir, "\n")
cat("BASE_DIR   :", BASE_DIR, "\n")
cat("PANEL_DIR  :", PANEL_DIR, "\n")
if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)

cat("\n== rendering panel e - cohort composition ... ")
ans <- tryCatch(render_panel_cohort(PANEL_DIR),
                error = function(e) paste("error:", conditionMessage(e)))
cat(if (is.null(ans)) "OK\n" else paste("SKIPPED:", ans, "\n", sep = " "))
cat("\nDone. Panel written to:", PANEL_DIR, "\n")
