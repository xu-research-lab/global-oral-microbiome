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
PROJECT_ROOT <- normalizePath(file.path(script_dir, ".."))
DATA_DIR     <- file.path(PROJECT_ROOT, "data", "metadata-and-biom-v2-qc")
BASE_DIR     <- PROJECT_ROOT

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

metadata_path <- Sys.getenv(
  "FIG1_METADATA",
  unset = file.path(DATA_DIR, "metadata_oral_qc.tsv")
)

panel_png <- c(
  DEPTH  = "seq_depth_histogram_0928.png"
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

render_panel_depth <- function(out_dir) {
  if (!is.null(err <- load_libs(c("tidyverse", "scales"), AUTO_INSTALL))) return(err)
  if (!file.exists(metadata_path))
    return(sprintf("metadata table not found: %s", metadata_path))
  base_sz <- fig1_base_sz("FIG1_DEPTH_BASE_SZ", 24)

  depth <- readr::read_tsv(metadata_path, show_col_types = FALSE) %>%
    filter(!is.na(total_reads), total_reads > 0)
  if (nrow(depth) == 0) return("no usable total_reads in metadata")
  med_reads <- median(depth$total_reads)

  probe <- ggplot(depth, aes(x = total_reads)) +
    geom_histogram(bins = 60) +
    scale_x_log10(limits = c(1e3, NA))
  y_top <- max(ggplot2::ggplot_build(probe)$data[[1]]$count, na.rm = TRUE)
  y_top <- ceiling(y_top / 1000) * 1000

  p <- ggplot(depth, aes(x = total_reads)) +
    geom_histogram(bins = 60, fill = "#3C5488", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = med_reads, linetype = "dashed", color = "#B75D69", linewidth = 1.0) +
    annotate("text", x = med_reads, y = Inf, vjust = 1.6, hjust = 0.5,
             label = paste0("median = ", scales::comma(round(med_reads))),
             color = "#B75D69", fontface = "bold", size = base_sz / 3.5) +
    scale_x_log10(limits = c(1e3, NA),
                  breaks = 10^(3:6),
                  labels = function(x) {
                    k <- round(log10(x))
                    sup <- c("⁰", "¹", "²", "³", "⁴",
                             "⁵", "⁶", "⁷", "⁸", "⁹")
                    paste0("10", sup[k + 1])
                  }) +
    scale_y_continuous(
      breaks = seq(0, y_top, by = 2000),
      limits = c(0, y_top),
      expand = expansion(mult = c(0, 0.02)),
      labels = function(x) ifelse(x == 0, "0", paste0(x / 1000, "K"))
    ) +
    labs(x = "Read counts", y = "Number of samples") +
    theme_classic(base_size = base_sz - 4, base_family = "sans") +
    theme(
      axis.text         = element_text(size = base_sz * 0.78, face = "bold", color = "#303030"),
      axis.title        = element_text(size = base_sz * 0.92, face = "bold", color = "#303030"),
      plot.background   = element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA)
    )

  save_plot_png_svg(p, file.path(out_dir, panel_png["DEPTH"]),
                    width = fig1_panel_in(5.2), height = fig1_panel_in(6.2))
  invisible(NULL)
}

cat("Script dir :", script_dir, "\n")
cat("BASE_DIR   :", BASE_DIR, "\n")
cat("PANEL_DIR  :", PANEL_DIR, "\n")
if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)

cat("\n== rendering panel c - sequencing-depth histogram ... ")
ans <- tryCatch(render_panel_depth(PANEL_DIR),
                error = function(e) paste("error:", conditionMessage(e)))
cat(if (is.null(ans)) "OK\n" else paste("SKIPPED:", ans, "\n", sep = " "))
cat("\nDone. Panel written to:", PANEL_DIR, "\n")
