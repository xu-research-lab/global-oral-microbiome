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
SERVER_METADATA <- file.path(
  SERVER_ROOT, "input", "metadata_oral_qc.tsv"
)

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

metadata_path <- Sys.getenv(
  "FIG1_METADATA",
  unset = if (on_server) {
    SERVER_METADATA
  } else {
    file.path(LOCAL_BASE, "input", "metadata_oral_qc.tsv")
  }
)

panel_png <- c(
  SHAN   = "shannon_diversity_0928.png"
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

render_panel_shannon <- function(out_dir) {
  if (!is.null(err <- load_libs(c("tidyverse", "scales"), AUTO_INSTALL))) return(err)
  if (!file.exists(metadata_path))
    return(sprintf("metadata table not found: %s", metadata_path))
  base_sz <- fig1_base_sz("FIG1_SHAN_BASE_SZ", 22)

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
    scale_y_continuous(labels = function(x) ifelse(x == 0, "0", paste0(x / 1000, "K"))) +
    labs(x = "Shannon diversity", y = "Number of samples") +
    theme_classic(base_size = base_sz - 4, base_family = "sans") +
    theme(
      axis.text         = element_text(size = base_sz * 0.72, face = "bold", color = "#303030"),
      axis.title        = element_text(size = base_sz * 0.82, face = "bold", color = "#303030"),
      plot.background   = element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA)
    )

  save_plot_png_svg(p, file.path(out_dir, panel_png["SHAN"]),
                    width = fig1_panel_in(7), height = fig1_panel_in(5.2))
  invisible(NULL)
}

cat("Script dir :", script_dir, "\n")
cat("BASE_DIR   :", BASE_DIR, "\n")
cat("PANEL_DIR  :", PANEL_DIR, "\n")
if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)

cat("\n== rendering panel f - alpha diversity (Shannon) ... ")
ans <- tryCatch(render_panel_shannon(PANEL_DIR),
                error = function(e) paste("error:", conditionMessage(e)))
cat(if (is.null(ans)) "OK\n" else paste("SKIPPED:", ans, "\n", sep = " "))
cat("\nDone. Panel written to:", PANEL_DIR, "\n")
