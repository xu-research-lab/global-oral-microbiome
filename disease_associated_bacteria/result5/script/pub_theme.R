# ============================================================================
# Shared Cell Press styling for publication panel exports.
# Sourced by each panel_*.R script under script/.
# ============================================================================

options(bitmapType = "cairo")
Sys.setenv(RGL_USE_NULL = "TRUE")

suppressPackageStartupMessages({
  library(ggplot2)
})

# Repository root: run scripts from it, or set DA_WORK_DIR.
WORK_DIR <- normalizePath(Sys.getenv("DA_WORK_DIR", "."), mustWork = TRUE)
COMP_DIR <- file.path(WORK_DIR, "method_comparison_three")
PUB_DIR <- file.path(WORK_DIR, "figure")
dir.create(PUB_DIR, showWarnings = FALSE, recursive = TRUE)

# Cell Press two-column research article widths (mm).
W_FULL_MM <- 174
W_HALF_MM <- 85
W_15_MM <- 114
MAX_H_MM <- 235

# Font sizes at final print size (pt). Cell Press recommends 6-8 pt.
BASE_PT <- 7
AXIS_TITLE_PT <- 8
STRIP_PT <- 7
LEGEND_PT <- 6.5
# PDF built-in Helvetica; replace with Arial in Illustrator before submission.
FONT <- "Helvetica"

mm_to_in <- function(mm) as.numeric(mm) / 25.4

theme_cell <- function(base_size = BASE_PT, base_family = FONT) {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text = element_text(family = base_family, size = base_size, colour = "black"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.title = element_text(size = AXIS_TITLE_PT, face = "plain"),
      axis.text = element_text(size = BASE_PT, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = "black",
                                     linewidth = 0.3),
      strip.text = element_text(size = STRIP_PT, face = "bold",
                               margin = margin(2, 2, 2, 2)),
      legend.title = element_text(size = LEGEND_PT, face = "bold"),
      legend.text = element_text(size = LEGEND_PT),
      legend.key.size = unit(3.2, "mm"),
      legend.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(2, 4, 2, 4),
      complete = TRUE
    )
}

#' Save a ggplot as editable Helvetica PDF + 600 dpi PNG preview.
#'
#' @param plot ggplot object
#' @param stem file stem without extension (e.g. "panel_B_overlap")
#' @param width_mm figure width in millimetres (default full Cell Press width)
#' @param height_mm figure height in millimetres (capped at MAX_H_MM)
save_panel <- function(plot, stem, width_mm = W_FULL_MM, height_mm) {
  if (missing(height_mm) || is.null(height_mm) || !is.finite(height_mm)) {
    stop("height_mm is required for ", stem)
  }
  if (height_mm > MAX_H_MM) {
    warning(sprintf(
      "%s: height %.1f mm exceeds Cell Press max %.0f mm; clamping",
      stem, height_mm, MAX_H_MM
    ))
    height_mm <- MAX_H_MM
  }
  w_in <- mm_to_in(width_mm)
  h_in <- mm_to_in(height_mm)
  pdf_path <- file.path(PUB_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(PUB_DIR, paste0(stem, ".png"))

  # grDevices::pdf keeps text editable (not outlined) and embeds Helvetica.
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = w_in,
    height = h_in,
    units = "in",
    device = grDevices::pdf,
    family = FONT,
    useDingbats = FALSE,
    bg = "white"
  )
  ggsave(
    filename = png_path,
    plot = plot,
    width = w_in,
    height = h_in,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  cat(sprintf(
    "Wrote %s  (%.0f x %.0f mm) and PNG preview\n",
    pdf_path, width_mm, height_mm
  ))
  invisible(list(pdf = pdf_path, png = png_path,
                width_mm = width_mm, height_mm = height_mm))
}

disease_short <- function(x) {
  map <- c(
    "Head and Neck Cancer" = "H&N Cancer",
    "Periodontal Disease" = "Periodontal",
    "Colorectal Cancer" = "CRC",
    "Gastroesophageal Reflux Disease" = "GERD",
    "Obstructive Sleep Apnea" = "OSA",
    "Oral Lichen Planus" = "OLP",
    "Esophageal Cancer" = "Esophageal Ca",
    "Breast Cancer" = "Breast Ca",
    "Lung Cancer" = "Lung Ca",
    "Thyroid Cancer" = "Thyroid Ca",
    "Dental Caries" = "Dental Caries",
    "Aphthous" = "RAS"
  )
  x <- as.character(x)
  ifelse(x %in% names(map), unname(map[x]), x)
}

# Two/three-letter cancer abbreviations for figures that opt in; accepts full
# or disease_short() names, everything else falls back to disease_short().
disease_abbrev <- function(x) {
  map <- c(
    "Head and Neck Cancer" = "HNC", "H&N Cancer" = "HNC",
    "Thyroid Cancer" = "TC", "Thyroid Ca" = "TC",
    "Esophageal Cancer" = "EC", "Esophageal Ca" = "EC",
    "Lung Cancer" = "LC", "Lung Ca" = "LC",
    "Breast Cancer" = "BC", "Breast Ca" = "BC"
  )
  x <- as.character(x)
  ifelse(x %in% names(map), unname(map[x]), disease_short(x))
}
