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
LOCAL_BASE  <- normalizePath(file.path(script_dir, ".."))
on_server   <- dir.exists(file.path(SERVER_ROOT, "input"))
if (on_server) {
  BASE_DIR  <- SERVER_ROOT
  PANEL_DIR <- file.path(SERVER_ROOT, "overview_figure_1", "results")
} else {
  BASE_DIR  <- LOCAL_BASE
  PANEL_DIR <- file.path(script_dir, "results")
}
BASE_DIR    <- Sys.getenv("FIG1_BASE_DIR",  unset = BASE_DIR)
PANEL_DIR   <- Sys.getenv("FIG1_PANEL_DIR", unset = PANEL_DIR)

panel_png <- c(
  A     = "pipeline_two_stage_funnel_HD.png",
  B     = "rarefaction_curve.png",
  C     = "global_sampling_world_map_final_cropped_260904.png",
  DEPTH = "seq_depth_histogram.png",
  PREV  = "prevalence_phyla_class_order.png",
  DENS  = "top5_phylum_density.png",
  SHAN  = "shannon_diversity.png",
  CONT  = "supp_phylum_continent_site_stacked.png",
  AGE   = "supp_phylum_age_site_stacked.png"
)

compo_base    <- "Figure1_combined_260907"
compo_w       <- 16L
compo_dpi     <- 300L
compo_mg      <- 0.14
compo_gap     <- 0.18
compo_tag_inset <- c(x = 0.12, y = 0.10)
compo_tag_size  <- 30
compo_tag_gap   <- 0.08
compo_tag_space <- 0.5

composite_rows <- list(
  c(A = "C"),
  c(B = "A",   C = "DEPTH"),
  c(D = "PREV"), 
  c(E = "DENS", F = "SHAN", G = "B")
)

# Composite assembly
combine_panels <- function(out_dir, png_map, rows, out_base, out_width, out_dpi,
                           mg = 0.14, gap = 0.18,
                           tag_inset = c(x = 0.12, y = 0.10), tag_size = 30,
                           tag_gap = 0.08, tag_space = 0.5) {
  letrs <- unlist(lapply(rows, names), use.names = FALSE)
  keys  <- unname(unlist(rows))
  if (length(letrs) != length(keys) || any(!nzchar(letrs)))
    stop("combine_panels: every row must be a NAMED vector of panel_png keys")

  paths <- file.path(out_dir, png_map[keys])
  names(paths) <- letrs
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    warning("Cannot build composite, missing panel file(s): ",
            paste(sprintf("%s (%s)", names(missing), basename(missing)),
                  collapse = ", "))
    return(invisible(FALSE))
  }

  imgs <- lapply(paths, png::readPNG)
  asp  <- vapply(imgs, function(im) ncol(im) / nrow(im), numeric(1))
  names(asp) <- letrs

  usable_w <- out_width - 2 * mg
  row_h <- vapply(rows, function(rw) {
    pn <- names(rw)
    (usable_w - gap * (length(pn) - 1)) / sum(asp[pn])
  }, numeric(1))
  if (!all(is.finite(row_h)) || any(row_h <= 0))
    stop("combine_panels: row heights are not finite - aspect ratios NA?")

  Htot <- mg + length(rows) * tag_space + sum(row_h)

  pos <- list()
  y0  <- mg

  for (i in rev(seq_along(rows))) {
    pn <- names(rows[[i]])
    h  <- row_h[i]
    x  <- mg
    for (p in pn) {
      w <- as.numeric(asp[p]) * h
      pos[[p]] <- c(x0 = x, x1 = x + w, y0 = y0, y1 = y0 + h)
      x <- x + w + gap
    }
    y0 <- y0 + h + tag_space
  }

  rect_ok <- function(pr) {
    !is.null(pr) && length(pr) >= 4 &&
      all(c("x0", "x1", "y0", "y1") %in% names(pr)) &&
      all(is.finite(pr[c("x0", "x1", "y0", "y1")])) &&
      pr[["x1"]] > pr[["x0"]] && pr[["y1"]] > pr[["y0"]]
  }

  draw_figure <- function() {
    grid::grid.newpage()
    tag_gp <- grid::gpar(fontface = "bold", fontsize = tag_size, col = "black")
    for (p in letrs) {
      pr <- pos[[p]]
      if (!rect_ok(pr)) {
        warning("combine_panels: skipping letter '", p, "' (invalid rect)")
        next
      }
      grid::pushViewport(grid::viewport(
        x = grid::unit((pr[["x0"]] + pr[["x1"]]) / 2, "in"),
        y = grid::unit((pr[["y0"]] + pr[["y1"]]) / 2, "in"),
        width  = grid::unit(pr[["x1"]] - pr[["x0"]], "in"),
        height = grid::unit(pr[["y1"]] - pr[["y0"]], "in"),
        just = "center"))
      grid::grid.raster(imgs[[p]], interpolate = TRUE)
      grid::popViewport()
      grid::grid.text(p,
        x = grid::unit(pr[["x0"]] + tag_inset["x"], "in"),
        y = grid::unit(pr[["y1"]] + tag_gap, "in"),
        just = c("left", "bottom"), gp = tag_gp)
    }
  }

  out_png <- file.path(out_dir, paste0(out_base, ".png"))

  grDevices::png(out_png, width = out_width, height = Htot, units = "in",
                 res = out_dpi, bg = "white")
  draw_figure()
  grDevices::dev.off()

  message("Saved: ", out_png)
  invisible(TRUE)
}

cat("Script dir :", script_dir, "\n")
cat("PANEL_DIR  :", PANEL_DIR, "\n")
cat("\n== Assembling composite figure ==\n")
cat("Composite layout:", length(composite_rows), "rows /",
    length(unlist(composite_rows)), "panels:",
    paste(names(unlist(composite_rows)), collapse = " "), "\n")
combine_panels(PANEL_DIR, panel_png, composite_rows, compo_base, compo_w,
               compo_dpi, mg = compo_mg, gap = compo_gap,
               tag_inset = compo_tag_inset, tag_size = compo_tag_size,
               tag_gap = compo_tag_gap, tag_space = compo_tag_space)
cat("Done.\n")
