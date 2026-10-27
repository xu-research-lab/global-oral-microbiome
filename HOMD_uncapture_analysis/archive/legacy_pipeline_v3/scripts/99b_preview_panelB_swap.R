#!/usr/bin/env Rscript
# Preview only: the master composite rebuilt with the second take on panel B
# (13y) in place of the first (13i).  Nothing here overwrites a finished figure --
# every output carries the _previewB suffix, and 99_combine_master_figs.R is left
# as it is.  The step is identical to that script except for TOP, so that the swap
# can be looked at before it is made.
FIGDIR <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/analysis_report/figs"
TOP    <- file.path(FIGDIR, "fig_v3_panelB_coverage_v2.png")
BOT    <- file.path(FIGDIR, "fig_v3_homd_species_quadrant.png")
RIGHT  <- file.path(FIGDIR, "fig_v3_supported_multi_asv_profile.png")
OUT_PNG <- file.path(FIGDIR, "fig_v3_master_combined_previewB.png")
LAYOUT_TSV <- file.path(FIGDIR, ".fig_v3_master_layout_previewB.tsv")

LEFT_W <- 1500
FUZZ   <- "1%"

t1 <- tempfile(fileext = ".png"); t2 <- tempfile(fileext = ".png")
t3 <- tempfile(fileext = ".png"); t4 <- tempfile(fileext = ".png")

convert <- function(args) if (system2("convert", args = args) != 0) stop("convert failed")
dims <- function(path, fmt) {
  as.integer(system2("identify", args = c("-format", fmt, path), stdout = TRUE))
}
trim_to <- function(inp, width, out) {
  convert(c(inp, "-fuzz", FUZZ, "-trim", "+repage", "-resize", width, out))
}

trim_to(TOP, paste0(LEFT_W, "x"), t1)
trim_to(BOT, paste0(LEFT_W, "x"), t2)
convert(c(t1, t2, "-append", t3))
h  <- dims(t3, "%h"); h1 <- dims(t1, "%h")
trim_to(RIGHT, paste0("x", h), t4)

convert(c(t3, t4, "+append", OUT_PNG))
writeLines(c(sprintf("split\t%d", h1), sprintf("left\t%d", LEFT_W),
             sprintf("right\t%d", dims(t4, "%w")), sprintf("height\t%d", h)),
           LAYOUT_TSV)
message("[99b] preview composite: ", OUT_PNG, "  (split=", h1, ", h=", h, ")")
