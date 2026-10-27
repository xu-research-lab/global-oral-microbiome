#!/usr/bin/env Rscript
# Final master composite using ImageMagick (keeps each finished figure intact):
#   left column = 13i phylum composite (top) + 13l species quadrant (bottom)
#   right      = 12j multi-member cluster profile (singletons are a supplementary figure)
#
# Each panel is trimmed to its ink before it is sized.  A panel rendered on its
# own device always carries the blank margin ggsave left around the plot, and
# that margin is scaled down together with the content, so a panel with generous
# margins ends up looking smaller than one without.  Trimming first makes the
# two columns comparable and lets each panel fill the space it is given.
#
# The geometry settled on here is written to LAYOUT_TSV: the compositor
# (14_schematic_overview.R) needs to know where the phylum and species panels
# meet in order to place the panel letters.
#
# The two columns carry different amounts of room under their plotting areas: the
# species panel's axis sits close to its own foot, while the cluster profile's
# 45-degree column labels take a third of its height.  Butted together at one
# height the two axes would end at visibly different levels, so the left column is
# given the bottom margin that brings its axis level with the other's.  The
# right-hand figure is resized to the column height, so its axis row is a fixed
# fraction of that height and the margin solves in one step rather than by
# iteration.

FIGDIR <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/analysis_report/figs"
TOP    <- file.path(FIGDIR, "fig_v3_panelB_coverage_v4.png")
BOT    <- file.path(FIGDIR, "fig_v3_homd_species_quadrant.png")
RIGHT  <- file.path(FIGDIR, "fig_v3_supported_multi_asv_profile.png")
OUT_PNG <- file.path(FIGDIR, "fig_v3_master_combined.png")
OUT_PDF <- file.path(FIGDIR, "fig_v3_master_combined.pdf")
LAYOUT_TSV <- file.path(FIGDIR, ".fig_v3_master_layout.tsv")

LEFT_W <- 1500
FUZZ   <- "1%"
GAP    <- 30                    # white between the two left-hand panels

t1 <- tempfile(fileext = ".png")
t2 <- tempfile(fileext = ".png")
t3 <- tempfile(fileext = ".png")
t4 <- tempfile(fileext = ".png")
t5 <- tempfile(fileext = ".png")
t6 <- tempfile(fileext = ".png")

convert <- function(args) {
  ok <- system2("convert", args = args)
  if (ok != 0) stop("ImageMagick convert failed")
}
dims <- function(path, fmt) {
  as.integer(system2("identify", args = c("-format", fmt, path), stdout = TRUE))
}

# trim to the ink, reset the canvas, then size to the requested width
trim_to <- function(inp, width, out) {
  convert(c(inp, "-fuzz", FUZZ, "-trim", "+repage", "-resize", width, out))
}

# The row inside a band of the image that carries the most ink.  A panel's axis
# line is the only thing that dark across that width, so this finds the axis
# whatever the panel happens to look like above it.  ImageMagick reduces the band
# to one column of pixels first: each value is the mean of its row, and the mean
# is lowest where the line runs.
axis_row <- function(path, x, y, w, h) {
  txt <- system2("convert",
                 c(path, "-crop", sprintf("%dx%d+%d+%d", w, h, x, y), "+repage",
                   "-colorspace", "Gray", "-filter", "box",
                   "-resize", sprintf("1x%d!", h), "-depth", "8", "txt:-"),
                 stdout = TRUE, stderr = FALSE)
  v <- as.integer(sub("^.*gray\\(([0-9]+).*$", "\\1", txt[-1]))
  stopifnot(length(v) == h)
  y + which.min(v) - 1L
}

trim_to(TOP, paste0(LEFT_W, "x"), t1)
trim_to(BOT, paste0(LEFT_W, "x"), t2)
# a blank strip between them: the upper panel ends in labels set at 45 degrees,
# so its trimmed edge is a diagonal, and butted straight against the panel below
# the two rows of text read as one crowded band
convert(c("-size", paste0(LEFT_W, "x", GAP), "xc:white", t5))
convert(c(t1, t5, t2, "-append", t3))

h  <- dims(t3, "%h")
h1 <- dims(t1, "%h") + GAP      # where the lower panel begins

# the right column, first at the left column's height: its axis row there gives
# the fraction of the height the axis sits at, which is the same at any height
trim_to(RIGHT, paste0("x", h), t4)
w4 <- dims(t4, "%w")

# the two axes: the species panel's, across the foot of the left column; the
# cluster profile's, across the foot of its boxplot column, which is the rightmost
# tenth of the figure
c_ax <- axis_row(t3, round(0.10 * LEFT_W), h - round(0.15 * h),
                 round(0.87 * LEFT_W), round(0.15 * h))
d_ax <- axis_row(t4, round(0.88 * w4), round(0.70 * h),
                 round(0.10 * w4), h - round(0.70 * h))
stopifnot(c_ax > h / 2, d_ax > h / 2)

h_new <- as.integer(round(c_ax / (d_ax / h)))
pad   <- h_new - h
if (pad > 0) {
  convert(c("-size", paste0(LEFT_W, "x", pad), "xc:white", t6))
  convert(c(t3, t6, "-append", t3))
  trim_to(RIGHT, paste0("x", h_new), t4)
} else {
  # the species panel's axis already sits below the other's: the left column
  # cannot be shortened without redrawing it, so leave the two as they are
  message("[99] axes already level (pad ", pad, "); left as is")
  h_new <- h
}

convert(c(t3, t4, "+append", OUT_PNG))
convert(c(OUT_PNG, OUT_PDF))

writeLines(c(sprintf("split\t%d",  h1),        # top of the species panel
             sprintf("left\t%d",   LEFT_W),
             sprintf("right\t%d",  w4),
             sprintf("height\t%d", h_new)),
           LAYOUT_TSV)

message(sprintf("[99] axes: species %d, profile %d of %d -> height %d (+%d)",
                c_ax, d_ax, h, h_new, pad))

message("[99] master composite written: ", paste(OUT_PNG, OUT_PDF, sep = " | "))
