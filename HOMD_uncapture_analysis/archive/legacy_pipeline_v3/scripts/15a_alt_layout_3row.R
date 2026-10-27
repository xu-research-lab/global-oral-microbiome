#!/usr/bin/env Rscript
# ALTERNATIVE page layout for Figure 1 -- three rows:
#
#   row 1   A  the data-set flow                       1/10 of the panel height
#   row 2   B (coverage)  |  C (species quadrant)      4/10
#   row 3   D  the multi-ASV profile                    5/10
#
# Row 2 is what fixes the size: B and C each take half of the figure width,
# minus the gutter, and their own aspects then decide how tall that row is.
# Rows 1 and 3 take their heights from it as 1/4 and 5/4 of row 2.
#
# Each panel is trimmed to its ink first and then scaled uniformly -- the panels
# keep the aspect they were drawn at, and only the space they are placed in is
# chosen here.  D is portrait (2148 x 2442), so at 5/10 of the height it cannot
# fill the width: it is centred, and the script reports how much of the row it
# actually covers.
#
# Run 13z (panel B), 13l (panel C), 12j (panel D) and -- with FLOW_ONLY set --
# 14_schematic_overview.R before this.  Then 99_combine_master_figs.R is not
# needed: this composites the three panels directly.

FIGDIR <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/analysis_report/figs"

PANELS <- c(A = file.path(FIGDIR, "fig_v3_flow_only.png"),
            B = file.path(FIGDIR, "fig_v3_panelB_coverage_v4.png"),
            C = file.path(FIGDIR, "fig_v3_homd_species_quadrant.png"),
            D = file.path(FIGDIR, "fig_v3_supported_multi_asv_profile.png"))
stopifnot(all(file.exists(PANELS)))

OUT_PNG <- file.path(FIGDIR, "fig_v3_master_alt_3row.png")
OUT_PDF <- file.path(FIGDIR, "fig_v3_master_alt_3row.pdf")

W_PX  <- 2700          # 9 in at 300 dpi: the figure's width before margins
GUT_X <- 80            # white between B and C
LET_H <- 52            # white band above rows 2 and 3, holding their letters
PAD_T <- 52            # band above row 1, holding A
PAD_L <- 52            # left margin, holding B and D
FUZZ  <- "1%"
LET_PT <- 46           # panel letters, px at 300 dpi (= 11 pt on the page)
LET_COL <- "#2B2B2B"

convert <- function(args) {
  if (system2("convert", args = args) != 0) stop("ImageMagick convert failed")
}
# no space in the -format string: system2 does not quote its arguments, so a
# space there would be split into a second argument
dims <- function(path) {
  as.integer(strsplit(system2("identify",
                              args = c("-format", "%w,%h", path),
                              stdout = TRUE), ",", fixed = TRUE)[[1]])
}

work <- tempfile("alt_layout_")
dir.create(work)

# trim to the ink and note the size every panel will be scaled from
ink <- vapply(names(PANELS), function(k) {
  out <- file.path(work, paste0(k, ".png"))
  convert(c(PANELS[[k]], "-fuzz", FUZZ, "-trim", "+repage", out))
  out
}, character(1))

size <- t(vapply(ink, dims, integer(2)))
rownames(size) <- names(ink)
colnames(size) <- c("w", "h")
message("[15] trimmed: ", paste(sprintf("%s %dx%d", rownames(size),
                                       size[, "w"], size[, "h"]), collapse = "  "))

# ---- geometry --------------------------------------------------------------
CELL_W <- (W_PX - GUT_X) %/% 2

# each of B and C is set to the cell width; the taller one sets the row
hB <- round(CELL_W * size["B", "h"] / size["B", "w"])
hC <- round(CELL_W * size["C", "h"] / size["C", "w"])
H2 <- max(hB, hC)
H1 <- round(H2 / 4)          # 1/10 of the total, since H2 is 4/10
H3 <- round(H2 * 5 / 4)      # 5/10

# fit a panel inside a box, keeping its aspect, then centre it on the row width
fit <- function(key, box_w, box_h, out) {
  w <- size[key, "w"]; h <- size[key, "h"]
  s <- min(box_w / w, box_h / h)
  convert(c(ink[[key]], "-resize", sprintf("%dx%d!", round(w * s), round(h * s)),
            "-background", "white", "-gravity", "center",
            "-extent", sprintf("%dx%d", box_w, box_h), out))
  c(w = round(w * s), h = round(h * s))
}
white_strip <- function(h, out) convert(c("-size", sprintf("%dx%d", W_PX, h),
                                          "xc:white", out))
# a panel in its cell: B flush left, C flush right, the gutter between them
cell <- function(key, h_cell, gravity, cell_w, out) {
  pad_w <- if (gravity == "West") CELL_W + GUT_X else CELL_W
  convert(c(ink[[key]], "-resize", sprintf("%dx", CELL_W),
            "-background", "white", "-gravity", "center",
            "-extent", sprintf("%dx%d", CELL_W, h_cell),
            "-background", "white", "-gravity", gravity,
            "-extent", sprintf("%dx%d", pad_w, h_cell), out))
}

row1 <- file.path(work, "row1.png")
a_dim <- fit("A", W_PX, H1, row1)

band1 <- file.path(work, "band1.png"); white_strip(LET_H, band1)

b_cell <- file.path(work, "b_cell.png")
c_cell <- file.path(work, "c_cell.png")
cell("B", H2, "West", CELL_W, b_cell)
cell("C", H2, "East", CELL_W, c_cell)
row2 <- file.path(work, "row2.png")
convert(c(b_cell, c_cell, "+append", row2))

band2 <- file.path(work, "band2.png"); white_strip(LET_H, band2)

row3 <- file.path(work, "row3.png")
d_dim <- fit("D", W_PX, H3, row3)

top <- file.path(work, "top.png"); white_strip(PAD_T, top)

body <- file.path(work, "body.png")
convert(c(top, row1, band1, row2, band2, row3, "-append", body))

h_body <- sum(dims(body)[2])
final  <- file.path(work, "final.png")
convert(c(body, "-background", "white", "-gravity", "West",
           "-extent", sprintf("%dx%d", W_PX + PAD_L, h_body), final))

# ---- panel letters ---------------------------------------------------------
# placed in the white bands, flush with the band's top and left, all four at the
# same size so they read as one grid
let_y <- c(A = 4,
           B = PAD_T + H1 + 4,
           C = PAD_T + H1 + 4,
           D = PAD_T + H1 + LET_H + H2 + 4)
let_x <- c(A = PAD_L + 4,
           B = PAD_L + 4,
           C = PAD_L + W_PX - CELL_W - GUT_X + 4,
           D = PAD_L + 4)

# -annotate under a northwest gravity takes x,y as the text box's corner, so the
# letter lands exactly in the band rather than sitting on an unseen baseline
# system2 runs its command through a shell, and an unquoted "#2B2B2B" is a shell
# comment: the rest of the line simply disappears
let_args <- c(final, "-font", "DejaVu-Sans-Bold", "-pointsize", LET_PT,
              "-fill", shQuote(LET_COL),
              as.vector(rbind("-gravity", "northwest",
                              "-annotate", sprintf("+%d+%d", let_x, let_y),
                              names(let_y))),
              OUT_PNG)
convert(let_args)
convert(c(OUT_PNG, OUT_PDF))

message(sprintf("[15] rows %d / %d / %d px (%.1f%% / %.1f%% / %.1f%% of the panel height)",
                H1, H2, H3, 100 * H1 / (H1 + H2 + H3),
                100 * H2 / (H1 + H2 + H3), 100 * H3 / (H1 + H2 + H3)))
message(sprintf("[15] A fills %d px of %d wide; D fills %d px of %d (%.0f%% of its row)",
                a_dim[["w"]], W_PX, d_dim[["w"]], W_PX,
                100 * d_dim[["w"]] / W_PX))
message("[15] alternative layout written: ", OUT_PNG, " | ", OUT_PDF)
