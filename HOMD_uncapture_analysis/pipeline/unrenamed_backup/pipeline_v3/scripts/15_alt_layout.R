#!/usr/bin/env Rscript
# Step 19f: the master figure -- two columns:
#
#   left  (2/3 of the width)   A  the data-set flow, above
#                              B  coverage  |  C  per-genus identity curves
#   right (1/3 of the width)   D  the multi-ASV profile
#
# D's width is a third of the figure, and D is portrait, so that width fixes its
# height and with it the height of the whole page.  The left column is set to
# the same height; A fills the space above B and C.
#
# The three lower panels are aligned on their *borders*, not on their trimmed
# edges: the amount of ink below a panel's axis differs (B's 45-degree tick
# labels, C's axis title, D's rotated site labels and "Prevalence"), so trimming
# and butting the images together leaves the three frames at visibly different
# levels.  The frame row is found in each panel first and the panels are then
# placed so those rows coincide.
#
# The panel letters sit in the margin to the left of each panel rather than in
# a band of their own, so they cost no height.
#
# Every panel is trimmed to its ink and scaled uniformly: only the space each
# one is placed in is decided here, never its aspect.  The PNG is composited by
# ImageMagick; the PDF is composed from the panels' own cairo_pdf pages by
# LaTeX, so it stays vector.
#
# Run 19a (panel B), 19c (panel C), 19d (panel D) and 19e (panel A) before this.
# Panels are read from 19_figures/panels/ and the master figure is written to
# 19_figures/final/.

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline"
PANELDIR <- file.path(PIPE, "19_figures", "panels")
OUTDIR   <- file.path(PIPE, "19_figures", "final")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

PANELS <- c(A = file.path(PANELDIR, "fig_v3_flow_only.png"),
            B = file.path(PANELDIR, "fig_v3_panelB_coverage_v4.png"),
            C = file.path(PANELDIR, "fig_v3_homd_genus_identity_curves_numbered.png"),
            D = file.path(PANELDIR, "fig_v3_supported_multi_asv_profile.png"))
PANEL_PDF <- sub("\\.png$", ".pdf", PANELS)
stopifnot(all(file.exists(PANELS)), all(file.exists(PANEL_PDF)))

OUT_PNG <- file.path(OUTDIR, "fig_v3_master_alt.png")
OUT_PDF <- file.path(OUTDIR, "fig_v3_master_alt.pdf")

# D's width is a third of the content width, so the content width is what sets
# everything: a narrower page means a smaller D and, because D is portrait, a
# proportionally shorter page.  The height therefore cannot be bought from the
# layout; it is taken from the margins instead, which is why they are more
# generous than the letter alone needs.
# The page is the 2744-px layout scaled by 0.9: content, gutters, margins and the
# panel letters all shrink together, so nothing moves relative to anything else
# and the figure is simply 10% smaller.  Every panel is a raster scaled by the
# same factor, so their type shrinks by 10% with them.
W_PX    <- 2304        # content width, 7.68 in at 300 dpi
GUT_X   <- 54          # white between B and C
GUT_COL <- 54          # white between the left column and D
PAD_T   <- 36          # top edge of the page
PAD_B   <- 54          # bottom edge: D's "Prevalence" label otherwise ends flush
PAD_L   <- 40          # left margin: wide enough for one panel letter
PAD_R   <- 40          # right margin: D's last axis label otherwise ends on the page edge
PAGE_W  <- W_PX + PAD_L + PAD_R
FUZZ    <- "1%"
LET_PT  <- 41          # panel letters, px at 300 dpi (= 10 pt on the page)
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

# The row inside a band that carries the most ink.  A panel's axis line is the
# only thing that dark across that width, so this finds the frame whatever the
# panel looks like just above it.  ImageMagick reduces the band to one column of
# pixels first: each value is the mean of its row, and the mean is lowest where
# the line runs.  x, y, w, h are fractions of the image.
axis_row <- function(path, x, y, w, h) {
  d  <- dims(path)
  W  <- d[1]; H <- d[2]
  x0 <- round(x * W); y0 <- round(y * H)
  ww <- round(w * W); hh <- round(h * H)
  txt <- system2("convert",
                 c(path, "-crop", sprintf("%dx%d+%d+%d", ww, hh, x0, y0),
                   "+repage", "-colorspace", "Gray", "-filter", "box",
                   "-resize", sprintf("1x%d!", hh), "-depth", "8", "txt:-"),
                 stdout = TRUE, stderr = FALSE)
  v <- as.integer(sub("^.*gray\\(([0-9]+).*$", "\\1", txt[-1]))
  stopifnot(length(v) == hh)
  y0 + which.min(v) - 1L
}

work <- tempfile("alt_layout_")
dir.create(work)

# trim to the ink, then set each panel to the width it will be shown at
ink <- vapply(names(PANELS), function(k) {
  out <- file.path(work, paste0(k, "_ink.png"))
  convert(c(PANELS[[k]], "-fuzz", FUZZ, "-trim", "+repage", out))
  out
}, character(1))

# ---- width of each column ---------------------------------------------------
D_W    <- W_PX %/% 3
LEFT_W <- W_PX - D_W - GUT_COL
CELL_W <- (LEFT_W - GUT_X) %/% 2

sized <- vapply(names(ink), function(k) {
  out <- file.path(work, paste0(k, "_w.png"))
  convert(c(ink[[k]], "-resize", paste0(if (k == "D") D_W else CELL_W, "x"), out))
  out
}, character(1))
# A takes the whole left column, not one cell
convert(c(ink[["A"]], "-resize", paste0(LEFT_W, "x"),
           file.path(work, "A_w.png")))
sized[["A"]] <- file.path(work, "A_w.png")

size <- t(vapply(sized, dims, integer(2)))
rownames(size) <- names(sized); colnames(size) <- c("w", "h")
message("[19f] shown at: ", paste(sprintf("%s %dx%d", rownames(size),
                                        size[, "w"], size[, "h"]), collapse = "  "))

# ---- the frames -------------------------------------------------------------
# B's dark facet strips sit higher up (~0.55 of the panel), so the search starts
# below them; its frame is at ~0.84.  D's frame is its rightmost boxplot axis,
# which is the only line across that stretch.
bB <- axis_row(sized[["B"]], 0.05, 0.72, 0.90, 0.20)
# C's key strip sits above the curves and its frame is the lowest line on the
# page, just above the x tick labels
bC <- axis_row(sized[["C"]], 0.10, 0.80, 0.80, 0.16)
bD <- axis_row(sized[["D"]], 0.93, 0.80, 0.06, 0.20)
message(sprintf("[19f] frame rows: B %d/%d  C %d/%d  D %d/%d",
                bB, size["B", "h"], bC, size["C", "h"], bD, size["D", "h"]))

# The frame line is the same in all three lower panels, and A has to sit above
# B | C at the full column width, which fixes its height.  The frame line is
# pushed down as far as A needs rather than A being scaled: the page is what
# grows, and every panel keeps the size it was drawn at.  (D's height sets the
# page, and D gets shorter as its label column gets wider, so the room above
# B | C is not something the layout can count on.)
GAP_AB <- 30                                # white between A and the B | C row
y_base <- max(max(bB, bC, bD),
              size["A", "h"] + GAP_AB + max(bB, bC))   # the common frame line
top    <- c(B = y_base - bB, C = y_base - bC, D = y_base - bD)
bottom <- top + c(B = size["B", "h"], C = size["C", "h"], D = size["D", "h"])
H_CONTENT <- max(bottom)

# A sits above the B | C row, which starts where the higher of the two begins
row_top <- min(top[["B"]], top[["C"]])
H_A     <- row_top
stopifnot(H_A >= size["A", "h"] + GAP_AB)

# ---- the left column --------------------------------------------------------
left <- file.path(work, "left.png")
convert(c("-size", sprintf("%dx%d", LEFT_W, H_CONTENT), "xc:white",
          sized[["A"]], "-geometry", sprintf("+0+%d", (H_A - size["A", "h"]) %/% 2),
          "-composite",
          sized[["B"]], "-geometry", sprintf("+0+%d", top[["B"]]), "-composite",
          sized[["C"]], "-geometry",
          sprintf("+%d+%d", CELL_W + GUT_X, top[["C"]]), "-composite",
          left))

# ---- the right column -------------------------------------------------------
right <- file.path(work, "right.png")
convert(c("-size", sprintf("%dx%d", D_W, H_CONTENT), "xc:white",
          sized[["D"]], "-geometry", sprintf("+0+%d", top[["D"]]),
          "-composite", right))

# ---- the page ---------------------------------------------------------------
gutter <- file.path(work, "gutter.png")
convert(c("-size", sprintf("%dx%d", GUT_COL, H_CONTENT), "xc:white", gutter))
cols <- file.path(work, "cols.png")
convert(c(left, gutter, right, "+append", cols))

body <- file.path(work, "body.png")
convert(c("-size", sprintf("%dx%d", PAGE_W, H_CONTENT), "xc:white",
          cols, "-geometry", sprintf("+%d+0", PAD_L), "-composite", body))

top_band <- file.path(work, "top_band.png")
convert(c("-size", sprintf("%dx%d", PAGE_W, PAD_T), "xc:white", top_band))
bot_band <- file.path(work, "bot_band.png")
convert(c("-size", sprintf("%dx%d", PAGE_W, PAD_B), "xc:white", bot_band))
final <- file.path(work, "final.png")
convert(c(top_band, body, bot_band, "-append", final))

# ---- panel letters ---------------------------------------------------------
# one letter per panel, in the margin to its left and level with its top edge;
# C and D take the gutter to their left, which is the only white they have
let_y <- c(A = PAD_T + (H_A - size["A", "h"]) %/% 2,
           B = PAD_T + top[["B"]],
           C = PAD_T + top[["C"]],
           D = PAD_T + top[["D"]])
let_x <- c(A = 4,
           B = 4,
           C = PAD_L + CELL_W + (GUT_X - 30) %/% 2,
           D = PAD_L + LEFT_W + (GUT_COL - 30) %/% 2)

# -annotate under a northwest gravity takes x,y as the text box's corner, so the
# letter lands exactly where it is put rather than on an unseen baseline.
# system2 runs its command through a shell, and an unquoted "#2B2B2B" is a shell
# comment: the rest of the line simply disappears.
let_args <- c(final, "-font", "DejaVu-Sans-Bold", "-pointsize", LET_PT,
              "-fill", shQuote(LET_COL),
              as.vector(rbind("-gravity", "northwest",
                              "-annotate", sprintf("+%d+%d", let_x, let_y),
                              names(let_y))),
              # the pixels are at 300 dpi; without this the PDF page is the
              # pixel count in points instead of inches
              "-units", "PixelsPerInch", "-density", 300,
              OUT_PNG)
convert(let_args)

# ---- the PDF: composed from the vector panels ------------------------------
# The PNG is an ImageMagick composite, so wrapping it in a PDF would make a
# raster page that blurs when enlarged.  Each of the four panels is also written
# as a cairo_pdf; this places those pages at the same coordinates, with the same
# trim, in a one-page LaTeX document, so the PDF stays vector.  The letters are
# drawn by R into a transparent vector page and laid over the top.
#
# Pixels at 300 dpi and PostScript points are related by 72/300; every coordinate
# below is the one used for the PNG, converted with PT.
PT <- 0.24

# "WxH+X+Y" of the ink inside a panel page, origin top left
trimbox <- function(path) {
  s <- system2("convert", c(path, "-fuzz", FUZZ, "-format", "%@", "info:"),
               stdout = TRUE, stderr = FALSE)
  m <- regmatches(s, regexec("([0-9]+)x([0-9]+)\\+([0-9]+)\\+([0-9]+)", s))[[1]]
  if (length(m) != 5L) stop("could not read the trim box of ", path)
  as.integer(m[2:5])                       # w, h, x, y
}

PAGE_H <- PAD_T + H_CONTENT + PAD_B
pos <- rbind(A = c(0,                (H_A - size["A", "h"]) %/% 2),
             B = c(0,                 top[["B"]]),
             C = c(CELL_W + GUT_X,    top[["C"]]),
             D = c(LEFT_W + GUT_COL,  top[["D"]]))
pos[, 1] <- pos[, 1] + PAD_L
pos[, 2] <- pos[, 2] + PAD_T

fig <- vapply(names(PANELS), function(k) {
  pb <- dims(PANELS[[k]])
  tb <- trimbox(PANELS[[k]])
  dw <- size[k, "w"]; dh <- size[k, "h"]
  sprintf(paste0("\\put(%.3f,%.3f){\\includegraphics",
                 "[trim=%.3f %.3f %.3f %.3f,clip,width=%.3fbp,height=%.3fbp]{%s}}%%"),
          pos[k, 1] * PT, (PAGE_H - pos[k, 2] - dh) * PT,
          tb[3] * PT, (pb[2] - tb[4] - tb[2]) * PT,
          (pb[1] - tb[3] - tb[1]) * PT, tb[4] * PT,
          dw * PT, dh * PT, PANEL_PDF[[k]])
}, character(1))

letters_pdf <- file.path(work, "letters.pdf")
grDevices::cairo_pdf(letters_pdf, width = PAGE_W * PT / 72, height = PAGE_H * PT / 72,
                     bg = "transparent")
grid::grid.newpage()
grid::grid.text(names(let_y),
                x = grid::unit(let_x * PT, "pt"),
                y = grid::unit(1, "npc") - grid::unit(let_y * PT, "pt"),
                just = c("left", "top"),
                gp = grid::gpar(fontsize = LET_PT * PT, fontface = "bold",
                                col = LET_COL, fontfamily = "DejaVu Sans"))
invisible(grDevices::dev.off())

tex <- file.path(work, "master.tex")
writeLines(c(
  "\\documentclass{article}",
  sprintf("\\usepackage[paperwidth=%.3fbp,paperheight=%.3fbp,margin=0pt]{geometry}",
          PAGE_W * PT, PAGE_H * PT),
  "\\usepackage{graphicx}",
  "\\pagestyle{empty}",
  # a paragraph's indent is added before the picture, which pushed every panel
  # 15 bp right and ran D's last axis label off the page
  "\\setlength{\\parindent}{0pt}%",
  "\\begin{document}%",
  "\\setlength{\\unitlength}{1bp}%",
  sprintf("\\begin{picture}(%.3f,%.3f)%%", PAGE_W * PT, PAGE_H * PT),
  fig,
  sprintf("\\put(0,0){\\includegraphics[width=%.3fbp,height=%.3fbp]{%s}}%%",
          PAGE_W * PT, PAGE_H * PT, letters_pdf),
  "\\end{picture}%",
  "\\end{document}"), tex)

st <- system2("pdflatex", c("-interaction=nonstopmode", "-halt-on-error",
                            "-output-directory", work, tex),
               stdout = FALSE, stderr = FALSE)
built <- file.path(work, "master.pdf")
if (st == 0 && file.exists(built)) {
  stopifnot(file.copy(built, OUT_PDF, overwrite = TRUE))
  message("[19f] PDF composed from the vector panels")
} else {
  warning("[15] pdflatex failed; the PDF is the raster composite")
  convert(c(OUT_PNG, "-units", "PixelsPerInch", "-density", 300, OUT_PDF))
}

d <- dims(OUT_PNG)
message(sprintf("[19f] page %d x %d px (%.2f x %.2f in at 300 dpi)",
                d[1], d[2], d[1] / 300, d[2] / 300))
message(sprintf("[19f] D %d of %d wide (%.1f%%); A row %d px; A shown %dx%d",
                D_W, W_PX, 100 * D_W / W_PX, H_A,
                size["A", "w"], size["A", "h"]))
message("[19f] alternative layout written: ", OUT_PNG, " | ", OUT_PDF)
