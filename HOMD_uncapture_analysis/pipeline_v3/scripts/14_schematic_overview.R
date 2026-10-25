#!/usr/bin/env Rscript
# pipeline_v3: conceptual schematic for the HOMD-comprehensiveness analysis.
#
# Four blocks, assembled with patchwork:
#   a  the data set, from BioProjects to the clean ASV set
#   b  the same ASVs counted three ways (ASVs / reads / subject detections)
#   c  what "represented by HOMD" means (97% identity + 99% query coverage)
#   d  the two complementary ways the question is asked, and where each is shown
#
# Every block is drawn on a fixed x/y canvas with theme_void, so the layout can
# be nudged by editing coordinates.  Requires ggplot2 + patchwork only.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

# ---- change these in one place if the figure is renumbered ----------------
FIG_A <- "Fig. 1B, 1C"
FIG_C <- "Fig. 1D"

# ---- palette: matches the quantitative figures ---------------------------
COL_GAP  <- "#B2182B"   # below the species boundary / not represented
COL_OK   <- "#1A9850"   # at or above the species boundary / represented
COL_LINE <- "#7F7F7F"
COL_BOX  <- "#F7F7F7"
COL_HEAD <- "#E9E9E9"
COL_TXT  <- "#2B2B2B"
COL_INK  <- "#1E293B"  # ink in the flow cards (a touch cooler than COL_TXT)
COL_SUB  <- "#475569"  # the state each card leaves the ASVs in
COL_FLOW <- "#94A3B8"  # the arrows that join the cards

SZ_TAG <- 3.4   # block letter
SZ_HD  <- 2.8   # block title
SZ_TX  <- 2.4   # body text
SZ_SM  <- 2.1   # small text

void <- function() {
  theme_void(base_size = 8) + theme(plot.margin = margin(0, 0, 0, 0))
}
# `expand = FALSE` gives exactly the limits asked for.  The default is to pad them
# by 5% at each end, which a block drawn to fill its canvas does not want and a
# block with hand-placed text can trip over.
canvas <- function(xlim, ylim, expand = TRUE) {
  list(coord_cartesian(xlim = xlim, ylim = ylim, clip = "off", expand = expand),
       void())
}
# block letter + title on one line; `sc` rescales the type when the same block
# is redrawn at a different physical width (see the master composite below).
# `title = NULL` prints the letter alone, for blocks whose title is carried by
# the figure legend instead.
heading <- function(letter, title, y, sc = 1, let_mm = SZ_TAG, title_x = 6.5) {
  out <- list(annotate("text", x = 0, y = y, label = letter, hjust = 0,
                       vjust = 0.5, size = let_mm * sc, fontface = "bold",
                       colour = COL_TXT))
  if (!is.null(title)) {
    out <- c(out, list(annotate("text", x = title_x, y = y, label = title,
                                hjust = 0, vjust = 0.5, size = SZ_HD * sc,
                                fontface = "bold", colour = COL_TXT)))
  }
  out
}
rect <- function(xmin, xmax, ymin, ymax, fill, colour = NA, lw = 0.3) {
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
           fill = fill, colour = colour, linewidth = lw)
}
htext <- function(x, y, label, size = SZ_TX, colour = COL_TXT, hjust = 0.5,
                  vjust = 0.5, face = "plain", lh = 1.0) {
  annotate("text", x = x, y = y, label = label, size = size, colour = colour,
           hjust = hjust, vjust = vjust, fontface = face, lineheight = lh)
}

# =========================================================================
# a | the data set
# =========================================================================
# The chain starts at the dereplicated catalogue and runs through the three
# screens that were applied to it.  Each screen is a card: it names the screen,
# the count it left behind, and the state the ASVs are in afterwards.  The upstream
# steps (BioProjects -> per-project DADA2 ASVs -> exact dereplication) are left to
# the figure legend rather than given cards of their own, and the last card is not
# a screen at all but the split of the clean set across the 97% species boundary.
#
# The block is sized from its own text rather than laid out on a canvas chosen in
# advance.  Lines are stated as the size they should read at on the page, the width
# each one needs is measured, each card is as wide as its widest line plus that
# line's margin, and the canvas is whatever the five cards, their gaps and the block
# margins come to.  Two things follow from that which a fixed canvas cannot give:
# the type is as large as the block can carry rather than as large as a guessed
# canvas happens to allow, and no card carries padding it has no text for.  Editing
# a count or renaming a screen re-lays the row instead of overflowing a card.
#
#
# Sizes below are in inches of the 7.09-in design width until `FLOW_K` converts
# them; `FLOW_K` is what the solve produces, not an input.

FLOW_W_TARGET <- 9.0                  # width the block reads at, inches

# The block spans the figure, but the panels below it are drawn at 9 in and then
# scaled down by the composite, so they read at a little under 7 pt.  The block
# cannot be shrunk the same way -- it is already laid out at its final width --
# so its type is set at that size directly and its margin constants brought in
# with it.  The vertical stack is pulled in further still: the block was two and
# a half card-heights deep, and most of what was in the cards was space.
BOX_GAP  <- 0.42                     # between cards
BORDER   <- 0.20                     # block margin, left and right
BOX_R    <- 0.06                     # corner radius
BOX_LW   <- 1.1                      # card outline, sketch points
BADGE_R  <- 0.09                     # number badge radius
HEAD_DX  <- 0.30                     # a card title starts clear of its badge
INSET_L  <- 0.35                     # ... so a title line carries this much margin
INSET_C  <- 0.22                     # a centred line carries this much instead
MAP_GAP  <- 0.15                     # between the two mapping counts
ARROW_LW <- 2.8                      # the arrows: heavier than a card outline
ARROW_HD <- 0.06                     # and their heads

# the vertical stack is set further down, once the canvas width is solved: its
# spacings are points on the finished block, so they have to be converted with
# FLOW_K rather than written as canvas inches.  The block is deliberately flat --
# the cards are read as a sequence of counts, not as the figure's subject -- and
# the letter band is only as deep as the letter.

# The size each kind of line reads at, in points, on the 9-in block.  The block
# is shown at 4.9 in in the composite, so these are ~1.8x what a reader finally
# sees; 1.4x the values used while the panels were still drawn at their design
# size, which brings the flow up to the same apparent size as B and C.
FIN <- c(head = 8.7, num = 11.5, sub = 8.0, badge = 8.3,
         mnum = 10.9, mlab = 8.0, note = 8.0)

# the five cards.  Each title is one line: at the size the block carries, the
# longest of them still fits the card the row can give it.
CARDS <- data.frame(
  head  = c("Input ASVs", "Chimera filtering", "Organelle filtering",
            "Contaminant filtering", "Mapping to eHOMD"),
  num   = c("2,844,174", "2,099,459", "2,083,093", "1,871,022", NA),
  sub   = c("non-redundant ASVs", "non-chimeric ASVs", "organelle-free ASVs",
            "contaminant-free ASVs", NA),
  # one palette for the whole block: the first card's blue.  A different hue per
  # step made the flow read as five unrelated stages when it is one chain, and
  # the hues competed with the red/blue the panels below use for meaning.
  face  = rep("#F0F7FF", 5),
  edge  = rep("#3B82F6", 5),
  stringsAsFactors = FALSE
)

# what the last card reports: every clean ASV falls on one side of the species
# boundary, and the two counts close the chain.  They are set side by side, so the
# pair costs the card width rather than its height.
MAP <- data.frame(
  num = c("965,488", "905,534"),
  lab = c("≥97% identity", "<97% identity"),
  # the covered count takes the block's blue; the gap count keeps the figure's
  # red, which is what panel c marks the <97% side with
  col = c("#3B82F6", COL_GAP),
  stringsAsFactors = FALSE
)

## (the contaminant criterion is given in the figure caption, not drawn here)

# ---- measure the type ------------------------------------------------------
# A scratch device gives each line's width in em, which does not depend on the size
# the line is finally set at -- so one measurement serves the solve below.
#
# It has to be the device the PNG goes out on, not the PDF one: cairo_pdf and ragg
# resolve `sans` to different faces, and the PNG's is around 8% the wider of the two
# ("Contaminant" is 6.11 em on one and 6.59 on the other).  Measuring the wider of
# the two is the safe direction -- the PDF then carries a little spare room rather
# than the PNG spilling over its card.  SLACK covers a third face on another machine.
em_width <- function(strings, bold) {
  f <- tempfile(fileext = ".png")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(f, width = 20, height = 10, units = "in", res = 72)
  } else {
    grDevices::cairo_pdf(f, width = 20, height = 10)
  }
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  vapply(seq_along(strings), function(i) {
    par(ps = 100, family = "sans", font = if (bold[i]) 2 else 1)
    strwidth(strings[i], units = "in") * 72 / 100
  }, numeric(1))
}
SLACK <- 1.03

LINES <- rbind(
  data.frame(card = 1:5, role = "head", s = CARDS$head),
  data.frame(card = 1:4, role = "num",  s = CARDS$num[1:4]),
  data.frame(card = 1:4, role = "sub",  s = CARDS$sub[1:4])
)
LINES$bold  <- LINES$role != "sub"
LINES$fin   <- c(head = "head", num = "num", sub = "sub")[LINES$role]
LINES$em    <- SLACK * em_width(LINES$s, LINES$bold)
# width of a line = alpha * FLOW_W + beta, both in canvas inches
LINES$alpha <- LINES$em * FIN[LINES$fin] / (72 * FLOW_W_TARGET)
LINES$beta  <- ifelse(LINES$fin == "head", INSET_L, INSET_C)

# the last card's two counts, side by side; whichever of a count and its label is
# the wider sets the column
MAP_COL <- max(SLACK * em_width(c(MAP$num[1], MAP$lab[1]), c(TRUE, FALSE)) *
                 FIN[c("mnum", "mlab")] / (72 * FLOW_W_TARGET))
LINES <- rbind(LINES, data.frame(card = 5, role = "map", s = NA_character_,
                                 bold = NA, fin = NA_character_, em = NA_real_,
                                 alpha = 2 * MAP_COL, beta = MAP_GAP + INSET_C))

# ---- solve the canvas ------------------------------------------------------
# total width - FLOW_W is decreasing in FLOW_W wherever the cards do not between
# them already need the whole canvas, so a bisection finds the fixed point.
by_card <- split(LINES, LINES$card)
card_w  <- function(W) vapply(by_card, function(d) max(d$alpha * W + d$beta),
                              numeric(1))
slack   <- function(W) sum(card_w(W)) + 4 * BOX_GAP + 2 * BORDER - W
lo <- 1; hi <- 200
for (i in seq_len(60)) {
  mid <- (lo + hi) / 2
  if (slack(mid) > 0) lo <- mid else hi <- mid
}
FLOW_W  <- (lo + hi) / 2
FLOW_K  <- FLOW_W_TARGET / FLOW_W    # a canvas inch reads as this many inches
PT      <- FIN / FLOW_K              # so this is the size to set a line at
BOX_W   <- card_w(FLOW_W)
BOX_X0  <- BORDER + c(0, cumsum(head(BOX_W, -1) + BOX_GAP))
BOX_MID <- BOX_X0 + BOX_W / 2

# ---- the vertical stack ----------------------------------------------------
# The block is always shown at FLOW_W_TARGET wide, so a spacing that should read
# as 13 pt on the page is 13/72 device inches, or 13/72/FLOW_K canvas inches.  Set
# as canvas inches instead, the stack would have tightened whenever the solve
# landed on a wider canvas -- which is what happened when the titles became one
# line and the canvas went from 9.5 to 11.7.
BOX_Y0   <-  9.8 / 72 / FLOW_K      # bottom of a card
BOX_T    <- 13.0 / 72 / FLOW_K      # badge centre, and the title line
BOX_NUM  <- 14.4 / 72 / FLOW_K      # count, below the title
BOX_SUB  <- 12.9 / 72 / FLOW_K      # state label, below the count
BOX_BOT  <-  9.1 / 72 / FLOW_K      # from the state label to the card foot
BOX_H    <- BOX_T + BOX_NUM + BOX_SUB + BOX_BOT
LET_BAND <- 12.1 / 72 / FLOW_K      # above the cards, for the block letter
FLOW_H   <- BOX_Y0 + BOX_H + LET_BAND

# the block has to fit its canvas; if a future edit breaks that, fail here rather
# than in a figure nobody measures
stopifnot(sum(BOX_W) + 4 * BOX_GAP + 2 * BORDER <= FLOW_W + 1e-6,
          PT[["sub"]] >= 3)

message(sprintf("[14] data-set flow solved: canvas %.2f x %.2f, cards %s",
                FLOW_W, FLOW_H, paste(sprintf("%.2f", BOX_W), collapse = " ")))

ARROWS <- data.frame(
  x    = BOX_X0[-5] + BOX_W[-5] + 0.03,
  xend = BOX_X0[-1] - 0.03,
  y    = BOX_Y0 + BOX_H / 2
)

# `sc` = physical width of the block relative to its 7.09-in design width;
# `title = NULL` drops the block title (the composite relies on the Figure 1
# legend for it).
build_flow <- function(sc = 1, letter = "a", title = "the data set",
                       let_mm = SZ_TAG) {
  k   <- (7.09 * sc) / FLOW_W           # canvas inches -> device inches
  sz  <- function(pt) pt * k / 2.845    # sketch points -> ggplot's mm `text size`
  umm <- function(inch) inch * k * 25.4 # sketch inches -> ggplot's mm `mark size`
  top <- BOX_Y0 + BOX_H                 # the top edge of a card
  y_t1  <- top - BOX_T
  y_num <- y_t1 - BOX_NUM
  y_sub <- y_num - BOX_SUB

  BADGE <- data.frame(x = BOX_X0 + 0.15, y = y_t1)
  TITLE <- data.frame(x = BOX_X0 + HEAD_DX, y1 = y_t1, l1 = CARDS$head)
  BODY  <- data.frame(x = BOX_MID, y = y_num, num = CARDS$num,
                      ys = y_sub, sub = CARDS$sub)
  BODY  <- BODY[!is.na(BODY$num), ]

  # the last card's two counts, side by side about its centre
  MAPS <- data.frame(
    x   = rep(BOX_MID[5], 2) + c(-1, 1) * (MAP_GAP / 2 + MAP_COL * FLOW_W / 2),
    y   = y_num, num = MAP$num, ys = y_sub, lab = MAP$lab, col = MAP$col,
    stringsAsFactors = FALSE)

  ggplot() +
    # the cards.  roundrectGrob rather than geom_rect: the corners are given in
    # absolute units, so they stay circular however the panel is stretched.
    lapply(seq_len(nrow(CARDS)), function(i) {
      annotation_custom(
        grid::roundrectGrob(
          r   = unit(BOX_R * k, "in"),
          gp  = grid::gpar(fill = CARDS$face[i], col = CARDS$edge[i],
                           lwd = BOX_LW * k)),
        xmin = BOX_X0[i], xmax = BOX_X0[i] + BOX_W[i],
        ymin = BOX_Y0, ymax = top)
    }) +
    # the number badge: a disc in the card's own edge colour, the number reversed out
    # of it in white.  The outline is "transparent", not NA: GeomPoint$handle_na treats
    # `colour` as a must-not-be-missing aesthetic, so a constant colour of NA makes
    # every row missing and the layer draws nothing -- at render time only, which
    # building the plot does not reveal.
    geom_point(data = BADGE, aes(x = x, y = y, fill = CARDS$edge),
               shape = 21, size = umm(2 * BADGE_R), colour = "transparent",
               stroke = 0) +
    scale_fill_identity() +
    geom_text(data = BADGE, aes(x = x, y = y, label = seq_len(nrow(BADGE))),
              colour = "white", fontface = "bold", size = sz(PT[["badge"]])) +
    geom_text(data = TITLE, aes(x = x, y = y1, label = l1),
              hjust = 0, fontface = "bold", size = sz(PT[["head"]]), colour = COL_INK) +
    geom_text(data = BODY, aes(x = x, y = y, label = num),
              fontface = "bold", size = sz(PT[["num"]]), colour = COL_INK) +
    geom_text(data = BODY, aes(x = x, y = ys, label = sub),
              size = sz(PT[["sub"]]), colour = COL_SUB) +
    geom_text(data = MAPS, aes(x = x, y = y, label = num), colour = MAPS$col,
              fontface = "bold", size = sz(PT[["mnum"]])) +
    geom_text(data = MAPS, aes(x = x, y = ys, label = lab), colour = MAPS$col,
              size = sz(PT[["mlab"]])) +
    # the flow, arrow by arrow between the cards
    geom_segment(data = ARROWS, aes(x = x, xend = xend, y = y, yend = y),
                 arrow = arrow(length = unit(ARROW_HD * k, "in"), type = "closed"),
                 linewidth = sz(ARROW_LW), colour = COL_FLOW, lineend = "butt") +
    heading(letter, title, top + LET_BAND - 0.10, sc, let_mm, title_x = 0.30) +
    # exactly the solved canvas: no scale expansion, so the row of cards fills the
    # width it is given instead of losing 5% at each end
    canvas(c(0, FLOW_W), c(0, FLOW_H), expand = FALSE)
}

p_a <- build_flow()

# ---- the flow on its own ---------------------------------------------------
# An alternative page layout gives the flow its own row above the panels.  Set
# FLOW_ONLY=<path without extension> to write that block alone and stop.  It is
# drawn at the composite's width and without its own letter: the compositor
# places A with the other three letters, in the same band.
if (nzchar(Sys.getenv("FLOW_ONLY"))) {
  fl    <- Sys.getenv("FLOW_ONLY")
  fl_w  <- 9.0                  # the master figure's width
  pfl <- build_flow(fl_w / 7.09, letter = "", title = NULL)
  hfl <- fl_w * (FLOW_H / FLOW_W)
  ggsave(paste0(fl, ".png"), pfl, width = fl_w, height = hfl,
         dpi = 300, bg = "white")
  ggsave(paste0(fl, ".pdf"), pfl, width = fl_w, height = hfl,
         device = grDevices::cairo_pdf)
  message("[14] flow-only block written: ", fl, ".png | .pdf")
  quit(save = "no")
}

# =========================================================================
# b | the same ASVs, counted three ways
# =========================================================================
METRICS <- data.frame(
  label = c("ASVs", "Reads", "Subjects"),
  gap   = c(0.473, 0.049, 0.163),
  y     = c(58, 36, 14),
  stringsAsFactors = FALSE
)
METRICS$h    <- 15
METRICS$x0   <- 26
METRICS$x1   <- 96
METRICS$n    <- (METRICS$x1 - METRICS$x0) * METRICS$gap
METRICS$gapn <- sprintf("%.1f%%", 100 * METRICS$gap)
METRICS$okn  <- sprintf("%.1f%%", 100 * (1 - METRICS$gap))

p_b <- ggplot(METRICS) +
  geom_rect(aes(xmin = x0, xmax = x0 + n, ymin = y, ymax = y + h),
            fill = COL_GAP) +
  geom_rect(aes(xmin = x0 + n, xmax = x1, ymin = y, ymax = y + h),
            fill = COL_OK) +
  geom_text(aes(x = x0 - 2, y = y + h / 2, label = label),
            hjust = 1, size = SZ_TX, colour = COL_TXT) +
  geom_text(aes(x = x0 + n / 2, y = y + h + 4, label = gapn),
            size = SZ_TX, colour = COL_GAP, fontface = "bold") +
  geom_text(aes(x = x1 - 2, y = y + h / 2, label = okn),
            hjust = 1, size = SZ_TX, colour = "white", fontface = "bold") +
  heading("b", "the same ASVs, counted three ways", 95) +
  htext(0, 85, "red: best eHOMD hit below the 97% species boundary",
        size = SZ_SM, hjust = 0, colour = "grey35") +
  htext(0, 3, "each ASV counts once in the ASV total, but in proportion to its reads in the read total",
        size = SZ_SM, hjust = 0, vjust = 0) +
  canvas(c(0, 100), c(0, 100))

# =========================================================================
# c | what "represented by HOMD" means
# =========================================================================
REF_X0 <- 8; REF_X1 <- 62
REF_Y  <- 74; REF_H <- 9

cases <- data.frame(
  y    = c(56, 34, 12),
  xend = c(REF_X1, REF_X1, 34),           # where the query alignment stops
  col  = c(COL_OK, COL_GAP, COL_OK),
  note = c("99.4% identity\n100% of the query aligned",
           "91.2% identity\nbelow the 97% boundary",
           "99.0% identity, but only\npart of the query aligned"),
  verdict = c("represented", "gap \u2014 not represented", "gap \u2014 not represented"),
  vcol    = c(COL_OK, COL_GAP, COL_GAP),
  stringsAsFactors = FALSE
)
cases$ymax <- cases$y + REF_H

p_c <- ggplot() +
  # guides at the two alignment limits, so the shared origin is visible
  annotate("segment", x = REF_X0, xend = REF_X0, y = 4, yend = REF_Y + REF_H,
           linetype = "dotted", colour = "grey60", linewidth = 0.25) +
  annotate("segment", x = REF_X1, xend = REF_X1, y = 4, yend = REF_Y + REF_H,
           linetype = "dotted", colour = "grey60", linewidth = 0.25) +
  rect(REF_X0, REF_X1, REF_Y, REF_Y + REF_H, "#D9D9D9", COL_LINE) +
  htext((REF_X0 + REF_X1) / 2, REF_Y + REF_H / 2, "eHOMD reference (16S rRNA)",
        size = SZ_SM) +
  geom_rect(data = cases, aes(xmin = REF_X0, xmax = REF_X1,
                              ymin = y, ymax = ymax),
            fill = "#EDEDED", colour = "grey85", linewidth = 0.25) +
  geom_rect(data = cases, aes(xmin = REF_X0, xmax = xend,
                              ymin = y, ymax = ymax, fill = col), colour = NA) +
  scale_fill_identity() +
  geom_text(data = cases, aes(x = REF_X1 + 5, y = y + REF_H / 2 + 1, label = note),
            hjust = 0, size = SZ_SM, colour = COL_TXT, lineheight = 0.95) +
  geom_text(data = cases, aes(x = REF_X1 + 5, y = y - 3, label = verdict),
            hjust = 0, size = SZ_SM, colour = cases$vcol, fontface = "bold") +
  heading("c", "what counts as represented by eHOMD", 96) +
  htext(REF_X0 - 2, REF_Y + REF_H / 2, "reference", size = SZ_SM, hjust = 1) +
  htext(REF_X0 - 2, 56, "query ASVs", size = SZ_SM, hjust = 1) +
  canvas(c(-16, 100), c(0, 100))

# =========================================================================
# d | two ways of asking the question
# =========================================================================
box_d <- function(x0, x1, header, question, bullets, footer) {
  list(
    rect(x0, x1, 6, 88, COL_BOX, COL_LINE),
    rect(x0, x1, 74, 88, COL_HEAD),
    htext(x0 + 2.5, 81, header, size = SZ_HD, hjust = 0, face = "bold"),
    htext(x0 + 2.5, 66, question, size = SZ_TX, hjust = 0, face = "italic"),
    htext(x0 + 2.5, 55, bullets, size = SZ_TX, hjust = 0, vjust = 1, lh = 1.25),
    htext(x0 + 2.5, 11, footer, size = SZ_TX, hjust = 0, face = "bold")
  )
}

p_d <- ggplot() +
  box_d(0, 48, "Direct comparison with eHOMD", "How large is the gap?",
        paste("best-hit identity of every ASV, scored against eHOMD",
              "97% identity = the conventional species boundary",
              "weighted three ways: ASVs, reads, subjects",
              sep = "\n"),
        paste0("\u2192  shown in ", FIG_A)) +
  box_d(52, 100, "Phylogenetic placement", "Which organisms are in the gap?",
        paste("ASVs in \u22652 BioProjects with \u226550 reads (n = 54,277)",
              "placed in the reference tree, clustered at 0.03 (\u224897%)",
              "clusters in >50 subjects, \u22652 BioProjects, \u22652 continents",
              "86.6% of the 486 recurrent clusters have no eHOMD species",
              sep = "\n"),
        paste0("\u2192  shown in ", FIG_C)) +
  htext(50, -3, "both read the same ~97% species boundary",
        size = SZ_SM, face = "italic") +
  heading("d", "two ways of asking the same question", 96) +
  canvas(c(0, 100), c(-8, 100))

# =========================================================================
# assemble
# =========================================================================
# The flow is drawn on a fixed 18 x 4.5 canvas, so at this width it has exactly one
# height.  The row is weighted to that height rather than to a round number, or the
# cards come out stretched wide and short; the figure is then as tall as its rows add
# up to.  The two rows below it keep the heights they were drawn at.
W_STAND <- 7.09
H_FLOW  <- W_STAND * (FLOW_H / FLOW_W)     # the canvas aspect, at this width
H_LOWER <- c(1.96, 2.03)

schematic <- p_a / (p_b | p_c) / p_d +
  plot_layout(heights = c(H_FLOW, H_LOWER), widths = c(1, 1))

fig <- file.path(FIGOUT, "fig_v3_schematic_overview")
dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
H_STAND <- H_FLOW + sum(H_LOWER)
ggsave(paste0(fig, ".pdf"), schematic, width = W_STAND, height = H_STAND, device = dev)
ggsave(paste0(fig, ".png"), schematic, width = W_STAND, height = H_STAND, dpi = 300)
message("[14] schematic written: ", fig, ".pdf | .png")

# =========================================================================
# composite: the data-set flow above the finished quantitative figure
# =========================================================================
# The flow is redrawn at the width of the master figure (rather than pasted in as
# a PNG) so that its text stays vector in the PDF; the master figure itself is
# already a raster composite, so it is embedded as an image at its own aspect.
# These two, and the output name below, can be pointed elsewhere so that a change
# to one panel can be looked at in place without overwriting the finished figure.
MASTER_PNG <- Sys.getenv("MASTER_PNG",
                         file.path(FIGOUT, "fig_v3_master_combined.png"))
LAYOUT_TSV <- Sys.getenv("MASTER_LAYOUT",
                         file.path(FIGOUT, ".fig_v3_master_layout.tsv"))
stopifnot(file.exists(MASTER_PNG), file.exists(LAYOUT_TSV))

W_TOTAL <- 9.0                  # final width, inches
sc      <- W_TOTAL / 7.09       # flow design width is 7.09 in

# The flow is drawn on a fixed canvas, so its height is set by that canvas's aspect
# rather than chosen: anything else stretches the cards.
h_flow <- W_TOTAL * (FLOW_H / FLOW_W)

# 99_combine_master_figs.R trims each panel to its ink and stitches them; the
# geometry it settled on comes back here, because where the phylum and species
# panels meet decides where the letter C goes.
lay <- read.delim(LAYOUT_TSV, header = FALSE, col.names = c("k", "v"),
                  stringsAsFactors = FALSE)
lay <- setNames(as.integer(lay$v), lay$k)
stopifnot(all(c("split", "left") %in% names(lay)))
LEFT_W  <- lay[["left"]]
SPLIT_Y <- lay[["split"]]

# The only white the composite carries is a band along the top and left of the
# stitched panels, and it exists solely to hold the panel letters: it is sized
# from the letters themselves, so it costs no more than they do.
SZ_LET <- 11                    # panel-letter size, pt at the final figure size
PAD_T <- 56; PAD_L <- 52

img <- png::readPNG(MASTER_PNG)
h <- dim(img)[1]; w <- dim(img)[2]
padded <- array(1, dim = c(h + PAD_T, w + PAD_L, dim(img)[3]))
padded[(PAD_T + 1):(PAD_T + h), (PAD_L + 1):(PAD_L + w), ] <- img
HP <- dim(padded)[1]; WP <- dim(padded)[2]
h_master <- W_TOTAL / (WP / HP)

# Panel letters, given in padded-image pixels and converted to panel npc.  Each
# one is flush with the top and left of the band, so B/C/D share an x and B/D
# share a y and the three read as one grid.
LETTERS <- data.frame(
  lab = c("B", "C", "D"),
  px  = c(4, 4, PAD_L + LEFT_W + 4),
  py  = c(4, PAD_T + SPLIT_Y + 4, 4)
)
LETTERS$npcx <- LETTERS$px / WP
LETTERS$npcy <- 1 - LETTERS$py / HP

p_img <- ggplot() +
  annotation_custom(grid::rasterGrob(padded, interpolate = TRUE),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  lapply(seq_len(nrow(LETTERS)), function(i) {
    annotation_custom(grid::textGrob(
      LETTERS$lab[i], x = unit(LETTERS$npcx[i], "npc"),
      y = unit(LETTERS$npcy[i], "npc"), hjust = 0, vjust = 1,
      gp = grid::gpar(fontsize = SZ_LET, fontface = "bold",
                                    col = COL_TXT)))
  }) +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

# The flow block's letter goes through ggplot's mm size scale and is then
# multiplied by sc, so it is handed SZ_LET converted to mm and divided back out.
combo <- build_flow(sc, letter = "A", title = NULL,
                    let_mm = SZ_LET / (2.845 * sc)) / p_img +
  plot_layout(heights = c(h_flow, h_master))

out <- Sys.getenv("MASTER_OUT", file.path(FIGOUT, "fig_v3_master_with_flow"))
ggsave(paste0(out, ".pdf"), combo, width = W_TOTAL, height = h_flow + h_master,
       device = dev)
ggsave(paste0(out, ".png"), combo, width = W_TOTAL, height = h_flow + h_master,
       dpi = 300)
message("[14] flow + master composite written: ", out, ".pdf | .png")
