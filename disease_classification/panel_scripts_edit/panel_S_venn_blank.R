#!/usr/bin/env Rscript
# ============================================================================
# Venn inset: pooled three-tool overlap (MaAsLin2 / ANCOM-BC2 / LinDA),
# sized to drop into the blank lower-right block of panel_B2_total
# (190 x 116 mm). The blank rectangle is ~160 mm wide x ~33 mm tall
# (rows HTN..Obesity, right of the last full-length bar). The S1 circle
# geometry is kept as TRUE circles (both axes share one mm-per-unit scale)
# and centred on the wide canvas; the set topology (7 regions) is unchanged.
# Tool captions run as two lines (name + n / consensus share) in the free
# corners and below the bottom circle.
# Canvas: 158 x 36 mm, fonts ~2x to match panel_B2_total / panel A abbrev.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggforce)
})

source("/beegfs/xuhuan/global_oral_tmp/panel_scripts_edit/pub_theme.R")

SUM_PATH <- file.path(WORK_DIR, "overlap_summary_3way.csv")
if (!file.exists(SUM_PATH)) stop("Missing: ", SUM_PATH)

ov <- read.csv(SUM_PATH, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(ov) == 20L)

totals <- list(
  all3    = sum(ov$n_all3),
  m2_ld   = sum(ov$n_maaslin2_linda),
  m2_ac   = sum(ov$n_maaslin2_ancombc2),
  ac_ld   = sum(ov$n_ancombc2_linda),
  only_m2 = sum(ov$n_only_maaslin2),
  only_ac = sum(ov$n_only_ancombc2),
  only_ld = sum(ov$n_only_linda),
  n_m2    = sum(ov$n_maaslin2),
  n_ac    = sum(ov$n_ancombc2),
  n_ld    = sum(ov$n_linda)
)
totals$consensus <- totals$all3 + totals$m2_ac + totals$m2_ld + totals$ac_ld

COL_M2 <- "#E6550D"
COL_AC <- "#31A354"
COL_LD <- "#756BB1"
COL_ALL3 <- "#08519C"

# Original S1 coordinate space: x in [-2, 3.2], y in [-2.15, 2.25].
# One single mm-per-unit scale for BOTH axes so the ellipses of the wide-canvas
# version become true circles again; the narrow circle cluster (3.3 units
# wide) is centred horizontally, and the top caption row keeps an 0.8 mm
# margin above the tallest circle.
W_MM <- 158
H_MM <- 36    # 33 mm slot height + 3 mm text margin so captions are not clipped
S <- 8.6                 # identical mm-per-unit scale in x and y
to_x <- function(x) (x - 0.60) * S + W_MM / 2
VY0 <- (H_MM - 0.8) - 1.75 * S  # circle top unit y = 0.70 + 1.05
to_y <- function(y) y * S + VY0

R <- 1.05
circ <- data.frame(
  x0 = to_x(c(0.00, 1.20, 0.60)),
  y0 = to_y(c(0.70, 0.70, -0.30)),
  a = R * S,
  b = R * S,
  tool = c("MaAsLin2", "ANCOM-BC2", "LinDA"),
  colour = c(COL_M2, COL_AC, COL_LD),
  stringsAsFactors = FALSE
)

# Region count labels; anchors re-centred for the un-stretched circle
# geometry (m2_ac sits mid-lens, all3 at the triple-intersection centre,
# m2_ld in the lower-left of its lens, ac_ld pushed right into its own moon).
lab_df <- data.frame(
  x = to_x(c(-0.55, 1.75, 0.60, 0.60, -0.05, 1.45, 0.60)),
  y = to_y(c(1.15, 1.15, -0.95, 1.15, -0.05, -0.02, 0.50)),
  label = c(
    as.character(totals$only_m2), as.character(totals$only_ac),
    as.character(totals$only_ld), as.character(totals$m2_ac),
    as.character(totals$m2_ld), as.character(totals$ac_ld),
    as.character(totals$all3)
  ),
  face = c("plain", "plain", "plain", "bold", "bold", "bold", "bold"),
  colour = c("grey45", "grey45", "grey45", "grey60",
             COL_ALL3, COL_ALL3, COL_ALL3),
  stringsAsFactors = FALSE
)
pct_m2 <- round(100 * (totals$all3 + totals$m2_ac + totals$m2_ld) / totals$n_m2, 0)
pct_ac <- round(100 * (totals$all3 + totals$m2_ac + totals$ac_ld) / totals$n_ac, 0)
pct_ld <- round(100 * (totals$all3 + totals$m2_ld + totals$ac_ld) / totals$n_ld, 0)
# Tool captions as two lines tucked into the free canvas margins
# (top-left, top-right and below the bottom circle):
#   line 1: tool name + sample size
#   line 2: share of that tool's calls that are 2+/3-tool consensus
name_df <- data.frame(
  x = c(1.8, W_MM - 1.8, to_x(0.60)),
  y = c(H_MM - 0.8, H_MM - 0.8, 4.3),
  hjust = c(0, 1, 0.5),
  vjust = c(1, 1, 0.5),
  label = c(
    sprintf("MaAsLin2 (n = %s)\n%d%% consensus",
            format(totals$n_m2, big.mark = ","), pct_m2),
    sprintf("ANCOM-BC2 (n = %s)\n%d%% consensus",
            format(totals$n_ac, big.mark = ","), pct_ac),
    sprintf("LinDA (n = %s)\n%d%% consensus",
            format(totals$n_ld, big.mark = ","), pct_ld)
  ),
  colour = c(COL_M2, COL_AC, COL_LD),
  stringsAsFactors = FALSE
)

p_venn <- ggplot() +
  geom_ellipse(
    data = circ, aes(x0 = x0, y0 = y0, a = a, b = b,
                     angle = 0, fill = tool),
    colour = NA, alpha = 0.10, inherit.aes = FALSE
  ) +
  geom_ellipse(
    data = circ, aes(x0 = x0, y0 = y0, a = a, b = b,
                     angle = 0, colour = tool),
    fill = NA, linewidth = 0.35, inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c("MaAsLin2" = COL_M2, "ANCOM-BC2" = COL_AC, "LinDA" = COL_LD),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c("MaAsLin2" = COL_M2, "ANCOM-BC2" = COL_AC, "LinDA" = COL_LD),
    guide = "none"
  ) +
  geom_text(
    data = lab_df,
    aes(x = x, y = y, label = label, colour = I(colour), fontface = I(face)),
    size = 10.5 / ggplot2::.pt, family = FONT,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = name_df,
    aes(x = x, y = y, label = label, colour = I(colour),
        hjust = I(hjust), vjust = I(vjust)),
    size = 10 / ggplot2::.pt, family = FONT, fontface = "bold",
    lineheight = 0.95, inherit.aes = FALSE
  ) +
  coord_fixed(xlim = c(0, W_MM), ylim = c(0, H_MM), expand = FALSE,
              clip = "off") +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

save_panel(p_venn, "panel_S_venn_blank", width_mm = W_MM, height_mm = H_MM)
cat(sprintf(
  "Venn inset: %g x %g mm | only_m2=%d only_ac=%d only_ld=%d | m2_ac=%d m2_ld=%d ac_ld=%d | all3=%d\n",
  W_MM, H_MM, totals$only_m2, totals$only_ac, totals$only_ld,
  totals$m2_ac, totals$m2_ld, totals$ac_ld, totals$all3
))