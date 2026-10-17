#!/usr/bin/env Rscript
# ============================================================================
# Panel B Venn: all significant (disease, genus, direction) calls from the
# three DA tools (MaAsLin2 / ANCOM-BC2 / LinDA), pooled across diseases,
# not the >=2-tool consensus. Output: panel_b_venn.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggforce)
})

source(file.path(Sys.getenv("DA_WORK_DIR", "."), "script", "pub_theme.R"))

SUM_PATH <- file.path(COMP_DIR, "overlap_summary_3way.csv")
if (!file.exists(SUM_PATH)) stop("Missing: ", SUM_PATH)

ov <- read.csv(SUM_PATH, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(ov) == 20L)

check_row <- function(r) {
  ok_m2 <- r$n_maaslin2 == r$n_all3 + r$n_maaslin2_ancombc2 +
    r$n_maaslin2_linda + r$n_only_maaslin2
  ok_ac <- r$n_ancombc2 == r$n_all3 + r$n_maaslin2_ancombc2 +
    r$n_ancombc2_linda + r$n_only_ancombc2
  ok_ld <- r$n_linda == r$n_all3 + r$n_maaslin2_linda +
    r$n_ancombc2_linda + r$n_only_linda
  ok_m2 && ok_ac && ok_ld
}
if (!all(vapply(seq_len(nrow(ov)), function(i) check_row(ov[i, ]), logical(1)))) {
  stop("Region arithmetic does not recover tool totals for some diseases")
}

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

r_circ <- 1.05
cx <- c(m2 = 0.00, ac = 1.20, ld = 0.60)
cy <- c(m2 = 0.70, ac = 0.70, ld = -0.30)

circles <- data.frame(
  x0 = unname(cx), y0 = unname(cy), r = r_circ,
  tool = c("MaAsLin2", "ANCOM-BC2", "LinDA"),
  colour = c(COL_M2, COL_AC, COL_LD),
  stringsAsFactors = FALSE
)
lab_df <- data.frame(
  x = c(-0.55, 1.75, 0.60, 0.60, -0.05, 1.25, 0.60),
  y = c(1.15, 1.15, -0.95, 1.35, 0.05, 0.05, 0.45),
  label = c(
    as.character(totals$only_m2), as.character(totals$only_ac),
    as.character(totals$only_ld), as.character(totals$m2_ac),
    as.character(totals$m2_ld), as.character(totals$ac_ld),
    as.character(totals$all3)
  ),
  face = c("plain", "plain", "plain", "plain", "bold", "bold", "bold"),
  size_pt = c(7, 7, 7, 6.5, 8.5, 7, 9),
  colour = c("grey45", "grey45", "grey45", "grey60",
             COL_ALL3, COL_ALL3, COL_ALL3),
  stringsAsFactors = FALSE
)
pct_m2 <- round(100 * (totals$all3 + totals$m2_ac + totals$m2_ld) / totals$n_m2, 0)
pct_ac <- round(100 * (totals$all3 + totals$m2_ac + totals$ac_ld) / totals$n_ac, 0)
pct_ld <- round(100 * (totals$all3 + totals$m2_ld + totals$ac_ld) / totals$n_ld, 0)
name_df <- data.frame(
  x = c(-1.35, 2.55, 0.60),
  y = c(1.95, 1.95, -1.48),
  vjust = c(0.5, 0.5, 1),
  label = c(
    sprintf("MaAsLin2\nn = %s\n%d%% in consensus",
            format(totals$n_m2, big.mark = ","), pct_m2),
    sprintf("ANCOM-BC2\nn = %s\n%d%% in consensus",
            format(totals$n_ac, big.mark = ","), pct_ac),
    sprintf("LinDA\nn = %s\n%d%% in consensus",
            format(totals$n_ld, big.mark = ","), pct_ld)
  ),
  colour = c(COL_M2, COL_AC, COL_LD),
  stringsAsFactors = FALSE
)

p_venn <- ggplot() +
  geom_circle(
    data = circles, aes(x0 = x0, y0 = y0, r = r, colour = tool),
    fill = NA, linewidth = 1.1, inherit.aes = FALSE
  ) +
  geom_circle(
    data = circles, aes(x0 = x0, y0 = y0, r = r, fill = tool),
    colour = NA, alpha = 0.10, inherit.aes = FALSE
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
    data = lab_df, aes(x = x, y = y, label = label),
    colour = lab_df$colour, fontface = lab_df$face,
    size = lab_df$size_pt / ggplot2::.pt, family = FONT, inherit.aes = FALSE
  ) +
  geom_text(
    data = name_df, aes(x = x, y = y, label = label, colour = I(colour)),
    vjust = name_df$vjust, size = 6.5 / ggplot2::.pt, family = FONT,
    fontface = "bold", lineheight = 0.95, inherit.aes = FALSE
  ) +
  coord_fixed(xlim = c(-2.0, 3.2), ylim = c(-2.15, 2.25), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_cell() +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.margin = margin(2, 8, 0, 8)
  )

save_panel(p_venn, "panel_b_venn", width_mm = W_FULL_MM, height_mm = 79.48)

cat(sprintf(
  "Venn: diseases=%d | MaAsLin2=%d ANCOM-BC2=%d LinDA=%d | consensus=%d\n",
  nrow(ov), totals$n_m2, totals$n_ac, totals$n_ld, totals$consensus
))
