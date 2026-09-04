library(ggplot2)
library(ggforce)
library(grid)

# ---- Locate script directory & output folder (mirrors figure1_composite_all_in_one.R) ----
find_script_dir <- function() {
  if (interactive()) return(normalizePath(getwd()))
  args <- commandArgs(trailingOnly = FALSE)
  f <- args[grep("^--file=", args)]
  if (length(f) == 1) return(dirname(normalizePath(sub("^--file=", "", f))))
  normalizePath(getwd())
}
script_dir <- find_script_dir()
SERVER_ROOT <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-map-2/04-patch-pipeline-map/03-patch-pipeline-map"
on_server <- dir.exists(file.path(SERVER_ROOT, "input"))
PANEL_DIR <- if (on_server) file.path(SERVER_ROOT, "overview_figure_1", "results") else file.path(script_dir, "results")
PANEL_DIR <- Sys.getenv("FIG1_PANEL_DIR", unset = PANEL_DIR)
if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)

# Palette (Morandi / low saturation)
palette <- c(blue = "#E9F1F8", orange = "#F9EDE0", yellow = "#FBF4E2", green = "#EAF3E9")
zone_title_pal <- c(blue = "#2F5D7E", orange = "#A5672A", yellow = "#8A6D1F", green = "#3E6B4F")

node_fill   <- "white"
node_border <- "#9AA3AB"
arrow_col   <- "#7E8A94"
bullet_col  <- "#8A94A0"
txt_main    <- "#262626"
txt_sub     <- "#3F3F3F"
txt_detail  <- "#6B6B6B"
txt_excl    <- "#B75D69"

# Master font scale (same knobs as the composite script)
A_TXT     <- as.numeric(Sys.getenv("FIG1_A_TXT_SCALE",     unset = "1.40"))
A_TXT_SUB <- as.numeric(Sys.getenv("FIG1_A_TXT_SCALE_SUB", unset = "1.25"))

# Geometry helpers
rect_xy <- function(cx, cy, w, h, grp = NULL) {
  d <- data.frame(
    x = c(cx - w/2, cx + w/2, cx + w/2, cx - w/2),
    y = c(cy - h/2, cy - h/2, cy + h/2, cy + h/2)
  )
  if (!is.null(grp)) d$grp <- grp
  d
}
arrowhead_down <- function(cx, cy, a = 3.2, b = 1.9)
  data.frame(x = c(cx - b, cx + b, cx), y = c(cy, cy, cy - a))

# Canvas & zones (x domain widened so the enlarged text fits inside its boxes)
zones <- data.frame(
  name = c("1 . Sample retrieval from NCBI SRA", "2 . Metadata harmonization",
           "3 . Filtering and sequencing-data linkage", "4 . Final dataset"),
  xmin = c(-40, 63, -40, -40), xmax = c(63, 160, 160, 160),
  ymin = c(82, 82, 17, -20), ymax = c(127.5, 127.5, 80, 15),
  key  = c("blue", "orange", "yellow", "green")
)
zones$fill <- palette[zones$key]
zones$tcol <- zone_title_pal[zones$key]
zones$ty <- zones$ymax - 2.2

zones_poly <- do.call(rbind, lapply(seq_len(nrow(zones)), function(i)
  rect_xy((zones$xmin[i] + zones$xmax[i])/2, (zones$ymin[i] + zones$ymax[i])/2,
          zones$xmax[i] - zones$xmin[i], zones$ymax[i] - zones$ymin[i], paste0("z", i))))
zones_poly$fill <- rep(zones$fill, each = 4)

# Top modules (parallel tasks)
bullet <- data.frame(
  x = c(rep(6.5, 4), rep(94.0, 4)),
  y = c(121, 118, 115, 112)
)
row_txt <- data.frame(
  x = c(rep(9.5, 4), rep(97.0, 4)),
  y = c(121, 118, 115, 112),
  label = c("116,643 samples from NCBI", "334 BioProjects from NCBI",
            "16S rRNA amplicons · DADA2", "All body sites (oral + extra-oral)",
            "Country (47 countries)", "Body site (7 oral categories)",
            "Health status (healthy + diseased)", "Country · site · platform · age · sex")
)

top_pill <- data.frame(id = c("tblue", "torange"), cx = c(11.5, 111), cy = c(95.5, 92),
                       w = c(74, 84), h = c(29, 21))
top_poly <- do.call(rbind, lapply(seq_len(nrow(top_pill)), function(i)
  rect_xy(top_pill$cx[i], top_pill$cy[i], top_pill$w[i], top_pill$h[i], top_pill$id[i])))

top_txt <- do.call(rbind, list(
  data.frame(x = 11.5,
             y = 95.5 + c(11.2, 8.0, 4.8, 1.6, -1.6, -4.8, -8.0, -11.2),
             label = c("116,643 samples", "after Step 1 · Sample retrieval",
                       "116,624 samples", "after Step 2 · SRA download validation",
                       "114,959 samples", "after Step 3 · Quality filtering & DADA2",
                       "Excluded: low-quality or denoising failures", "−1,665 samples"),
             size = c(5.2, 3.8, 5.2, 3.8, 5.2, 3.8, 3.5, 3.5),
             fontface = c("bold", "plain", "bold", "plain", "bold", "plain", "italic", "italic"),
             color = c(txt_main, txt_sub, txt_main, txt_sub, txt_main, txt_sub, txt_detail, txt_excl)),
  data.frame(x = 111.5,
             y = 92 + c(3.6, 0, -3.6),
             label = c("Harmonized attribute table", "country · site · status · age · sex",
                       "standardized metadata for linkage"),
             size = c(5.2, 3.8, 3.5), fontface = c("bold", "plain", "italic"),
             color = c(txt_main, txt_sub, txt_detail))
))

row_arrow <- data.frame(x = 111, y = 108, xend = 111, yend = 101.0)

merge_seg <- data.frame(x = c(11.5, 111, 11.5, 60), y = c(82, 83.25, 80, 80),
                        xend = c(11.5, 111, 111, 60), yend = c(80, 80, 80, 73.0))
merge_head <- arrowhead_down(60, 73.0, a = 2.2, b = 1.3)
merge_dot  <- data.frame(x = c(11.5, 111), y = c(80, 80))

# Funnel nodes
step_pill <- data.frame(
  id = paste0("s", 1:4),
  cx = rep(60, 4), cy = c(64, 44, 24, -3),
  w = rep(120, 4), h = c(17.5, 17.5, 17.5, 34)
)
step_poly <- do.call(rbind, lapply(seq_len(nrow(step_pill)), function(i)
  rect_xy(step_pill$cx[i], step_pill$cy[i], step_pill$w[i], step_pill$h[i], step_pill$id[i])))

step_line <- do.call(rbind, list(
  data.frame(x = 60, y = 64 + c(6.0, 2.0, -2.0, -6.0),
             label = c("113,233 samples", "after Step 4 · Chimera Filtering",
                       "Excluded: >30% of oral samples with a chimera rate >25%",
                       "−1,726 samples · −22 studies"),
             size = c(5.2, 4.2, 3.8, 3.8), fontface = c("bold", "plain", "italic", "italic"),
             color = c(txt_main, txt_sub, txt_detail, txt_excl)),
  data.frame(x = 60, y = 44 + c(6.0, 2.0, -2.0, -6.0),
             label = c("108,691 samples", "after Step 5 · Closed-Reference Mapping",
                       "Excluded: oral aggregate GG2 closed-reference mapping rate <70%",
                       "−4,542 samples · −24 studies"),
             size = c(5.2, 4.2, 3.8, 3.8), fontface = c("bold", "plain", "italic", "italic"),
             color = c(txt_main, txt_sub, txt_detail, txt_excl)),
  data.frame(x = 60, y = 24 + c(6.0, 2.0, -2.0, -6.0),
             label = c("80,788 samples", "after Step 6 · Site Filtering",
                       "Excluded: non-oral samples (site_group ≠ oral)",
                       "−27,903 samples"),
             size = c(5.2, 4.2, 3.8, 3.8), fontface = c("bold", "plain", "italic", "italic"),
             color = c(txt_main, txt_sub, txt_detail, txt_excl)),
  data.frame(x = 60, y = -3 + c(14.375, 9.75, 5.125, 0.5, -5.125, -9.75, -14.375),
             label = c("51,186 samples", "after Step 7 · Deduplication",
                       "Excluded: >1 sample per subject × site",
                       "−29,602 samples", "final curated dataset · 288 studies",
                       "42,224 subjects · 61,235 features",
                       "47 countries · 7 oral site categories"),
             size = c(5.2, 4.2, 3.8, 3.8, 4.2, 3.8, 3.8),
             fontface = c("bold", "plain", "italic", "italic", "plain", "italic", "italic"),
             color = c(txt_main, txt_sub, txt_detail, txt_excl, txt_sub, txt_detail, txt_detail))
))

final_div <- data.frame(x = 36, y = -5.25, xend = 84, yend = -5.25)

# Apply the text-scale factor to the in-panel text
top_txt$size   <- top_txt$size  * A_TXT
step_line$size <- step_line$size * A_TXT
step_arrow <- data.frame(x = rep(60, 3), y = c(55.25, 35.25, 15.25),
                         xend = rep(60, 3), yend = c(52.75, 32.75, 14.20))

# Left-align zone bullets with the "." of each "N . Title" step label,
# measured at render time so the offset adapts to device/fonts.
zone_font_pt <- 4.0 * A_TXT_SUB * ggplot2::.pt
grid::pushViewport(grid::viewport(gp = grid::gpar(fontsize = zone_font_pt, fontface = "bold")))
pref_blue_mm   <- as.numeric(grid::convertWidth(grid::stringWidth("1 "), "mm", valueOnly = TRUE))
pref_orange_mm <- as.numeric(grid::convertWidth(grid::stringWidth("2 "), "mm", valueOnly = TRUE))
grid::popViewport()
mg_mm       <- 6 / ggplot2::.pt
dev_w_mm    <- 13 * 25.4 - 2 * mg_mm
dev_h_mm    <- 10.6 * 25.4 - 2 * mg_mm
xrange      <- (160 - (-40)) * 1.1
yrange      <- (134 - (-22)) * 1.1
panel_w_mm  <- min(dev_w_mm, dev_h_mm * xrange / yrange)
units_per_mm <- xrange / panel_w_mm
dot_blue   <- (-40 + 2.5) + pref_blue_mm   * units_per_mm
dot_orange <- (63  + 2.5) + pref_orange_mm * units_per_mm
row_txt$x  <- c(rep(dot_blue, 4), rep(dot_orange, 4))
bullet$x   <- row_txt$x - 3.0

# Assemble plot
p2 <- ggplot() +
  geom_shape(data = zones_poly, aes(x = x, y = y, group = grp),
             fill = zones_poly$fill, radius = unit(0.35, "cm"), color = NA) +
  geom_text(data = zones, aes(x = xmin + 2.5, y = ty, label = name, color = tcol),
            hjust = 0, vjust = 1, fontface = "bold", size = 4.0 * A_TXT_SUB) +
  scale_color_identity() +
  geom_point(data = bullet, aes(x = x, y = y), shape = 16, size = 1.6, color = bullet_col) +
  geom_text(data = row_txt, aes(x = x, y = y, label = label), hjust = 0, size = 4.0 * A_TXT_SUB, color = txt_sub) +
  geom_segment(data = row_arrow, aes(x = x, y = y, xend = xend, yend = yend),
               color = arrow_col, linewidth = 1.0,
               arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_segment(data = merge_seg, aes(x = x, y = y, xend = xend, yend = yend),
               color = arrow_col, linewidth = 1.0) +
  geom_polygon(data = merge_head, aes(x = x, y = y), fill = arrow_col) +
  geom_point(data = merge_dot, aes(x = x, y = y), shape = 21, size = 1.8,
             fill = arrow_col, color = arrow_col) +
  geom_shape(data = top_poly, aes(x = x, y = y, group = grp),
             fill = node_fill, color = node_border, radius = unit(0.35, "cm"), linewidth = 0.8) +
  geom_text(data = top_txt, aes(x = x, y = y, label = label),
            size = top_txt$size, fontface = top_txt$fontface, color = top_txt$color) +
  geom_shape(data = step_poly, aes(x = x, y = y, group = grp),
             fill = node_fill, color = node_border, radius = unit(0.35, "cm"), linewidth = 0.8) +
  geom_text(data = step_line, aes(x = x, y = y, label = label),
            size = step_line$size, fontface = step_line$fontface, color = step_line$color) +
  geom_segment(data = final_div, aes(x = x, y = y, xend = xend, yend = yend),
               color = node_border, linewidth = 0.5) +
  geom_segment(data = step_arrow, aes(x = x, y = y, xend = xend, yend = yend),
               color = arrow_col, linewidth = 1.0,
               arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  coord_fixed(ratio = 1, xlim = c(-40, 160), ylim = c(-22, 134)) +
  theme_void() + theme(plot.margin = margin(6, 6, 6, 6, "pt"))

ggsave(file.path(PANEL_DIR, "pipeline_two_stage_funnel_HD.png"), plot = p2,
       width = 13, height = 10.6, device = "png", dpi = 600, bg = "white")
message("Saved: ", file.path(PANEL_DIR, "pipeline_two_stage_funnel_HD.png"))
