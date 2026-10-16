#!/usr/bin/env Rscript
# ============================================================================
# Panel B: consensus (genus, direction) pairs per disease.
# Consensus = >= 2 DA tools, q <= 0.05, same direction.
# Unit = (GTDB genus, direction). 400 genera -> 845 associations -> 506 pairs
# (340 disease-specific + 166 shared). Obesity is included (0 consensus).
#
# One bar per disease sorted by total; within each direction,
# disease-specific pairs are hatched. Output: panel_B_total.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
})

source(file.path(Sys.getenv("DA_WORK_DIR", "."), "script", "pub_theme.R"))

CONS_PATH <- file.path(COMP_DIR, "consensus_genera.csv")
CAT_PATH <- file.path(COMP_DIR, "cohort_overview", "disease_category_map.csv")
if (!file.exists(CONS_PATH)) stop("Missing: ", CONS_PATH)
if (!file.exists(CAT_PATH)) stop("Missing: ", CAT_PATH)

DIR_COLS <- c("Health-enriched" = "#2166AC", "Disease-enriched" = "#B2182B")

# Disease names abbreviated as in panel_A_sankey_by_disease_abbrev.
SANKEY_ABBREV_EXTRA <- c("Periodontal" = "PD", "Dental Caries" = "DC",
                         "Hypertension" = "HTN")
sankey_disease_label <- function(x) {
  out <- disease_abbrev(x)
  hit <- out %in% names(SANKEY_ABBREV_EXTRA)
  out[hit] <- unname(SANKEY_ABBREV_EXTRA[out[hit]])
  out
}

cons <- read.csv(CONS_PATH, stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(!is.na(Genus), Genus != "", !is.na(Disease), Disease != "") %>%
  distinct(Disease, Genus, Direction)
stopifnot(nrow(cons) == 845L)
stopifnot(dplyr::n_distinct(cons$Genus) == 400L)

cat_map <- read.csv(CAT_PATH, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(cat_map) == 20L)
if (!all(unique(cons$Disease) %in% cat_map$Disease)) {
  stop("Consensus disease missing from the category map")
}

cons$DirLabel <- ifelse(cons$Direction == "+", "Disease-enriched", "Health-enriched")
cons$Element <- paste(cons$Genus, cons$Direction)

elem_meta <- cons %>% count(Element, name = "degree")
stopifnot(nrow(elem_meta) == 506L)
stopifnot(sum(elem_meta$degree == 1L) == 340L)
stopifnot(sum(elem_meta$degree >= 2L) == 166L)

# One corner-to-corner diagonal per square tile repeats into continuous 45°
# stripes. Needs a device with TilingPattern support (pdf, cairo png).
hatch_fill <- function(base, spacing_mm = 1.1, lwd = 0.4) {
  tile <- grid::grobTree(
    grid::rectGrob(gp = grid::gpar(fill = base, col = NA)),
    grid::linesGrob(
      c(0, 1), c(0, 1),
      gp = grid::gpar(col = "white", lwd = lwd, lineend = "square")
    )
  )
  grid::pattern(
    tile,
    width = unit(spacing_mm, "mm"), height = unit(spacing_mm, "mm"),
    extend = "repeat"
  )
}

d_total <- expand.grid(
  Disease = cat_map$Disease,
  DirLabel = c("Health-enriched", "Disease-enriched"),
  Sharing = c("Shared (>=2 diseases)", "Disease-specific"),
  stringsAsFactors = FALSE
) %>%
  left_join(
    cons %>%
      left_join(elem_meta, by = "Element") %>%
      mutate(Sharing = ifelse(degree >= 2L,
                              "Shared (>=2 diseases)", "Disease-specific")) %>%
      count(Disease, DirLabel, Sharing, name = "n"),
    by = c("Disease", "DirLabel", "Sharing")
  ) %>%
  mutate(n = ifelse(is.na(n), 0L, as.integer(n)))

n_by <- function(dir, sharing) {
  sum(d_total$n[d_total$DirLabel == dir & d_total$Sharing == sharing])
}
n_hc_shared <- n_by("Health-enriched", "Shared (>=2 diseases)")
n_hc_spec <- n_by("Health-enriched", "Disease-specific")
n_ds_shared <- n_by("Disease-enriched", "Shared (>=2 diseases)")
n_ds_spec <- n_by("Disease-enriched", "Disease-specific")
stopifnot(n_distinct(d_total$Disease) == 20L)
stopifnot(sum(d_total$n) == 845L)
stopifnot(n_hc_shared == 355L, n_hc_spec == 174L)
stopifnot(n_ds_shared == 150L, n_ds_spec == 166L)
stopifnot(n_hc_shared + n_ds_shared == 505L, n_hc_spec + n_ds_spec == 340L)

# First y level is drawn at the bottom: ascending total puts the largest on
# top; ties go to more disease-enriched, then alphabetical from the top.
total_order <- d_total %>%
  group_by(Disease) %>%
  summarise(
    n_total = sum(n),
    n_disease = sum(n[DirLabel == "Disease-enriched"]),
    .groups = "drop"
  ) %>%
  mutate(Label = sankey_disease_label(Disease)) %>%
  arrange(n_total, n_disease, desc(Label))

# position_stack(reverse = TRUE) draws the first level leftmost:
# blue solid, blue hatched, red solid, red hatched.
FILL_LEVELS <- c(
  "Health-enriched | shared", "Health-enriched | disease-specific",
  "Disease-enriched | shared", "Disease-enriched | disease-specific"
)
fill_values <- list(
  "Health-enriched | shared" = DIR_COLS[["Health-enriched"]],
  "Health-enriched | disease-specific" = hatch_fill(DIR_COLS[["Health-enriched"]]),
  "Disease-enriched | shared" = DIR_COLS[["Disease-enriched"]],
  "Disease-enriched | disease-specific" = hatch_fill(DIR_COLS[["Disease-enriched"]])
)
fill_labels <- c(
  "Health-enriched | shared" =
    sprintf("Health-enriched, shared (n = %d)", n_hc_shared),
  "Health-enriched | disease-specific" =
    sprintf("Health-enriched, disease-specific (n = %d)", n_hc_spec),
  "Disease-enriched | shared" =
    sprintf("Disease-enriched, shared (n = %d)", n_ds_shared),
  "Disease-enriched | disease-specific" =
    sprintf("Disease-enriched, disease-specific (n = %d)", n_ds_spec)
)

d_total <- d_total %>%
  mutate(
    Label = factor(sankey_disease_label(Disease), levels = total_order$Label),
    FillKey = factor(
      paste(DirLabel, ifelse(Sharing == "Disease-specific",
                             "disease-specific", "shared"), sep = " | "),
      levels = FILL_LEVELS
    )
  )
stopifnot(!anyNA(d_total$FillKey))
zero_total <- total_order %>%
  filter(n_total == 0L) %>%
  mutate(Label = factor(Label, levels = total_order$Label))

p_total <- ggplot(d_total, aes(x = n, y = Label, fill = FillKey)) +
  geom_col(
    width = 0.78, colour = "white", linewidth = 0.12,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    data = zero_total,
    aes(x = 0, y = Label, label = "0"),
    inherit.aes = FALSE, hjust = -0.35, colour = "grey45",
    size = 5.2 / ggplot2::.pt, family = FONT
  ) +
  scale_fill_manual(
    values = fill_values, breaks = FILL_LEVELS, labels = fill_labels,
    name = NULL,
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(x = "Consensus (genus, direction) pairs", y = NULL) +
  theme_cell() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "left",
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.35, colour = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 5.5, colour = "black"),
    plot.margin = margin(4, 6, 2, 2)
  )

save_panel(p_total, "panel_B_total", width_mm = W_FULL_MM, height_mm = 90)
top_row <- tail(total_order, 1)
bottom_row <- head(total_order, 1)
cat(sprintf(
  "B total: top = %s (%d) | bottom = %s (%d)\n",
  top_row$Label, top_row$n_total, bottom_row$Label, bottom_row$n_total
))
cat("B total order (top → bottom):",
    paste(rev(total_order$Label), collapse = ", "), "\n")
