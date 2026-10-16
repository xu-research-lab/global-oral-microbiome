#!/usr/bin/env Rscript
# ============================================================================
# Figure 5 supplement B: UpSet of every shared consensus (genus, direction)
# pair (degree >= 2; 166 pairs, 101 intersections), sorted by size L→R and
# bars stacked by direction.
# Consensus = >= 2 DA tools, q <= 0.05, same direction.
# 20 diseases ordered by the five Panel A categories, then short name.
# Output: panel_s5b_upset_all.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ComplexUpset)
})

source(file.path(Sys.getenv("DA_WORK_DIR", "."), "script", "pub_theme.R"))

CONS_PATH <- file.path(COMP_DIR, "consensus_genera.csv")
CAT_PATH <- file.path(COMP_DIR, "cohort_overview", "disease_category_map.csv")
if (!file.exists(CONS_PATH)) stop("Missing: ", CONS_PATH)
if (!file.exists(CAT_PATH)) stop("Missing: ", CAT_PATH)

DIR_COLS <- c("Health-enriched" = "#2166AC", "Disease-enriched" = "#B2182B")
CAT_ORDER <- c(
  "Oral and dental diseases",
  "Malignancies",
  "Infectious and immune-mediated diseases",
  "Cardiometabolic and sleep-related disorders",
  "Gastrointestinal diseases"
)

upset_label <- function(x) {
  x <- disease_short(x)
  gsub("&", "-", x, fixed = TRUE)
}

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

# Top → bottom disease order.
dz <- read.csv(CAT_PATH, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    Category = factor(Category, levels = CAT_ORDER),
    SetLabel = sankey_disease_label(Disease)
  ) %>%
  arrange(Category, tolower(upset_label(Disease)))
stopifnot(nrow(dz) == 20L, !anyDuplicated(dz$SetLabel))
if (!all(unique(cons$Disease) %in% dz$Disease)) {
  stop("Consensus disease missing from the category map")
}

cons <- cons %>% left_join(dz %>% select(Disease, SetLabel), by = "Disease")
cons$Element <- paste(cons$Genus, cons$Direction)

set_labels <- dz$SetLabel
disease_sets <- split(cons$Element, cons$SetLabel)
for (lab in set_labels) {
  if (is.null(disease_sets[[lab]])) disease_sets[[lab]] <- character(0)
}
disease_sets <- disease_sets[set_labels]

all_elements <- sort(unique(unlist(disease_sets)))
mat <- as.data.frame(
  lapply(disease_sets, function(g) as.integer(all_elements %in% g)),
  check.names = FALSE
)
rownames(mat) <- all_elements
colnames(mat) <- set_labels

element_degree <- as.integer(rowSums(mat))
element_dir <- ifelse(
  grepl(" \\+$", rownames(mat)),
  "Disease-enriched", "Health-enriched"
)
stopifnot(length(all_elements) == 506L)
stopifnot(sum(element_degree == 1L) == 340L)
stopifnot(sum(element_degree >= 2L) == 166L)

shared_idx <- element_degree >= 2L
mat_shared <- mat[shared_idx, , drop = FALSE]
mat_shared$DirLabel <- factor(
  element_dir[shared_idx],
  levels = c("Health-enriched", "Disease-enriched")
)
# ComplexUpset draws the first set at the bottom.
set_bottom_to_top <- rev(set_labels)

p_upset <- upset(
  mat_shared, set_bottom_to_top,
  name = "Disease combination (all shared intersections)",
  width_ratio = 0.16,
  height_ratio = 0.68,
  min_size = 1,
  n_intersections = 200,
  sort_intersections = "descending",
  sort_intersections_by = "cardinality",
  sort_sets = FALSE,
  encode_sets = FALSE,
  keep_empty_groups = TRUE,
  wrap = FALSE,
  base_annotations = list(
    "Shared pairs" = intersection_size(
      counts = FALSE,
      mapping = aes(fill = DirLabel)
    ) +
      scale_fill_manual(values = DIR_COLS, name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      theme(
        # Right-hand bars are all size 1, leaving the upper right empty.
        legend.position = "inside",
        legend.position.inside = c(0.99, 0.97),
        legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 5.5, family = FONT)
      )
  ),
  set_sizes = FALSE,
  matrix = intersection_matrix(
    geom = geom_point(size = 0.7),
    segment = geom_segment(linewidth = 0.2)
  ),
  themes = upset_modify_themes(list(
    default = theme(text = element_text(size = 6.0, family = FONT)),
    intersections_matrix = theme(
      axis.text.y = element_text(size = 4.8, family = FONT, colour = "black"),
      axis.title.x = element_text(size = 6.2, family = FONT)
    )
  ))
)

# Wider than Cell Press full width, so not via save_panel().
w_mm <- 350
h_mm <- 140
stem <- "panel_s5b_upset_all"
grDevices::png(
  file.path(PUB_DIR, paste0(stem, ".png")),
  width = mm_to_in(w_mm), height = mm_to_in(h_mm),
  units = "in", res = 200, type = "cairo", bg = "white"
)
print(p_upset)
invisible(dev.off())
grDevices::pdf(
  file.path(PUB_DIR, paste0(stem, ".pdf")),
  width = mm_to_in(w_mm), height = mm_to_in(h_mm),
  family = FONT, useDingbats = FALSE
)
print(p_upset)
invisible(dev.off())
cat(sprintf("Wrote %s (%.0f x %.0f mm) | sets=%d shared pairs=%d\n",
            stem, w_mm, h_mm, length(set_labels), sum(shared_idx)))
