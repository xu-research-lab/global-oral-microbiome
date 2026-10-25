#!/usr/bin/env Rscript
# ============================================================================
# Panel B (main): cross-disease sharing of consensus (genus, direction) pairs.
# Consensus = >= 2 DA tools, q <= 0.05, same direction.
# Unit = (GTDB genus, direction). 400 genera -> 845 associations -> 506 pairs
# (340 disease-specific + 166 shared). Obesity is included (0 consensus).
#
#   B1 — UpSet of shared pairs only (degree >= 2), bars stacked by direction.
#        20 diseases ordered and coloured by the five Panel A categories.
#   B2 — per-disease disease-specific vs shared counts, stacked by direction.
#
# The three-tool Venn of all significant calls is the supplement
# panel_S_method_overlap.R. Panel D remains a table exporter.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggtext)
  library(ComplexUpset)
  library(grid)
})

source("/beegfs/xuhuan/global_oral_tmp/panel_scripts_edit/pub_theme.R")

CONS_PATH <- file.path(WORK_DIR, "consensus_genera.csv")
CAT_PATH <- file.path(WORK_DIR, "cohort_overview", "disease_category_map.csv")
if (!file.exists(CONS_PATH)) stop("Missing: ", CONS_PATH)
if (!file.exists(CAT_PATH)) stop("Missing: ", CAT_PATH)

N_SHARED_INTERSECTIONS <- 40L
DIR_COLS <- c("Health-enriched" = "#2166AC", "Disease-enriched" = "#B2182B")
CAT_ORDER <- c(
  "Oral and dental diseases",
  "Malignancies",
  "Infectious and immune-mediated diseases",
  "Cardiometabolic and sleep-related disorders",
  "Gastrointestinal diseases"
)
CAT_COLOURS <- c(
  "Oral and dental diseases" = "#2C7FB8",
  "Malignancies" = "#B2182B",
  "Infectious and immune-mediated diseases" = "#7B3294",
  "Cardiometabolic and sleep-related disorders" = "#E08214",
  "Gastrointestinal diseases" = "#1B7837"
)
CAT_SHORT <- c(
  "Oral and dental diseases" = "Oral/dental",
  "Malignancies" = "Malignancies",
  "Infectious and immune-mediated diseases" = "Infectious/immune",
  "Cardiometabolic and sleep-related disorders" = "Cardio/sleep",
  "Gastrointestinal diseases" = "GI"
)

# Anatomical site of disease, same groups, order and colours as
# panel_A_sankey_no_category. Used only by panel_B2_total_anatomical.
ANATOMICAL_MEMBERS <- list(
  "Oral cavity" = c(
    "Aphthous", "Dental Caries", "Oral Lichen Planus",
    "Periodontal Disease", "OSCC"
  ),
  "Head/neck & upper airway" = c(
    "Head and Neck Cancer", "Obstructive Sleep Apnea", "Thyroid Cancer"
  ),
  "Gastrointestinal tract" = c(
    "Esophageal Cancer", "Gastroesophageal Reflux Disease",
    "Colorectal Cancer", "IBD"
  ),
  "Respiratory / thoracic" = c("Lung Cancer", "COVID-19", "Breast Cancer"),
  "Systemic" = c("Diabetes", "Hypertension", "Obesity", "HIV", "RA")
)
ANATOMICAL_COLOURS <- c(
  "Oral cavity" = "#2C7FB8",
  "Head/neck & upper airway" = "#7B3294",
  "Gastrointestinal tract" = "#1B7837",
  "Respiratory / thoracic" = "#B2182B",
  "Systemic" = "#E08214"
)

upset_label <- function(x) {
  x <- disease_short(x)
  gsub("&", "-", x, fixed = TRUE)
}

cons <- read.csv(CONS_PATH, stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(!is.na(Genus), Genus != "", !is.na(Disease), Disease != "") %>%
  distinct(Disease, Genus, Direction)
stopifnot(nrow(cons) == 845L)
stopifnot(dplyr::n_distinct(cons$Genus) == 400L)

cat_map <- read.csv(CAT_PATH, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    Category = factor(Category, levels = CAT_ORDER),
    DiseaseLabel = upset_label(Disease)
  ) %>%
  arrange(Category, tolower(DiseaseLabel))
stopifnot(nrow(cat_map) == 20L)
if (!all(unique(cons$Disease) %in% cat_map$Disease)) {
  stop("Consensus disease missing from the category map")
}

# Top → bottom. ComplexUpset draws the first set at the bottom, so the
# vector passed to upset() is reversed after labels are built.
dz <- cat_map %>%
  group_by(Category) %>%
  mutate(is_first = row_number() == 1L) %>%
  ungroup() %>%
  mutate(
    cat_col = unname(CAT_COLOURS[as.character(Category)]),
    cat_lab = unname(CAT_SHORT[as.character(Category)]),
    axis_lab = ifelse(
      is_first,
      sprintf("<span style='color:%s'><b>%s</b>  %s</span>", cat_col, cat_lab, DiseaseLabel),
      sprintf("<span style='color:%s'>%s</span>", cat_col, DiseaseLabel)
    )
  )

cons <- cons %>%
  left_join(dz %>% select(Disease, DiseaseLabel), by = "Disease")
cons$DirLabel <- ifelse(cons$Direction == "+", "Disease-enriched", "Health-enriched")
cons$Element <- paste(cons$Genus, cons$Direction)

set_labels <- dz$DiseaseLabel
lab_map <- setNames(dz$axis_lab, dz$DiseaseLabel)

disease_sets <- split(cons$Element, cons$DiseaseLabel)
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
stopifnot(sum(element_dir == "Health-enriched") +
            sum(element_dir == "Disease-enriched") == 506L)

# ---------------------------------------------------------------------------
# B1: shared UpSet. Set order is bottom → top inside ComplexUpset, so pass
# the visual top-to-bottom order reversed.
# ---------------------------------------------------------------------------
shared_idx <- element_degree >= 2L
mat_shared <- mat[shared_idx, , drop = FALSE]
mat_shared$DirLabel <- factor(
  element_dir[shared_idx],
  levels = c("Health-enriched", "Disease-enriched")
)
set_bottom_to_top <- rev(set_labels)

# Column names are the coloured markdown labels (encode_sets = FALSE).
colnames(mat_shared)[match(set_bottom_to_top, colnames(mat_shared))] <-
  unname(lab_map[set_bottom_to_top])
set_md <- unname(lab_map[set_bottom_to_top])

p_b1 <- upset(
  mat_shared, set_md,
  name = "Disease combination (shared only)",
  width_ratio = 0.30,
  height_ratio = 0.70,
  min_size = 1,
  n_intersections = N_SHARED_INTERSECTIONS,
  sort_intersections_by = c("degree", "cardinality"),
  sort_sets = FALSE,
  encode_sets = FALSE,
  keep_empty_groups = TRUE,
  wrap = FALSE,
  base_annotations = list(
    "Shared pairs" = intersection_size(
      text = list(size = 1.6, family = FONT),
      mapping = aes(fill = DirLabel)
    ) +
      scale_fill_manual(values = DIR_COLS, name = NULL) +
      theme(
        legend.position = "top",
        legend.justification = "left",
        legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 5.5, family = FONT)
      )
  ),
  set_sizes = (
    upset_set_size(
      geom = geom_bar(aes(fill = DirLabel), width = 0.55)
    ) +
      scale_fill_manual(values = DIR_COLS, name = NULL, guide = "none") +
      ylab("shared pairs") +
      theme(
        axis.text.x = element_text(size = 5.0, family = FONT),
        axis.title.x = element_text(size = 6.0, family = FONT)
      )
  ),
  matrix = intersection_matrix(
    geom = geom_point(size = 1.05),
    segment = geom_segment(linewidth = 0.28)
  ),
  themes = upset_modify_themes(list(
    default = theme(text = element_text(size = 6.0, family = FONT)),
    intersections_matrix = theme(
      axis.text.y = element_markdown(size = 4.8, family = FONT),
      axis.title.x = element_text(size = 6.2, family = FONT)
    )
  ))
)

# ggplot stores left-axis text separately; only plots that already show
# y labels get markdown, so category colours render and set-size stays blank.
color_set_labels <- function(p) {
  if (inherits(p, "patchwork") && length(p$patches$plots)) {
    for (i in seq_along(p$patches$plots)) {
      p$patches$plots[[i]] <- color_set_labels(p$patches$plots[[i]])
    }
  }
  if (inherits(p, "ggplot") && !inherits(p$theme$axis.text.y, "element_blank")) {
    md <- element_markdown(size = 4.7, family = FONT, lineheight = 0.8)
    p$theme$axis.text.y <- md
    p$theme$axis.text.y.left <- md
  }
  p
}
p_b1 <- color_set_labels(p_b1)

# Inspection only: every shared intersection (101), sorted by size L→R.
# Does not replace B1. Plain black set names abbreviated as in
# panel_A_sankey_by_disease_abbrev.
SANKEY_ABBREV_EXTRA <- c("Periodontal" = "PD", "Dental Caries" = "DC",
                         "Hypertension" = "HTN")
sankey_disease_label <- function(x) {
  out <- disease_abbrev(x)
  hit <- out %in% names(SANKEY_ABBREV_EXTRA)
  out[hit] <- unname(SANKEY_ABBREV_EXTRA[out[hit]])
  out
}
plain_map <- setNames(sankey_disease_label(dz$Disease), dz$DiseaseLabel)
stopifnot(!anyDuplicated(plain_map))
mat_shared_plain <- mat[shared_idx, , drop = FALSE]
mat_shared_plain$DirLabel <- mat_shared$DirLabel
colnames(mat_shared_plain)[match(set_labels, colnames(mat_shared_plain))] <-
  unname(plain_map[set_labels])
set_plain <- unname(plain_map[set_bottom_to_top])

p_b1_all <- upset(
  mat_shared_plain, set_plain,
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
w_all_mm <- 350
h_all_mm <- 140
grDevices::png(
  file.path(PUB_DIR, "panel_B1_upset_all.png"),
  width = mm_to_in(w_all_mm), height = mm_to_in(h_all_mm),
  units = "in", res = 200, type = "cairo", bg = "white"
)
print(p_b1_all)
invisible(dev.off())
grDevices::pdf(
  file.path(PUB_DIR, "panel_B1_upset_all.pdf"),
  width = mm_to_in(w_all_mm), height = mm_to_in(h_all_mm),
  family = FONT, useDingbats = FALSE
)
print(p_b1_all)
invisible(dev.off())
cat(sprintf("Wrote panel_B1_upset_all (%.0f x %.0f mm), not a replacement for B1\n",
            w_all_mm, h_all_mm))

# ---------------------------------------------------------------------------
# B2: disease-specific vs shared, same top-to-bottom disease order
# ---------------------------------------------------------------------------
elem_meta <- data.frame(
  Element = rownames(mat),
  degree = element_degree,
  stringsAsFactors = FALSE
)
y_levels <- rev(set_labels)
d2 <- cons %>%
  left_join(elem_meta, by = "Element") %>%
  mutate(
    Sharing = ifelse(degree >= 2L, "Shared (>=2 diseases)", "Disease-specific"),
    DiseaseLabel = factor(DiseaseLabel, levels = y_levels)
  ) %>%
  count(DiseaseLabel, Sharing, DirLabel, name = "n")

grid_d2 <- expand.grid(
  DiseaseLabel = factor(set_labels, levels = y_levels),
  Sharing = c("Disease-specific", "Shared (>=2 diseases)"),
  DirLabel = c("Health-enriched", "Disease-enriched"),
  stringsAsFactors = FALSE
)
d2 <- grid_d2 %>%
  left_join(d2, by = c("DiseaseLabel", "Sharing", "DirLabel")) %>%
  mutate(n = ifelse(is.na(n), 0L, as.integer(n)))
d2$Sharing <- factor(
  d2$Sharing,
  levels = c("Disease-specific", "Shared (>=2 diseases)")
)
d2$DirLabel <- factor(d2$DirLabel, levels = c("Health-enriched", "Disease-enriched"))

n_hc <- sum(cons$DirLabel == "Health-enriched")
n_ds <- sum(cons$DirLabel == "Disease-enriched")
stopifnot(n_hc == 529L, n_ds == 316L)
dir_labels <- c(
  "Health-enriched"  = sprintf("Health-enriched (n = %s)", format(n_hc, big.mark = ",")),
  "Disease-enriched" = sprintf("Disease-enriched (n = %s)", format(n_ds, big.mark = ","))
)

zero_rows <- d2 %>%
  group_by(DiseaseLabel, Sharing) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  filter(n == 0)

p_b2 <- ggplot(d2, aes(x = n, y = DiseaseLabel, fill = DirLabel)) +
  geom_col(
    width = 0.78, colour = "white", linewidth = 0.12,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    data = zero_rows,
    aes(x = 0, y = DiseaseLabel, label = "0"),
    inherit.aes = FALSE, hjust = -0.35, colour = "grey45",
    size = 5.2 / ggplot2::.pt, family = FONT
  ) +
  facet_wrap(~ Sharing, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = DIR_COLS, labels = dir_labels, name = NULL) +
  scale_y_discrete(labels = lab_map) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.08)),
    breaks = scales::pretty_breaks(n = 4)
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
    axis.text.y = element_markdown(size = 4.8, lineheight = 0.75),
    axis.text.y.left = element_markdown(size = 4.8, lineheight = 0.75),
    strip.text = element_text(size = 5.0, face = "bold"),
    plot.margin = margin(1.5, 3, 1, 1)
  )

# ---------------------------------------------------------------------------
# B2 total (standalone): one bar per disease sorted by total; within each
# direction, disease-specific pairs are hatched. Not part of the stacked
# panel_B_consensus.
# ---------------------------------------------------------------------------
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

p_b2_total <- ggplot(d_total, aes(x = n, y = Label, fill = FillKey)) +
  geom_col(
    width = 0.78, colour = "white", linewidth = 0.24,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    data = zero_total,
    aes(x = 0, y = Label, label = "0"),
    inherit.aes = FALSE, hjust = -0.35, colour = "grey45",
    size = 10.4 / ggplot2::.pt, family = FONT
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
    axis.line.x = element_line(linewidth = 0.7, colour = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    legend.text = element_text(size = 10.5),
    legend.key.size = unit(5.2, "mm"),
    legend.spacing.x = unit(0.5, "mm"),
    plot.margin = margin(8, 12, 4, 4)
  )

# 2026-10-10: matched to panel_A_sankey_by_disease_abbrev (190 x 116 mm,
# font_scale 2 there), so the fonts are doubled as well.
save_panel(p_b2_total, "panel_B2_total", width_mm = 190, height_mm = 116)
top_row <- tail(total_order, 1)
bottom_row <- head(total_order, 1)
cat(sprintf(
  "B2 total: top = %s (%d) | bottom = %s (%d)\n",
  top_row$Label, top_row$n_total, bottom_row$Label, bottom_row$n_total
))
cat("B2 total order (top → bottom):",
    paste(rev(total_order$Label), collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# B2 total by anatomical site: identical bars, but grouped top-to-bottom by
# anatomical site (panel_A_sankey_no_category order); within a site still
# sorted by total. panel_B2_total keeps the pure count order.
# ---------------------------------------------------------------------------
anat_map <- data.frame(
  Disease = unlist(ANATOMICAL_MEMBERS, use.names = FALSE),
  Group = rep(names(ANATOMICAL_MEMBERS), lengths(ANATOMICAL_MEMBERS)),
  stringsAsFactors = FALSE
)
if (anyDuplicated(anat_map$Disease)) {
  stop("Disease in more than one anatomical group: ",
       paste(anat_map$Disease[duplicated(anat_map$Disease)], collapse = ", "))
}
stopifnot(setequal(anat_map$Disease, cat_map$Disease))

# First y level draws at the bottom, so build the reverse of the intended
# top-to-bottom order: last site first, smallest total first.
anat_order <- d_total %>%
  group_by(Disease) %>%
  summarise(
    n_total = sum(n),
    n_disease = sum(n[DirLabel == "Disease-enriched"]),
    .groups = "drop"
  ) %>%
  left_join(anat_map, by = "Disease") %>%
  mutate(
    Group = factor(Group, levels = names(ANATOMICAL_MEMBERS)),
    Label = disease_abbrev(Disease)
  ) %>%
  arrange(desc(as.integer(Group)), n_total, n_disease, desc(Label))
stopifnot(!anyNA(anat_order$Group), n_distinct(anat_order$Label) == 20L)

anat_lab_map <- stats::setNames(
  sprintf(
    "<span style='color:%s'>%s</span>",
    ANATOMICAL_COLOURS[as.character(anat_order$Group)], anat_order$Label
  ),
  anat_order$Label
)

d_anat <- d_total %>%
  mutate(Label = factor(disease_abbrev(Disease), levels = anat_order$Label))
zero_anat <- anat_order %>%
  filter(n_total == 0L) %>%
  mutate(Label = factor(Label, levels = anat_order$Label))

# Invisible points carry the anatomical colour legend; fill is already spoken
# for by direction x sharing.
legend_anat <- data.frame(
  n = 0,
  Label = factor(anat_order$Label[1], levels = anat_order$Label),
  Group = factor(names(ANATOMICAL_MEMBERS), levels = names(ANATOMICAL_MEMBERS))
)

p_b2_anat <- ggplot(d_anat, aes(x = n, y = Label, fill = FillKey)) +
  geom_col(
    width = 0.78, colour = "white", linewidth = 0.12,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    data = zero_anat,
    aes(x = 0, y = Label, label = "0"),
    inherit.aes = FALSE, hjust = -0.35, colour = "grey45",
    size = 5.2 / ggplot2::.pt, family = FONT
  ) +
  geom_point(
    data = legend_anat,
    aes(x = n, y = Label, colour = Group),
    inherit.aes = FALSE, alpha = 0, size = 0
  ) +
  scale_fill_manual(
    values = fill_values, breaks = FILL_LEVELS, labels = fill_labels,
    name = NULL,
    guide = guide_legend(order = 1, nrow = 2, byrow = TRUE)
  ) +
  scale_colour_manual(
    values = ANATOMICAL_COLOURS,
    name = "Anatomical site",
    guide = guide_legend(
      order = 2, nrow = 1, title.position = "left",
      override.aes = list(alpha = 1, size = 2.6, shape = 15)
    )
  ) +
  scale_y_discrete(labels = anat_lab_map) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(x = "Consensus (genus, direction) pairs", y = NULL) +
  theme_cell() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "left",
    legend.box.just = "left",
    legend.spacing.y = unit(1.0, "mm"),
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.35, colour = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_markdown(size = 5.5),
    axis.text.y.left = element_markdown(size = 5.5),
    plot.margin = margin(4, 6, 2, 2)
  )

save_panel(p_b2_anat, "panel_B2_total_anatomical",
           width_mm = W_FULL_MM, height_mm = 96)
cat("B2 anatomical order (top → bottom):\n")
anat_order %>%
  arrange(desc(row_number())) %>%
  transmute(Group = as.character(Group), Label, n_total) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ---------------------------------------------------------------------------
# Save B1 and B2, then stack. ComplexUpset + patchwork double-draws labels.
# ---------------------------------------------------------------------------
h_b1 <- 148
h_b2 <- 84
height_mm <- h_b1 + h_b2
if (height_mm > MAX_H_MM) stop("Panel B exceeds Cell Press height")
w_in <- mm_to_in(W_FULL_MM)

save_one <- function(plot_obj, stem, h_mm) {
  pdf_path <- file.path(PUB_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(PUB_DIR, paste0(stem, ".png"))
  h_in <- mm_to_in(h_mm)
  grDevices::pdf(pdf_path, width = w_in, height = h_in,
                 family = FONT, useDingbats = FALSE)
  print(plot_obj)
  invisible(dev.off())
  grDevices::png(png_path, width = w_in, height = h_in, units = "in",
                 res = 600, type = "cairo", bg = "white")
  print(plot_obj)
  invisible(dev.off())
  png_path
}

png_b1 <- save_one(p_b1, "panel_B1_upset_shared", h_b1)
png_b2 <- save_one(p_b2, "panel_B2_sharing", h_b2)

img1 <- png::readPNG(png_b1)
img2 <- png::readPNG(png_b2)
png_path <- file.path(PUB_DIR, "panel_B_consensus.png")
pdf_path <- file.path(PUB_DIR, "panel_B_consensus.pdf")
h_in <- mm_to_in(height_mm)
draw_stack <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    2, 1, heights = unit(c(h_b1, h_b2), "mm")
  )))
  grid.raster(img1, vp = viewport(layout.pos.row = 1), interpolate = TRUE)
  grid.text("B1", x = unit(1.5, "mm"), y = unit(1, "npc") - unit(2, "mm"),
            just = c("left", "top"),
            gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT),
            vp = viewport(layout.pos.row = 1))
  grid.raster(img2, vp = viewport(layout.pos.row = 2), interpolate = TRUE)
  grid.text("B2", x = unit(1.5, "mm"), y = unit(1, "npc") - unit(1.5, "mm"),
            just = c("left", "top"),
            gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT),
            vp = viewport(layout.pos.row = 2))
}
grDevices::png(png_path, width = w_in, height = h_in, units = "in",
               res = 600, type = "cairo", bg = "white")
draw_stack()
invisible(dev.off())
grDevices::pdf(pdf_path, width = w_in, height = h_in,
               family = FONT, useDingbats = FALSE)
draw_stack()
invisible(dev.off())

cat(sprintf("Wrote %s  (%.0f x %.0f mm)\n", png_path, W_FULL_MM, height_mm))
cat(sprintf(
  "B1 sets=%d shared pairs=%d intersections<=%d | B2 specific=%d shared=%d associations=%d (%d health / %d disease)\n",
  length(set_labels), sum(shared_idx), N_SHARED_INTERSECTIONS,
  sum(element_degree == 1L), sum(element_degree >= 2L),
  nrow(cons), n_hc, n_ds
))
