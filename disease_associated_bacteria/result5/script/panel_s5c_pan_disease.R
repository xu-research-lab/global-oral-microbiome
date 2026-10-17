#!/usr/bin/env Rscript
# ============================================================================
# Panel E: pan-disease signals separated by disease breadth and directional
# coherence (any-two consensus; GTDB genus names kept intact).
#
#   E1 — curated genus x disease signed heatmap; row names coloured by class
#   E2 — GTDB lineage-resolution diverging bars
#
# The former breadth x balance scatter is exported as a supplement
# (panel_S_breadth_balance / panel_S_breadth_consistency).
#
# Consensus unit = (Disease, Genus, Direction) from consensus_genera.csv.
# Do NOT collapse GTDB suffixes (Fusobacterium_A ≠ Fusobacterium_C).
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(patchwork)
  library(scales)
  library(grid)
})

source(file.path(Sys.getenv("DA_WORK_DIR", "."), "script", "pub_theme.R"))

CONS_PATH <- file.path(COMP_DIR, "consensus_genera.csv")
CAT_PATH <- file.path(COMP_DIR, "cohort_overview", "disease_category_map.csv")
if (!file.exists(CONS_PATH)) stop("Missing: ", CONS_PATH)
if (!file.exists(CAT_PATH)) stop("Missing: ", CAT_PATH)

DIR_COLS <- c("Health-enriched" = "#2166AC", "Disease-enriched" = "#B2182B")
CLASS_COLS <- c(
  "Bidirectional"                = "#636363",
  "Consistent disease-enriched"  = "#B2182B",
  "Consistent health-enriched"   = "#2166AC",
  "Few diseases"                 = "#D9D9D9"
)

# Curated labels for E1 (and core of E2 rows)
LABEL_BIDIR <- c(
  "Dialister_36337", "Dialister_37762", "Bacteroides_H_857956",
  "Parvimonas", "Selenomonas_A_42762", "Treponema_B_986141",
  "Veillonella_A", "Fusobacterium_C", "Leptotrichia_A_993758"
)
LABEL_DS <- c("Limosilactobacillus", "Bifidobacterium_388775", "Scardovia")
LABEL_HC <- c(
  "Haemophilus_D_735815", "Cardiobacterium", "Corynebacterium",
  "Gemella", "Leptotrichia_A_993641", "Rothia"
)
LABEL_ALL <- c(LABEL_BIDIR, LABEL_DS, LABEL_HC)

# Extra heatmap rows (context-skewed + remaining core bidirectional from plan)
HEATMAP_EXTRA_BIDIR <- c(
  "Treponema_B_986142", "Peptostreptococcus", "Tannerella",
  "Treponema_C", "Treponema_D", "Prevotella", "Anaeroglobus_38550"
)
HEATMAP_EXTRA_HC <- c("Lachnoanaerobaculum", "Actinomyces")

# E3 lineage complexes
LINEAGE_MAP <- list(
  Fusobacterium = c("Fusobacterium_A", "Fusobacterium_C"),
  Dialister     = c("Dialister_36337", "Dialister_37762"),
  Treponema     = c(
    "Treponema_B_986141", "Treponema_B_986142", "Treponema_C", "Treponema_D"
  ),
  Leptotrichia  = c("Leptotrichia_A_993758", "Leptotrichia_A_993641"),
  Haemophilus   = c(
    "Haemophilus_D_735815", "Haemophilus_A",
    "Haemophilus_D_734546", "Haemophilus_D_736121"
  )
)

# Anatomical groups as in panel_A_sankey_AnatomicalGroup; used to validate the
# disease set and to order the E1 source-table columns.
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

# Expected heatmap size per minimum breadth (rows = genera, cols = diseases).
E1_EXPECTED_ROWS <- c("3" = 109L, "4" = 74L, "5" = 52L)
E1_EXPECTED_COLS <- c("3" = 18L, "4" = 17L, "5" = 17L)

# ---------------------------------------------------------------------------
# Load + metrics
# ---------------------------------------------------------------------------
cons <- read.csv(CONS_PATH, stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(!is.na(Genus), Genus != "", !is.na(Disease), Disease != "") %>%
  mutate(
    Direction = as.character(Direction),
    n_methods = as.integer(n_methods)
  )
stopifnot(nrow(cons) == 845L)
stopifnot(n_distinct(cons$Disease) == 19L)
stopifnot(n_distinct(cons$Genus) == 400L)

classify_signal <- function(n_diseases, n_plus, n_minus, min_breadth) {
  cls <- case_when(
    n_diseases < min_breadth ~ "Few diseases",
    n_minus == 0L & n_plus == n_diseases ~ "Consistent disease-enriched",
    n_plus == 0L & n_minus == n_diseases ~ "Consistent health-enriched",
    n_plus > 0L & n_minus > 0L ~ "Bidirectional",
    TRUE ~ "Few diseases"
  )
  factor(cls, levels = names(CLASS_COLS))
}

metrics <- cons %>%
  group_by(Genus) %>%
  summarise(
    n_diseases = n(),
    n_plus = sum(Direction == "+"),
    n_minus = sum(Direction == "-"),
    n_2of3 = sum(n_methods == 2L),
    n_3of3 = sum(n_methods == 3L),
    diseases_plus = paste(sort(Disease[Direction == "+"]), collapse = "; "),
    diseases_minus = paste(sort(Disease[Direction == "-"]), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    Directional_balance = (n_plus - n_minus) / n_diseases,
    Directional_consistency = pmax(n_plus, n_minus) / n_diseases,
    signal_class = classify_signal(n_diseases, n_plus, n_minus, 5L)
  )

# Assertions against plan examples
assert_genus <- function(g, n, np, nm) {
  row <- metrics %>% filter(Genus == g)
  if (nrow(row) != 1L) stop("Missing genus: ", g)
  if (row$n_diseases != n || row$n_plus != np || row$n_minus != nm) {
    stop(sprintf(
      "%s expected %d/%d+/%d- got %d/%d+/%d-",
      g, n, np, nm, row$n_diseases, row$n_plus, row$n_minus
    ))
  }
}
assert_genus("Dialister_36337", 8, 4, 4)
assert_genus("Limosilactobacillus", 7, 7, 0)
assert_genus("Haemophilus_D_735815", 7, 0, 7)
assert_genus("Fusobacterium_C", 6, 2, 4)
assert_genus("Fusobacterium_A", 1, 1, 0)
assert_genus("Dialister_37762", 6, 3, 3)
assert_genus("Leptotrichia_A_993758", 7, 2, 5)
assert_genus("Leptotrichia_A_993641", 5, 0, 5)

n_ge5 <- sum(metrics$n_diseases >= 5L)
if (n_ge5 != 52L) {
  stop(sprintf("Expected 52 genera with breadth>=5, got %d", n_ge5))
}

metrics_out <- file.path(COMP_DIR, "consensus_bubble", "genus_directional_metrics.csv")
dir.create(dirname(metrics_out), showWarnings = FALSE, recursive = TRUE)
write.csv(
  metrics %>%
    arrange(desc(n_diseases), Directional_consistency, Genus) %>%
    select(
      Genus, n_diseases, n_plus, n_minus,
      Directional_consistency, Directional_balance, signal_class,
      diseases_plus, diseases_minus, n_2of3, n_3of3
    ),
  metrics_out,
  row.names = FALSE
)
cat("Wrote:", metrics_out, "\n")
cat(sprintf(
  "Genera: %d | breadth>=5: %d | disease-all+: %d | health-all-: %d | bidirectional: %d\n",
  nrow(metrics), n_ge5,
  sum(metrics$signal_class == "Consistent disease-enriched"),
  sum(metrics$signal_class == "Consistent health-enriched"),
  sum(metrics$signal_class == "Bidirectional")
))

# ---------------------------------------------------------------------------
# E1: breadth x signed directional balance
# ---------------------------------------------------------------------------
metrics$label_flag <- metrics$Genus %in% LABEL_ALL
label_df <- metrics %>% filter(label_flag)

# Slight jitter only on integer x for overlapping points at same breadth
set.seed(42)
plot_e1 <- metrics %>%
  mutate(
    x_plot = n_diseases + runif(n(), -0.12, 0.12),
    point_alpha = ifelse(n_diseases >= 5L, 0.85, 0.25),
    point_col = as.character(signal_class)
  )

p_e1 <- ggplot(plot_e1, aes(x = x_plot, y = Directional_balance)) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf, ymin = -0.34, ymax = 0.34,
    fill = "grey90", alpha = 0.55
  ) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey40") +
  geom_hline(
    yintercept = c(-0.34, 0.34),
    linetype = "dashed", linewidth = 0.3, colour = "grey45"
  ) +
  geom_hline(
    yintercept = c(-0.60, 0.60),
    linetype = "dotted", linewidth = 0.3, colour = "grey55"
  ) +
  geom_point(
    aes(size = n_3of3, colour = point_col, alpha = point_alpha),
    stroke = 0
  ) +
  geom_point(
    data = label_df,
    aes(x = n_diseases, y = Directional_balance, size = n_3of3),
    shape = 21, fill = NA, colour = "black", stroke = 0.45, alpha = 1
  ) +
  geom_text_repel(
    data = label_df,
    aes(x = n_diseases, y = Directional_balance, label = Genus),
    size = 2.0, family = FONT, fontface = "italic",
    max.overlaps = 40, box.padding = 0.28, point.padding = 0.15,
    segment.size = 0.2, segment.colour = "grey50",
    min.segment.length = 0, seed = 42, force = 0.7
  ) +
  scale_colour_manual(values = CLASS_COLS, name = NULL) +
  scale_size_continuous(
    name = "No. 3/3 associations",
    range = c(1.2, 4.8),
    breaks = c(0, 1, 2, 3, 5)
  ) +
  scale_alpha_identity() +
  scale_x_continuous(
    breaks = 1:8,
    limits = c(0.5, 8.7),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(-1.08, 1.08),
    breaks = c(-1, -0.6, -0.34, 0, 0.34, 0.6, 1),
    labels = c("-1", "-0.60", "-0.34", "0", "+0.34", "+0.60", "+1")
  ) +
  annotate(
    "text",
    x = 8.55, y = 0, label = "Strong bidirectional\nzone",
    hjust = 1, vjust = 0.5, size = 2.0, colour = "grey40", family = FONT,
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = 0.7, y = 1.0, label = "Disease-enriched",
    hjust = 0, size = 2.1, colour = DIR_COLS[["Disease-enriched"]],
    family = FONT, fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.7, y = -1.0, label = "Health-enriched",
    hjust = 0, size = 2.1, colour = DIR_COLS[["Health-enriched"]],
    family = FONT, fontface = "bold"
  ) +
  labs(
    x = "Disease breadth (n diseases)",
    y = "Directional balance  (n+ - n-) / n"
  ) +
  theme_cell() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(2.8, "mm"),
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 6),
    plot.margin = margin(2, 4, 2, 2)
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 2.2, alpha = 1)),
    size = guide_legend(order = 2)
  )

# ---------------------------------------------------------------------------
# E1 (heatmap): genera with breadth >= min_breadth; names coloured by class
# relative to that threshold; columns ordered by number of cells shown.
# ---------------------------------------------------------------------------
cat_map <- read.csv(CAT_PATH, stringsAsFactors = FALSE, check.names = FALSE)

anat_map <- data.frame(
  Disease = unlist(ANATOMICAL_MEMBERS, use.names = FALSE),
  Group = rep(names(ANATOMICAL_MEMBERS), lengths(ANATOMICAL_MEMBERS)),
  stringsAsFactors = FALSE
)
if (anyDuplicated(anat_map$Disease)) {
  stop("Disease in more than one anatomical group: ",
       paste(anat_map$Disease[duplicated(anat_map$Disease)], collapse = ", "))
}
if (!setequal(anat_map$Disease, cat_map$Disease)) {
  stop("Anatomical groups disagree with disease_category_map.csv: ",
       paste(sort(union(setdiff(anat_map$Disease, cat_map$Disease),
                        setdiff(cat_map$Disease, anat_map$Disease))),
             collapse = ", "))
}
stopifnot(all(unique(cons$Disease) %in% anat_map$Disease))
anat_map <- anat_map %>%
  mutate(
    Group = factor(Group, levels = names(ANATOMICAL_MEMBERS)),
    DiseaseLabel = disease_short(Disease)
  )

# Column labels as in panel_A_sankey_by_disease_abbrev.
SANKEY_ABBREV_EXTRA <- c("Periodontal" = "PD", "Dental Caries" = "DC",
                         "Hypertension" = "HTN")
sankey_disease_label <- function(x) {
  out <- disease_abbrev(x)
  hit <- out %in% names(SANKEY_ABBREV_EXTRA)
  out[hit] <- unname(SANKEY_ABBREV_EXTRA[out[hit]])
  out
}

# Shared top chrome so the name column and heatmap rows line up:
# angled x labels of the same size.
strip_theme <- theme(
  axis.text.x.top = element_text(
    size = 4.8, angle = 55, hjust = 0, vjust = 0, family = FONT
  )
)

# Tile geometry fixed in millimetres: a wide null panel would stretch the
# cells to the full 174 mm and leave the heatmap looking sparse. The width
# left over carries the legends.
E1_CELL_W_MM <- 5.0
E1_HEAT_PAD_MM <- 4

build_e1 <- function(min_breadth, table_suffix = "") {
  key <- as.character(min_breadth)

  row_meta <- metrics %>%
    filter(n_diseases >= min_breadth) %>%
    mutate(
      signal_class = classify_signal(n_diseases, n_plus, n_minus, min_breadth),
      block = case_when(
        signal_class == "Consistent disease-enriched" ~ 1L,
        signal_class == "Bidirectional" ~ 2L,
        TRUE ~ 3L
      ),
      frac_plus = n_plus / n_diseases
    ) %>%
    arrange(
      block,
      # Bidirectional: higher disease-enriched fraction toward the top
      if_else(block == 2L, -frac_plus, 0),
      if_else(block == 2L, -as.numeric(n_plus), 0),
      desc(n_diseases),
      Genus
    )
  if (nrow(row_meta) != E1_EXPECTED_ROWS[[key]]) {
    stop(sprintf("Expected %d genera with breadth>=%d, got %d",
                 E1_EXPECTED_ROWS[[key]], min_breadth, nrow(row_meta)))
  }
  stopifnot(!any(row_meta$signal_class == "Few diseases"))
  bidir_ord <- row_meta %>%
    filter(block == 2L) %>%
    transmute(Genus, n_plus, n_minus, n_diseases, frac_plus)
  cat(sprintf("[breadth>=%d] Bidirectional order (top → bottom):\n", min_breadth))
  print(bidir_ord, row.names = FALSE)
  row_levels <- rev(row_meta$Genus)

  # E1 source table: genera (plot order, top → bottom) × 19 consensus diseases.
  # Cells: + disease-enriched; − health-enriched; blank = no consensus.
  # Obesity is omitted (0 consensus); some diseases may be all-blank.
  dz_tbl <- anat_map %>%
    filter(Disease %in% unique(cons$Disease)) %>%
    arrange(Group, DiseaseLabel)
  e1_order <- as.character(row_meta$Genus)
  dir_mat <- matrix(
    "",
    nrow = length(e1_order),
    ncol = nrow(dz_tbl),
    dimnames = list(e1_order, dz_tbl$DiseaseLabel)
  )
  heat_e1 <- cons %>% filter(Genus %in% e1_order)
  for (i in seq_len(nrow(heat_e1))) {
    dir_mat[heat_e1$Genus[[i]], disease_short(heat_e1$Disease[[i]])] <-
      heat_e1$Direction[[i]]
  }
  e1_tbl <- row_meta %>%
    transmute(
      Genus,
      signal_class = as.character(signal_class),
      n_diseases, n_plus, n_minus,
      frac_plus = round(frac_plus, 3)
    )
  e1_tbl <- cbind(
    e1_tbl,
    as.data.frame(dir_mat[e1_tbl$Genus, , drop = FALSE], check.names = FALSE)
  )
  if (ncol(dir_mat) != 19L) {
    stop(sprintf("E1 table expected 19 disease columns, got %d", ncol(dir_mat)))
  }
  if (sum(dir_mat != "") != nrow(heat_e1)) {
    stop("E1 matrix cell count disagrees with consensus_genera.csv")
  }
  tbl_stem <- file.path(
    COMP_DIR, "consensus_bubble",
    paste0("panel_E1_genus_disease_direction_matrix", table_suffix)
  )
  write.csv(e1_tbl, paste0(tbl_stem, ".csv"), row.names = FALSE, na = "")
  write.table(
    e1_tbl, paste0(tbl_stem, ".tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE, na = ""
  )
  cat("Wrote:", paste0(tbl_stem, ".csv"), "\n")

  # Disease columns: only those with a cell among the shown genera.
  heat_cells <- cons %>%
    filter(Genus %in% row_levels) %>%
    mutate(DiseaseLabel = disease_short(Disease))

  dz_keep <- unique(cons$Disease[cons$Genus %in% row_levels])
  dz_drop <- setdiff(unique(cons$Disease), dz_keep)
  if (length(dz_drop)) {
    cat(sprintf("[breadth>=%d] omitting diseases with no cells: %s\n",
                min_breadth, paste(sort(dz_drop), collapse = ", ")))
  }

  # Columns ordered by the number of shown genera with a call, most first.
  dz_meta <- anat_map %>%
    filter(Disease %in% dz_keep) %>%
    left_join(count(heat_cells, Disease, name = "n_shown"), by = "Disease") %>%
    arrange(desc(n_shown), DiseaseLabel)
  cat(sprintf("[breadth>=%d] column order: %s\n", min_breadth,
              paste(sprintf("%s (%d)", dz_meta$DiseaseLabel, dz_meta$n_shown),
                    collapse = ", ")))

  col_levels <- dz_meta$DiseaseLabel
  if (length(col_levels) != E1_EXPECTED_COLS[[key]]) {
    stop(sprintf("Expected %d disease columns at breadth>=%d, got %d",
                 E1_EXPECTED_COLS[[key]], min_breadth, length(col_levels)))
  }
  stopifnot(setequal(as.character(unique(heat_cells$DiseaseLabel)), col_levels))
  cat(sprintf(
    "[breadth>=%d] E1 matrix: %d genera x %d diseases | non-empty cells=%d\n",
    min_breadth, length(row_levels), length(col_levels), sum(dir_mat != "")
  ))

  # Complete grid for blanks
  grid_df <- expand.grid(
    Genus = row_levels,
    DiseaseLabel = col_levels,
    stringsAsFactors = FALSE
  ) %>%
    left_join(heat_cells %>% select(Genus, DiseaseLabel, Direction),
              by = c("Genus", "DiseaseLabel")) %>%
    left_join(row_meta %>% select(Genus, signal_class), by = "Genus") %>%
    mutate(
      Genus = factor(Genus, levels = row_levels),
      DiseaseLabel = factor(as.character(DiseaseLabel), levels = col_levels),
      DirFill = factor(
        ifelse(is.na(Direction), NA_character_,
               ifelse(Direction == "+", "Disease-enriched", "Health-enriched")),
        levels = c("Health-enriched", "Disease-enriched")
      )
    )

  block_breaks <- row_meta %>%
    mutate(Genus = factor(Genus, levels = row_levels), y = as.integer(Genus)) %>%
    group_by(block) %>%
    summarise(ymax = max(y), .groups = "drop")
  sep_y <- head(sort(unique(block_breaks$ymax)), -1) + 0.5

  ann_present <- intersect(
    names(CLASS_COLS),
    as.character(unique(row_meta$signal_class))
  )
  lab_df <- row_meta %>%
    mutate(
      Genus = factor(Genus, levels = row_levels),
      Class = factor(as.character(signal_class), levels = names(CLASS_COLS))
    )

  # Invisible points so the class legend uses dots, not the "a" text glyph.
  legend_pts <- lab_df %>%
    distinct(Class) %>%
    mutate(
      DiseaseLabel = factor(col_levels[1], levels = col_levels),
      Genus = factor(row_levels[1], levels = row_levels)
    )

  # Name column: flush right, fixed millimetre width, no HTML axis.
  # Its white dummy x label reserves the same top chrome as the heatmap's
  # angled disease labels.
  p_labs <- ggplot(lab_df, aes(x = 1, y = Genus, colour = Class, label = Genus)) +
    geom_text(
      hjust = 1, vjust = 0.5, size = 1.65, fontface = "italic", family = FONT,
      show.legend = FALSE
    ) +
    geom_point(aes(x = 0.5, y = Genus), alpha = 0, size = 1.6) +
    scale_colour_manual(
      values = CLASS_COLS[ann_present],
      breaks = ann_present,
      name = "Signal class",
      guide = guide_legend(
        order = 1, ncol = 1, title.position = "top",
        override.aes = list(size = 2.0, alpha = 1)
      )
    ) +
    scale_x_continuous(
      limits = c(0, 1), expand = c(0, 0), position = "top",
      breaks = 1, labels = "COVID-19"
    ) +
    scale_y_discrete(expand = expansion(mult = 0, add = 0.55)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_cell() +
    strip_theme +
    theme(
      axis.text.x.top = element_text(colour = "white"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      axis.line = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.justification = "top",
      legend.key.size = unit(2.6, "mm"),
      legend.text = element_text(size = 5.4),
      legend.title = element_text(size = 5.8, face = "bold"),
      plot.margin = margin(10, 0.4, 4, 1.5)
    )

  p_heat <- ggplot(grid_df, aes(x = DiseaseLabel, y = Genus, fill = DirFill)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    geom_point(
      data = legend_pts,
      aes(x = DiseaseLabel, y = Genus, colour = Class),
      inherit.aes = FALSE,
      size = 1.6, alpha = 0
    ) +
    geom_hline(yintercept = sep_y, colour = "grey30", linewidth = 0.3) +
    scale_fill_manual(
      values = DIR_COLS,
      na.value = "grey95",
      na.translate = FALSE,
      name = "Direction",
      guide = guide_legend(order = 2, ncol = 1, title.position = "top")
    ) +
    scale_colour_manual(
      values = CLASS_COLS[ann_present],
      breaks = ann_present,
      name = "Signal class",
      guide = "none"
    ) +
    scale_x_discrete(position = "top", expand = c(0, 0), drop = TRUE,
                     labels = sankey_disease_label) +
    scale_y_discrete(expand = expansion(mult = 0, add = 0.55)) +
    labs(x = NULL, y = NULL) +
    theme_cell() +
    strip_theme +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.length.y = unit(0, "pt"),
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.justification = "top",
      legend.box.just = "left",
      legend.spacing.y = unit(2.0, "mm"),
      legend.key.size = unit(2.6, "mm"),
      legend.text = element_text(size = 5.4),
      legend.title = element_text(size = 5.8, face = "bold"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.35),
      plot.margin = margin(10, 4, 4, 0)
    )

  # Longest GTDB name at the label size, plus a small pad.
  grDevices::pdf(NULL)
  name_w_mm <- max(vapply(row_levels, function(s) {
    as.numeric(grid::convertX(grid::grobWidth(grid::textGrob(
      s, gp = grid::gpar(fontsize = 4.8, fontface = "italic", fontfamily = FONT)
    )), "mm"))
  }, numeric(1))) + 3.2
  grDevices::dev.off()

  heat_w_mm <- length(col_levels) * E1_CELL_W_MM + E1_HEAT_PAD_MM
  cat(sprintf("[breadth>=%d] E1 name column: %.1f mm | heatmap: %.1f mm\n",
              min_breadth, name_w_mm, heat_w_mm))
  p <- p_labs + p_heat + patchwork::guide_area() +
    plot_layout(
      design = "ABC",
      widths = unit(c(name_w_mm, heat_w_mm, 1), c("mm", "mm", "null")),
      guides = "collect"
    ) +
    # guide_area() centres the guide box vertically; the patchwork-level theme
    # is what pins it to the top of the tall variants.
    plot_annotation(theme = theme(legend.justification = "top"))
  list(plot = p, n_row = length(row_levels), genera = e1_order)
}

# Row pitch as before; the legend no longer takes a row of its own, and the
# anatomical strips (~5 mm) are gone.
e1_height <- function(n_row) 17.35 + 2.55 * n_row

save_e1 <- function(plot, stem, height_mm) {
  if (height_mm <= MAX_H_MM) {
    return(save_panel(plot, stem, width_mm = W_FULL_MM, height_mm = height_mm))
  }
  w_in <- mm_to_in(W_FULL_MM)
  h_in <- mm_to_in(height_mm)
  pdf_path <- file.path(PUB_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(PUB_DIR, paste0(stem, ".png"))
  ggsave(
    filename = pdf_path, plot = plot, width = w_in, height = h_in,
    units = "in", device = grDevices::pdf, family = FONT,
    useDingbats = FALSE, bg = "white", limitsize = FALSE
  )
  ggsave(
    filename = png_path, plot = plot, width = w_in, height = h_in,
    units = "in", dpi = 600, bg = "white", limitsize = FALSE
  )
  cat(sprintf("Wrote %s  (%.0f x %.0f mm, above the %.0f mm journal max) and PNG preview\n",
              pdf_path, W_FULL_MM, height_mm, MAX_H_MM))
}

e1_variants <- list(
  list(min_breadth = 5L, stem = "panel_s5c_pan_disease", suffix = ""),
  list(min_breadth = 4L, stem = "panel_s5c_pan_disease_min4", suffix = "_min4"),
  list(min_breadth = 3L, stem = "panel_s5c_pan_disease_min3", suffix = "_min3")
)
e1_built <- list()
for (v in e1_variants) {
  h <- e1_height(E1_EXPECTED_ROWS[[as.character(v$min_breadth)]])
  res <- build_e1(v$min_breadth, v$suffix)
  save_e1(patchwork::wrap_elements(full = res$plot), v$stem, h)
  e1_built[[as.character(v$min_breadth)]] <- res
}
p_e2 <- e1_built[["5"]]$plot

# ---------------------------------------------------------------------------
# E2 (lineage): GTDB lineage-resolution diverging bars
# ---------------------------------------------------------------------------
lineage_rows <- list()
for (complex_name in names(LINEAGE_MAP)) {
  for (g in LINEAGE_MAP[[complex_name]]) {
    row <- metrics %>% filter(Genus == g)
    if (!nrow(row)) {
      warning("E3 missing genus: ", g)
      next
    }
    lineage_rows[[length(lineage_rows) + 1L]] <- data.frame(
      Complex = complex_name,
      Genus = g,
      n_plus = row$n_plus,
      n_minus = row$n_minus,
      n_diseases = row$n_diseases,
      stringsAsFactors = FALSE
    )
  }
}
lineage_df <- bind_rows(lineage_rows) %>%
  mutate(
    Complex = factor(Complex, levels = names(LINEAGE_MAP)),
    Genus = factor(Genus, levels = rev(unlist(LINEAGE_MAP, use.names = FALSE)))
  )

lineage_long <- bind_rows(
  lineage_df %>% transmute(
    Complex, Genus, side = "Health-enriched", n = -n_minus, label = n_minus
  ),
  lineage_df %>% transmute(
    Complex, Genus, side = "Disease-enriched", n = n_plus, label = n_plus
  )
) %>%
  mutate(side = factor(side, levels = c("Health-enriched", "Disease-enriched")))

p_e3 <- ggplot(lineage_long, aes(x = n, y = Genus, fill = side)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "black") +
  geom_text(
    data = lineage_long %>% filter(label > 0, side == "Health-enriched"),
    aes(label = label),
    hjust = 1.25, size = 2.0, family = FONT, colour = "grey20"
  ) +
  geom_text(
    data = lineage_long %>% filter(label > 0, side == "Disease-enriched"),
    aes(label = label),
    hjust = -0.25, size = 2.0, family = FONT, colour = "grey20"
  ) +
  facet_grid(Complex ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = DIR_COLS, name = NULL) +
  scale_x_continuous(
    limits = c(-8, 8),
    breaks = c(-6, -4, -2, 0, 2, 4, 6),
    labels = c("6", "4", "2", "0", "2", "4", "6"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "No. diseases (health-enriched  |  disease-enriched)",
    y = NULL
  ) +
  theme_cell() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = unit(2.6, "mm"),
    strip.text.y.left = element_text(
      size = 5.5, face = "bold", angle = 0, family = FONT
    ),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey92", colour = NA),
    axis.text.y = element_text(size = 5.0, face = "italic"),
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.35, colour = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(4, 6, 2, 8)
  )

# Verify both Leptotrichia lineages present
stopifnot(all(c("Leptotrichia_A_993758", "Leptotrichia_A_993641") %in%
                as.character(lineage_df$Genus)))
stopifnot(all(c("Fusobacterium_A", "Fusobacterium_C") %in%
                as.character(lineage_df$Genus)))

# ---------------------------------------------------------------------------
# Assemble + save (heatmap + lineage only)
# ---------------------------------------------------------------------------
h_e2 <- e1_height(E1_EXPECTED_ROWS[["5"]])
h_e3 <- 63
height_mm <- min(MAX_H_MM, h_e2 + h_e3)

# Freeze E1 so stacking with E2 cannot stretch the name column.
p_e1_alone <- patchwork::wrap_elements(full = p_e2)
combined <- (p_e1_alone / p_e3) +
  plot_layout(heights = c(h_e2, h_e3)) +
  plot_annotation(tag_levels = list(c("E1", "E2")))

save_panel(combined, "panel_E_pan_disease",
           width_mm = W_FULL_MM, height_mm = height_mm)

save_panel(p_e1, "panel_S_breadth_balance",
           width_mm = W_FULL_MM, height_mm = 95)

# ---------------------------------------------------------------------------
# Supplementary: breadth x consistency (cannot distinguish + from −)
# ---------------------------------------------------------------------------
p_s <- ggplot(
  metrics %>% filter(n_diseases >= 2L),
  aes(x = n_diseases, y = Directional_consistency)
) +
  geom_point(
    aes(colour = signal_class, size = n_3of3),
    alpha = 0.7, stroke = 0
  ) +
  geom_hline(yintercept = 0.67, linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = 0.80, linetype = "dotted", linewidth = 0.3) +
  geom_text_repel(
    data = label_df,
    aes(label = Genus),
    size = 1.9, family = FONT, fontface = "italic",
    max.overlaps = 30, seed = 42, segment.size = 0.2
  ) +
  scale_colour_manual(values = CLASS_COLS, name = NULL) +
  scale_size_continuous(name = "No. 3/3", range = c(1.0, 4.0)) +
  scale_y_continuous(limits = c(0.48, 1.02), breaks = seq(0.5, 1, 0.1)) +
  labs(
    x = "Disease breadth (n diseases)",
    y = "Directional consistency  max(n+, n-) / n"
  ) +
  theme_cell() +
  theme(
    legend.position = "right",
    legend.key.size = unit(2.6, "mm"),
    legend.text = element_text(size = 5.5)
  )

save_panel(p_s, "panel_S_breadth_consistency",
           width_mm = W_FULL_MM, height_mm = 95)

cat(sprintf(
  "Done. Heatmap rows (breadth>=3/4/5) = %d/%d/%d | lineage bars=%d\n",
  e1_built[["3"]]$n_row, e1_built[["4"]]$n_row, e1_built[["5"]]$n_row,
  nrow(lineage_df)
))
