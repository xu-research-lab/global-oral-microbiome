#!/usr/bin/env Rscript
# ============================================================================
# Figure 5 supplement A: per-disease significant genera (q <= 0.05) split
# into the 7 mutually exclusive overlap regions of MaAsLin2 / ANCOM-BC2 /
# LinDA. Disease order: Panel A categories, then short name.
# Output: panel_s5a_by_methods.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)
})

source(file.path(Sys.getenv("DA_WORK_DIR", "."), "script", "pub_theme.R"))

SUM_PATH <- file.path(COMP_DIR, "overlap_summary_3way.csv")
CAT_PATH <- file.path(COMP_DIR, "cohort_overview", "disease_category_map.csv")
if (!file.exists(SUM_PATH)) stop("Missing: ", SUM_PATH)
if (!file.exists(CAT_PATH)) stop("Missing: ", CAT_PATH)

CAT_ORDER <- c(
  "Oral and dental diseases",
  "Malignancies",
  "Infectious and immune-mediated diseases",
  "Cardiometabolic and sleep-related disorders",
  "Gastrointestinal diseases"
)

# Disease names abbreviated as in panel_A_sankey_by_disease_abbrev.
SANKEY_ABBREV_EXTRA <- c("Periodontal" = "PD", "Dental Caries" = "DC",
                         "Hypertension" = "HTN")
sankey_disease_label <- function(x) {
  out <- disease_abbrev(x)
  hit <- out %in% names(SANKEY_ABBREV_EXTRA)
  out[hit] <- unname(SANKEY_ABBREV_EXTRA[out[hit]])
  out
}

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

cat_map <- read.csv(CAT_PATH, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(Category = factor(Category, levels = CAT_ORDER))
if (nrow(cat_map) != 20L) stop("Expected 20 diseases in the category map")
if (!setequal(cat_map$Disease, ov$Disease)) {
  stop("Category map diseases do not match overlap_summary_3way.csv")
}

region_levels <- c(
  "MaAsLin2 only", "ANCOM-BC2 only", "LinDA only",
  "MaAsLin2 & ANCOM-BC2", "MaAsLin2 & LinDA", "ANCOM-BC2 & LinDA",
  "All three"
)
fill_cols <- c(
  "MaAsLin2 only"        = "#FDD0A2",
  "ANCOM-BC2 only"       = "#C7E9C0",
  "LinDA only"           = "#DADAEB",
  "MaAsLin2 & ANCOM-BC2" = "#BCBD22",
  "MaAsLin2 & LinDA"     = "#1F9E89",
  "ANCOM-BC2 & LinDA"    = "#E7298A",
  "All three"            = "#08519C"
)

plot_df <- ov %>%
  transmute(
    Disease,
    DiseaseLabel = disease_short(Disease),
    `MaAsLin2 only` = as.integer(n_only_maaslin2),
    `ANCOM-BC2 only` = as.integer(n_only_ancombc2),
    `LinDA only` = as.integer(n_only_linda),
    `MaAsLin2 & ANCOM-BC2` = as.integer(n_maaslin2_ancombc2),
    `MaAsLin2 & LinDA` = as.integer(n_maaslin2_linda),
    `ANCOM-BC2 & LinDA` = as.integer(n_ancombc2_linda),
    `All three` = as.integer(n_all3)
  ) %>%
  left_join(cat_map %>% select(Disease, Category), by = "Disease") %>%
  mutate(
    n_union = `All three` + `MaAsLin2 & LinDA` + `ANCOM-BC2 & LinDA` +
      `MaAsLin2 & ANCOM-BC2` + `MaAsLin2 only` + `ANCOM-BC2 only` + `LinDA only`
  )
# First level at the bottom; reverse so Oral/dental is at the top.
plot_df <- plot_df %>% arrange(Category, tolower(DiseaseLabel))
y_levels <- rev(plot_df$DiseaseLabel)
plot_df$DiseaseLabel <- factor(plot_df$DiseaseLabel, levels = y_levels)

long <- plot_df %>%
  select(DiseaseLabel, all_of(region_levels)) %>%
  pivot_longer(all_of(region_levels), names_to = "Region", values_to = "n") %>%
  mutate(Region = factor(Region, levels = region_levels))

zero_df <- dplyr::filter(plot_df, n_union == 0)
legend_breaks <- setdiff(region_levels, "MaAsLin2 & ANCOM-BC2")

p_methods <- ggplot(long, aes(x = n, y = DiseaseLabel, fill = Region)) +
  geom_col(
    width = 0.78, colour = "white", linewidth = 0.12,
    position = position_stack(reverse = TRUE)
  ) +
  geom_point(
    data = zero_df, aes(x = 0, y = DiseaseLabel), inherit.aes = FALSE,
    shape = 21, fill = "white", colour = "grey55", size = 0.9, stroke = 0.3
  ) +
  geom_text(
    data = zero_df, aes(x = 0, y = DiseaseLabel, label = "0"),
    inherit.aes = FALSE, hjust = -0.7, colour = "grey45",
    size = BASE_PT / ggplot2::.pt, family = FONT
  ) +
  scale_fill_manual(values = fill_cols, breaks = legend_breaks, name = NULL) +
  scale_y_discrete(labels = sankey_disease_label) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.04)),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  labs(x = "Significant genera (q <= 0.05)", y = NULL) +
  theme_cell() +
  theme(
    axis.text.y = element_markdown(size = 5.8, colour = "black"),
    # Legend sits in the empty lower right.
    legend.position = "inside",
    legend.position.inside = c(0.99, 0.02),
    legend.justification = c(1, 0),
    legend.direction = "horizontal",
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.35, colour = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(2, 4, 2, 2)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

save_panel(p_methods, "panel_s5a_by_methods", width_mm = W_FULL_MM,
           height_mm = 89.52)

cat(sprintf(
  "By methods: diseases=%d | zero-union: %s\n",
  nrow(ov),
  if (nrow(zero_df)) paste(zero_df$Disease, collapse = ", ") else "(none)"
))
