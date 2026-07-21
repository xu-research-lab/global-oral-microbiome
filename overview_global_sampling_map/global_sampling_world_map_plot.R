#!/usr/bin/env Rscript
# Global Sampling World Map
# Single output figure: world map + site×continent bubble matrix inset
#
# Output: global_sampling_world_map.pdf / .png

library(tidyverse)
library(countrycode)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(scales)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(scales)

sf_use_s2(FALSE)

# ---------------------------------------------------------------------------
# Global theme & font sizing
# ---------------------------------------------------------------------------
base_sz <- 17
text_sz <- base_sz / .pt

update_geom_defaults("text", list(family = "sans"))
update_geom_defaults("label", list(family = "sans"))

theme_set(theme_minimal(base_size = base_sz, base_family = "sans"))

# ---------------------------------------------------------------------------
# Color palettes
# ---------------------------------------------------------------------------
continent_colors <- c(
  "NAM" = "#D32F2F", "SA" = "#AB47BC", "EU" = "#4CAF50",
  "AS"  = "#FF9800", "AF" = "#795548", "OA" = "#29B6F6"
)
continent_order <- c("NAM", "SA", "EU", "AS", "AF", "OA")

site_cols_ordered <- c("Saliva", "Tooth", "Oral Mucosa", "Oral Swab",
                       "Tongue", "Gingiva", "Others")
site_colors <- c(
  "Saliva"      = "#E31A1C",
  "Tooth"       = "#FFB300",
  "Oral Mucosa" = "#1565C0",
  "Oral Swab"   = "#8E24AA",
  "Tongue"      = "#33A02C",
  "Gingiva"     = "#B15928",
  "Others"      = "#A6CEE3"
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
metadata_path <- paste0(
  "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-pipeline-map/",
  "input/metadata_oral_qc.tsv"
)
output_dir <- paste0(
  "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-pipeline-map/",
  "output"
)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Data extraction and country summarization
# ---------------------------------------------------------------------------
metadata <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  rename(
    BioProject      = bioproject_id,
    Country_Cleaned = country_cleaned,
    Site_Cleaned    = site_cleaned
  ) %>%
  mutate(
    Site_Cleaned = str_to_title(Site_Cleaned),
    Site_Cleaned = case_when(
      Site_Cleaned == "Other" ~ "Others",
      TRUE ~ Site_Cleaned
    )
  ) %>%
  filter(
    !(Site_Cleaned %in% c("Nasal", "Tonsils")),
    !is.na(Country_Cleaned),
    Country_Cleaned != ""
  ) %>%
  mutate(iso3c = countrycode(Country_Cleaned,
                              origin = "country.name",
                              destination = "iso3c"))

south_america_iso3c <- c(
  "ARG", "BOL", "BRA", "CHL", "COL", "ECU", "GUY",
  "PRY", "PER", "SUR", "URY", "VEN", "FLK", "GUF"
)

meta_clean <- metadata %>%
  mutate(
    Temp_Continent = countrycode(iso3c, "iso3c", "continent"),
    Continent = case_when(
      Temp_Continent == "Americas" &
        (iso3c %in% south_america_iso3c)                  ~ "SA",
      Temp_Continent == "Americas" &
        !(iso3c %in% south_america_iso3c)                 ~ "NAM",
      Temp_Continent == "Europe"                           ~ "EU",
      Temp_Continent == "Asia"                             ~ "AS",
      Temp_Continent == "Africa"                           ~ "AF",
      Temp_Continent == "Oceania"                          ~ "OA",
      iso3c == "ATA"                                       ~ "Antarctica",
      TRUE                                                 ~ Temp_Continent
    )
  ) %>%
  filter(!Continent %in% c("Other", "Antarctica")) %>%
  mutate(Continent = factor(Continent, levels = continent_order))

# World map base
world_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") %>%
  filter(name != "Antarctica")

# Country-level summary
country_summary <- meta_clean %>%
  group_by(Country_Cleaned, Continent) %>%
  summarise(
    n_samples = n(),
    n_studies = n_distinct(BioProject),
    .groups = "drop"
  ) %>%
  mutate(iso3c = countrycode(Country_Cleaned, "country.name", "iso3c"))

centroids <- world_sf %>%
  mutate(iso3c = iso_a3) %>%
  filter(!is.na(iso3c)) %>%
  st_point_on_surface() %>%
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(iso3c, lon, lat) %>%
  filter(!is.na(iso3c))

country_map_data <- country_summary %>%
  left_join(centroids, by = "iso3c") %>%
  mutate(
    lon = ifelse(is.na(lon) & iso3c == "BRB", -59.5, lon),
    lat = ifelse(is.na(lat) & iso3c == "BRB", 13.2, lat)
  ) %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  mutate(Continent = factor(Continent, levels = continent_order))

# Labels for countries with >= 100 samples
country_map_labels <- country_map_data %>%
  filter(n_samples >= 100) %>%
  mutate(
    country_label = paste0(
      Country_Cleaned, "\n(",
      n_studies, " bp; n=", scales::comma(n_samples), ")"
    )
  )

# ---------------------------------------------------------------------------
# Demographic summary table
# ---------------------------------------------------------------------------
col_names <- colnames(metadata)
col_lower <- tolower(col_names)

age_col <- if ("age_cleaned" %in% col_lower) {
  col_names[which(col_lower == "age_cleaned")]
} else if ("age" %in% col_lower) {
  col_names[which(col_lower == "age")]
} else NULL

sex_col <- if ("sex_cleaned" %in% col_lower) {
  col_names[which(col_lower == "sex_cleaned")]
} else if ("sex" %in% col_lower) {
  col_names[which(col_lower == "sex")]
} else NULL

bmi_col <- if ("bmi_cleaned" %in% col_lower) {
  col_names[which(col_lower == "bmi_cleaned")]
} else if ("bmi" %in% col_lower) {
  col_names[which(col_lower == "bmi")]
} else NULL

demo_df <- meta_clean
if (!is.null(age_col)) demo_df[[age_col]] <- as.numeric(demo_df[[age_col]])
if (!is.null(bmi_col)) demo_df[[bmi_col]] <- as.numeric(demo_df[[bmi_col]])

tbl_base <- demo_df %>%
  group_by(Continent) %>%
  summarise(
    Samples   = n(),
    Countries = n_distinct(Country_Cleaned),
    Studies   = n_distinct(BioProject),
    .groups   = "drop"
  )

if (!is.null(age_col)) {
  age_summ <- demo_df %>%
    filter(!is.na(!!sym(age_col))) %>%
    group_by(Continent) %>%
    summarise(
      Age = paste0(
        sprintf("%.1f", mean(!!sym(age_col), na.rm = TRUE)),
        "\u00B1",
        sprintf("%.1f", sd(!!sym(age_col), na.rm = TRUE))
      )
    )
  tbl_base <- tbl_base %>% left_join(age_summ, by = "Continent")
}

if (!is.null(bmi_col)) {
  bmi_summ <- demo_df %>%
    filter(!is.na(!!sym(bmi_col)),
           !!sym(bmi_col) > 10,
           !!sym(bmi_col) < 60) %>%
    group_by(Continent) %>%
    summarise(
      BMI = paste0(
        sprintf("%.1f", mean(!!sym(bmi_col), na.rm = TRUE)),
        "\u00B1",
        sprintf("%.1f", sd(!!sym(bmi_col), na.rm = TRUE))
      )
    )
  tbl_base <- tbl_base %>% left_join(bmi_summ, by = "Continent")
}

if (!is.null(sex_col)) {
  sex_summ <- demo_df %>%
    filter(!is.na(!!sym(sex_col))) %>%
    mutate(
      sex_label = case_when(
        tolower(!!sym(sex_col)) %in% c("male", "m")   ~ "M",
        tolower(!!sym(sex_col)) %in% c("female", "f") ~ "F",
        TRUE ~ "O"
      )
    ) %>%
    filter(sex_label %in% c("M", "F")) %>%
    count(Continent, sex_label) %>%
    group_by(Continent) %>%
    mutate(prop = round(n / sum(n) * 100, 0)) %>%
    pivot_wider(
      id_cols     = Continent,
      names_from  = sex_label,
      values_from = prop,
      values_fill = 0
    ) %>%
    mutate(`Sex (F/M)` = paste0(F, "% / ", M, "%")) %>%
    select(Continent, `Sex (F/M)`)
  tbl_base <- tbl_base %>% left_join(sex_summ, by = "Continent")
}

tbl_final <- data.frame(Continent = factor(continent_order,
                                           levels = continent_order)) %>%
  left_join(tbl_base, by = "Continent") %>%
  arrange(desc(Samples)) %>%
  mutate(
    Samples    = format(Samples, big.mark = ","),
    across(everything(), as.character),
    across(everything(), ~ replace_na(.x, "-")),
    Continent  = ifelse(Continent == "NAM", "NA", as.character(Continent))
  )

my_alignments <- c(0, 0, 0, 0, 0, 0, 0)

my_theme <- ttheme_minimal(
  core = list(
    fg_params = list(
      fontsize   = base_sz - 8,
      fontfamily = "sans",
      col        = "#303030",
      hjust      = my_alignments,
      x          = my_alignments
    ),
    bg_params = list(
      fill = c("#F9F9F9", "#FFFFFF"),
      col  = "grey90"
    )
  ),
  colhead = list(
    fg_params  = list(
      fontsize   = base_sz - 8,
      fontface   = "bold",
      fontfamily = "sans",
      col        = "#303030",
      hjust      = my_alignments,
      x          = my_alignments
    ),
    bg_params = list(fill = "#ECEFF1", col = "grey80")
  ),
  padding = unit(c(2.5, 2.5), "mm")
)

summary_table_grob <- tableGrob(tbl_final, rows = NULL, theme = my_theme)

# =============================================================================
# FIGURE 1: World map + site×continent bubble matrix inset
# =============================================================================

# --- World map base ---
p_b_base <- ggplot() +
  geom_sf(data = world_sf,
          fill = "#E8E0D4", color = "white", linewidth = 0.2) +
  geom_point(
    data = country_map_data,
    aes(x = lon, y = lat, size = n_samples, fill = Continent),
    shape = 21, color = "white", stroke = 0.35, alpha = 0.8
  ) +
  scale_size_continuous(
    name   = "Number of samples",
    range  = c(1, 25),
    limits = c(1, 50000),
    breaks = c(10, 100, 1000, 10000, 40000),
    labels = c("10", "100", "1K", "10K", "40K"),
    guide  = guide_legend(
      override.aes = list(
        shape = 21, fill = "grey70",
        color = "grey30", stroke = 0.5
      )
    )
  ) +
  scale_fill_manual(values = continent_colors, guide = "none") +
  geom_text_repel(
    data          = country_map_labels,
    aes(x = lon, y = lat, label = country_label, color = Continent),
    size          = text_sz * 0.60,
    fontface      = "bold",
    lineheight    = 0.85,
    segment.color = "#999999",
    segment.size  = 0.2,
    box.padding   = 0.4,
    max.overlaps  = Inf
  ) +
  scale_color_manual(values = continent_colors, guide = "none") +
  annotation_custom(
    summary_table_grob,
    xmin = -40, xmax = 120,
    ymin = -60, ymax = -30
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  theme_void(base_size = base_sz, base_family = "sans") +
  theme(
    legend.position = "top",
    legend.text     = element_text(size = base_sz * 0.8, color = "#303030"),
    legend.title    = element_text(size = base_sz * 0.9, face = "bold",
                                   color = "#303030", hjust = 0.5),
    legend.margin   = margin(t = 0, b = -20),
    plot.margin      = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")
  )

# --- Site x Continent bubble matrix (inset, no size legend) ---
site_continent_data <- meta_clean %>%
  filter(!is.na(Site_Cleaned)) %>%
  count(Continent, Site_Cleaned, name = "count") %>%
  group_by(Continent) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  mutate(
    Site_Cleaned    = factor(Site_Cleaned, levels = site_cols_ordered),
    Site_X          = as.numeric(Site_Cleaned) * 0.52,
    Continent_Label = ifelse(Continent == "NAM", "NA",
                             as.character(Continent)),
    Continent_Label = factor(Continent_Label,
                             levels = c("NA", "AS", "EU", "AF", "OA", "SA"))
  )

p_d <- ggplot(site_continent_data,
              aes(x = Site_X, y = Continent_Label,
                  size = count, fill = Site_Cleaned)) +
  geom_point(shape = 21, color = "white", stroke = 0.6, alpha = 0.90) +
  scale_size_continuous(
    name   = "Number of samples",
    range  = c(1.5, 8),
    breaks = c(50, 500, 5000, 15000),
    labels = c("50", "0.5K", "5K", "15K"),
    guide  = "none"
  ) +
  scale_fill_manual(name = "Site", values = site_colors,
                    limits = site_cols_ordered) +
  scale_x_continuous(
    breaks  = seq_along(site_cols_ordered) * 0.52,
    labels  = site_cols_ordered,
    expand  = expansion(add = c(0.12, 0.08))
  ) +
  scale_y_discrete(
    limits  = rev(levels(site_continent_data$Continent_Label)),
    expand  = expansion(add = 0.30)
  ) +
  labs(x = NULL, y = "Continent") +
  guides(fill = "none") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = base_sz - 4, base_family = "sans") +
  theme(
    axis.text.x   = element_text(
      size   = (base_sz - 4) * 0.70,
      face   = "bold", color = "#303030",
      angle  = 30, hjust = 1, vjust = 1
    ),
    axis.text.y   = element_text(
      size = (base_sz - 4) * 0.70,
      face = "bold", color = "#303030"
    ),
    axis.title.y  = element_text(
      size = (base_sz - 4) * 0.75,
      face = "bold", color = "#303030"
    ),
    axis.line     = element_line(color = "#303030", linewidth = 0.4),
    axis.ticks    = element_line(color = "#303030", linewidth = 0.3),
    legend.position   = "none",
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(4, 4, 4, 8, unit = "pt")
  )

# --- Composite ---
p_world_map <- p_b_base +
  inset_element(
    p_d,
    left   = 0.02,
    bottom = 0.05,
    right  = 0.18,
    top    = 0.40
  )

ggsave(file.path(output_dir, "global_sampling_world_map.pdf"),
       plot = p_world_map, width = 18, height = 10, device = "pdf")
ggsave(file.path(output_dir, "global_sampling_world_map.png"),
       plot = p_world_map, width = 18, height = 10,
       device = "png", dpi = 600, bg = "white")

message("Saved: global_sampling_world_map.pdf / .png")
