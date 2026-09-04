suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(countrycode)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(ggtext)
})

sf_use_s2(FALSE)

# ---- Locate script directory & input/output folders (mirrors figure1_composite_all_in_one.R) ----
find_script_dir <- function() {
  if (interactive()) return(normalizePath(getwd()))
  args <- commandArgs(trailingOnly = FALSE)
  f <- args[grep("^--file=", args)]
  if (length(f) == 1) return(dirname(normalizePath(sub("^--file=", "", f))))
  normalizePath(getwd())
}
script_dir <- find_script_dir()
SERVER_ROOT <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-map-2/04-patch-pipeline-map/03-patch-pipeline-map"
LOCAL_BASE  <- normalizePath(file.path(script_dir, ".."))
on_server   <- dir.exists(file.path(SERVER_ROOT, "input"))
BASE_DIR  <- if (on_server) SERVER_ROOT else LOCAL_BASE
BASE_DIR  <- Sys.getenv("FIG1_BASE_DIR",  unset = BASE_DIR)
metadata_path <- Sys.getenv("FIG1_METADATA",
                            unset = file.path(BASE_DIR, "input", "metadata_oral_qc.tsv"))
PANEL_DIR <- if (on_server) file.path(SERVER_ROOT, "overview_figure_1", "results") else file.path(script_dir, "results")
PANEL_DIR <- Sys.getenv("FIG1_PANEL_DIR", unset = PANEL_DIR)
if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)

if (!file.exists(metadata_path))
  stop("metadata table not found: ", metadata_path)

# base_sz = master font-size knob for panel C.
base_sz <- as.numeric(Sys.getenv("FIG1_C_BASE_SZ", unset = "22"))
text_sz <- base_sz / .pt

theme_map <- theme_void(base_size = base_sz, base_family = "sans") +
  theme(
    legend.title   = element_text(size = base_sz * 0.8, face = "bold", color = "#303030"),
    legend.text    = element_text(size = base_sz * 0.7, color = "#303030"),
    legend.position = "top",
    legend.box      = "horizontal",
    legend.spacing.x = unit(3, "mm"),
    legend.key      = element_blank(),
    plot.background  = element_rect(fill = "#F8F5F0", color = NA),
    panel.background = element_rect(fill = "#F8F5F0", color = NA),
    plot.margin      = margin(5, 5, 5, 5, "pt")
  )

continent_colors <- c(
  "NAM" = "#FDD835", "SA" = "#AB47BC", "EU" = "#4CAF50",
  "AS"  = "#FF9800", "AF" = "#795548", "OA" = "#F48FB1"
)
continent_order <- c("NAM", "SA", "EU", "AS", "AF", "OA")

site_cols_ordered <- c("Saliva", "Tooth", "Oral Mucosa", "Oral Swab", "Tongue", "Gingiva", "Others")
site_colors <- c(
  "Saliva"      = "#08306B",
  "Tooth"       = "#E7298A",
  "Oral Mucosa" = "#666666",
  "Oral Swab"   = "#E6AB02",
  "Tongue"      = "#A6761D",
  "Gingiva"     = "#CAB2D6",
  "Others"      = "#117733"
)

# Read and clean metadata
metadata <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  rename(
    BioProject      = bioproject_id,
    Country_Cleaned = country_cleaned,
    Site_Cleaned    = site_cleaned
  ) %>%
  mutate(
    Site_Cleaned = str_to_title(Site_Cleaned),
    Site_Cleaned = case_when(Site_Cleaned == "Other" ~ "Others", TRUE ~ Site_Cleaned)
  ) %>%
  filter(!(Site_Cleaned %in% c("Nasal", "Tonsils")), !is.na(Country_Cleaned),
         Country_Cleaned != "") %>%
  mutate(iso3c = countrycode(Country_Cleaned, origin = "country.name", destination = "iso3c"))

south_america_iso3c <- c("ARG", "BOL", "BRA", "CHL", "COL", "ECU", "GUY", "PRY",
                         "PER", "SUR", "URY", "VEN", "FLK", "GUF")

meta_clean <- metadata %>%
  mutate(
    Temp_Continent = countrycode(iso3c, "iso3c", "continent"),
    Continent = case_when(
      Temp_Continent == "Americas" & (iso3c %in% south_america_iso3c) ~ "SA",
      Temp_Continent == "Americas" & !(iso3c %in% south_america_iso3c) ~ "NAM",
      Temp_Continent == "Europe"   ~ "EU",
      Temp_Continent == "Asia"     ~ "AS",
      Temp_Continent == "Africa"   ~ "AF",
      Temp_Continent == "Oceania"  ~ "OA",
      iso3c == "ATA"               ~ "Antarctica",
      TRUE                         ~ Temp_Continent
    ),
    age_num = suppressWarnings(as.numeric(age_cleaned)),
    bmi_num = suppressWarnings(as.numeric(bmi_cleaned)),
    sex_label = case_when(
      tolower(sex_cleaned) %in% c("f", "female") ~ "F",
      tolower(sex_cleaned) %in% c("m", "male")   ~ "M",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!Continent %in% c("Other", "Antarctica")) %>%
  mutate(Continent = factor(Continent, levels = continent_order))

# Continent-level summary
continent_summary <- meta_clean %>%
  group_by(Continent) %>%
  summarise(
    n_countries = n_distinct(Country_Cleaned),
    n_bp   = n_distinct(BioProject),
    n_subj = n_distinct(subject_id),
    n_samp = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Full_Name = case_when(
      Continent == "NAM" ~ "North America (NA)",
      Continent == "SA"  ~ "South America (SA)",
      Continent == "EU"  ~ "Europe (EU)",
      Continent == "AS"  ~ "Asia (AS)",
      Continent == "AF"  ~ "Africa (AF)",
      Continent == "OA"  ~ "Oceania (OA)"
    ),
    card_label = paste0(
      "<span style='font-size: ", round(base_sz * 0.70), "pt;'><b>", Full_Name, "</b></span><br>",
      "<span style='font-size: ", round(base_sz * 0.52), "pt;'>",
      "Countries = ", n_countries, "<br>",
      "Bioprojects = ", n_bp, "<br>",
      "Subjects = ", scales::comma(n_subj), "<br>",
      "Samples = ", scales::comma(n_samp),
      "</span>"
    )
  )

# Absolute safe-zone coordinates for country info boxes
continent_coords <- data.frame(
  Continent = c("NAM", "SA", "EU", "AS", "AF", "OA"),
  lon_pie   = c(-175, -135,  -40,  145,   -40,   75),
  lat_pie   = c(  30,  -30,   55,   10,  -40,  -35)
)
continent_summary <- continent_summary %>%
  left_join(continent_coords, by = "Continent") %>%
  mutate(
    lon_lab = if_else(
      Continent == "EU",
      lon_pie[Continent == "AF"] + 9.5,
      lon_pie + 9.5
    )
  )

# Country-level summary (BioProject count per country, for gray points)
country_summary <- meta_clean %>%
  group_by(Country_Cleaned, Continent) %>%
  summarise(
    n_samples = n(),
    n_studies = n_distinct(BioProject),
    .groups = "drop"
  ) %>%
  mutate(iso3c = countrycode(Country_Cleaned, "country.name", "iso3c"))

world_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") %>%
  filter(name != "Antarctica") %>%
  mutate(iso3c = iso_a3) %>%
  mutate(
    Temp_Continent = countrycode(iso3c, "iso3c", "continent"),
    Continent = case_when(
      Temp_Continent == "Americas" & (iso3c %in% south_america_iso3c) ~ "SA",
      Temp_Continent == "Americas" & !(iso3c %in% south_america_iso3c) ~ "NAM",
      Temp_Continent == "Europe"   ~ "EU",
      Temp_Continent == "Asia"     ~ "AS",
      Temp_Continent == "Africa"   ~ "AF",
      Temp_Continent == "Oceania"  ~ "OA",
      name == "Somaliland"         ~ "AF",
      name == "N. Cyprus"          ~ "AS",
      name == "Kosovo"             ~ "EU",
      TRUE                         ~ NA_character_
    )
  ) %>%
  select(-Temp_Continent)

# Randomly distributed gray points: one per BioProject within each country
set.seed(42)
bp_points <- lapply(seq_len(nrow(country_summary)), function(i) {
  cs <- country_summary[i, ]
  idx <- which(world_sf$iso3c == cs$iso3c)
  if (length(idx) == 0) return(NULL)
  pts <- tryCatch(st_sample(world_sf$geometry[idx], size = cs$n_studies), error = function(e) NULL)
  if (is.null(pts) || length(pts) == 0) {
    if (cs$iso3c == "BRB") {
      return(data.frame(
        iso3c = "BRB", lon = rep(-59.5, cs$n_studies),
        lat = rep(13.2, cs$n_studies), n_studies = cs$n_studies
      ))
    }
    return(NULL)
  }
  coords <- st_coordinates(pts)
  data.frame(
    iso3c = cs$iso3c,
    lon   = coords[, 1],
    lat   = coords[, 2],
    n_studies = cs$n_studies
  )
}) %>%
  bind_rows()

# Pie charts per continent as grobs
site_continent_data <- meta_clean %>%
  filter(!is.na(Site_Cleaned)) %>%
  count(Continent, Site_Cleaned) %>%
  group_by(Continent) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(Site_Cleaned = factor(Site_Cleaned, levels = site_cols_ordered))

pie_grobs <- lapply(continent_order, function(cont) {
  df <- site_continent_data %>% filter(Continent == cont)
  if (nrow(df) == 0) return(grid::nullGrob())
  p <- ggplot(df, aes(x = "", y = prop, fill = Site_Cleaned)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = site_colors, limits = site_cols_ordered) +
    theme_void() +
    theme(legend.position = "none", plot.background = element_blank(),
          panel.background = element_blank())
  ggplotGrob(p)
})

r_pie <- 7.5
pie_layers <- lapply(seq_along(continent_order), function(i) {
  cont <- continent_order[i]
  coords <- continent_coords %>% filter(Continent == cont)
  annotation_custom(
    grob = pie_grobs[[i]],
    xmin = coords$lon_pie - r_pie, xmax = coords$lon_pie + r_pie,
    ymin = coords$lat_pie - r_pie, ymax = coords$lat_pie + r_pie
  )
})

# Core map rendering
p_map_base <- ggplot() +
  geom_sf(data = world_sf, aes(fill = Continent), color = "#D5D0C5", linewidth = 0.3) +
  scale_fill_manual(values = scales::alpha(continent_colors, 0.6),
                    na.value = "#EFEBDD", guide = "none") +
  geom_point(data = bp_points, aes(x = lon, y = lat, size = n_studies),
             shape = 16, color = "#0D47A1", alpha = 1) +
  scale_size_continuous(
    name   = "Number of samples",
    range  = c(2, 6),
    limits = c(0, 100),
    breaks = c(25, 50, 100),
    labels = c("25", "50", "100"),
    guide  = guide_legend(override.aes = list(color = "#0D47A1", alpha = 1,
                                              fill = "transparent", shape = 16))
  ) +
  pie_layers +
  geom_richtext(
    data = continent_summary,
    aes(x = lon_lab, y = lat_pie, label = card_label),
    hjust = 0, vjust = 0.5,
    fill = "#F8F5F0", color = "#303030",
    label.color = "grey75",
    label.padding = unit(0.4, "lines"),
    label.r = unit(0.15, "lines"),
    lineheight = 1.3,
    alpha = 0.95
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-55, 75), expand = FALSE, clip = "off") +
  theme_map

# Transparent legend
p_site_legend <- ggplot(data.frame(Site = factor(site_cols_ordered, levels = site_cols_ordered)),
                        aes(x = 1, y = 1, fill = Site)) +
  geom_point(shape = 21, size = 0, color = NA, stroke = 0) +
  scale_fill_manual(
    name = "Sampling Site",
    values = site_colors,
    limits = site_cols_ordered,
    guide = guide_legend(override.aes = list(size = 5, color = "white",
                                             stroke = 0.5, alpha = 1))
  ) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(size = base_sz * 0.85, face = "bold", color = "#303030"),
    legend.text = element_text(size = base_sz * 0.75, color = "#303030"),
    legend.key.size = unit(5, "mm"),
    legend.background = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# Composite map
p_final_map <- p_map_base +
  inset_element(
    p_site_legend,
    left = 0.01, bottom = 0.05, right = 0.20, top = 0.35,
    align_to = "full", clip = FALSE
  ) +
  plot_annotation(
    theme = theme(
      plot.background  = element_rect(fill = "#F8F5F0", color = NA),
      panel.background = element_rect(fill = NA, color = NA)
    )
  )

ggsave(file.path(PANEL_DIR, "global_sampling_world_map_final_cropped_260904.png"),
       plot = p_final_map, width = 18, height = 7.2,
       device = "png", dpi = 600, bg = "#F8F5F0")
message("Saved: ", file.path(PANEL_DIR, "global_sampling_world_map_final_cropped_260904.png"))
