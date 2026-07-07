library(tidyverse)
library(countrycode)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(scatterpie)
library(ggrepel)
library(ggsci)

# 数据加载与清洗
metadata_path <- "/home/hexincheng/projects/hxc/oral_project/oral_plot_by_yl_20260602/oral_code_results_new/input/oral_input/metadata_oral_qc_dedup.tsv"

metadata <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  rename(
    BioProject = bioproject_id,
    Country_Cleaned = country_cleaned,
    Site_Cleaned = site_cleaned
  ) %>%
  mutate(
    Site_Cleaned = case_when(
      Site_Cleaned == "saliva" ~ "Saliva",
      Site_Cleaned == "oral mucosa" ~ "Oral Mucosa",
      Site_Cleaned == "tooth" ~ "Tooth",
      Site_Cleaned == "tongue" ~ "Tongue",
      Site_Cleaned == "oral swab" ~ "Oral_swab",
      Site_Cleaned == "other" ~ "Others",
      Site_Cleaned == "gingiva" ~ "Gingiva",
      Site_Cleaned == "nasal" ~ "Nasal",
      Site_Cleaned == "tonsils" ~ "Tonsils",
      TRUE ~ Site_Cleaned
    )
  ) %>%
  filter(!(Site_Cleaned %in% c("Nasal", "Tonsils"))) %>%
  filter(!is.na(Country_Cleaned), Country_Cleaned != "")

metadata <- metadata %>%
  mutate(Country_Cleaned = case_when(
    Country_Cleaned == "Norwegian" ~ "Norway",
    Country_Cleaned %in% c("Hong Kong", "Macau", "Taiwan") ~ "China",
    TRUE ~ Country_Cleaned
  )) %>%
  mutate(iso3c = countrycode(Country_Cleaned, origin = "country.name", destination = "iso3c"))

south_america_iso3c <- c("ARG", "BOL", "BRA", "CHL", "COL", "ECU", "GUY", "PRY", "PER", "SUR", "URY", "VEN", "FLK", "GUF")

metadata <- metadata %>%
  mutate(
    Temp_Continent = countrycode(iso3c, origin = "iso3c", destination = "continent")
  ) %>%
  mutate(Continent = case_when(
    Temp_Continent == "Americas" & (iso3c %in% south_america_iso3c) ~ "South America",
    Temp_Continent == "Americas" & !(iso3c %in% south_america_iso3c) ~ "North America",
    Temp_Continent == "Oceania" ~ "Oceania",
    Temp_Continent == "Asia" ~ "Asia",
    Temp_Continent == "Europe" ~ "Europe",
    Temp_Continent == "Africa" ~ "Africa",
    iso3c == "ATA" ~ "Antarctica",
    TRUE ~ "Other"
  )) %>%
  select(-Temp_Continent)

other_check_fixed <- metadata %>%
  filter(Continent == "Other") %>%
  count(Country_Cleaned, sort = TRUE)

print(other_check_fixed)

# 01 采样分布世界地图
world_map <- map_data("world") %>%
  mutate(
    country_for_code = case_when(
      region == "USA" ~ "United States",
      region == "UK" ~ "United Kingdom",
      region == "South Korea" ~ "Republic of Korea",
      region == "North Korea" ~ "Democratic People's Republic of Korea",
      TRUE ~ region
    ),
    iso3c = countrycode(country_for_code, origin = "country.name", destination = "iso3c"),
    Temp_Continent = countrycode(iso3c, origin = "iso3c", destination = "continent"),
    Continent = case_when(
      Temp_Continent == "Americas" & (iso3c %in% south_america_iso3c) ~ "South America",
      Temp_Continent == "Americas" & !(iso3c %in% south_america_iso3c) ~ "North America",
      Temp_Continent == "Oceania" ~ "Oceania",
      Temp_Continent == "Asia" ~ "Asia",
      Temp_Continent == "Europe" ~ "Europe",
      Temp_Continent == "Africa" ~ "Africa",
      TRUE ~ "Other"
    )
  )

world_filtered <- world_map %>%
  filter(Continent != "Other")

region_coords <- data.frame(
  Continent = c("Africa", "North America", "South America", "Asia", "Europe", "Oceania"),
  label_lon = c(20, -100, -60, 85, 15, 135),
  label_lat = c(0, 45, -15, 30, 50, -25),
  text_lon = c(-15, -130, -90, 85, -20, 150),
  text_lat = c(-20, 45, -15, 55, 65, -45),
  hjust_val = c(1, 1, 1, 0.5, 1, 0.5)
)

region_labels <- metadata %>%
  filter(Continent != "Other", Continent != "Antarctica") %>%
  group_by(Continent) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  mutate(label_text = paste0(Continent, "\n(n=", scales::comma(n_samples), ")")) %>%
  left_join(region_coords, by = "Continent")

country_centers <- world_filtered %>%
  group_by(region, Continent) %>%
  summarise(
    center_lon = (min(long, na.rm = TRUE) + max(long, na.rm = TRUE)) / 2,
    center_lat = (min(lat, na.rm = TRUE) + max(lat, na.rm = TRUE)) / 2,
    .groups = "drop"
  )

bioproject_counts <- metadata %>%
  filter(Continent != "Other", Continent != "Antarctica", !is.na(BioProject)) %>%
  distinct(BioProject, Country_Cleaned, Continent) %>%
  mutate(
    region = case_when(
      Country_Cleaned == "United States" ~ "USA",
      Country_Cleaned == "United Kingdom" ~ "UK",
      Country_Cleaned == "Republic of Korea" ~ "South Korea",
      TRUE ~ Country_Cleaned
    )
  ) %>%
  count(region, Continent, name = "bp_count") %>%
  left_join(country_centers, by = c("region", "Continent")) %>%
  filter(!is.na(center_lon), !is.na(center_lat))

set.seed(123)
bioproject_points_scatter <- bioproject_counts %>%
  tidyr::uncount(bp_count) %>%
  mutate(
    lon = pmin(pmax(center_lon + stats::runif(dplyr::n(), -4, 4), -175), 175),
    lat = pmin(pmax(center_lat + stats::runif(dplyr::n(), -3, 3), -55), 80)
  )

continent_bg_colors <- c(
  "North America" = "#E64B35",
  "Oceania" = "#4DBBD5",
  "Europe" = "#4CAF50",
  "Asia" = "#F5D241",
  "South America" = "#F39B7F",
  "Africa" = "#F781BF"
)

p1 <- ggplot() +
  geom_polygon(
    data = world_filtered,
    aes(x = long, y = lat, group = group, fill = Continent),
    color = "white",
    linewidth = 0.2,
    alpha = 0.4
  ) +
  scale_fill_manual(
    name = "Region",
    values = continent_bg_colors,
    breaks = names(continent_bg_colors)
  ) +
  geom_point(
    data = bioproject_points_scatter,
    aes(x = lon, y = lat, color = "BioProject"),
    shape = 16,
    size = 2.0,
    alpha = 0.8
  ) +
  scale_color_manual(name = NULL, values = c("BioProject" = "#377EB8")) +
  geom_segment(
    data = region_labels,
    aes(x = text_lon, y = text_lat, xend = label_lon, yend = label_lat),
    color = "grey40",
    linewidth = 0.6
  ) +
  geom_text(
    data = region_labels,
    aes(x = text_lon, y = text_lat, label = label_text, hjust = hjust_val),
    size = 4.5,
    fontface = "bold",
    color = "black",
    lineheight = 1.1
  ) +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-60, 85), expand = TRUE, clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(t = 20, r = 60, b = 10, l = 60, unit = "pt"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, "cm")
  )

print(p1)

# 02 采样堆积条形图
site_cols_ordered <- c(
  "Saliva",
  "Oral Mucosa",
  "Tooth",
  "Tongue",
  "Oral_swab",
  "Others",
  "Gingiva"
)

my_custom_colors <- c(
  "#3CE390",
  "#5D75DF",
  "#F5D241",
  "#47346A",
  "#F58731",
  "#841B14",
  "#D7421D"
)

bar_data_final <- metadata %>%
  filter(Continent != "Other", Continent != "Antarctica") %>%
  filter(!is.na(Site_Cleaned)) %>%
  group_by(Continent, Site_Cleaned) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(Continent, desc(count)) %>%
  group_by(Continent) %>%
  mutate(order_id = row_number()) %>%
  ungroup()

p2 <- ggplot(bar_data_final, aes(x = reorder(Continent, -count, sum),
                                 y = count,
                                 group = order_id,
                                 fill = Site_Cleaned)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(
    name = "Site",
    limits = site_cols_ordered,
    values = my_custom_colors
  ) +
  scale_y_continuous(
    limits = c(0, 20000),
    breaks = c(0, 5000, 10000, 15000, 20000),
    labels = c("0", "5,000", "10,000", "15,000", "20,000"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 20, l = 30, unit = "pt"),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black", size = 15),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.x = element_text(size = 15, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 15, face = "bold", margin = margin(r = 15)),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"),
    aspect.ratio = 1
  ) +
  labs(title = "", x = "Region", y = "Number of Samples")

print(p2)

# 保存
output_dir <- "/home/hexincheng/projects/hxc/oral_project/oral_plot_by_yl_20260602/oral_code_results_new/output/oral_output/sampling_distribution"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("已自动创建输出文件夹！")
}

p1_pdf_path <- file.path(output_dir, "sampling_distribution_world_map.pdf")
p1_png_path <- file.path(output_dir, "sampling_distribution_world_map.png")

ggsave(filename = p1_pdf_path, plot = p1, width = 10, height = 6, device = "pdf")
ggsave(filename = p1_png_path, plot = p1, width = 10, height = 6, device = "png", dpi = 300, bg = "white")

if (file.exists(p1_pdf_path) && file.exists(p1_png_path)) {
  message("图 1 已成功保存为 PDF 和 PNG 格式！")
} else {
  warning("图 1 保存失败，请检查绘图依赖环境。")
}

p2_pdf_path <- file.path(output_dir, "sampling_distribution_stacked_bar.pdf")
p2_png_path <- file.path(output_dir, "sampling_distribution_stacked_bar.png")

ggsave(filename = p2_pdf_path, plot = p2, width = 8, height = 8, device = "pdf")
ggsave(filename = p2_png_path, plot = p2, width = 8, height = 8, device = "png", dpi = 300, bg = "white")

message("图 2 已成功保存为 PDF 和 PNG 格式！")
