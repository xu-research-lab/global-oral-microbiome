#!/usr/bin/env Rscript
# Figure 1 v2: Composite figure for global oral microbiome dataset overview
# Changes: phylum stacked bar chart (Panel c) is now embedded as an inset into the left side (Pacific ocean) of the map (Panel b).

library(tidyverse)
library(mgcv)
library(countrycode)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(patchwork)
library(cowplot)
library(gridExtra)
library(data.table)
library(scales)

sf_use_s2(FALSE)

# GLOBAL THEME & FONT SIZING
base_sz <- 17
text_sz <- base_sz / .pt

update_geom_defaults("text", list(family = "sans"))
update_geom_defaults("label", list(family = "sans"))

tag_style <- theme(
  plot.tag = element_text(size = base_sz * 1.5, face = "bold", color = "#303030", hjust = 0, vjust = 1)
)

theme_set(theme_minimal(base_size = base_sz, base_family = "sans"))

# COLOR PALETTE
continent_colors <- c(
  "NAM" = "#D32F2F", "SA" = "#AB47BC", "EU" = "#4CAF50",
  "AS" = "#FF9800", "AF" = "#795548", "OA" = "#29B6F6"
)
continent_order <- c("NAM", "SA", "EU", "AS", "AF", "OA")

site_cols_ordered <- c("Saliva", "Tooth", "Oral Mucosa", "Oral Swab", "Tongue",  "Gingiva", "Others")
site_colors <- c(
  "Saliva"      = "#E31A1C",
  "Tooth"       = "#FFB300",
  "Oral Mucosa" = "#1565C0",
  "Oral Swab"   = "#8E24AA",
  "Tongue"      = "#33A02C",
  "Gingiva"     = "#B15928",
  "Others"      = "#A6CEE3"
)

# =============================================================================
# Panel a: 2-Phase Pipeline funnel
# =============================================================================
funnel_data <- data.frame(
  step = 1:7,
  n_inside = c("116,643", "86,255", "84,564", "82,336", "82,336", "81,227", "73,340"),
  s_inside = c("334 studies", "334 studies", "312 studies", "293 studies", "293 studies", "293 studies", "293 studies"),
  n_num = c(116643, 86255, 84564, 82336, 82336, 81227, 73340),
  y = c(17.5, 15.5, 13.5, 11.5, 4.5, 2.5, 0.5)
)

right_anno <- data.frame(
  y_mid = c(16.5, 14.5, 12.5, 3.5, 1.5),
  line1 = c("Step 1: Dataset Metadata Filter", "Step 2: Dataset Chimera Filter", "Step 3: Dataset ASV & CR Filter",
            "Step 4: Sample Chimera Filter", "Step 5: Sample ASV & CR Filter"),
  line2 = c("(remove non-oral samples)", "(datasets >30% oral samples with chimera >25%)", "(datasets >50% oral samples with ASV ≤30 and CR <30%)",
            "(samples with chimera >25%)", "(samples with ASV ≤30 and CR <30%)"),
  line3 = c("−30,388 samples", "−22 studies", "−19 studies",
            "−1,109 samples", "−7,887 samples")
)

max_w <- 13.5
min_w <- 8.5
bar_h <- 0.60

funnel_data$half_w <- min_w + (funnel_data$n_num - min(funnel_data$n_num)) / (max(funnel_data$n_num - min(funnel_data$n_num))) * (max_w - min_w)

poly_blocks <- do.call(rbind, lapply(1:nrow(funnel_data), function(i) {
  b <- funnel_data[i, ]; hw <- b$half_w; yc <- b$y
  data.frame(x = c(-hw, hw, hw, -hw), y = c(yc - bar_h, yc - bar_h, yc + bar_h, yc + bar_h), step = b$step, grp = paste0("bar_", i))
}))

conn_indices <- c(1, 2, 3, 5, 6)
conn_blocks <- do.call(rbind, lapply(conn_indices, function(i) {
  top_hw <- funnel_data$half_w[i]; bot_hw <- funnel_data$half_w[i + 1]
  top_y_bottom <- funnel_data$y[i] - bar_h; bot_y_top <- funnel_data$y[i + 1] + bar_h
  data.frame(x = c(-top_hw, top_hw, bot_hw, -bot_hw), y = c(top_y_bottom, top_y_bottom, bot_y_top, bot_y_top), step = i, grp = paste0("conn_", i))
}))

arrow_data <- data.frame(
  x = 0,
  y_start = funnel_data$y[conn_indices] - bar_h,
  y_end = funnel_data$y[conn_indices + 1] + bar_h
)

bar_fills <- c("#263238", "#455A64", "#78909C", "#90A4AE",  "#1976D2", "#43A047", "#F57C00")
bar_text_colors <- rep("white", 7)

anno_x <- max_w + 2.0
anno_w <- 75.0
anno_h <- 1.5

right_anno$rect_xmin <- anno_x - 0.2
right_anno$rect_xmax <- anno_x + anno_w + 0.5
right_anno$rect_ymin <- right_anno$y_mid - anno_h / 2
right_anno$rect_ymax <- right_anno$y_mid + anno_h / 2

dataset_rect_xmax <- anno_x + anno_w + 0.5
dataset_rect_ymax <- 19.5
sample_rect_xmax <- anno_x + anno_w + 0.5
sample_rect_ymax <- 7.25

p_a <- ggplot() +
  annotate("rect", xmin = -max_w - 0.5, xmax = anno_x + anno_w + 0.5,
           ymin = 10.25, ymax = 19.5, fill = "#F8F9FA", color = "#B0BEC5", linetype = "dashed", alpha = 0.5) +
  annotate("text", x = dataset_rect_xmax - 1.0, y = dataset_rect_ymax - 0.4,
           label = "Dataset-Level Filters", fontface = "bold", hjust = 1, size = text_sz * 1.1, color = "#37474F") +
  annotate("rect", xmin = -max_w - 0.5, xmax = anno_x + anno_w + 0.5,
           ymin = -1.5, ymax = 7.25, fill = "#F3E5F5", color = "#CE93D8", linetype = "dashed", alpha = 0.3) +
  annotate("text", x = sample_rect_xmax - 1.0, y = sample_rect_ymax - 0.4,
           label = "Sample-Level Filters", fontface = "bold", hjust = 1, size = text_sz * 1.1, color = "#6A1B9A") +
  geom_polygon(data = conn_blocks, aes(x = x, y = y, group = grp), fill = "#CFD8DC", alpha = 0.45) +
  geom_polygon(data = poly_blocks, aes(x = x, y = y, group = grp), fill = rep(bar_fills, each = 4), color = NA) +
  geom_segment(data = funnel_data, aes(x = -half_w, xend = half_w, y = y + bar_h, yend = y + bar_h), color = "white", linewidth = 0.4) +
  geom_segment(data = funnel_data, aes(x = -half_w, xend = half_w, y = y - bar_h, yend = y - bar_h), color = "white", linewidth = 0.4) +
  geom_segment(data = arrow_data, aes(x = x, y = y_start, xend = x, yend = y_end), arrow = arrow(length = unit(0.3, "cm"), type = "closed"), linewidth = 0.8, color = "#546E7A") +
  annotate("segment", x = 0, xend = 0, y = 9.95, yend = 7.55,
           arrow = arrow(length = unit(0.4, "cm"), type = "closed"), linewidth = 2.0, color = "#37474F") +
  geom_text(data = funnel_data, aes(x = 0, y = y + 0.25, label = n_inside), size = text_sz * 1.0, color = bar_text_colors, fontface = "bold", hjust = 0.5) +
  geom_text(data = funnel_data, aes(x = 0, y = y - 0.30, label = s_inside), size = text_sz * 0.75, color = bar_text_colors, fontface = "plain", hjust = 0.5) +
  geom_rect(data = right_anno,
            aes(xmin = rect_xmin, xmax = rect_xmax, ymin = rect_ymin, ymax = rect_ymax),
            fill = "#FFFFFF", color = "#B0BEC5", linewidth = 0.4, alpha = 0.95) +
  geom_text(data = right_anno, aes(x = anno_x + 0.4, y = y_mid + 0.45, label = line1),
            size = text_sz * 0.85, fontface = "bold", color = "#212121", hjust = 0) +
  geom_text(data = right_anno, aes(x = anno_x + 0.4, y = y_mid, label = line2),
            size = text_sz * 0.65, fontface = "italic", color = "#616161", hjust = 0) +
  geom_text(data = right_anno, aes(x = anno_x + 0.4, y = y_mid - 0.45, label = line3),
            size = text_sz * 0.80, fontface = "bold", color = "#D32F2F", hjust = 0) +
  geom_segment(data = right_anno, aes(x = max_w, xend = anno_x - 0.15, y = y_mid, yend = y_mid), color = "#78909C", linewidth = 0.5, alpha = 0.8, linetype = "dashed") +
  coord_cartesian(xlim = c(-14.5, 95.0), ylim = c(-0.5, 20.5), clip = "off") +
  theme_void(base_size = base_sz, base_family = "sans") +
  theme(plot.margin = margin(10, 5, 10, 5, unit = "pt")) + labs(tag = "a")

# =============================================================================
# Panel b Prep: Data Extraction
# =============================================================================
metadata_path <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/01-left-pipeline-right-map/input/metadata_oral_qc_update.tsv"
metadata <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  rename(BioProject = bioproject_id, Country_Cleaned = country_cleaned, Site_Cleaned = site_cleaned) %>%
  mutate(Site_Cleaned = str_to_title(Site_Cleaned), Site_Cleaned = case_when(Site_Cleaned == "Other" ~ "Others", TRUE ~ Site_Cleaned)) %>%
  filter(!(Site_Cleaned %in% c("Nasal", "Tonsils")), !is.na(Country_Cleaned), Country_Cleaned != "") %>%
  mutate(iso3c = countrycode(Country_Cleaned, origin = "country.name", destination = "iso3c"))

south_america_iso3c <- c("ARG","BOL","BRA","CHL","COL","ECU","GUY","PRY","PER","SUR","URY","VEN","FLK","GUF")
meta_clean <- metadata %>%
  mutate(Temp_Continent = countrycode(iso3c, "iso3c", "continent"),
         Continent = case_when(
           Temp_Continent == "Americas" & (iso3c %in% south_america_iso3c) ~ "SA",
           Temp_Continent == "Americas" & !(iso3c %in% south_america_iso3c) ~ "NAM",
           Temp_Continent == "Europe"  ~ "EU",
           Temp_Continent == "Asia"    ~ "AS",
           Temp_Continent == "Africa"  ~ "AF",
           Temp_Continent == "Oceania" ~ "OA",
           iso3c == "ATA" ~ "Antarctica",
           TRUE ~ Temp_Continent
         )) %>%
  filter(!Continent %in% c("Other", "Antarctica")) %>%
  mutate(Continent = factor(Continent, levels = continent_order))

world_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") %>% filter(name != "Antarctica")

country_summary <- meta_clean %>%
  group_by(Country_Cleaned, Continent) %>%
  summarise(n_samples = n(), n_studies = n_distinct(BioProject), .groups = "drop") %>%
  mutate(iso3c = countrycode(Country_Cleaned, "country.name", "iso3c"))

centroids <- world_sf %>% mutate(iso3c = iso_a3) %>% filter(!is.na(iso3c)) %>% st_point_on_surface() %>%
  mutate(lon = st_coordinates(geometry)[,1], lat = st_coordinates(geometry)[,2]) %>% st_drop_geometry() %>% select(iso3c, lon, lat) %>% filter(!is.na(iso3c))

country_map_data <- country_summary %>% left_join(centroids, by = "iso3c") %>%
  mutate(lon = ifelse(is.na(lon) & iso3c == "BRB", -59.5, lon), lat = ifelse(is.na(lat) & iso3c == "BRB", 13.2, lat)) %>%
  filter(!is.na(lon), !is.na(lat)) %>% mutate(Continent = factor(Continent, levels = continent_order))

country_map_labels <- country_map_data %>%
  filter(n_samples >= 100) %>%
  mutate(
    country_label = paste0(Country_Cleaned, "\n(", n_studies, " bp; n=", scales::comma(n_samples), ")")
  )

col_names <- colnames(metadata); col_lower <- tolower(col_names)
age_col <- if ("age_cleaned" %in% col_lower) col_names[which(col_lower == "age_cleaned")] else if ("age" %in% col_lower) col_names[which(col_lower == "age")] else NULL
sex_col <- if ("sex_cleaned" %in% col_lower) col_names[which(col_lower == "sex_cleaned")] else if ("sex" %in% col_lower) col_names[which(col_lower == "sex")] else NULL
bmi_col <- if ("bmi_cleaned" %in% col_lower) col_names[which(col_lower == "bmi_cleaned")] else if ("bmi" %in% col_lower) col_names[which(col_lower == "bmi")] else NULL

demo_df <- meta_clean
if (!is.null(age_col)) demo_df[[age_col]] <- as.numeric(demo_df[[age_col]])
if (!is.null(bmi_col)) demo_df[[bmi_col]] <- as.numeric(demo_df[[bmi_col]])

tbl_base <- demo_df %>% group_by(Continent) %>%
  summarise(Samples = n(), Countries = n_distinct(Country_Cleaned), Studies = n_distinct(BioProject), .groups = "drop")
if (!is.null(age_col)) {
  age_summ <- demo_df %>% filter(!is.na(!!sym(age_col))) %>% group_by(Continent) %>%
    summarise(Age = paste0(sprintf("%.1f", mean(!!sym(age_col), na.rm=T)), "±", sprintf("%.1f", sd(!!sym(age_col), na.rm=T))))
  tbl_base <- tbl_base %>% left_join(age_summ, by="Continent")
}
if (!is.null(bmi_col)) {
  bmi_summ <- demo_df %>% filter(!is.na(!!sym(bmi_col)), !!sym(bmi_col)>10, !!sym(bmi_col)<60) %>% group_by(Continent) %>%
    summarise(BMI = paste0(sprintf("%.1f", mean(!!sym(bmi_col), na.rm=T)), "±", sprintf("%.1f", sd(!!sym(bmi_col), na.rm=T))))
  tbl_base <- tbl_base %>% left_join(bmi_summ, by="Continent")
}
if (!is.null(sex_col)) {
  sex_summ <- demo_df %>% filter(!is.na(!!sym(sex_col))) %>%
    mutate(sex_label = case_when(tolower(!!sym(sex_col)) %in% c("male","m") ~ "M", tolower(!!sym(sex_col)) %in% c("female","f") ~ "F", TRUE ~ "O")) %>%
    filter(sex_label %in% c("M","F")) %>% count(Continent, sex_label) %>% group_by(Continent) %>% mutate(prop = round(n/sum(n)*100,0)) %>%
    pivot_wider(id_cols = Continent, names_from = sex_label, values_from = prop, values_fill = 0) %>%
    mutate(`Sex (F/M)` = paste0(F, "% / ", M, "%")) %>% select(Continent, `Sex (F/M)`)
  tbl_base <- tbl_base %>% left_join(sex_summ, by="Continent")
}

tbl_final <- data.frame(Continent = factor(continent_order, levels = continent_order)) %>%
  left_join(tbl_base, by="Continent") %>%
  arrange(desc(Samples)) %>%
  mutate(
    Samples = format(Samples, big.mark=","),
    across(everything(), as.character),
    across(everything(), ~replace_na(.x, "-")),
    Continent = ifelse(Continent == "NAM", "NA", as.character(Continent))
  )

my_alignments <- c(0, 0, 0, 0, 0, 0, 0)

my_theme <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize = base_sz - 8, fontfamily = "sans", col = "#303030", hjust = my_alignments, x = my_alignments),
    bg_params = list(fill = c("#F9F9F9", "#FFFFFF"), col = "grey90")
  ),
  colhead = list(
    fg_params = list(fontsize = base_sz - 8, fontface = "bold", fontfamily = "sans", col = "#303030", hjust = my_alignments, x = my_alignments),
    bg_params = list(fill = "#ECEFF1", col = "grey80")
  ),
  padding = unit(c(2.5, 2.5), "mm")
)

summary_table_grob <- tableGrob(tbl_final, rows = NULL, theme = my_theme)

# =============================================================================
# Panel c: Phylum-level relative abundance stacked bar chart (INSET for Panel b)
# =============================================================================
metadata_raw <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  rename(BioProject = bioproject_id, SampleID = sample_id)

abundance_file <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/01-left-pipeline-right-map/input/phylum-relative-table.tsv"
rel_dt <- fread(abundance_file, header = TRUE, sep = "\t", skip = 1)
setnames(rel_dt, 1, "Phylum")

phylum_cols <- setdiff(names(rel_dt), "Phylum")
rel_long <- melt(rel_dt, id.vars = "Phylum", measure.vars = phylum_cols, 
                 variable.name = "sample_id", value.name = "Abundance")
rel_long[, sample_id := as.character(sample_id)]
rel_long[, Phylum_clean := gsub("^d__Bacteria;\\s*p__", "", Phylum)]
rel_long[, Phylum_clean := gsub("^p__", "", Phylum_clean)]
rel_long <- rel_long[!is.na(Phylum_clean) & Phylum_clean != ""]

meta_sub <- as.data.table(metadata_raw)[, .(sample_id = SampleID, site_cleaned, continent_group)]
meta_sub[, sample_id := as.character(sample_id)]

merged_dt <- merge(rel_long, meta_sub, by = "sample_id", all.x = FALSE)

target_sites_lower   <- c("saliva", "tooth", "tongue", "oral mucosa")
target_sites_display <- c("Saliva", "Tooth", "Tongue", "Oral Mucosa")
site_label_map       <- setNames(target_sites_display, target_sites_lower)
target_continents    <- c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")

merged_dt <- merged_dt[!is.na(site_cleaned) & site_cleaned != ""]
merged_dt <- merged_dt[!is.na(continent_group) & continent_group != "" & continent_group != "NA"]
merged_dt[, site_cleaned := tolower(trimws(as.character(site_cleaned)))]
merged_dt[, continent_group := trimws(as.character(continent_group))]

merged_dt <- merged_dt[site_cleaned %in% target_sites_lower]
merged_dt <- merged_dt[continent_group %in% target_continents]
merged_dt[, site_display := site_label_map[site_cleaned]]
merged_dt[, site_display := factor(site_display, levels = target_sites_display)]

global_top9 <- merged_dt[, .(Total = sum(Abundance, na.rm = TRUE)), by = Phylum_clean][
                         order(-Total)][seq_len(min(9, .N))]$Phylum_clean
merged_dt[, Phylum_group := ifelse(Phylum_clean %in% global_top9, Phylum_clean, "Others")]
plot_data_c <- merged_dt[, .(Mean_Abundance = mean(Abundance, na.rm = TRUE)), 
                         by = .(site_display, continent_group, Phylum_group)]
plot_data_c[, Mean_Abundance := Mean_Abundance / sum(Mean_Abundance), by = .(site_display, continent_group)]

dominant_phylum <- global_top9[1]
combined_levels <- character(0)
for (s in target_sites_display) {
  sub_c <- plot_data_c[site_display == s & Phylum_group == dominant_phylum]
  if (nrow(sub_c) == 0) {
    ordered_conts <- intersect(target_continents, unique(plot_data_c[site_display == s]$continent_group))
  } else {
    ordered_conts <- sub_c[order(-Mean_Abundance)]$continent_group
  }
  combined_levels <- c(combined_levels, paste(s, ordered_conts, sep = "|"))
}

plot_data_c[, facet_var := factor(paste(site_display, continent_group, sep = "|"), levels = combined_levels)]
strip_facet_prefix <- function(x) {
  raw_names <- sapply(x, function(label) gsub("^[^|]*\\|", "", as.character(label)))
  abbr_map <- c("North America" = "NA", "South America" = "SA", 
                "Europe" = "EU", "Asia" = "AS", 
                "Africa" = "AF", "Oceania" = "OA")
  return(unname(abbr_map[raw_names]))
}

color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
unified_colors <- setNames(color_palette[seq_len(length(global_top9))], global_top9)
unified_colors["Others"] <- "#CCCCCC"

other_phyla <- setdiff(global_top9, dominant_phylum)
legend_order <- c("Others", rev(other_phyla), dominant_phylum)
plot_data_c[, Phylum_group := factor(Phylum_group, levels = legend_order)]

p_c <- ggplot(plot_data_c, aes(x = facet_var, y = Mean_Abundance, fill = Phylum_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.90) +
  facet_grid(. ~ site_display, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = unified_colors, name = "Phylum", breaks = legend_order, labels = legend_order) +
  scale_x_discrete(labels = strip_facet_prefix) +
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"), 
                     breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Relative Abundance (%)", tag = "c") +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) + 
  theme_classic(base_size = base_sz - 4, base_family = "sans") +
  theme(
    axis.text.x = element_text(size = base_sz - 6, face = "bold", color = "#303030", angle = 30, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = base_sz - 5, face = "bold", color = "#303030"),
    axis.title.y = element_text(size = base_sz - 4, face = "bold", color = "#303030", margin = margin(r = 4)),
    legend.position = "top",
    legend.title = element_text(size = base_sz - 6, face = "bold", color = "#303030"),
    legend.text = element_text(size = base_sz - 7, color = "#303030"),
    legend.key.size = unit(3.0, "mm"),
    legend.margin = margin(b = -5),
    strip.text = element_text(size = base_sz - 8, face = "bold", color = "#303030", margin = margin(b = 4, t = 4)),
    strip.background = element_rect(fill = "#ECEFF1", color = NA),
    panel.spacing = unit(0.5, "lines"),
    axis.line = element_line(color = "#303030", linewidth = 0.3),
    axis.ticks = element_line(color = "#303030", linewidth = 0.2),
    plot.margin = margin(2, 2, 2, 2, unit = "pt"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# =============================================================================
# Panel b: Rendering Map (with p_c as inset)
# =============================================================================
p_b_base <- ggplot() +
  geom_sf(data = world_sf, fill = "#E8E0D4", color = "white", linewidth = 0.2) +
  geom_point(data = country_map_data, aes(x = lon, y = lat, size = n_samples, fill = Continent), shape = 21, color = "white", stroke = 0.35, alpha = 0.8) +
  scale_size_continuous(name = "Number of samples",
                        range = c(1, 25),
                        limits = c(1, 50000),
                        breaks = c(10, 100, 1000, 10000, 40000),
                        labels = c("10", "100", "1K", "10K", "40K"),
                        guide = guide_legend(override.aes = list(shape = 21, fill = "grey70", color = "grey30", stroke = 0.5))) +
  scale_fill_manual(values = continent_colors, guide = "none") +
  geom_text_repel(data = country_map_labels, aes(x = lon, y = lat, label = country_label, color = Continent),
                  size = text_sz * 0.60, fontface = "bold", lineheight = 0.85, segment.color = "#999999", segment.size = 0.2, box.padding = 0.4, max.overlaps = Inf) +
  scale_color_manual(values = continent_colors, guide = "none") +
  # Table positioned below Africa and Australia
  annotation_custom(summary_table_grob,
                    xmin = -40, xmax = 120,
                    ymin = -60, ymax = -30) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  theme_void(base_size = base_sz, base_family = "sans") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = base_sz * 0.8, color = "#303030"),
    legend.title = element_text(size = base_sz * 0.9, face = "bold", color = "#303030", hjust = 0.5),
    legend.margin = margin(t = 0, b = -20),
    plot.margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")
  ) + labs(tag = "b")

p_b <- p_b_base + inset_element(p_c + tag_style, 
                                left = 0.12,
                                bottom = 0.15,  
                                right = 0.24,
                                top = 0.45)

# =============================================================================
# Panel d: Site x Continent Bubble Matrix
# =============================================================================
site_continent_data <- meta_clean %>% filter(!is.na(Site_Cleaned)) %>%
  count(Continent, Site_Cleaned, name = "count")

site_continent_data <- site_continent_data %>%
  group_by(Continent) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  mutate(
    Site_Cleaned = factor(Site_Cleaned, levels = site_cols_ordered),
    Site_X = as.numeric(Site_Cleaned) * 0.52,
    Continent_Label = ifelse(Continent == "NAM", "NA", as.character(Continent)),
    Continent_Label = factor(Continent_Label, levels = c("NA", "AS", "EU", "AF", "OA", "SA"))
  )

p_d <- ggplot(site_continent_data, aes(x = Site_X, y = Continent_Label, size = count, fill = Site_Cleaned)) +
  geom_point(shape = 21, color = "white", stroke = 0.6, alpha = 0.90) +
  scale_size_continuous(name = "Number of samples", range = c(2, 12),
                         breaks = c(50, 500, 5000, 15000),
                         labels = c("50", "0.5K", "5K", "15K"),
                         guide = guide_legend(override.aes = list(fill = "grey70", color = "grey30", stroke = 0.5),
                                              label.position = "right", title.position = "left", nrow = 1, byrow = TRUE)) +
  scale_fill_manual(name = "Site", values = site_colors, limits = site_cols_ordered) +
  scale_x_continuous(breaks = seq_along(site_cols_ordered) * 0.52,
                     labels = site_cols_ordered,
                     expand = expansion(add = c(0.16, 0.10))) +
  scale_y_discrete(limits = rev(levels(site_continent_data$Continent_Label)), expand = expansion(add = 0.45)) +
  labs(x = NULL, y = "Continent", tag = "d") +
  guides(fill = "none",
         size = guide_legend(override.aes = list(fill = "grey70", color = "grey30", stroke = 0.5),
                             label.position = "right", title.position = "left", order = 1, nrow = 1, byrow = TRUE)) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = base_sz - 2, base_family = "sans") +
  theme(
    axis.text.x = element_text(size = (base_sz - 2) * 0.85, face = "bold", color = "#303030",
                                angle = 30, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = (base_sz - 2) * 0.85, face = "bold", color = "#303030"),
    axis.title.y = element_text(size = (base_sz - 2) * 1.0, face = "bold", color = "#303030"),
    axis.line = element_line(color = "#303030", linewidth = 0.4),
    axis.ticks = element_line(color = "#303030", linewidth = 0.3),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = (base_sz - 2) * 0.68, face = "bold", color = "#303030",
                                margin = margin(r = 3, unit = "pt")),
    legend.text = element_text(size = (base_sz - 2) * 0.62, color = "#303030"),
    legend.key.size = unit(2.5, "mm"),
    legend.spacing.x = unit(2, "pt"),
    legend.spacing.y = unit(1, "mm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 18, unit = "pt")
  )

# =============================================================================
# Panel e: Saturation curve inset
# =============================================================================
rarefaction_file <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/01-left-pipeline-right-map/output/rarefaction_pcofg.tsv"

df_rare <- read.table(rarefaction_file, header = TRUE, sep = "\t") %>%
  pivot_longer(cols = c("Phylum","Class","Order","Family","Genus"), names_to = "Taxonomy", values_to = "Unique_taxa") %>%
  mutate(Taxonomy = factor(Taxonomy, levels = c("Genus","Family","Order","Class","Phylum")))

p_e <- ggplot(df_rare, aes(x = Samples, y = Unique_taxa, color = Taxonomy)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = 1.0) +
  scale_color_manual(values = c("Genus"="#5C6BC0", "Family"="#FFB300", "Order"="#EF5350", "Class"="#26C6DA", "Phylum"="#66BB6A")) +
  scale_x_continuous(breaks = seq(0, 80000, by = 20000),
                     limits = c(0, 98000),
                     labels = c("0", "20K", "40K", "60K", "80K"), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000),
                     labels = c("0", "1K", "2K", "3K", "4K", "5K"),
                     expand = expansion(mult = c(0.1, 0.06))) +
  annotate("text", x = 81000, y = 4900, label = "Genus", color = "#5C6BC0", fontface = "bold", size = base_sz/4.5, hjust = 0) +
  annotate("text", x = 81000, y = 1400, label = "Family", color = "#FFB300", fontface = "bold", size = base_sz/4.5, hjust = 0) +
  annotate("text", x = 81000, y = 800, label = "Order", color = "#EF5350", fontface = "bold", size = base_sz/4.5, hjust = 0) +
  annotate("text", x = 81000, y = 350, label = "Class", color = "#26C6DA", fontface = "bold", size = base_sz/4.5, hjust = 0) +
  annotate("text", x = 81000, y = 5, label = "Phylum", color = "#66BB6A", fontface = "bold", size = base_sz/4.5, hjust = 0) +
  labs(x = "Samples", y = "Observed taxa", tag = "e") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = base_sz - 6, base_family = "sans") +
  theme(
    axis.text.x = element_text(size = (base_sz - 2) * 0.85, face = "bold", color = "#303030"),
    axis.text.y = element_text(size = (base_sz - 2) * 0.85, face = "bold", color = "#303030"),
    axis.title.y = element_text(size = (base_sz - 2) * 1.0, face = "bold", color = "#303030", margin = margin(r = 6)),
    axis.title.x = element_text(size = (base_sz - 2) * 1.0, face = "bold", color = "#303030", margin = margin(t = 10)),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5, unit = "pt"),
    plot.background = element_blank(),
    panel.background = element_blank()
  )

# =============================================================================
# Panel f: Health/Disease Stacked Bar + Bubble Overlay
# =============================================================================
clean_status_value <- function(x) {
  x <- str_squish(as.character(x))
  na_if(x, "")
}

df_status <- metadata_raw %>%
  mutate(
    status_cleaned = clean_status_value(status_cleaned),
    status_cleaned = if_else(str_to_lower(status_cleaned) == "hc", "HC", status_cleaned),
    disease_study_cleaned = str_to_title(clean_status_value(disease_study_cleaned))
  ) %>%
  filter(!is.na(SampleID), !is.na(status_cleaned), !is.na(disease_study_cleaned),
         !str_to_lower(status_cleaned) %in% c("na", "n/a", "nan", "null", "none", "unknown", "unk"),
         !str_to_lower(disease_study_cleaned) %in% c("na", "n/a", "nan", "null", "none", "unknown", "unk"),
         disease_study_cleaned %in% c("Yes", "No"))

df_disease <- df_status %>% filter(disease_study_cleaned == "Yes")

project_disease <- df_disease %>%
  filter(status_cleaned != "HC") %>%
  group_by(BioProject, status_cleaned) %>%
  tally() %>%
  arrange(BioProject, desc(n)) %>%
  slice(1) %>%
  select(BioProject, Target_Disease = status_cleaned)

hc_only_in_disease <- df_disease %>%
  filter(!(BioProject %in% project_disease$BioProject)) %>%
  distinct(BioProject) %>%
  mutate(Target_Disease = "Unclassified disease")

project_disease <- bind_rows(project_disease, hc_only_in_disease)

df_disease_annotated <- df_disease %>%
  inner_join(project_disease, by = "BioProject") %>%
  mutate(Sample_Group = if_else(status_cleaned == "HC", "HC", "Disease"))

df_plot_e <- df_disease_annotated %>%
  mutate(Sample_Group = factor(Sample_Group, levels = c("Disease", "HC")))

sample_counts_e <- df_plot_e %>%
  group_by(Target_Disease, Sample_Group) %>%
  summarise(Count = n_distinct(SampleID), .groups = "drop")

project_counts_e <- df_plot_e %>%
  group_by(Target_Disease) %>%
  summarise(
    Project_Num = n_distinct(BioProject),
    Total_Samples = n_distinct(SampleID),
    .groups = "drop"
  )

top_diseases <- project_counts_e %>%
  filter(!str_to_lower(Target_Disease) %in% c("unclassified disease", "other", "others")) %>%
  arrange(desc(Total_Samples), Target_Disease) %>%
  slice_head(n = 15) %>%
  pull(Target_Disease)

disease_levels <- top_diseases
disease_label_map <- c(
  "Gastroesophageal Reflux Disease" = "GERD",
  "Acute Myeloid Leukemia" = "AML",
  "Head and Neck Cancer" = "HNC"
)
disease_axis_labels <- setNames(
  if_else(disease_levels %in% names(disease_label_map), disease_label_map[disease_levels], disease_levels),
  disease_levels
)

sample_counts_e <- sample_counts_e %>%
  filter(Target_Disease %in% disease_levels) %>%
  mutate(Target_Disease = factor(Target_Disease, levels = rev(disease_levels)))

project_counts_e_filtered <- project_counts_e %>%
  filter(Target_Disease %in% disease_levels) %>%
  mutate(
    Target_Disease = factor(Target_Disease, levels = rev(disease_levels)),
    Label_X = pmin(Total_Samples, 8000)
  )

status_colors <- c("Disease" = "#E64B35", "HC" = "#3C5488")

p_f <- ggplot() +
  geom_bar(data = sample_counts_e, aes(x = Count, y = Target_Disease, fill = Sample_Group),
           stat = "identity", position = "stack", width = 0.75) +
  geom_text(data = project_counts_e_filtered,
             aes(x = Label_X, y = Target_Disease,
                 label = paste0(Project_Num, " bp; n=", scales::comma(Total_Samples))),
             inherit.aes = FALSE, hjust = -0.04, size = base_sz / 5.8, color = "#303030", fontface = "plain") +
  scale_fill_manual(name = "Status", values = status_colors) +
  scale_x_continuous(limits = c(0, 8000),
                     breaks = c(0, 2000, 4000, 6000, 8000),
                     labels = c("0", "2K", "4K", "6K", "8K"),
                     expand = expansion(mult = c(0, 0.15))) +
  scale_y_discrete(labels = disease_axis_labels) +
  labs(x = "Samples", y = NULL, tag = "f") +
  guides(fill = guide_legend(order = 2, title.position = "left", nrow = 1, byrow = TRUE)) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = base_sz - 2, base_family = "sans") +
  theme(
    axis.text.y = element_text(size = (base_sz - 2) * 0.80, face = "bold", color = "#303030",
                                margin = margin(r = 10)),
    axis.text.x = element_text(size = (base_sz - 2) * 0.80, face = "bold", color = "#303030",
                               margin = margin(t = 2)),
    axis.title.x = element_text(size = (base_sz - 2) * 1.0, face = "bold", color = "#303030", margin = margin(t = 10)),
    axis.line = element_line(color = "#303030", linewidth = 0.4),
    axis.ticks = element_line(color = "#303030", linewidth = 0.3),
    axis.ticks.length.x = unit(3, "pt"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = (base_sz - 2) * 0.8, face = "bold", color = "#303030"),
    legend.text = element_text(size = (base_sz - 2) * 0.7, color = "#303030"),
    legend.key.size = unit(3.5, "mm"),
    legend.margin = margin(t = -8, r = 0, b = 0, l = 0, unit = "pt"),
    legend.box.margin = margin(t = -8, r = 0, b = 0, l = 0, unit = "pt"),
    legend.spacing.y = unit(1, "mm"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 0, 0, 10, unit = "pt")
  )

# =============================================================================
# Compose Figure 1 v2: a | (b_with_c_inset / (d | e | f))
# =============================================================================
output_dir <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/01-left-pipeline-right-map/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

p_a <- p_a + tag_style
p_d <- p_d + tag_style
p_e <- p_e + tag_style
p_f <- p_f + tag_style 

p_b_base <- p_b_base + tag_style

p_b <- p_b_base + inset_element(p_c + tag_style, 
                                left = -0.16,
                                bottom = 0.05,  
                                right = 0.25,
                                top = 0.55)

right_col_bottom <- (p_d | p_e | p_f) +
  plot_layout(widths = c(1.0, 1.0, 1.4), guides = "keep")

right_col <- p_b / right_col_bottom +
  plot_layout(heights = c(4.6, 1.4))

p_figure1 <- p_a + right_col + plot_layout(widths = c(1.2, 3.6))

ggsave(file.path(output_dir, "Figure1_composite_v2.pdf"), plot = p_figure1, width = 26, height = 13, device = "pdf")
ggsave(file.path(output_dir, "Figure1_composite_v2.png"), plot = p_figure1, width = 26, height = 13, device = "png", dpi = 600, bg = "white")