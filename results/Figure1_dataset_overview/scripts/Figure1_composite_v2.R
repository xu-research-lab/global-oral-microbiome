#!/usr/bin/env Rscript
# Figure 1 v2: Composite figure for global oral microbiome dataset overview
# Changes from v1: saturation curve inset in panel b (bottom-left), panel c replaced by phylum stacked bar chart

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
anno_w <- 48.0
anno_h <- 1.5

right_anno$rect_xmin <- anno_x - 0.2
right_anno$rect_xmax <- anno_x + anno_w + 0.5
right_anno$rect_ymin <- right_anno$y_mid - anno_h / 2
right_anno$rect_ymax <- right_anno$y_mid + anno_h / 2

dataset_rect_xmax <- anno_x + anno_w + 0.5
dataset_rect_ymax <- 19.5
sample_rect_xmax <- anno_x + anno_w + 0.5
sample_rect_ymax <- 6.5

p_a <- ggplot() +
  annotate("rect", xmin = -max_w - 0.5, xmax = anno_x + anno_w + 0.5,
           ymin = 10.8, ymax = 19.5, fill = "#F8F9FA", color = "#B0BEC5", linetype = "dashed", alpha = 0.5) +
  annotate("text", x = dataset_rect_xmax - 1.0, y = dataset_rect_ymax - 0.4,
           label = "Dataset-Level Filters", fontface = "bold", hjust = 1, size = text_sz * 1.1, color = "#37474F") +
  annotate("rect", xmin = -max_w - 0.5, xmax = anno_x + anno_w + 0.5,
           ymin = -1.5, ymax = 7.5, fill = "#F3E5F5", color = "#CE93D8", linetype = "dashed", alpha = 0.3) +
  annotate("text", x = sample_rect_xmax - 1.0, y = sample_rect_ymax - 0.4,
           label = "Sample-Level Filters", fontface = "bold", hjust = 1, size = text_sz * 1.1, color = "#6A1B9A") +
  geom_polygon(data = conn_blocks, aes(x = x, y = y, group = grp), fill = "#CFD8DC", alpha = 0.45) +
  geom_polygon(data = poly_blocks, aes(x = x, y = y, group = grp), fill = rep(bar_fills, each = 4), color = NA) +
  geom_segment(data = funnel_data, aes(x = -half_w, xend = half_w, y = y + bar_h, yend = y + bar_h), color = "white", linewidth = 0.4) +
  geom_segment(data = funnel_data, aes(x = -half_w, xend = half_w, y = y - bar_h, yend = y - bar_h), color = "white", linewidth = 0.4) +
  geom_segment(data = arrow_data, aes(x = x, y = y_start, xend = x, yend = y_end), arrow = arrow(length = unit(0.3, "cm"), type = "closed"), linewidth = 0.8, color = "#546E7A") +
  annotate("segment", x = 0, xend = 0, y = 10.5, yend = 7.8,
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
  coord_cartesian(xlim = c(-14.5, 65.0), ylim = c(-0.5, 20.5), clip = "off") +
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
    country_label = paste0(Country_Cleaned, "\n(bp=", n_studies, "; n=", scales::comma(n_samples), ")")
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
# Saturation curve inset (to be embedded in panel b)
# =============================================================================
rarefaction_file <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/01-left-pipeline-right-map/output/rarefaction_pcofg.tsv"

df_rare <- read.table(rarefaction_file, header = TRUE, sep = "\t") %>%
  pivot_longer(cols = c("Phylum","Class","Order","Family","Genus"), names_to = "Taxonomy", values_to = "Unique_taxa") %>%
  mutate(Taxonomy = factor(Taxonomy, levels = c("Genus","Family","Order","Class","Phylum")))

p_c_inset <- ggplot(df_rare, aes(x = Samples, y = Unique_taxa, color = Taxonomy)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = 1.0) +
  scale_color_manual(values = c("Genus"="#5C6BC0", "Family"="#FFB300", "Order"="#EF5350", "Class"="#26C6DA", "Phylum"="#66BB6A")) +
  scale_x_continuous(breaks = seq(0, 80000, by = 20000),
                     limits = c(0, 82000),
                     labels = c("0", "20K", "40K", "60K", "80K"), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000),
                     labels = c("0", "1K", "2K", "3K", "4K", "5K"),
                     expand = expansion(mult = c(0.1, 0.06))) +
  annotate("text", x = 80000, y = 5200, label = "Genus", color = "#5C6BC0", fontface = "bold", size = base_sz/4.5) +
  annotate("text", x = 80000, y = 3700, label = "Family", color = "#FFB300", fontface = "bold", size = base_sz/4.5) +
  annotate("text", x = 80000, y = 2200, label = "Order", color = "#EF5350", fontface = "bold", size = base_sz/4.5) +
  annotate("text", x = 80000, y = 800, label = "Class", color = "#26C6DA", fontface = "bold", size = base_sz/4.5) +
  annotate("text", x = 80000, y = 200, label = "Phylum", color = "#66BB6A", fontface = "bold", size = base_sz/4.5) +
  labs(x = "Samples", y = "Observed taxa") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = base_sz - 6, base_family = "sans") +
  theme(
    axis.text = element_text(size = (base_sz - 6) * 0.9),
    axis.title = element_text(size = (base_sz - 6) * 1.0, face = "bold"),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5, unit = "pt"),
    plot.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA)
  )

# =============================================================================
# Panel b: Rendering Map (with saturation curve inset in bottom-left)
# =============================================================================
p_b <- ggplot() +
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
                    xmin = -5, xmax = 130,
                    ymin = -78, ymax = -48) +
  # Saturation curve inset in bottom-left (South Pacific, moved up)
  annotation_custom(ggplotGrob(p_c_inset),
                    xmin = -175, xmax = -100,
                    ymin = -75, ymax = -20) +
  coord_sf(xlim = c(-170, 180), ylim = c(-85, 75), expand = FALSE) +
  theme_void(base_size = base_sz, base_family = "sans") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = base_sz * 0.8, color = "#303030"),
    legend.title = element_text(size = base_sz * 0.9, face = "bold", color = "#303030", hjust = 0.5),
    legend.margin = margin(t = 0, b = -20),
    plot.margin = margin(t = 10, r = -30, b = 0, l = 20, unit = "pt")
  ) + labs(tag = "b")

# =============================================================================
# Panel c: Phylum relative abundance stacked bar chart (NEW - replaces old saturation curve)
# =============================================================================
abundance_file <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/01-left-pipeline-right-map/input/phylum-relative-table.tsv"
metadata_file_phylum <- metadata_path

target_sites_lower   <- c("saliva", "tooth", "tongue", "oral mucosa")
target_sites_display <- c("Saliva", "Tooth", "Tongue", "Oral Mucosa")
site_label_map       <- setNames(target_sites_display, target_sites_lower)
target_continents    <- c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")

# Read phylum relative abundance data
rel <- fread(abundance_file, header = TRUE, sep = "\t", skip = 1)
setnames(rel, 1, "Phylum")

# Read metadata for phylum chart
meta_phylum <- fread(metadata_file_phylum, header = TRUE, sep = "\t")
required_cols <- c("sample_id", "site_cleaned", "continent_group")
missing_cols  <- setdiff(required_cols, names(meta_phylum))
if (length(missing_cols) > 0) {
  # Try alternate column names
  if ("bioproject_id" %in% names(meta_phylum)) {
    meta_phylum[, sample_id := sample_id]
  }
  if ("country_cleaned" %in% names(meta_phylum) & !("continent_group" %in% names(meta_phylum))) {
    meta_phylum[, continent_group := countrycode(country_cleaned, "country.name", "continent")]
  }
}

# Melt to long format
phylum_cols <- setdiff(names(rel), "Phylum")
rel_long <- data.table::melt(rel, id.vars = "Phylum", measure.vars = phylum_cols,
                              variable.name = "sample_id", value.name = "Abundance")
rel_long[, sample_id := as.character(sample_id)]
rel_long[, Phylum_clean := gsub("^d__Bacteria;\\s*p__", "", Phylum)]
rel_long[, Phylum_clean := gsub("^p__", "", Phylum_clean)]
rel_long <- rel_long[!is.na(Phylum_clean) & Phylum_clean != ""]

meta_sub <- meta_phylum[, .(sample_id, site_cleaned, continent_group)]
meta_sub[, sample_id := as.character(sample_id)]

merged <- merge(rel_long, meta_sub, by = "sample_id", all.x = FALSE)

# Filter target sites and continents
merged <- merged[!is.na(site_cleaned) & site_cleaned != ""]
merged <- merged[!is.na(continent_group) & continent_group != "" & continent_group != "NA"]
merged[, site_cleaned := tolower(trimws(as.character(site_cleaned)))]
merged[, continent_group := trimws(as.character(continent_group))]
merged <- merged[site_cleaned %in% target_sites_lower]
merged <- merged[continent_group %in% target_continents]
merged[, site_display := site_label_map[site_cleaned]]
merged[, site_display := factor(site_display, levels = target_sites_display)]

# Top 9 phyla + Others
global_top9 <- merged[, .(Total = sum(Abundance, na.rm = TRUE)), by = Phylum_clean][order(-Total)][seq_len(min(9, .N))]$Phylum_clean
merged[, Phylum_group := ifelse(Phylum_clean %in% global_top9, Phylum_clean, "Others")]

plot_data <- merged[, .(Mean_Abundance = mean(Abundance, na.rm = TRUE)), by = .(site_display, continent_group, Phylum_group)]
plot_data[, Mean_Abundance := Mean_Abundance / sum(Mean_Abundance), by = .(site_display, continent_group)]

# Per-facet x-axis ordering by dominant phylum
dominant_phylum <- global_top9[1]
combined_levels <- character(0)
for (s in target_sites_display) {
  sub <- plot_data[site_display == s & Phylum_group == dominant_phylum]
  if (nrow(sub) == 0) {
    ordered_conts <- intersect(target_continents, unique(plot_data[site_display == s]$continent_group))
  } else {
    ordered_conts <- sub[order(-Mean_Abundance)]$continent_group
  }
  combined_levels <- c(combined_levels, paste(s, ordered_conts, sep = "|"))
}
plot_data[, facet_var := factor(paste(site_display, continent_group, sep = "|"), levels = combined_levels)]

strip_facet_prefix <- function(x) {
  sapply(x, function(label) gsub("^[^|]*\\|", "", as.character(label)))
}

# Color palette and stacking order
color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
unified_colors <- setNames(color_palette[seq_along(global_top9)], global_top9)
unified_colors["Others"] <- "#CCCCCC"

other_phyla  <- setdiff(global_top9, dominant_phylum)
legend_order <- c("Others", rev(other_phyla), dominant_phylum)
plot_data[, Phylum_group := factor(Phylum_group, levels = legend_order)]

# Short continent names for x-axis
continent_short <- c("Africa" = "AF", "Asia" = "AS", "Europe" = "EU", "North America" = "NAM", "Oceania" = "OA", "South America" = "SA")

p_c <- ggplot(plot_data, aes(x = facet_var, y = Mean_Abundance, fill = Phylum_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.95) +
  facet_grid(. ~ site_display, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = unified_colors, name = "Phylum", breaks = legend_order, labels = legend_order) +
  scale_x_discrete(labels = function(x) sapply(x, function(lab) {
    cont <- gsub("^[^|]*\\|", "", as.character(lab))
    if (cont %in% names(continent_short)) continent_short[cont] else cont
  })) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Continent", y = "Relative Abundance", tag = "c") +
  theme_minimal(base_size = base_sz - 3, base_family = "sans") +
  theme(
    axis.text.x = element_text(size = (base_sz - 3) * 0.85, face = "bold", color = "#303030", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = (base_sz - 3) * 0.85, face = "bold", color = "#303030"),
    axis.title = element_text(size = (base_sz - 3) * 1.0, face = "bold", color = "#303030"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = (base_sz - 3) * 1.0, face = "bold", color = "#303030"),
    legend.text = element_text(size = (base_sz - 3) * 0.8, color = "#303030"),
    legend.key.size = unit(3, "mm"),
    strip.text = element_text(size = (base_sz - 3) * 1.0, face = "bold", color = "#303030"),
    strip.background = element_rect(fill = "#E3F2FD", color = "#90CAF9", linewidth = 0.4),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    plot.margin = margin(10, 10, 10, 10, unit = "pt")
  )

# =============================================================================
# Panel d: Donut charts (same as v1)
# =============================================================================
donut_data <- meta_clean %>% filter(!is.na(Site_Cleaned)) %>%
  count(Continent, Site_Cleaned, name = "count") %>% group_by(Continent) %>%
  mutate(prop = count / sum(count), Site_Cleaned = factor(Site_Cleaned, levels = site_cols_ordered)) %>%
  ungroup() %>%
  mutate(
    Continent_Label = ifelse(Continent == "NAM", "NA", as.character(Continent)),
    Continent_Label = factor(Continent_Label, levels = c("NA", "AS", "EU", "AF", "OA", "SA"))
  )

p_d <- ggplot(donut_data, aes(x = 2, y = prop, fill = Site_Cleaned)) +
  geom_col(width = 1, color = "white", linewidth = 0.15) + coord_polar(theta = "y", direction = -1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ Continent_Label, nrow = 2) + scale_fill_manual(name = "Site", values = site_colors, limits = site_cols_ordered) +
  theme_void(base_size = base_sz, base_family = "sans") +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(face = "bold", size = base_sz * 0.7),
    legend.text = element_text(size = base_sz * 0.65, color = "#303030"),
    legend.key.size = unit(2.5, "mm"),
    legend.box.spacing = unit(-5, "mm"),
    legend.spacing.y = unit(0, "mm"),
    legend.margin = margin(l = -30, unit = "pt"),
    legend.box.margin = margin(l = -20, unit = "pt"),
    strip.text = element_text(size = base_sz * 0.8, face = "bold", color = "#303030"),
    plot.margin = margin(t = -15, r = 0, b = 0, l = -15, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    plot.tag = element_text(size = base_sz * 1.5, face = "bold", color = "#303030", hjust = 0, vjust = 1, margin = margin(b = 0, unit = "pt"))
  ) + labs(tag = "d")

# =============================================================================
# Compose Figure 1 v2: a | (b / (c | d))
# =============================================================================
output_dir <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/01-left-pipeline-right-map/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

p_a <- p_a + tag_style
p_b <- p_b + tag_style
p_c <- p_c + tag_style
p_d <- p_d + tag_style

bottom_row <- (p_c | p_d) + plot_layout(widths = c(0.8, 1.3))
right_col <- (p_b / bottom_row) + plot_layout(heights = c(2.6, 0.9))
p_figure1 <- p_a + right_col + plot_layout(widths = c(1.2, 3.0))

ggsave(file.path(output_dir, "Figure1_composite_v2.pdf"), plot = p_figure1, width = 22, height = 12, device = "pdf")
ggsave(file.path(output_dir, "Figure1_composite_v2.png"), plot = p_figure1, width = 22, height = 12, device = "png", dpi = 600, bg = "white")