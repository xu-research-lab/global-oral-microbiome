# Global Sampling World Map

![Global Sampling World Map](global_sampling_world_map.png)

---

## Figure Legend

**Global distribution of curated oral microbiome samples with site composition inset.**

**(Main panel) World map of sample distribution.** Country centroids (from Natural Earth, scale 110) are shown as circles colored by continent (NA = red, SA = purple, EU = green, AS = orange, AF = brown, OA = cyan). Circle sizes are scaled proportionally to the number of samples per country (legend breaks: 10, 100, 1K, 10K, 40K). Country labels display the BioProject count (bp) and sample count (n) in the format "Country (bp=X; n=Y)"; countries with fewer than 100 samples are shown as circles without text labels to reduce visual clutter. A demographic summary table embedded in the Southern Ocean area reports, per continent: total sample count, number of countries, number of studies, mean age ± SD, mean BMI ± SD (filtered to values between 10 and 60), and sex ratio (F%/M%, excluding non-binary entries). North America is abbreviated as "NA" in the table. Antarctica is excluded from the map. Barbados is manually positioned at (−59.5°, 13.2°) due to missing centroid data.

**(Inset, bottom-left) Oral site composition bubble matrix.** A bubble matrix displays the distribution of seven oral sampling sites (Saliva, Tooth, Oral Mucosa, Oral Swab, Tongue, Gingiva, Others) across six continents (NA, AS, EU, AF, OA, SA). Each bubble's size is proportional to the number of samples for that (site, continent) combination. Bubbles are color-coded by site: Saliva (red, #E31A1C), Tooth (amber, #FFB300), Oral Mucosa (blue, #1565C0), Oral Swab (purple, #8E24AA), Tongue (green, #33A02C), Gingiva (brown, #B15928), Others (light blue, #A6CEE3). The x-axis lists oral sites and the y-axis lists continents (abbreviated). The inset shares the world map's "Number of samples" legend and does not display its own size legend. The inset is positioned at the lower-left of the world map (left = 0.02, bottom = 0.05, right = 0.18, top = 0.40) to avoid overlapping South America and the left boundary.

---

## Methods

### World Map: Geographic Distribution and Demographic Summary

Country-level sample and study counts were computed from the curated metadata by grouping on `Country_Cleaned` and `Continent`. ISO 3166-1 alpha-3 country codes were derived using the `countrycode` R package. Country centroids were obtained from Natural Earth (scale 110) via `rnaturalearth::ne_countries()` and `sf::st_point_on_surface()`. Barbados was manually assigned coordinates (−59.5°, 13.2°) due to missing centroid data in the Natural Earth dataset. Antarctica was excluded from the map. The map projection uses `coord_sf` with `xlim = c(−180, 180)` and `ylim = c(−60, 85)`.

Bubble sizes were mapped to sample counts (range = c(1, 25); breaks at 10, 100, 1K, 10K, 40K). Country labels show study count (bp) and sample count (n) in the format "Country\n(bp=X; n=Y)", placed using `ggrepel::geom_text_repel` (segment.color = "#999999", segment.size = 0.2, box.padding = 0.4, max.overlaps = Inf). Only countries with ≥100 samples are labeled; countries with fewer samples are represented by circles without text.

A demographic summary table was computed per continent, sorted by descending sample count. Columns: Samples, Countries, Studies, Age (mean ± SD), BMI (mean ± SD, filtered to values >10 and <60), and Sex (F%/M%, excluding non-binary entries and computed only from male/female categories). The table was rendered using `gridExtra::tableGrob` with `ttheme_minimal` (alternating white/grey rows, bold header on grey-blue background) and embedded in the Southern Ocean area via `annotation_custom` at `xmin = −40, xmax = 120, ymin = −60, ymax = −30`. North America is displayed as "NA" in the table.

### Inset: Site × Continent Bubble Matrix

Oral site composition was computed per continent from the curated metadata. Sites "Nasal" and "Tonsils" were excluded; "Other" was renamed to "Others". For each (Continent, Site_Cleaned) pair, the sample count was tabulated. The bubble matrix was generated using `ggplot2::geom_point(shape = 21, ...)` with site as fill color and sample count as point size (range = c(1.5, 8)). Sites are ordered along the x-axis as: Saliva, Tooth, Oral Mucosa, Oral Swab, Tongue, Gingiva, Others. Continent labels on the y-axis (NA, AS, EU, AF, OA, SA) are displayed in reverse order for visual alignment. The size legend is suppressed (`guide = "none"`); the inset shares the world map's "Number of samples" legend. The inset uses `theme_classic(base_size = base_sz − 4)` with reduced axis text sizes (`(base_sz − 4) × 0.70`) and compact margins (`margin(4, 4, 4, 8, "pt")`). The inset is embedded into the lower-left corner of the world map via `patchwork::inset_element` (left = 0.02, bottom = 0.05, right = 0.18, top = 0.40), positioned to avoid overlapping South America and the left plot boundary.

---

## Output

- **Script**: `global_sampling_world_map_plot.R`
- **Output**: `global_sampling_world_map.pdf`, `global_sampling_world_map.png` (18 × 10 inches, 600 DPI)

## Abbreviations

- **bp**: number of BioProjects (studies) for a given country
- **n**: number of samples for a given country
- **NA**: North America
- **SA**: South America
- **EU**: Europe
- **AS**: Asia
- **AF**: Africa
- **OA**: Oceania

## Input Data

- `input/metadata_oral_qc.tsv` — Final curated metadata with `bioproject_id`, `sample_id`, `site_cleaned`, `continent_group`, `country_cleaned`, `age_cleaned`, `sex_cleaned`, `bmi_cleaned` columns