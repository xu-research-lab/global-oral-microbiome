# Figure 1. Overview of the Global Oral Microbiome Dataset

![Figure 1](results/Figure1_composite_v1.png)

---

## Figure Legend

**Figure 1. Overview of the global oral microbiome dataset: two-phase quality control pipeline, geographic distribution with taxonomic saturation, phylum composition by oral site, and site composition across continents.**

**(a) Two-phase quality control funnel.** Seven sequential steps are grouped into two phases. The Dataset-Level Filters phase (Steps 1–3, enclosed by a grey dashed rectangle) progressively removed non-oral samples (Step 1, −30,388 samples), datasets with >30% oral samples and chimera rate >25% (Step 2, −22 studies), and datasets with >50% oral samples and ASV count ≤30 or closed-reference mapping rate <30% (Step 3, −19 studies). The Sample-Level Filters phase (Steps 4–5, enclosed by a purple dashed rectangle) then removed individual samples with chimera rate >25% (Step 4, −1,109 samples) and samples with ASV count ≤30 and closed-reference mapping rate <30% (Step 5, −7,887 samples). Block widths are proportional to retained sample count. Each block displays the sample count (bold) and study count (plain). Right-side annotation boxes detail the step name (bold), filtering criterion (italic), and number of removed samples or studies (bold red). Inter-block arrows indicate the direction of filtering; a thick arrow separates the two phases. The final curated dataset comprises 73,340 samples from 293 studies.

**(b) Global distribution of curated oral microbiome samples.** Country centroids (from Natural Earth, scale 110) are shown as circles colored by continent (NAM = red, SA = purple, EU = green, AS = orange, AF = brown, OA = cyan). Circle sizes are scaled on a log10 axis proportional to the number of samples per country (legend breaks: 10, 100, 1K, 10K, 40K). Country labels display the BioProject count (bp) and sample count (n) in the format "Country (bp=X; n=Y)"; **countries with fewer than 100 samples are shown as circles without text labels** to reduce visual clutter. A demographic summary table embedded in the Southern Ocean area reports, per continent: total sample count, number of countries, number of studies, mean age ± SD, mean BMI ± SD (filtered to values between 10 and 60), and sex ratio (F%/M%, excluding non-binary entries). North America is abbreviated as "NA" in the table. Antarctica is excluded from the map. Barbados is manually positioned at (−59.5°, 13.2°) due to missing centroid data. An inset in the bottom-left (South Pacific) displays taxonomic saturation curves (see panel b inset description below).

**(b inset) Taxonomic saturation curves.** Accumulation of unique taxa across five taxonomic ranks (Genus, Family, Order, Class, Phylum) is shown as a function of the number of samples (up to 73,340). Curves were GAM-smoothed (mgcv, cubic spline, y ~ s(x, bs = "cs")) and color-coded: Genus (indigo, #5C6BC0), Family (amber, #FFB300), Order (red, #EF5350), Class (cyan, #26C6DA), Phylum (green, #66BB6A). Inline text labels identify each rank at x = 80,000; y-coordinates are spread to avoid overlap. The inset has a white background with a grey border and no legend. For details on rarefaction computation, see the Methods section.

**(c) Phylum-level relative abundance across oral sites and continents.** Stacked bar charts show the mean relative abundance of the top 9 most abundant bacterial phyla (by global total abundance) plus an "Others" category for all remaining phyla. Four oral sites are displayed as separate facets (Saliva, Tooth, Tongue, Oral Mucosa). Within each facet, continents on the x-axis are ordered by descending abundance of the dominant phylum (Firmicutes). Continent names are abbreviated: AF = Africa, AS = Asia, EU = Europe, NAM = North America, OA = Oceania, SA = South America. The dominant phylum is stacked at the bar bottom; the legend order places "Others" at the top and the dominant phylum at the bottom. Mean relative abundances are renormalized to 100% per (site, continent) combination. Sites other than the four target sites (e.g., Oral Swab, Gingiva) and continent group "NA" are excluded from this panel. Colors follow the Set1 palette with grey (#CCCCCC) for Others.

**(d) Oral site composition by continent.** Donut charts display the proportional distribution of seven oral sampling sites (Saliva, Tooth, Oral Mucosa, Oral Swab, Tongue, Gingiva, Others) for each of six continents (NA, AS, EU, AF, OA, SA), arranged in a 2 × 3 grid. The remaining part of the oral is named "Others". Proportions are calculated as site count divided by total count within each continent. North America is abbreviated as "NA" in the facet labels. Colors: Saliva (red, #E31A1C), Tooth (amber, #FFB300), Oral Mucosa (blue, #1565C0), Oral Swab (purple, #8E24AA), Tongue (green, #33A02C), Gingiva (brown, #B15928), Others (light blue, #A6CEE3). The Site legend is positioned to the right of the donut grid.

---

## Methods

### Panel a: Quality Control Pipeline

The merged metadata table was sequentially filtered through five steps organized into two phases. In the Dataset-Level phase, Step 1 removed non-oral samples (−30,388 samples, 70.6% retained), Step 2 excluded datasets with >30% oral samples and chimera rates >25% (−22 studies), and Step 3 excluded datasets with >50% oral samples and ASV counts ≤30 or closed-reference mapping rates <30% (−19 studies). In the Sample-Level phase, Step 4 removed individual samples with chimera rates >25% (−1,109 samples, 89.7% retained), and Step 5 removed samples with ASV counts ≤30 and closed-reference mapping rates <30% (−7,887 samples, 99.3% retained). The final dataset comprised 73,340 samples from 293 studies.

In the funnel visualization, block widths are proportional to retained sample count (half_width = min_w + (n − n_min) / (n_max − n_min) × (max_w − min_w), where min_w = 8.5, max_w = 13.5). Inter-block connectors (trapezoids) and directional arrows link consecutive steps. A thick downward arrow separates the two phases. Right-side annotation boxes (width = 48.0, height = 1.5) detail each step name (bold), filtering criterion (italic), and removed count (bold red), connected to the inter-block arrows via dashed lines. The right edges of the annotation boxes align with the right edges of the dashed phase-grouping rectangles (rect_xmax = anno_x + anno_w + 0.5).

### Panel b: Geographic Distribution and Demographic Summary

Country-level sample and study counts were computed from the curated metadata by grouping on `Country_Cleaned` and `Continent`. ISO 3166-1 alpha-3 country codes were derived using the `countrycode` R package. Country centroids were obtained from Natural Earth (scale 110) via `rnaturalearth::ne_countries()` and `sf::st_point_on_surface()`. Barbados was manually assigned coordinates (−59.5°, 13.2°) due to missing centroid data in the Natural Earth dataset. Antarctica was excluded from the map. The map projection uses `coord_sf` with `xlim = c(−170, 180)` and `ylim = c(−85, 75)`.

Bubble sizes were mapped to sample counts on a log10 scale (range = c(1, 25); breaks at 10, 100, 1K, 10K, 40K). Country labels show study count (bp) and sample count (n) in the format "Country\n(bp=X; n=Y)", placed using `ggrepel::geom_text_repel` (segment.color = "#999999", segment.size = 0.2, box.padding = 0.4, max.overlaps = Inf). **Only countries with ≥100 samples are labeled**; countries with fewer samples are represented by circles without text.

A demographic summary table was computed per continent, sorted by descending sample count. Columns: Samples, Countries, Studies, Age (mean ± SD), BMI (mean ± SD, filtered to values >10 and <60), and Sex (F%/M%, excluding non-binary entries and computed only from male/female categories). The table was rendered using `gridExtra::tableGrob` with `ttheme_minimal` (alternating white/grey rows, bold header on grey-blue background) and embedded in the Southern Ocean area via `annotation_custom` at `xmin = −5, xmax = 130, ymin = −78, ymax = −48`. North America is displayed as "NA" in the table.

### Panel b inset: Taxonomic Saturation Curves

To assess whether sampling depth was sufficient to capture microbial diversity, rarefaction curves were computed across five taxonomic ranks (Phylum, Class, Order, Family, Genus). The BIOM feature table (`input/biom_oral_qc.biom`) and taxonomy annotations (`input/taxonomy.tsv`) were loaded. For each ASV/OTU, taxonomic levels were parsed from semicolon-delimited strings using the prefixes p__, c__, o__, f__, and g__. Incompletely annotated taxa (empty strings after prefix removal) were excluded from counting. Samples were randomly shuffled (seed = 42) and accumulated at increasing checkpoints (1, 5, 10, 50, 100, 200, 500, 800, every 1K from 1K–10K, every 5K from 10K to the full 73,340). At each checkpoint, the sparse matrix was subset to the corresponding sample indices, and present features were identified via `scipy.sparse.csc_matrix.getnnz(axis=1)`. Unique taxa per rank were counted in parallel across CPU cores using Python `multiprocessing.Pool`. The resulting rarefaction table (`results/rarefaction_pcofg.tsv`) was plotted in R using GAM-smoothed curves (`mgcv::gam`, formula y ~ s(x, bs = "cs")) with `theme_classic` (base_size = base_sz − 6). Five ranks were color-coded: Genus (#5C6BC0), Family (#FFB300), Order (#EF5350), Class (#26C6DA), Phylum (#66BB6A). Inline text labels are positioned at x = 80,000 with evenly spread y-coordinates (5200, 3700, 2200, 800, 200) to avoid overlap. The inset was rendered with a white background and grey border, and embedded in panel b via `annotation_custom(ggplotGrob(...))` at South Pacific coordinates (xmin = −175, xmax = −100, ymin = −75, ymax = −20).

### Panel c: Phylum Relative Abundance

Phylum-level relative abundance data were read from `input/phylum-relative-table.tsv` (skipping one BIOM comment line). The curated metadata (`input/metadata_oral_qc_update.tsv`) provided `sample_id`, `site_cleaned`, and `continent_group` columns for merging. Data were filtered to four oral sites (saliva, tooth, tongue, oral mucosa) and six continents (Africa, Asia, Europe, North America, Oceania, South America); records with missing or "NA" continent_group were excluded. Taxonomy prefixes (d__Bacteria;p__ or p__) were stripped from phylum names. The top 9 phyla by global total abundance were identified; all remaining phyla were grouped as "Others". For each (site, continent) combination, mean relative abundance was computed and renormalized to 100%. Within each facet, continents were ordered by descending abundance of the dominant phylum (the phylum with the highest global total abundance). Stacked bars were generated using `ggplot2::geom_bar(stat = "identity", position = "stack")` with `facet_grid(. ~ site_display, scales = "free_x", space = "free_x")`. The y-axis displays percentage via `scales::percent_format(accuracy = 1)` with breaks at 0%, 25%, 50%, 75%, 100%. Continent names on the x-axis were abbreviated (AF, AS, EU, NAM, OA, SA). The dominant phylum is stacked at the bar bottom; the legend order places "Others" at the top and the dominant phylum at the bottom. Colors follow the Set1 palette with grey (#CCCCCC) for Others. The facet header background is light blue (#E3F2FD) with a blue border (#90CAF9).

### Panel d: Site Composition

Proportional composition of seven oral sampling sites (Saliva, Tooth, Oral Mucosa, Oral Swab, Tongue, Gingiva, Others) was computed per continent from the curated metadata. Sites "Nasal" and "Tonsils" were excluded from the dataset; "Other" was renamed to "Others". Proportions were calculated as count / sum(count) within each continent group. Donut charts were generated using `ggplot2::geom_col(aes(x = 2, y = prop, fill = Site_Cleaned))` with `coord_polar(theta = "y", direction = -1)`, displayed as a 2 × 3 grid via `facet_wrap(~ Continent_Label, nrow = 2)`. The x = 2 aesthetic produces larger donuts compared to x = 1. North America is abbreviated as "NA" in the facet labels. Continent order in the grid: NA, AS, EU, AF, OA, SA. The Site legend is positioned to the right with minimized spacing (legend.box.spacing = unit(−5, "mm"), legend.margin = margin(l = −30, unit = "pt")). Colors: Saliva = #E31A1C, Tooth = #FFB300, Oral Mucosa = #1565C0, Oral Swab = #8E24AA, Tongue = #33A02C, Gingiva = #B15928, Others = #A6CEE3.

---

## Composite Layout

- Layout: `p_a + (p_b / (p_c | p_d))` with `plot_layout(widths = c(1.2, 3.0))` and inner `heights = c(2.6, 0.9)`, `widths = c(0.8, 1.3)`
- Output: `results/Figure1_composite_v1.png`, 22 × 12 inches, 600 DPI

## Abbreviations

- **bp**: number of BioProjects (studies) for a given country
- **n**: number of samples for a given country
- **NA**: North America
- **SA**: South America
- **EU**: Europe
- **AS**: Asia
- **AF**: Africa
- **OA**: Oceania
- **CR**: Closed-reference mapping rate
- **ASV**: Amplicon sequence variant

## Input Data

All input files are located in the `input/` directory, scripts in `scripts/`, and outputs in `results/`:

- `input/metadata_oral_qc_update.tsv` — Final curated metadata with `bioproject_id`, `sample_id`, `site_cleaned`, `continent_group`, `country_cleaned`, `age_cleaned`, `sex_cleaned`, `bmi_cleaned` columns
- `input/phylum-relative-table.tsv` — Phylum-level relative abundance table (BIOM format with one comment line to skip)
- `input/taxonomy.tsv` — Taxonomy annotations with `p__/c__/o__/f__/g__` prefixes for rarefaction computation (Python script)
- `input/biom_oral_qc.biom` — BIOM feature table used for rarefaction computation (Python script)
- `results/rarefaction_pcofg.tsv` — Pre-computed rarefaction/saturation data across five taxonomic ranks (generated by `scripts/oral_sampling_saturation_curve_analysis.py`)
- `results/Figure1_composite_v1.png` — Composite figure output (22 × 12 inches, 600 DPI)
---

## Version 2 (Figure1_composite_v2)

![Figure 1 v2](results/Figure1_composite_v2.png)

### Changes from Version 1

Version 2 reorganizes the figure layout: (a) the phylum stacked bar chart (formerly panel c) is now embedded as an inset in the Pacific Ocean area of the world map (panel b); (b) the saturation curves (formerly panel b inset) are now a standalone panel e; (c) the donut charts (formerly panel d) are replaced by a site × continent bubble matrix (panel d) and a health/disease stacked bar chart (panel f). The overall composition is `a | (b_with_c_inset / (d | e | f))`.

---

## Figure Legend (v2)

**(a) Two-phase quality control funnel.** Same structure as v1: seven sequential steps organized into Dataset-Level Filters (Steps 1–3, grey dashed rectangle) and Sample-Level Filters (Steps 4–5, purple dashed rectangle). Block widths are proportional to retained sample count. Each block displays sample count (bold) and study count (plain). Right-side annotation boxes are wider (width = 75.0 vs. 48.0 in v1) and positioned at `x = max_w + 2.0`, with the "Dataset-Level Filters" and "Sample-Level Filters" labels right-aligned to the dashed rectangle right edges. A thick separator arrow divides the two phases. Final curated dataset: 73,340 samples from 293 studies.

**(b) Global distribution of curated oral microbiome samples.** Country centroids from Natural Earth (scale 110) shown as circles colored by continent (NAM = red, SA = purple, EU = green, AS = orange, AF = brown, OA = cyan). Circle sizes scaled on a log10 axis (legend breaks: 10, 100, 1K, 10K, 40K). Country labels display the BioProject count (bp) and sample count (n) in the format "Country\n(X bp; n=Y)"; countries with fewer than 100 samples are shown as circles without text labels. The map uses `coord_sf(xlim = c(−180, 180), ylim = c(−60, 85))`. A demographic summary table embedded in the Southern Ocean area (xmin = −40, xmax = 120, ymin = −60, ymax = −30) reports per continent: total sample count, number of countries, number of studies, mean age ± SD, mean BMI ± SD (filtered to 10–60), and sex ratio (F%/M%). North America is abbreviated as "NA" in the table. Barbados is manually positioned at (−59.5°, 13.2°). Panel c (phylum stacked bar chart) is embedded as an inset in the Pacific Ocean area of panel b.

**(c) Phylum-level relative abundance across oral sites and continents (inset in panel b).** Stacked bar charts of mean relative abundance for the top 9 most abundant phyla plus "Others". Four oral sites (Saliva, Tooth, Tongue, Oral Mucosa) are displayed as separate facets. Continents on x-axis are ordered by descending abundance of the dominant phylum. Continent abbreviations: AF, AS, EU, NAM, OA, SA. Colors follow the Set1 palette with grey (#CCCCCC) for Others. The inset is positioned at approximate Pacific Ocean coordinates (left = −0.16, bottom = 0.05, right = 0.25, top = 0.55 in relative NPC units) using `patchwork::inset_element`. The inset has a transparent background and uses `theme_classic(base_size = 13)` for compact rendering.

**(d) Site × continent bubble matrix.** Bubbles (points with `shape = 21`) display the sample count for each combination of 7 oral sites (x-axis: Saliva, Tooth, Oral Mucosa, Oral Swab, Tongue, Gingiva, Others) and 6 continents (y-axis: NA, AS, EU, AF, OA, SA from top to bottom). Bubble fill corresponds to site; bubble size reflects sample count. The "Number of samples" legend uses breaks at 50, 500, 5,000, and 15,000 (labeled "50", "0.5K", "5K", "15K") in a single horizontal row. X-axis labels are angled at 30° for readability. Continents are displayed with `scale_y_discrete(limits = rev(levels))` to maintain top-to-bottom geographic order. The plot uses `coord_cartesian(clip = "off")` to prevent label clipping.

**(e) Taxonomic saturation curves (standalone).** Accumulation of unique taxa across five taxonomic ranks (Genus, Family, Order, Class, Phylum) as a function of sample count (up to 73,340). Same rarefaction data as v1, plotted with GAM-smoothed curves (`mgcv::gam`, formula y ~ s(x, bs = "cs")). Curve colors: Genus (#5C6BC0), Family (#FFB300), Order (#EF5350), Class (#26C6DA), Phylum (#66BB6A). Inline text labels at x = 81,000 identify each rank. X-axis limits extended to 98,000. The plot uses `theme_classic(base_size = 11)` with axis breaks at 0, 20K, 40K, 60K, 80K. Legend is suppressed; rank identification is via inline labels.

**(f) Health/disease stacked bar chart.** For BioProjects with `disease_study_cleaned == "Yes"`, the primary disease per BioProject is determined as the most common non-"HC" `status_cleaned` value. BioProjects containing only HC samples within disease studies are excluded from display. The top 15 diseases by total sample count are shown as horizontal stacked bars, where each bar comprises HC (healthy control, #3C5488) and Disease (#E64B35) sample segments. Bar length represents total sample count per disease. End-of-bar labels display "X bp; n=Y" (BioProject count and total sample count). Disease names are abbreviated where applicable: Gastroesophageal Reflux Disease → GERD, Acute Myeloid Leukemia → AML, Head and Neck Cancer → HNC. X-axis range: 0–8,000, with breaks at 0, 2K, 4K, 6K, 8K. The Status legend is placed below the x-axis. The plot uses `coord_cartesian(clip = "off")` and `theme_classic(base_size = 15)` with minimal right margin.

---

## Methods (v2)

### Panel a: Quality Control Pipeline

The pipeline funnel uses the same underlying data and filtering logic as v1. The annotation boxes are widened to `anno_w = 75.0` (from 48.0 in v1) positioned at `anno_x = max_w + 2.0`, with `coord_cartesian(xlim = c(−14.5, 95.0))`. The "Dataset-Level Filters" and "Sample-Level Filters" labels are right-aligned to `rect_xmax − 1.0` of their respective dashed rectangles. The thick separator arrow spans from y = 9.95 to y = 7.55 with `linewidth = 2.0`.

### Panel b: Geographic Distribution and Demographic Summary

Country-level sample and study counts were computed from the curated metadata by grouping on `Country_Cleaned` and `Continent`. ISO 3166-1 alpha-3 codes were derived via `countrycode`. South America is manually assigned for countries in the Americas with ISO codes matching the `south_america_iso3c` list. Country centroids were obtained from Natural Earth (scale 110) via `sf::st_point_on_surface()`. Barbados manually assigned coordinates (−59.5°, 13.2°). Antarctica is excluded.

The map projection uses `coord_sf(xlim = c(−180, 180), ylim = c(−60, 85), expand = FALSE)`. Bubble sizes mapped to sample counts on a log10 scale (`range = c(1, 25)`; breaks at 10, 100, 1K, 10K, 40K). Country labels use the format "Country\n(X bp; n=Y)" via `ggrepel::geom_text_repel` (segment.color = "#999999", segment.size = 0.2, box.padding = 0.4, max.overlaps = Inf). Only countries with ≥100 samples are labeled.

The demographic summary table (computed identically to v1) is embedded via `annotation_custom` at `xmin = −40, xmax = 120, ymin = −60, ymax = −30`, positioned below Africa and Australia. The table uses `ttheme_minimal` with alternating white/grey rows and a bold header on grey-blue background.

### Panel c: Phylum Relative Abundance (Inset in Panel b)

The phylum-level relative abundance analysis uses the same data processing pipeline as v1. The key change is that the resulting plot is rendered at a reduced size (`base_size = base_sz − 4 = 13`) and embedded as an inset within panel b via `patchwork::inset_element`. The inset is positioned at `left = −0.16, bottom = 0.05, right = 0.25, top = 0.55` (NPC coordinates relative to the map panel), placing it over the Pacific Ocean area. The inset background is transparent (`plot.background = element_rect(fill = "transparent")`) and the legend uses a compact 4-row layout.

### Panel d: Site × Continent Bubble Matrix

Sample counts for each combination of oral site (7 categories: Saliva, Tooth, Oral Mucosa, Oral Swab, Tongue, Gingiva, Others) and continent (6 categories: NA, AS, EU, AF, OA, SA) were tallied from the curated metadata. The bubble matrix uses `ggplot2::geom_point(shape = 21)` with:

- **x-axis**: Site mapped to numeric positions (`Site_X = as.numeric(Site_Cleaned) * 0.52`), with breaks at `seq_along(site_cols_ordered) * 0.52` labeled by site name. X-axis labels are rotated 30° for readability.
- **y-axis**: Continent displayed in reverse order via `scale_y_discrete(limits = rev(levels))` so that NA appears at the top. North America is abbreviated as "NA" in the continent labels.
- **Size scale**: Range 2–12, with breaks at 50, 500, 5,000, 15,000 and labels "50", "0.5K", "5K", "15K". The legend is presented as a single horizontal row (`guide_legend(nrow = 1)`).
- **Fill color**: Site colors from the global `site_colors` palette. The fill legend is suppressed (`guides(fill = "none")`) to avoid redundancy.

`coord_cartesian(clip = "off")` is used to prevent label clipping at the plot boundaries.

### Panel e: Taxonomic Saturation Curves

The same rarefaction data (`results/rarefaction_pcofg.tsv`, generated by `scripts/oral_sampling_saturation_curve_analysis.py`) is used as in v1. Five taxonomic ranks (Phylum, Class, Order, Family, Genus) are plotted with GAM-smoothed curves via `geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))`. The x-axis extends to 98,000 with breaks every 20,000 (labeled "0", "20K", "40K", "60K", "80K"). Inline text labels identify each rank at x = 81,000; y-positions are manually set to avoid overlap (4900, 1400, 800, 350, 5 for Genus through Phylum). The legend is suppressed. The plot uses `theme_classic(base_size = 11)` with transparent background.

### Panel f: Health/Disease Stacked Bar Chart

**Data cleaning.** The raw metadata was filtered as follows:
1. `status_cleaned` values were normalized: "hc" → "HC", empty strings → `NA`, and NA-like strings ("na", "n/a", "nan", "null", "none", "unknown", "unk") were filtered out.
2. `disease_study_cleaned` values were converted to title case and filtered to only "Yes" or "No".
3. Only BioProjects with `disease_study_cleaned == "Yes"` were retained for analysis.

**Disease assignment.** For each BioProject, the most common non-"HC" `status_cleaned` value was designated as the `Target_Disease`. BioProjects containing only HC samples were labeled "Unclassified disease".

**Top-15 selection.** Disease labels matching "Unclassified disease", "Other", or "Others" (case-insensitive) were excluded. The remaining diseases were ranked by total sample count (descending), and the top 15 were selected.

**Disease name abbreviation.** Three long disease names were abbreviated for display:
- Gastroesophageal Reflux Disease → GERD
- Acute Myeloid Leukemia → AML
- Head and Neck Cancer → HNC

**Visualization.** A horizontal stacked bar chart was generated using `ggplot2::geom_bar(stat = "identity", position = "stack")`:
- **Fill**: `Sample_Group` (Disease = #E64B35, HC = #3C5488), with legend named "Status".
- **x-axis**: Sample count, with `scale_x_continuous(limits = c(0, 8000), breaks = c(0, 2000, 4000, 6000, 8000), labels = c("0", "2K", "4K", "6K", "8K"))`.
- **y-axis**: Disease names with abbreviated labels via `scale_y_discrete(labels = disease_axis_labels)`. Y-axis text uses `margin(r = 10)` for spacing.
- **End-of-bar labels**: `geom_text` displaying "X bp; n=Y" (BioProject count and total sample count), positioned at `pmax(Total_Samples, 8000)` with `hjust = −0.04`.
- **Status legend**: Placed at the bottom, single row (`guide_legend(nrow = 1)`).
- The plot uses `coord_cartesian(clip = "off")` and `theme_classic(base_size = 15)` with `plot.margin = margin(10, 0, 0, 10)` (right margin = 0 for edge alignment).

---

## Composite Layout (v2)

- Layout: `p_a + ((p_b_with_inset) / (p_d | p_e | p_f))` with `plot_layout(widths = c(1.2, 3.6))`
  - Inner row: `heights = c(4.6, 1.4)`
  - Bottom-row panels: `widths = c(1.0, 1.0, 1.4)`, `guides = "keep"`
- Panel c (phylum bars) embedded in panel b via `patchwork::inset_element(p_c, left = −0.16, bottom = 0.05, right = 0.25, top = 0.55)`
- Output: `results/Figure1_composite_v2.png`, 26 × 13 inches, 600 DPI
- Output also saved as `results/Figure1_composite_v2.pdf` (same dimensions)

### Updated Output Files

- `results/Figure1_composite_v1.png` — v1 composite figure (22 × 12 inches, 600 DPI)
- `results/Figure1_composite_v2.png` — v2 composite figure (26 × 13 inches, 600 DPI)
- `results/Figure1_composite_v2.pdf` — v2 composite figure, vector format (26 × 13 inches)

### Updated Input Data

The v2 R script (`scripts/Figure1_composite_v2.R`) reads the same input files as v1:

- `input/metadata_oral_qc_update.tsv` — Curated metadata (used for panels b, d, f; reads `bioproject_id`, `sample_id`, `site_cleaned`, `continent_group`, `country_cleaned`, `age_cleaned`, `sex_cleaned`, `bmi_cleaned`, `status_cleaned`, `disease_study_cleaned`)
- `input/phylum-relative-table.tsv` — Phylum-level relative abundance table (panel c inset, loaded via `data.table::fread`)
- `input/taxonomy.tsv` — Taxonomy annotations (panel e, via Python rarefaction script)
- `input/biom_oral_qc.biom` — BIOM feature table (panel e, via Python rarefaction script)
- `results/rarefaction_pcofg.tsv` — Pre-computed rarefaction/saturation data (panel e, generated by `scripts/oral_sampling_saturation_curve_analysis.py`)
