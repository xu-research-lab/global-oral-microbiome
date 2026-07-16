# 4. Beta Diversity — PCoA 可视化方法对比

## 4.1 Method

Beta diversity analyses were based on UniFrac ordinations. UniFrac PCoA coordinates and pairwise PERMANOVA tests were computed as in the original report (Section 4). Four visualization methods were compared for the same PCoA data and PERMANOVA results:

**Method A — Scatter (2-color pairwise).** Pairwise matrix panels showing Weighted (upper triangle) and Unweighted (lower triangle) UniFrac PCoA. Each off-diagonal panel displays only two groups: the row group in red (`#CC3311`) and the column group in blue (`#0072B2`). Points are rasterized via `ggrastr` for performance. Panels are annotated with PERMANOVA R² and p values. Script: `beta_pcoa_20260715/03-plot-unweight-weight-pcoa.R`.

**Method B — Density contour.** Same matrix layout, but instead of individual points, two-dimensional kernel density estimate (KDE) contours are drawn as filled polygons (`stat_density_2d`, bins=6) with semi-transparent fills (α=0.15 per band, accumulating to ~0.9 at the center). Row group in red, column group in blue. Panels with fewer than 10 total samples or 5 per group are replaced by spacers. Script: `beta_pcoa_20260715/04-plot-density-contour-pcoa.R`.

**Method C — Hexbin.** Same matrix layout, but points are replaced by hexagonal bins (`geom_hex`, bins=15). Opacity encodes log₁₀(count+1) scaled to [0.05, 0.70]. Row group in red, column group in blue. Minimum data thresholds are the same as Method B. Script: `beta_pcoa_20260715/05-plot-hexbin-pcoa.R`.

**Method D — Differential Hexbin.** Same matrix layout as Method C, but each panel shows the *difference* between the row group and column group hexbin densities. Positive differences (row > column) are shown in red and negative differences (column > row) are shown in blue. This highlights where one group concentrates relative to the other, making group-specific distributional shifts immediately visible. Script: `beta_pcoa_20260715/07-plot-diff-hexbin-pcoa.R`.

All four methods use the same input data (UniFrac PCoA coordinates dated 20260715, PERMANOVA results dated 20260715, metadata `metadata_oral_qc.tsv`) and the same seven metadata variables: `site_cleaned`, `continent_group`, `sequencing_region_cleaned`, `age_cleaned_group`, `sex_cleaned`, `platform_cleaned`, `bmi_cleaned_group`. For `site_cleaned`, the "other" group is excluded. For `sequencing_region_cleaned`, groups with fewer than 500 samples are dropped.

Pairwise PERMANOVA was recomputed on updated metadata (20260715), producing `02-unifrac-unweight-distance-matrix-permanova-20260715.tsv` and `02-unifrac-weight-distance-matrix-permanova-20260715.tsv`. The key difference from the 0604 version is that the new results include `continent_group` pairwise comparisons that were absent previously.

## 4.2 Main Figures — Site

### Method A: Scatter (2-color pairwise)

![PCoA Scatter by site](figures/04_beta_diversity/04_beta_diversity_update_v1/scatter/PCoA-Matrix-site_cleaned-20260715.png)

**Figure 9A. Beta-diversity structure by oral sampling site — Scatter matrix.**
Weighted and unweighted UniFrac PCoA coordinates plotted in pairwise 2-color scatter matrix form. Each panel shows only the row group (red) and column group (blue); all other samples are hidden. Panels are annotated with pairwise PERMANOVA R² and p values.

### Method B: Density contour

![PCoA Density by site](figures/04_beta_diversity/04_beta_diversity_update_v1/density/PCoA-Density-site_cleaned-20260715.png)

**Figure 9B. Beta-diversity structure by oral sampling site — Density contour matrix.**
Same layout as Method A, but individual points are replaced by filled KDE contour polygons (6 contour levels). The row group density is shown in red and the column group in blue; overlapping regions show additive opacity. This visualization reveals the shape and extent of each group's distribution in PCoA space without overplotting.

### Method C: Hexbin

![PCoA Hexbin by site](figures/04_beta_diversity/04_beta_diversity_update_v1/hexbin/PCoA-Hexbin-site_cleaned-20260715.png)

**Figure 9C. Beta-diversity structure by oral sampling site — Hexbin matrix.**
Same layout as Method A, but points are replaced by hexagonal bins whose opacity encodes log₁₀-transformed sample density. The row group is shown in red and the column group in blue. Hexbin eliminates overplotting entirely and allows quantitative assessment of within-group density at the cost of losing individual outlier visibility.

### Method D: Differential Hexbin

![PCoA DiffHexbin by site](figures/04_beta_diversity/04_beta_diversity_update_v1/diff_hexbin/PCoA-DiffHexbin-site_cleaned-20260715.png)

**Figure 9D. Beta-diversity structure by oral sampling site — Differential Hexbin matrix.**
Same layout as Method C, but each panel displays the difference between row-group and column-group hexbin densities. Red regions indicate where the row group has higher density; blue regions indicate where the column group dominates. This differential view directly highlights distributional shifts between groups that may be subtle in the overlaid hexbin.

## 4.3 Main Figures — Continent

### Method A: Scatter

![PCoA Scatter by continent](figures/04_beta_diversity/04_beta_diversity_update_v1/scatter/PCoA-Matrix-continent_group-20260715.png)

**Figure 10A. Beta-diversity structure by geographic region — Scatter matrix.**

### Method B: Density contour

![PCoA Density by continent](figures/04_beta_diversity/04_beta_diversity_update_v1/density/PCoA-Density-continent_group-20260715.png)

**Figure 10B. Beta-diversity structure by geographic region — Density contour matrix.**

### Method C: Hexbin

![PCoA Hexbin by continent](figures/04_beta_diversity/04_beta_diversity_update_v1/hexbin/PCoA-Hexbin-continent_group-20260715.png)

**Figure 10C. Beta-diversity structure by geographic region — Hexbin matrix.**

### Method D: Differential Hexbin

![PCoA DiffHexbin by continent](figures/04_beta_diversity/04_beta_diversity_update_v1/diff_hexbin/PCoA-DiffHexbin-continent_group-20260715.png)

**Figure 10D. Beta-diversity structure by geographic region — Differential Hexbin matrix.**

## 4.4 Main Figures — Sequencing Region

### Method A: Scatter

![PCoA Scatter by sequencing region](figures/04_beta_diversity/04_beta_diversity_update_v1/scatter/PCoA-Matrix-sequencing_region_cleaned-20260715.png)

**Figure 11A. Beta-diversity structure by amplified 16S rRNA region — Scatter matrix.**

### Method B: Density contour

![PCoA Density by sequencing region](figures/04_beta_diversity/04_beta_diversity_update_v1/density/PCoA-Density-sequencing_region_cleaned-20260715.png)

**Figure 11B. Beta-diversity structure by amplified 16S rRNA region — Density contour matrix.**

### Method C: Hexbin

![PCoA Hexbin by sequencing region](figures/04_beta_diversity/04_beta_diversity_update_v1/hexbin/PCoA-Hexbin-sequencing_region_cleaned-20260715.png)

**Figure 11C. Beta-diversity structure by amplified 16S rRNA region — Hexbin matrix.**

### Method D: Differential Hexbin

![PCoA DiffHexbin by sequencing region](figures/04_beta_diversity/04_beta_diversity_update_v1/diff_hexbin/PCoA-DiffHexbin-sequencing_region_cleaned-20260715.png)

**Figure 11D. Beta-diversity structure by amplified 16S rRNA region — Differential Hexbin matrix.**

## 4.5 Supplementary Figures

Supplementary beta-diversity figures (age group, sex, platform, BMI group) for all four methods are available at the following paths (replace `{variable}` and `{Method}` accordingly):

- Scatter: `figures/04_beta_diversity/04_beta_diversity_update_v1/scatter/PCoA-Matrix-{variable}-20260715.png`
- Density: `figures/04_beta_diversity/04_beta_diversity_update_v1/density/PCoA-Density-{variable}-20260715.png`
- Hexbin: `figures/04_beta_diversity/04_beta_diversity_update_v1/hexbin/PCoA-Hexbin-{variable}-20260715.png`
- DiffHexbin: `figures/04_beta_diversity/04_beta_diversity_update_v1/diff_hexbin/PCoA-DiffHexbin-{variable}-20260715.png`

Where `{variable}` ∈ {`age_cleaned_group`, `sex_cleaned`, `platform_cleaned`, `bmi_cleaned_group`}.

## 4.6 Reproducibility

| Item | Path |
|---|---|
| Scatter script | `beta_pcoa_20260715/03-plot-unweight-weight-pcoa.R` |
| Density script | `beta_pcoa_20260715/04-plot-density-contour-pcoa.R` |
| Hexbin script | `beta_pcoa_20260715/05-plot-hexbin-pcoa.R` |
| DiffHexbin script | `beta_pcoa_20260715/07-plot-diff-hexbin-pcoa.R` |
| PERMANOVA script (unweighted) | `beta_pcoa_20260715/02-unweight-distance-matrix-permanova.py` |
| PERMANOVA script (weighted) | `beta_pcoa_20260715/02-weight-distance-matrix-permanova.py` |
| PCoA export (unweighted) | `beta_pcoa_20260715/01-export-unweight-pcoa-hdf5-to-csv.py` |
| PCoA export (weighted) | `beta_pcoa_20260715/01-export-weight-pcoa-hdf5-to-csv.py` |
| PERMANOVA results (unweighted) | `beta_pcoa_20260715/beta_pcoa_out_20260715/unifrac_unweight/02-unifrac-unweight-distance-matrix-permanova-20260715.tsv` |
| PERMANOVA results (weighted) | `beta_pcoa_20260715/beta_pcoa_out_20260715/unifrac_weight/02-unifrac-weight-distance-matrix-permanova-20260715.tsv` |