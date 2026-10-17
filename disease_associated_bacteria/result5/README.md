# Disease-associated oral bacteria: per-disease differential abundance

Genus-level differential abundance (DA) of each oral/systemic disease versus
healthy controls, pooled across BioProjects, using three tools (MaAsLin2,
ANCOM-BC2, LinDA), plus the publication figures built from the results.

All scripts use paths relative to this folder (`result5/`). Run them from here,
or set `DA_WORK_DIR` to this folder.

## Input data (not included here)

The raw inputs are the QC'd oral dataset of this repository, in
`global_oral/data/metadata-and-biom-v2-qc/` (16S rRNA amplicon data, features
classified with Greengenes2, id99). They are not duplicated here:

| File | Source in this repository |
|---|---|
| `metadata_oral_qc_dedup.tsv` | `data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv` |
| `biom_oral_qc_dedup.biom` | `data/metadata-and-biom-v2-qc/biom_oral_qc_dedup.biom` |
| `taxonomy.tsv` | inside `data/metadata-and-biom-v2-qc/gg2_taxonomy.qza` |

See `data/metadata-and-biom-v2-qc/README.md` for how these were produced.
To rerun the analysis, put the three files in `result5/data/`:

```bash
cd global_oral/disease_associated_bacteria/result5
mkdir -p data
cp ../../data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv data/
cp ../../data/metadata-and-biom-v2-qc/biom_oral_qc_dedup.biom data/
unzip -j ../../data/metadata-and-biom-v2-qc/gg2_taxonomy.qza '*/data/taxonomy.tsv' -d data/
```

Alternatively point the scripts elsewhere: `DA_DATA_DIR` (Python scripts,
folder holding all three files) and `DA_META_PATH` (R scripts, metadata file).

The figures do **not** need the raw data; they only read the tables in
`method_comparison_three/`.

## Analysis design

- **Eligible diseases**: a BioProject is *matched* for a disease when it has
  both cases and healthy controls (HC). A disease is analysed when, over its
  matched BioProjects, cases ≥ 30, HC ≥ 30 and matched BioProjects ≥ 2.
- **Genus tables**: samples with < 1,000 reads are dropped, features are
  aggregated to genus (`g__`) level, and a genus is kept if its prevalence is
  ≥ 1% in at least one eligible disease cohort. GTDB-style genus suffixes
  (e.g. `Prevotella_A_…`) are kept as separate genera.
- **Models** (disease vs HC, reference = HC, BioProject as random effect).
  Covariates age group, sex, sampling site, sequencing region and continent
  are added when they have ≥ 2 levels in the cohort (missing values coded as
  `Unknown`). If a model fails to fit, LinDA and ANCOM-BC2 refit with
  `Disease_Status` alone (recorded in `fix_formula`).
  - MaAsLin2: TSS + LOG transform, LM, random effect `bioproject_id`.
  - ANCOM-BC2: counts, `rand_formula = (1|bioproject_id)`.
  - LinDA: counts, `(1|bioproject_id)` in the formula.
  - Significance: BH-adjusted q ≤ 0.05.
- **Consensus genus**: significant (q ≤ 0.05) in at least two of the three
  tools with the same direction (845 disease–genus calls). `+` /
  disease-enriched (red, `#B2182B`) = higher in disease; `-` /
  health-enriched (blue, `#2166AC`) = higher in HC.

## Reproducing

```bash
# 1. preprocessing (Python: pandas, numpy, scipy, biom-format)
python3 script/compute_eligible_diseases.py      # eligible_diseases.txt / .csv
python3 script/prepare_genus_tables.py           # genus_table_counts.tsv, genus_table_relabund.tsv
python3 script/write_metadata_after_cleaning.py  # metadata_after_cleaning.tsv (cleaned metadata record)

# 2. DA per tool: without arguments every eligible disease is run and merged
#    into <tool>/<tool>_all_results.tsv
Rscript script/run_maaslin2_disease.R
Rscript script/run_linda_disease.R
Rscript script/run_ancombc2_disease.R
# or one disease per call (e.g. as cluster jobs), then merge
Rscript script/run_maaslin2_disease.R "Periodontal Disease"
Rscript script/run_maaslin2_disease.R --merge-only
# 3. figures (written to figure/)
Rscript script/panel_A_sankey.R
Rscript script/panel_b_total.R
Rscript script/panel_b_venn.R
Rscript script/panel_s5a_by_methods.R
Rscript script/panel_s5b_upset.R
Rscript script/panel_s5c_pan_disease.R
```

Step 2 writes per-disease folders to `<tool>/disease_results/`. The figure
scripts also write some extra figure variants and tables
(`method_comparison_three/consensus_bubble/`) that are not part of this
release.

## Files

```
script/compute_eligible_diseases.py      eligible diseases (cases, HC, matched BioProjects)
script/prepare_genus_tables.py           genus count and relative-abundance tables
script/write_metadata_after_cleaning.py  cleaned metadata record
script/run_maaslin2_disease.R            MaAsLin2 per disease -> maaslin2/
script/run_linda_disease.R               LinDA per disease -> linda/
script/run_ancombc2_disease.R            ANCOM-BC2 per disease -> ancombc2/
script/pub_theme.R                       shared paths, fonts, sizes, disease abbreviations, save_panel()
script/panel_A_sankey.R                  cohort Sankey -> panel_A_sankey_by_disease_abbrev
script/panel_b_total.R                   consensus pairs per disease, shared vs disease-specific -> panel_B_total
script/panel_b_venn.R                    three-tool Venn of all significant calls -> panel_b_venn
script/panel_s5a_by_methods.R            per-disease tool-overlap counts -> panel_s5a_by_methods
script/panel_s5b_upset.R                 UpSet of all shared consensus intersections -> panel_s5b_upset_all
script/panel_s5c_pan_disease.R           pan-disease genus heatmap -> panel_s5c_pan_disease

maaslin2/maaslin2_all_results.tsv.gz   all MaAsLin2 coefficients (all diseases)
linda/linda_all_results.tsv.gz         all LinDA results
ancombc2/ancombc2_all_results.tsv.gz   all ANCOM-BC2 results

method_comparison_three/
  consensus_genera.csv                 consensus calls: Disease, Genus, Direction, n_methods, methods
  overlap_summary_3way.csv             per-disease hit counts, overlaps and Jaccard indices of the three tools
  cohort_overview/cohort_overview_data.csv    per (disease, BioProject): category, continent, country, 16S region, n_case, n_hc
  cohort_overview/disease_category_map.csv    disease -> category and display name
  cohort_overview/site_sample_counts.csv      per (disease, BioProject, sampling site): number of case + HC samples

figure/                              the six figures, PNG (600 dpi) and PDF
```

Result-table notes:

- MaAsLin2 rows include covariates; the disease effect is
  `metadata == "Disease_Status"` (`coef` > 0 means higher in disease).
- LinDA / ANCOM-BC2: one row per disease × genus; `lfc` is disease vs HC,
  `Direction` is its sign, `fix_formula` / `rand_formula` give the fitted model.
- `Genus` is the cleaned genus name used for matching across tools.

## Software

R 4.3.1 with Maaslin2 1.16.0, ANCOMBC 2.4.0, MicrobiomeStat 1.4 (LinDA),
phyloseq 1.46.0, dplyr 1.1.4, tidyr 1.3.1, ggplot2 4.0.1, ggforce 0.4.2,
ggtext 0.1.2, gridtext 0.1.5, ComplexUpset 1.3.3, patchwork 1.3.0,
ggrepel 0.9.6, scales 1.4.0, png 0.1.8.

Python 3.9 with pandas 2.2.3, numpy 1.22.4, scipy 1.8.1, biom-format 2.1.10.
