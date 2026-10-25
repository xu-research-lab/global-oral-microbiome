# HOMD comprehensiveness of the global oral 16S ASV catalogue

This directory contains the analysis behind Figure 1, the supplementary figure,
the disease figures and the two supplementary tables. It asks two questions of a
global collection of oral 16S ASV tables:

1. how much of the catalogue eHOMD represents at species level, and
2. which lineages fall below that boundary and how they are distributed.

The final phylotype catalogue, the sample-level models and the figures are in
[`pipeline_v3/`](pipeline_v3/); that directory has its own
[README](pipeline_v3/README.md) and a per-file
[sync manifest](pipeline_v3/analysis_report/sync_manifest.md) that records the
exact script → table → figure chain and the excluded superseded versions.

## Terminology

- **ASV** — an amplicon sequence variant called by DADA2 within one project.
- **derep sequence** — an exact-dereplicated sequence pooled over all projects;
  `derep_id` (`DEREP_<n>`) names it.
- **phylotype** — a TreeCluster `max_clade` unit at t = 0.03 on the GG2 tree,
  made of placed ASVs and the reference tips that fall inside the same clade.
  The final id is `c<number>` for a phylotype of two or more ASVs, `s<number>`
  for a single ASV, and a bare number for a single-ASV phylotype that contains
  reference sequences. The pipeline still calls these units *clusters*.
- **eHOMD-captured / uncaptured** — captured means the ASV has an eHOMD best hit
  at ≥ 97 % identity over ≥ 99 % of the query; uncaptured means it does not.
  At phylotype level, a phylotype is captured when at least one member is
  annotated to a species that is in eHOMD.

## Layout

```
HOMD_uncapture_analysis/
├── README.md                        this file
├── process_data.sh                  per-project DADA2, batch template
├── process_forwards.R               per-project DADA2, single-end R body
├── 01_collect_asvs_and_abundance.py pool the per-project ASV tables, exact dereplication
├── 02_homd_mapping.sbatch           vsearch against eHOMD 16.03, HOMD-captured flag
├── 03_chimera_removal.sbatch        uchime_ref against eHOMD + GG2 + SILVA
├── C01 … C09                        control samples → contaminant catalogue
├── pipeline_v2/                     candidates, contamination firewall, support floor, SEPP
├── pipeline_v3/                     phylotypes, HOMD coverage, models, figures, tables
├── exact_derep_summary.tsv          derep catalogue
├── exact_derep_representatives.fa   derep sequences
├── asv_abundance_metadata.tsv       per-ASV abundance and metadata
├── homd_mapping/                    per-derep eHOMD annotation
├── chimera_removal/                 per-derep chimera calls
├── contaminant_catalog/             control-derived contaminant ASVs
└── control_quality_viz/             control QC outputs
```

## Order of execution

The two chains run in parallel until the contamination firewall: the ASV chain
(01–03) builds the pooled catalogue, and the control chain (C01–C09) builds the
contaminant catalogue that the firewall uses.

```
# per project, on the raw reads (one job per BioProject)
sbatch process_data.sh                       # calls process_forwards.R

# pooled ASV catalogue
python 01_collect_asvs_and_abundance.py
sbatch 02_homd_mapping.sbatch
sbatch 03_chimera_removal.sbatch

# control samples and contaminant catalogue
python C01_find_control_samples.py
python C02_extract_control_sample_ids.py
sbatch C04_convert_sra_to_fastq.sbatch
sbatch C05_run_dada2_controls.sbatch
sbatch C05b_dada2_PRJEB11419.sbatch
sbatch C06_dada2_extra_controls.sbatch
sbatch C07_control_quality_viz.sbatch
sbatch C08_filter_genuine_controls.sbatch
python C08b_filter_genuine_controls_v2.py    # the thresholds actually used
sbatch C09_build_contaminant_catalog.sbatch

# candidate set: uncaptured, non-chimeric, clean, supported
python pipeline_v2/01_chimera_filter.py
sbatch pipeline_v2/02_contamination_firewall.sbatch
python pipeline_v2/04_support_floor.py
sbatch pipeline_v2/04b_sepp_place.sbatch
bash   pipeline_v2/04c_sepp_missing.sh
python pipeline_v2/01b_subject_prevalence.py

# phylotypes, coverage, models, tables, figures (pipeline_v3)
# exact order in pipeline_v3/analysis_report/sync_manifest.md
```

`01b_subject_prevalence.py` reads the support-floor table, so it runs after
step 04 even though it keeps its original `01b` name.

## Per-project DADA2

| File | What it does |
|---|---|
| `process_data.sh` | SLURM array template: one BioProject per task, downloads the SRA runs, converts to FASTQ, calls `process_forwards.R`, writes `results/ASVs.fa` and `results/ASVs_counts.tsv`. The array index and the project list are environment-specific and must be edited per batch. |
| `process_forwards.R` | DADA2 body for one project: quality filtering, error learning, denoising, merger, chimera removal, ASV table. Takes the project directory as its first argument and `trim` as an optional second. |

## Pooled ASV catalogue

| File | What it does | Writes |
|---|---|---|
| `01_collect_asvs_and_abundance.py` | Reads every project's `ASVs.fa` and `ASVs_counts.tsv`, records per-ASV abundance, prevalence and amplicon region, and exact-dereplicates the sequences across projects. | `exact_derep_summary.tsv`, `exact_derep_representatives.fa`, `asv_abundance_metadata.tsv`, `project_input_summary.tsv`, `collect_stats.json` |
| `02_homd_mapping.sbatch` | `vsearch --usearch_global` of every derep sequence against eHOMD 16.03 (`--iddef 1`, best hit). Flags a sequence as HOMD-captured only when identity ≥ 97 % **and** query coverage ≥ 99 %. | `homd_mapping/exact_derep_homd_annotated.tsv`, `homd_mapping/exact_derep_vs_homd.hits.tsv`, `homd_mapping/homd_mapping_stats.json` |
| `03_chimera_removal.sbatch` | Builds the combined eHOMD + GG2 + SILVA 138.2 reference and runs `vsearch --uchime_ref` on all derep sequences; merges the chimera calls into the eHOMD annotation. | `chimera_removal/exact_derep_chimera_annotated.tsv`, `chimera_removal/chimeras.fa`, `chimera_removal/nonchimeras.fa`, `chimera_removal/chimera_stats.json` |

Counts: 2,844,174 derep sequences, 2,099,459 non-chimeric.

## Control samples and the contaminant catalogue

The `C` prefix separates this chain from the ASV chain; both originally used the
numbers 01, 02, 04… There is no `C03`: the gap is kept so that the script numbers
still match the order in which the intermediates were produced.

| File | What it does | Writes |
|---|---|---|
| `C01_find_control_samples.py` | Scans the `SraRunTable.csv` of every project for blank, negative-control, NTC and mock samples. | `control_samples_detailed.tsv`, `control_samples_summary.tsv` |
| `C02_extract_control_sample_ids.py` | Extracts the actual Run/sample ids of those controls. | `control_sample_ids.tsv`, `control_samples_by_category.tsv` |
| `C04_convert_sra_to_fastq.sbatch` | Downloads the control runs and converts them to FASTQ with `fastq-dump`. Mock communities are excluded. | `control_samples/fastq/…` |
| `C05_run_dada2_controls.sbatch` | DADA2 per BioProject for the control samples. | `control_dada2/<project>/` |
| `C05b_dada2_PRJEB11419.sbatch` | DADA2 for PRJEB11419, which contributes 921 controls and is run on its own. | `control_dada2/PRJEB11419/` |
| `C06_dada2_extra_controls.sbatch` | DADA2 for the negative-control and reagent-water samples that sit outside the main projects. | `control_dada2_extra/` |
| `C07_control_quality_viz.sbatch` | Per-control ASV richness and genus composition; assigns the SILVA genus of each control ASV. | `control_quality_viz/` |
| `C08_filter_genuine_controls.sbatch` | First authenticity screen (v1): keeps controls with `< 30` ASVs and `< 20 %` oral reads. | `control_quality_viz/genuine_controls_keep.csv` |
| `C08b_filter_genuine_controls_v2.py` | The screen actually used (v2): read-count, ASV-count and oral-fraction thresholds, per-study requirements, and a per-sample exclusion reason. | `control_quality_viz/control_classification_v2.csv`, `control_quality_viz/genuine_controls_keep_v2.csv`, `control_quality_viz/05_control_qc_v2.png` |
| `C09_build_contaminant_catalog.sbatch` | Builds the cross-lab contaminant ASV catalogue from the genuine controls, with taxonomy and prevalence. | `contaminant_catalog/contaminant_asv_catalog.tsv`, `contaminant_catalog/contaminant_catalog.fasta`, `contaminant_catalog/catalog_summary.json` |

## Candidate set (`pipeline_v2/`)

| File | What it does | Writes |
|---|---|---|
| `01_chimera_filter.py` | Keeps the sequences that are both HOMD-uncaptured and non-chimeric. | `01_chimera_filtered.tsv` (1,136,142) |
| `02_contamination_firewall.py` (+ `.sbatch`) | Stage 1 clusters the candidates with every control ASV at 97 % identity and flags candidates that share a cluster with a control; Stage 2 annotates the survivors with a GG2 genus and flags genera in `contaminant_genera.tsv`. | `02_contamination_flagged.tsv`, `02_summary.json` |
| `contaminant_genera.tsv` | The genus list used by Stage 2: 560 control-derived genera plus 101 genera from Salter 2014 and Weyrich 2019. Curated by hand. | — |
| `04_support_floor.py` | Applies the support floor (≥ 2 BioProjects, ≥ 50 reads) without dropping contamination-flagged ASVs, so the whole phylotype can be adjudicated later. Writes the SEPP query and reuses the placements already available. | `04_support_floor_passed.tsv` (54,290), `candidates_for_sepp.fasta`, `candidates_missing_sepp.fasta`, `sepp_placements_reused.tsv`, `04_summary.json` |
| `04b_sepp_place.sbatch` | Records the exact QIIME2 `fragment-insertion sepp` command used for the placement run. The insertion was run by hand, so this file is the record of that command rather than a job that was submitted unchanged. | `sepp/` |
| `04c_sepp_missing.sh` | SEPP for the candidates whose placements could not be reused. | `sepp_missing/placements_export/placements.json` |
| `01b_subject_prevalence.py` | Distinct-subject prevalence for each support-floor sequence, using `(bioproject, subject_id)` as the subject identity. | `01b_derep_subject_prevalence.tsv`, `01b_subject_stats.json` |

The placed ASVs were grafted onto the GG2 backbone with `gappa graft`; that
command was run by hand and is not recorded in a script. Its output is the tree
read by `pipeline_v3/01_treecluster_cluster.py`.

## Phylotypes, coverage and models (`pipeline_v3/`)

See [`pipeline_v3/README.md`](pipeline_v3/README.md) for the analysis itself and
[`pipeline_v3/analysis_report/sync_manifest.md`](pipeline_v3/analysis_report/sync_manifest.md)
for the script → table → figure map. In brief:

Scripts in `pipeline_v3/` itself:

| File | What it does | Writes |
|---|---|---|
| `01_treecluster_cluster.py` (+ `.sbatch`) | TreeCluster `max_clade` at 0.03; per-ASV and per-phylotype taxonomy from the reference tips in the clade; eHOMD relationship. | `clusters.tsv`, `asv_clusters.tsv`, `summary.json` |
| `02_cluster_aggregate.py` | Aggregates members to phylotypes, adds support and contamination status, defines `clean_cluster`. | `02_clusters_aggregated.tsv`, `02_summary.json` |
| `05_cluster_geography.py` | Adds bioproject, continent and phylum summaries; builds the candidate table. | `05_candidate_novel_with_phylum.tsv`, `05_summary.json` |
| `07_master_table.py` | Per-ASV master table: sequence, phylotype, eHOMD hit, distances, support. | `07_master_asv_table.tsv`, `07_summary.json` |
| `07b_add_organelle_flags.py` | Adds mitochondrion / chloroplast flags from `organelle_mapping.tsv`. | organelle columns in `02_*` and `07_*` |
| `08_supported_clusters.py` | The support gate: > 50 subjects, ≥ 2 BioProjects, ≥ 2 continents, no flagged and no organelle member. | `08_supported_clusters.tsv`, `08_supported_summary.json` |
| `12d_build_anchors_v3.py` | Plot input for the supported-catalogue ASV profile. | `12d_anchors_supported_v3.tsv` |
| `12d_build_headline_profile_v3.py` | One row per member ASV of a supported phylotype, with phylotype aggregates. | `12d_headline_profile_v3.tsv` |
| `14_cluster_feature_profile.py` | Phylotype × sample count and relative-abundance matrices from the projects' `DADA2.biom`. | `14_cluster_counts_dedup.tsv`, `14_cluster_relabund_dedup.tsv`, `14_cluster_sample_metadata.tsv`, `14_supported_cluster_annotation.tsv` |
| `15_cluster_cohort_detection.py` | Detection of every supported phylotype in each disease cohort. | `15_cluster_cohort_detection.tsv` |
| `16_cluster_differential_wilcox.R` | Wilcoxon screen that fixes the cohorts and features carried into MaAsLin3. | `16_cluster_differential_wilcox.tsv` |
| `17_maaslin3_disease.R` | MaAsLin3 abundance and prevalence models per cohort. | `17_maaslin3_all_results.tsv`, `17_maaslin3_significant.tsv`, `17_maaslin3_status.tsv`, `maaslin3/` |
| `20_rename_cluster_ids.py` | Turns the numeric TreeCluster ids into `c`/`s` phylotype ids in every table; run after any regeneration. | `cluster_id_rename_map.tsv` and rewritten ids |
| `15_plot_cluster_cohort_detection.R` | Detection-across-cohorts figure. | `fig_cluster_disease_cohort_prevalence*` |

Scripts in `pipeline_v3/scripts/`:

| File | What it does |
|---|---|
| `13b_extract_fasta.py` | Extracts the uncaptured ASV sequences for the organelle search. |
| `13b_organelle_map.sbatch` | `vsearch` of those sequences against the GG2 organelle reference. |
| `13c_organelle_annotate.py` | Adds the organelle flags to `pipeline_v2/02_contamination_flagged.tsv` in place. |
| `13d_homd_coverage_clean.py` | Per-phylum and per-species coverage from the clean ASVs. |
| `13r_homd_lt_fractions.py` | Fraction of ASVs / reads / subjects below each identity threshold, per phylum. |
| `13z2_homd_phylum_band_composition.py` | Phylum composition of the < 97 % and ≥ 97 % bands. |
| `13z3_homd_genus_identity_curves.py` | Per-genus identity curves. |
| `12j_plot_supported_asv_profiles_v3.R` | Figure 1 panel D and supplementary Figure S1. |
| `13z_plot_homd_phylum_coverage_panel_v3.R` | Figure 1 panel B; sources `13z2_panel_homd_band_composition.R`. |
| `13z2_panel_homd_band_composition.R` | Lower half of panel B. |
| `13z3_plot_homd_genus_identity_curves_numbered.R` | Figure 1 panel C. |
| `14_schematic_overview.R` | Figure 1 panel A when run with `FLOW_ONLY`; the panel A counts are hard-coded in its `CARDS` table. |
| `15_alt_layout.R` | Assembles Figure 1 from panels A–D. |
| `17b_plot_cluster_prevalence.R` | One-phylotype prevalence figure. |
| `17c_plot_significant_clusters.R` | Prevalence and abundance figures for the significant phylotypes. |
| `17d_capnocytophaga_control.py`, `17d_plot_capnocytophaga_control.R` | Control figure for the covered part of *Capnocytophaga*. |
| `17e_21907_cohort_sensitivity_fit.R`, `17e_plot_21907_cohort_sensitivity.R` | Cohort-definition sensitivity for `c21907`. |
| `17f_ctrl_only_fit.R`, `17f_plot_ctrl_only.R` | Re-fit of the significant phylotypes on cohorts with both arms. |
| `18_build_disease_feature_table.py` | Supplementary table of the phylotypes differential in head and neck cancer and periodontal disease. |
| `19_build_supported_phylotype_table.py` | Supplementary table of the 486 supported phylotypes. |

Figures land in `pipeline_v3/analysis_report/figs/` (PDF and PNG of each
panel) and the two supplementary tables in `pipeline_v3/analysis_report/tables/`.
`pipeline_v3/analysis_report/` also holds the Methods/Results draft
(`HOMD_comprehensiveness_methods_results.md`) and the manifest.

## Inputs not in this repository

| Input | Location | Used by |
|---|---|---|
| SEPP tree (GG2 backbone + 54,277 placed ASVs) | `pipeline/05_sepp/tree_export/tree.nwk` | `pipeline_v3/01` |
| GG2 2022.10 tree taxonomy | `pipeline_v2/gappa/gg2_2022.10_tree_taxonomy.tsv` | `pipeline_v3/01`, `05` |
| eHOMD 16.03 reference and taxonomy | `/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03*` | `02`, `03`, `13d` |
| GG2 and SILVA reference databases | `/beegfs/dongbiao/greengene2/…`, `/beegfs/dongbiao/silva_138.2/…` | `03`, `C07`–`C09` |
| per-project `DADA2.biom` | `/beegfs/dongbiao/oral/<project>/results/` | `14` |
| sample metadata | `data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv` | `14` |
| eHOMD annotation and chimera tables (1.3 GB and 1.5 GB) | `homd_mapping/`, `chimera_removal/` | `pipeline_v3/01`, `13*` |
| contamination-flagged table (630 MB) | `pipeline_v2/02_contamination_flagged.tsv` | `pipeline_v3/13*` |

Files over GitHub's 100 MB limit are not tracked; the scripts that generate them
are.

## Not part of the final analysis

These remain on disk for provenance but are superseded and should not be read as
part of the pipeline: `10_decontam_within_studies.sbatch`, `cluster97.sbatch`,
`flag_contaminants.sbatch`, `genus_decontam.sbatch`, `ggtree.sbatch`,
`ggtree_rect.sbatch`, `plot_figures.py`, `plot_phylotype_tree.R`,
`update_integrate.sbatch`; the July `pipeline/` chain; and in `pipeline_v2/` the
alternative phylotype route (`03_self_cluster*`, `05*`, `06*`, `07*`,
`08_phylotype_geography.py`, `09*`, `10*`, `11*`, `12*`, `13_build_master_table.py`,
`14_supported_phylotypes.py`, `run_02_contam.sh`) plus its `analysis_report/`.

Superseded in `pipeline_v3/`: `06_build_cluster_tree.py`,
`09_make_summary_figs.py`, `09b_local_tree_cluster_examples.py`, `12e`–`12i`,
`15a_alt_layout_3row.R`, `16_plot_cluster_differential_boxplot.R`,
`99_combine_master_figs.R`, `99b_preview_panelB_swap.R`,
`build_isolation_source_table.py`, `build_report_notebook_v3.py`, and in
`pipeline_v3/scripts/` every HOMD-coverage variant not listed above
(`13a`, `13e`–`13q`, `13s`–`13z4` except `13z2`/`13z3`/`13z`).
`pipeline_v3/analysis_report/sync_manifest.md` lists the same set in full.

## Environment

R at `/home/software/anaconda3/envs/renv/bin/Rscript`; Python at the matching
`python` (pandas, numpy, openpyxl, biom). DADA2 and SEPP run in QIIME 2
(`qiime2-2022.2`, `rachis-qiime2-2026.4`), VSEARCH 2.7.0 for the chimera screen
and 2.30.6 for the mapping and clustering steps, TreeCluster 1.0.5, MaAsLin3.
`pipeline_v3/scripts/15_alt_layout.R` also needs ImageMagick and `pdflatex`.
