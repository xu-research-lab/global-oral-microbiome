# Sync manifest — `HOMD_uncapture_analysis/pipeline_v3`

Target: <https://github.com/sean11sudo/global-oral-microbiome/tree/main/HOMD_uncapture_analysis>

This is the final version of everything behind the figures and the two
supplementary tables. Exploratory scripts and superseded figure versions are
listed at the end so that the exclusion is deliberate, not accidental.

All paths are relative to `pipeline_v3/`.

---

## 1. Figures

### Figure 1 — composite

| Panel | Files | Built by |
|---|---|---|
| A, data-set flow | `fig_v3_flow_only.{pdf,png}` | `scripts/14_schematic_overview.R`, run with `FLOW_ONLY=<path>/fig_v3_flow_only` (the script writes only this block and stops) |
| B, HOMD phylum coverage | `fig_v3_panelB_coverage_v4.{pdf,png}` | `scripts/13z_plot_homd_phylum_coverage_panel_v3.R`, which sources `scripts/13z2_panel_homd_band_composition.R` |
| C, genus identity curves | `fig_v3_homd_genus_identity_curves_numbered.{pdf,png}` | `scripts/13z3_plot_homd_genus_identity_curves_numbered.R` |
| D, multi-ASV profile | `fig_v3_supported_multi_asv_profile.{pdf,png}` | `scripts/12j_plot_supported_asv_profiles_v3.R` |
| **composite** | `fig_v3_master_alt.{pdf,png}` | `scripts/15_alt_layout.R` (needs ImageMagick `convert` and `pdflatex`) |

### Supplementary figure

| Figure | Files | Built by |
|---|---|---|
| S1, singleton top-40 profile | `fig_v3_supported_singleton_top40_profile.{pdf,png}` | `scripts/12j_plot_supported_asv_profiles_v3.R` |

`12j_plot_supported_asv_profiles_v3.R` also writes
`fig_v3_supported_asv_profiles_row.{pdf,png}` (panel D and S1 side by side).
It is a by-product, not used in the manuscript — decide whether to keep it in
the repository.

### Disease figures

| Figure | Files | Built by |
|---|---|---|
| detection across cohorts | `fig_cluster_disease_cohort_prevalence.{pdf,png}`, `..._top30.{pdf,png}` | `15_plot_cluster_cohort_detection.R` |
| significant phylotypes, prevalence | `fig_{colorectal_cancer,head_and_neck_cancer,lung_cancer,periodontal_disease}_prevalence_sig.{pdf,png}` | `scripts/17c_plot_significant_clusters.R` |
| significant phylotypes, abundance | `fig_head_and_neck_cancer_abundance_sig.{pdf,png}` | `scripts/17c_plot_significant_clusters.R` |
| single phylotype, OSCC | `fig_oscc_21907_prevalence.{pdf,png}` | `scripts/17b_plot_cluster_prevalence.R` |
| Capnocytophaga coverage control | `fig_oscc_capnocytophaga_covered_prevalence.{pdf,png}` | `scripts/17d_plot_capnocytophaga_control.R` |
| cohort-definition sensitivity | `fig_oscc_21907_cohort_sensitivity.{pdf,png}` | `scripts/17e_plot_21907_cohort_sensitivity.R` |
| control-only check | `fig_{head_and_neck_cancer,periodontal_disease}_ctrl_only.{pdf,png}` | `scripts/17f_plot_ctrl_only.R` |

---

## 2. Scripts

### Upstream chain (produces the tables the figures read)

| Script | Reads | Writes |
|---|---|---|
| `20_rename_cluster_ids.py` | `clusters.tsv`, `asv_clusters.tsv`, `02_clusters_aggregated.tsv`, `05_candidate_novel_with_phylum.tsv`, `07_master_asv_table.tsv`, `08_supported_clusters.tsv`, `12d_*`, `14_*`, `15_*`, `16_*` | rename map and rewritten ids |
| `01_treecluster_cluster.py` | `pipeline/05_sepp/tree_export/tree.nwk`, `pipeline_v2/gappa/gg2_2022.10_tree_taxonomy.tsv`, eHOMD fasta | `clusters.tsv`, `asv_clusters.tsv`, `summary.json` |
| `02_cluster_aggregate.py` | `asv_clusters.tsv`, `pipeline_v2/04_support_floor_passed.tsv`, `pipeline_v2/01b_derep_subject_prevalence.tsv`, `pipeline_v2/01_chimera_filtered.tsv` | `02_clusters_aggregated.tsv`, `02_summary.json` |
| `05_cluster_geography.py` | `02_clusters_aggregated.tsv`, `asv_clusters.tsv`, GG2 taxonomy | `05_candidate_novel_with_phylum.tsv`, `05_summary.json` |
| `07_master_table.py` | `asv_clusters.tsv`, `02_clusters_aggregated.tsv`, `pipeline_v2/04_support_floor_passed.tsv`, `pipeline_v2/01b_derep_subject_prevalence.tsv`, `pipeline_v2/01_chimera_filtered.tsv` | `07_master_asv_table.tsv`, `07_summary.json` |
| `07b_add_organelle_flags.py` | `organelle_mapping.tsv`, `02_clusters_aggregated.tsv`, `07_master_asv_table.tsv` | organelle columns in the latter two |
| `08_supported_clusters.py` | `02_clusters_aggregated.tsv`, `asv_clusters.tsv` | `08_supported_clusters.tsv`, `08_supported_summary.json` |
| `12d_build_anchors_v3.py` | `08_supported_clusters.tsv`, `05_candidate_novel_with_phylum.tsv`, eHOMD fasta | `12d_anchors_supported_v3.tsv` |
| `12d_build_headline_profile_v3.py` | `08_supported_clusters.tsv`, `05_candidate_novel_with_phylum.tsv` | `12d_headline_profile_v3.tsv` |
| `14_cluster_feature_profile.py` | `07_master_asv_table.tsv`, `08_supported_clusters.tsv`, `data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv`, each project's `DADA2.biom` | `14_cluster_counts_dedup.tsv`, `14_cluster_relabund_dedup.tsv`, `14_cluster_sample_metadata.tsv`, `14_supported_cluster_annotation.tsv`, `14_summary.json` |
| `15_cluster_cohort_detection.py` | `14_cluster_counts_dedup.tsv`, `14_cluster_sample_metadata.tsv` | `15_cluster_cohort_detection.tsv` |
| `16_cluster_differential_wilcox.R` | `14_cluster_relabund_dedup.tsv`, `14_cluster_sample_metadata.tsv` | `16_cluster_differential_wilcox.tsv` |
| `17_maaslin3_disease.R` | `14_cluster_relabund_dedup.tsv`, `14_cluster_sample_metadata.tsv`, `16_cluster_differential_wilcox.tsv` | `17_maaslin3_all_results.tsv`, `17_maaslin3_significant.tsv` (`qval_individual < 0.05`, sorted by it, with `direction`), `17_maaslin3_status.tsv`, `maaslin3/` |

### HOMD coverage chain (panels B and C)

| Script | Reads | Writes |
|---|---|---|
| `scripts/13b_extract_fasta.py` | `homd_mapping/exact_derep_homd_annotated.tsv` | `13b_uncaptured_asvs.fasta` |
| `scripts/13b_organelle_map.sbatch` | that fasta, `gg2_organelle_ref.fasta` | `13b_organelle.hits.tsv` |
| `scripts/13c_organelle_annotate.py` | hits, GG2 taxonomy, `pipeline_v2/02_contamination_flagged.tsv` | adds `organelle_flag` / `organelle_target` columns to that table **in place** |
| `scripts/13d_homd_coverage_clean.py` | eHOMD annotated table, chimera table, contamination table, HOMD taxonomy | `13d_homd_phylum_coverage.tsv` |
| `scripts/13r_homd_lt_fractions.py` | same inputs | `13r_homd_phylum_lt_fractions.tsv`, `13r_homd_total_lt_fractions.tsv` |
| `scripts/13z2_homd_phylum_band_composition.py` | same inputs | `13z2_homd_phylum_band_composition.tsv` |
| `scripts/13z3_homd_genus_identity_curves.py` | same inputs, HOMD taxonomy | `13z3_homd_genus_identity_curves.tsv` |

### Figure and table scripts

`scripts/12j_plot_supported_asv_profiles_v3.R`,
`scripts/13z_plot_homd_phylum_coverage_panel_v3.R`,
`scripts/13z2_panel_homd_band_composition.R`,
`scripts/13z3_plot_homd_genus_identity_curves_numbered.R`,
`scripts/14_schematic_overview.R`,
`scripts/15_alt_layout.R`,
`15_plot_cluster_cohort_detection.R`,
`scripts/17b_plot_cluster_prevalence.R`,
`scripts/17c_plot_significant_clusters.R`,
`scripts/17d_capnocytophaga_control.py`,
`scripts/17d_plot_capnocytophaga_control.R`,
`scripts/17e_21907_cohort_sensitivity_fit.R`,
`scripts/17e_plot_21907_cohort_sensitivity.R`,
`scripts/17f_ctrl_only_fit.R`,
`scripts/17f_plot_ctrl_only.R`,
`scripts/18_build_disease_feature_table.py`,
`scripts/19_build_supported_phylotype_table.py`.

---

## 3. Generated data to push

All are below GitHub's 100 MB per-file limit; total 207 MB.

| File | Size |
|---|---|
| `clusters.tsv` | 1.8 MB |
| `asv_clusters.tsv` | 18.2 MB |
| `cluster_id_rename_map.tsv` | 1.8 MB |
| `02_clusters_aggregated.tsv` | 13.8 MB |
| `05_candidate_novel_with_phylum.tsv` | 3.1 MB |
| `07_master_asv_table.tsv` | 42.5 MB |
| `08_supported_clusters.tsv` | 0.2 MB |
| `12d_anchors_supported_v3.tsv` | 0.04 MB |
| `12d_headline_profile_v3.tsv` | 0.4 MB |
| `13d_homd_phylum_coverage.tsv` | 0.01 MB |
| `13r_homd_phylum_lt_fractions.tsv` | 0.7 MB |
| `13r_homd_total_lt_fractions.tsv` | 0.04 MB |
| `13z2_homd_phylum_band_composition.tsv` | 0.01 MB |
| `13z3_homd_genus_identity_curves.tsv` | 0.5 MB |
| `14_cluster_counts_dedup.tsv` | 48.6 MB |
| `14_cluster_relabund_dedup.tsv` | 49.3 MB |
| `14_cluster_sample_metadata.tsv` | 22.2 MB |
| `14_supported_cluster_annotation.tsv` | 0.1 MB |
| `15_cluster_cohort_detection.tsv` | 0.7 MB |
| `16_cluster_differential_wilcox.tsv` | 0.6 MB |
| `17_maaslin3_all_results.tsv` | 0.1 MB |
| `17_maaslin3_significant.tsv` | 0.01 MB |
| `17_maaslin3_status.tsv` | 0.01 MB |
| `17d_capnocytophaga_control_{samples,persample,asvs}.tsv` | 1.2 MB |
| `17e_21907_cohort_sensitivity_{results,samples}.tsv` | 0.06 MB |
| `17f_ctrl_only_{results,samples}.tsv` | 0.3 MB |
| `summary.json`, `02_summary.json`, `05_summary.json`, `07_summary.json`, `08_supported_summary.json`, `14_summary.json` | < 0.01 MB |

Also push `gg2_organelle_ref.fasta` (9.4 MB) and `organelle_mapping.tsv`
(0.1 MB): the first is an input of the organelle step, the second is what
`07b` reads.

Not pushed: `13b_uncaptured_asvs.fasta` (296 MB, over the limit). It is
regenerated by `scripts/13b_extract_fasta.py`.

---

## 4. Supplementary tables

| File | Rows | Built by |
|---|---|---|
| `analysis_report/tables/TableS_supported_phylotypes.tsv` / `.xlsx` | 486 | `scripts/19_build_supported_phylotype_table.py` |
| `analysis_report/tables/TableS_differential_phylotypes_HNC_PD.tsv` | 49 | `scripts/18_build_disease_feature_table.py` |
| `analysis_report/tables/TableS_differential_phylotypes_HNC_PD_members.tsv` | 200 | same |
| `analysis_report/tables/TableS_differential_phylotypes_HNC_PD.xlsx` | — | same |
| `analysis_report/tables/README.md` | — | — |

---

## 5. Documentation

| File | Status |
|---|---|
| `analysis_report/sync_manifest.md` | this file |
| `analysis_report/tables/README.md` | current |
| `analysis_report/HOMD_comprehensiveness_methods_results.md` | **needs a decision before sync**: the panel C results subsection and legend still describe the species quadrant figure, which panel C no longer is. The panel B and D text matches the current figures. |
| `README.md` (pipeline_v3 root) | rewritten: layout, inputs, steps, outputs, environment |

---

## 6. Excluded — superseded scripts

Delete or leave behind; do not push.

- Backup copies: `01_treecluster_cluster.py.bak_*`, `02_cluster_aggregate.py.bak*`,
  `12d_build_headline_profile_v3.py.bak*`,
  `12f_plot_headline_prevalence_profile_v3.R.bak*`.
- Earlier summary/report figures: `09_make_summary_figs.py`,
  `09b_local_tree_cluster_examples.py`, `build_report_notebook_v3.py`,
  `build_isolation_source_table.py`, `12e_ggtree_circles_supported_v3.R`,
  `12f_plot_headline_prevalence_profile_v3.R`, `12g_plot_master_asv_overview_v3.R`,
  `12h_prep_background_dist.py`, `12h_prep_genus_density.py`,
  `12h_prep_genus_meddist.py`, `12h_plot_insertion_density_v3.R`,
  `12i_plot_supported_profile_v3.R`.
- Earlier HOMD coverage exploration: `13a_homd_genus_coverage.py`,
  `13e_plot_homd_coverage.R`, `13f_homd_phylum_ecdf.py`,
  `13g_plot_homd_phylum_ecdf.R`, `13h_plot_homd_phylum_identity_2panels.R`,
  `13i_plot_homd_phylum_identity_facet.R`, `13j_homd_species_coverage_fullname.py`,
  `13k_plot_homd_species_scatter.R`, `13l_plot_homd_species_quadrant.R`,
  `13m_homd_species_metric_coverage.py`, `13m_plot_homd_phylum_identity_dist.R`,
  `13m_plot_homd_species_metric_dotline.R`, `13n_homd_phylum_identity_bins.py`,
  `13o_homd_species_lt79.py`, `13o2_homd_species_lt97.py`,
  `13o_plot_homd_phylum_identity_bins.R`, `13p_homd_ratio_curves.py`,
  `13p_plot_homd_species_lt97_vs_refseq.R`, `13q_plot_homd_ratio_curves.R`,
  `13s_plot_homd_lt_fractions.R`, `13t_plot_homd_lt_fractions_facet.R`,
  `13u_homd_lt97_ge70_phylum.py`, `13v_homd_phylum_identity_density.py`,
  `13w_plot_homd_phylum_identity_density.R`,
  `13y_plot_homd_phylum_coverage_panel.R`,
  `13z2_plot_homd_phylum_band_composition.R`,
  `13z3_plot_homd_genus_identity_curves.R` (standalone version),
  `13z4_homd_species_identity_curves.py`, `13z4_plot_homd_species_identity_curves.R`.
- Earlier composites: `15a_alt_layout_3row.R`, `99b_preview_panelB_swap.R`,
  `99_combine_master_figs.R`.
- `06_build_cluster_tree.py`: only feeds the tree figures, which are not in the paper.
- `16_plot_cluster_differential_boxplot.R`: the figure it draws is not in the
  manuscript; `16_cluster_differential_wilcox_sig.tsv` is therefore not needed
  either (`17_maaslin3_disease.R` reads the complete `16_..._wilcox.tsv`).

## 7. Excluded — superseded figures and intermediates

- All `fig_v3_homd_*` variants except `fig_v3_homd_genus_identity_curves_numbered`.
- `fig_v3_master_combined*`, `fig_v3_master_with_flow`, `fig_v3_schematic_overview`.
- `fig_v3_panelB_coverage_v3`, `fig_v3_panelB_band_composition_preview`.
- `fig_v3_supported_profile`, `fig_v3_asv_homd_distance_distribution`,
  `fig_v3_insertion_density_tree`, `fig_v3_master_asv_overview`.
- `fig_supported`, `fig_funnel`, `fig_novelty`, `fig_cluster_phylum`,
  `fig_tree_circles_supported_v3`, `fig_v3_headline_prevalence_profile`.
- `fig_v3_homd_species_identity_curves`, `fig_v3_homd_species_quadrant`,
  `fig_v3_homd_species_share_ge_lt` — the species-level curve set, not in the
  composite.
- `.bak_*` PNGs, `.fig_v3_master_layout*.tsv`, `Schematic diagram.pdf`, `Rplots.pdf`.
- Data: `13a_*`, `13f_*`, `13j_*`, `13m_*`, `13n_*`, `13o*_*`, `13p_*`, `13u_*`,
  `13v_*`, `13l_species_xy_schematic.png`, `13z4_*`, `12d_annotation_isolation_sources.tsv`,
  `12f_headline_clusters_v3.tsv`, `12h_genus_density/`, `logs/`, `logs_02_run.log`,
  `01_run.log`, `14_run.log`, `17_maaslin3_four.log`, `treecluster/`, `maaslin3*/`
  (unless you want the raw MaAsLin3 model outputs).

## 8. Inputs that cannot be pushed

| Path | Size | Why |
|---|---|---|
| `HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv` | 1.3 GB | over the limit, untracked |
| `HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv` | 1.5 GB | over the limit, untracked |
| `HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv` | 630 MB | over the limit, untracked |
| `HOMD_uncapture_analysis/pipeline_v2/01_chimera_filtered.tsv` | 553 MB | over the limit, untracked |
| `pipeline/05_sepp/tree_export/tree.nwk` | 20 MB | under the limit, but it is an input from the earlier pipeline; decide whether to push it or just document it |
| `pipeline_v2/gappa/gg2_2022.10_tree_taxonomy.tsv` | 52 MB | same |
| `data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv` | 22 MB | already tracked in the repository |
| eHOMD taxonomy and fasta | — | outside the repository (`/beegfs/huangxiaochang/global_oral/`) |
| each project's `DADA2.biom` | — | outside the repository (`/beegfs/dongbiao/oral/`) |

The README should say where each of these comes from and which script consumes
it. Everything else in sections 1–4 can be pushed as is.

---

## 9. Running order

`Rscript` below is `/home/software/anaconda3/envs/renv/bin/Rscript`; Python is
the matching `python` in the same environment.

```
# phylotype catalog and profiles
python 01_treecluster_cluster.py
python 02_cluster_aggregate.py
python 05_cluster_geography.py
python 07_master_table.py
python 07b_add_organelle_flags.py
python 08_supported_clusters.py
python 12d_build_anchors_v3.py
python 12d_build_headline_profile_v3.py

# HOMD coverage (panels B and C)
python scripts/13b_extract_fasta.py
sbatch scripts/13b_organelle_map.sbatch
python scripts/13d_homd_coverage_clean.py
python scripts/13r_homd_lt_fractions.py
python scripts/13z2_homd_phylum_band_composition.py
python scripts/13z3_homd_genus_identity_curves.py

# sample-level models
python 14_cluster_feature_profile.py
python 15_cluster_cohort_detection.py
Rscript 16_cluster_differential_wilcox.R
Rscript 17_maaslin3_disease.R
python scripts/17d_capnocytophaga_control.py
Rscript scripts/17e_21907_cohort_sensitivity_fit.R
Rscript scripts/17f_ctrl_only_fit.R

# ids: run after every regeneration that writes cluster ids
python 20_rename_cluster_ids.py

# figures
FLOW_ONLY=$PWD/analysis_report/figs/fig_v3_flow_only \
  Rscript scripts/14_schematic_overview.R
Rscript scripts/13z_plot_homd_phylum_coverage_panel_v3.R
Rscript scripts/13z3_plot_homd_genus_identity_curves_numbered.R
Rscript scripts/12j_plot_supported_asv_profiles_v3.R
Rscript scripts/15_alt_layout.R
Rscript 15_plot_cluster_cohort_detection.R
Rscript scripts/17b_plot_cluster_prevalence.R
Rscript scripts/17c_plot_significant_clusters.R
Rscript scripts/17d_plot_capnocytophaga_control.R
Rscript scripts/17e_plot_21907_cohort_sensitivity.R
Rscript scripts/17f_plot_ctrl_only.R

# supplementary tables
python scripts/18_build_disease_feature_table.py
python scripts/19_build_supported_phylotype_table.py
```

Scripts that write more than the figures listed above: `05_cluster_geography.py`
also writes `fig_cluster_phylum.png`; `14_schematic_overview.R` without
`FLOW_ONLY` writes the old schematic composites; `12j` also writes
`fig_v3_supported_asv_profiles_row`. These by-products are not needed and can be
deleted after a run.

## 10. Checks before pushing

1. Figure 1 opens from `fig_v3_master_alt.pdf` and shows panels A–D in order.
2. Panel A is the flow only (no block b/c/d), produced with `FLOW_ONLY`.
3. Decide what the panel C results paragraph should say now that the figure is
   the genus identity curves rather than the species quadrant.
4. `TableS_supported_phylotypes.tsv` has 486 data rows, unique phylotype ids.
5. `TableS_differential_phylotypes_HNC_PD.tsv` has 49 rows, 48 distinct
   phylotypes (`c21907` appears in both cohorts).
6. `17_maaslin3_significant.tsv` has 57 rows and is exactly the
   `qval_individual < 0.05` subset of `17_maaslin3_all_results.tsv`, sorted by
   that column (the step is now part of `17_maaslin3_disease.R`).
7. No file in the push list is ≥ 100 MB.
8. No `.bak*` file or superseded figure is in the push list.

---

## 11. Upstream audit — Methods text against code

Scope: the two Methods sections "Assessment of global oral ASV diversity
against eHOMD" and "Identification of phylotypes not represented at species
level in eHOMD". Paths are relative to `HOMD_uncapture_analysis/`.

### 11.1 Steps that have final code

| Methods step | Code | Numbers |
|---|---|---|
| pooling, exact dereplication | `01_collect_asvs_and_abundance.py` | 2,844,174 |
| chimera screen, `uchime_ref`, VSEARCH 2.7.0, eHOMD 16.03 + GG2 + SILVA 138.2 | `03_chimera_removal.sbatch` | 2,099,459 non-chimeric |
| control samples → authenticity → catalogue | `C01_find_control_samples.py`, `C02_extract_control_sample_ids.py`, `C04_convert_sra_to_fastq.sbatch`, `C05_run_dada2_controls.sbatch`, `C05b_dada2_PRJEB11419.sbatch`, `C06_dada2_extra_controls.sbatch`, `C07_control_quality_viz.sbatch`, `C08_filter_genuine_controls.sbatch`, `C08b_filter_genuine_controls_v2.py`, `C09_build_contaminant_catalog.sbatch` | `contaminant_catalog/` |
| ≥97% identity to control ASVs | Stage 1 of `pipeline_v2/02_contamination_firewall.py` | 160,675 flagged |
| same genus as controls + literature contaminants | Stage 2 of the same script, `pipeline_v2/contaminant_genera.tsv` | 61,245 more, 217,970 in total |
| eHOMD mapping, best hit | `02_homd_mapping.sbatch` | 2,844,174 queries |
| <97% identity → placement | inside `02_homd_mapping.sbatch` (`homd_uncaptured`) | 1,136,142 (uncaptured and non-chimeric) |
| ≥2 BioProjects, ≥50 reads | `pipeline_v2/04_support_floor.py` | 54,290 |
| SEPP into the GG2 tree | `pipeline_v2/04b_sepp_place.sbatch`, `pipeline_v2/04c_sepp_missing.sh` | 54,132 reused + 145 of 158 placed = 54,277 |
| TreeCluster `max_clade` 0.03 | `pipeline_v3/01_treecluster_cluster.py`, `01_treecluster_cluster.sbatch` | |
| drop any phylotype with a flagged member | `pipeline_v3/02_cluster_aggregate.py` (`clean_cluster`), `08_supported_clusters.py` | |
| taxonomy from reference tips in the clade | `pipeline_v3/01_treecluster_cluster.py` | |
| support gate | `pipeline_v3/08_supported_clusters.py` | 486 |

### 11.2 Steps with no code on disk

1. **Assembly of the final tree.** `pipeline_v3/01_treecluster_cluster.py` reads
   `pipeline/05_sepp/tree_export/tree.nwk`, which has 385,546 tips and contains
   the 54,277 placed queries. That tip set is identical to
   `pipeline/05_sepp/gappa_graft/placements.newick`, i.e. the v2 placements were
   grafted onto the backbone with `gappa graft`. The SEPP run and the graft were
   both executed by hand; `pipeline_v2/04b_sepp_place.sbatch` records the
   `fragment-insertion sepp` command but was not submitted unchanged, and no
   script in the repository runs the graft. The two newicks differ only in
   serialisation (same tips, same topology).
2. **`pipeline_v3/organelle_mapping.tsv`.** 1,618 rows, 1,562 distinct ASVs;
   every one of the 1,562 is among the 54,290 support-floor candidates. Four
   columns (`derep_id, ref, identity, coverage`). `scripts/13b_organelle_map.sbatch`
   is a different run — query `13b_uncaptured_asvs.fasta` (1,136,142), seven
   columns, output `13b_organelle.hits.tsv`. Nothing writes `organelle_mapping.tsv`.
3. **`pipeline_v3/gg2_organelle_ref.fasta`.** 6,238 sequences: GG2 2022.10
   full-length sequences filtered on `f__Mitochondria` / `c__Chloroplast` (the
   release taxonomy labels 6,237). Only described in the `07b` docstring; no
   builder.
4. **`pipeline_v2/contaminant_genera.tsv`.** 661 genera, of which 560 carry
   source `control_HOMD_uncaptured` and 101 come from Salter 2014 / Weyrich
   2019. Read by Stage 2; no builder, and the literature half is a hand list
   that needs its citations recorded.
5. **`pipeline_v2/gappa/gg2_2022.10_tree_taxonomy.tsv`.** Exactly the 331,269
   backbone tips. Only 29,940 of its ids occur in the local GG2 exported
   taxonomy, so it is not derivable from files in the repository; it is a GG2
   2022.10 release file. It needs a provenance line (release and source path),
   not a script.

### 11.3 Where the code and the Methods text disagree

1. **Organelle scope.** The screen (`13b`/`13c`) runs only on the 1,136,142
   HOMD-uncaptured non-chimeric ASVs. Panel A subtracts the 16,366 organelle
   ASVs from all 2,099,459 non-chimeric ASVs, so the figure assumes no
   eHOMD-captured ASV is an organelle sequence. Either narrow the Methods
   sentence to the uncaptured set or run the screen on all non-chimeric ASVs.
2. **Order.** Methods: chimera → organelle → contaminant → eHOMD mapping. Code:
   eHOMD mapping on all 2,844,174 → chimera → organelle and contaminant within
   the uncaptured subset. The counts are unaffected; the stated order is not
   what ran.
3. **Contaminant-flagged ASVs are not removed before placement.**
   `04_support_floor.py` keeps them when they meet the floor (21,662 of the
   54,290) so that the whole phylotype can be adjudicated; `08` then drops any
   phylotype with a flagged member. That matches the later sentence but not
   "potential contaminants were stringently removed" before mapping.
4. **Species-level criterion.** Code: no eHOMD hit **or** identity <97% **or**
   query coverage <99%. Methods gives identity ≥97% only. The 99% coverage
   condition is what removes short partial alignments and belongs in the text.
5. **Identity definition.** `--iddef 1` (VSEARCH "all" definition) is used;
   confirm that "identity was calculated across the full query length" describes
   that definition and not `iddef 4`.
6. **VSEARCH version.** The chimera step runs 2.7.0 as stated; the eHOMD mapping
   step runs 2.30.6 (`/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch`),
   which the Methods does not state.
