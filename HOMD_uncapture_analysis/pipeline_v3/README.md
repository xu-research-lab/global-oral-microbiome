# HOMD comprehensiveness of the oral 16S ASV catalogue

The analysis behind Figure 1, the disease figures and the two supplementary
tables: cluster the SEPP-inserted ASVs into phylotypes, annotate each phylotype
against the reference tree, ask which eHOMD species it represents, and test the
supported phylotypes for differential prevalence and abundance across cohorts.

A *phylotype* is what the pipeline code calls a cluster: `c<number>` for two or
more ASVs, `s<number>` for a single ASV, and a plain number for a single-ASV
phylotype that contains reference sequences.

## Layout

```
pipeline_v3/
├── 01–17*.{py,R}          numbered analysis steps
├── 13*_*.tsv              generated tables (see "Outputs")
├── scripts/               figure scripts, HOMD coverage steps, table builders
└── analysis_report/
    ├── figs/              figures
    ├── tables/            supplementary tables + README
    ├── sync_manifest.md   what is synced and why
    └── HOMD_comprehensiveness_methods_results.md
```

## Inputs

These are not in this directory and are not in the repository (too large, or
produced by the earlier pipelines):

| Input | Location |
|---|---|
| SEPP tree | `HOMD_uncapture_analysis/pipeline/05_sepp/tree_export/tree.nwk` |
| GG2 tree taxonomy | `HOMD_uncapture_analysis/pipeline_v2/gappa/gg2_2022.10_tree_taxonomy.tsv` |
| chimera-filtered derep table | `HOMD_uncapture_analysis/pipeline_v2/01_chimera_filtered.tsv` |
| contamination / organelle flags | `HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv` |
| support floor | `HOMD_uncapture_analysis/pipeline_v2/04_support_floor_passed.tsv` |
| per-subject prevalence | `HOMD_uncapture_analysis/pipeline_v2/01b_derep_subject_prevalence.tsv` |
| eHOMD annotation per ASV | `HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv` |
| chimera calls per ASV | `HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv` |
| sample metadata | `data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv` |
| per-project counts | `<project>/results/DADA2.biom` under `/beegfs/dongbiao/oral` |
| eHOMD reference | `/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.*` |
| GG2 organelle taxonomy | `/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv` |

`13b_uncaptured_asvs.fasta` (296 MB) is generated from the eHOMD annotation
table by `scripts/13b_extract_fasta.py` and is not versioned.

## Steps

**Phylotype catalogue.** `01_treecluster_cluster.py` runs TreeCluster
`max_clade` at t = 0.03 (every phylotype is a clade whose maximum pairwise tree
distance is ≤ 0.03) and annotates every ASV from the reference tree:

- `nearest_homd_ref` — within 0.03 of a reference whose species is in eHOMD;
  annotated to that species.
- `cluster_refs` — otherwise, if the phylotype contains references, the most
  specific shared taxonomy of those reference labels.
- `nearest_ref` — single ASVs and ASV-only phylotypes take the nearest
  reference's taxonomy capped at genus; domain as the final fallback.

`02`, `05`, `07`, `07b` and `08` then aggregate to phylotypes, add geography,
build the master ASV table, flag organelles and apply the support gate
(detected in > 50 subjects, ≥ 2 bioprojects, ≥ 2 continents, no
contamination-flagged member, no organelle member). `20_rename_cluster_ids.py`
puts the final ids on every table and must run after any regeneration.

**HOMD coverage.** `scripts/13b`–`13z3` compare the ASVs with eHOMD at the
sequence and tree level and produce the per-phylum coverage, the identity
fractions and the per-genus identity curves used in panels B and C.

**Sample-level models.** `14_cluster_feature_profile.py` builds the
phylotype × sample count and relative-abundance matrices from the projects'
`DADA2.biom` tables. `15` reports detection across cohorts; `16` is the
Wilcoxon screen that defines the cohorts and features modelled; `17` fits
MaAsLin3 per cohort (disease arm plus the healthy controls of the same
bioprojects; covariates age, sex, site, sequencing region, continent; bioproject
as random intercept). `17b`–`17f` draw the disease figures and the sensitivity
checks.

**Tables.** `scripts/18` and `scripts/19` build the two supplementary tables.

The exact running order is in [`analysis_report/sync_manifest.md`](analysis_report/sync_manifest.md).

## Outputs

- `analysis_report/figs/fig_v3_master_alt.{pdf,png}` — Figure 1, panels A–D
- `analysis_report/figs/fig_v3_supported_singleton_top40_profile.*` — Fig. S1
- `analysis_report/figs/fig_*_prevalence_sig.*`, `fig_*_ctrl_only.*`,
  `fig_oscc_*` — disease figures
- `analysis_report/tables/TableS_supported_phylotypes.{tsv,xlsx}` — 486
  supported phylotypes
- `analysis_report/tables/TableS_differential_phylotypes_HNC_PD.{tsv,xlsx}` and
  `..._members.tsv` — phylotypes called differential in head and neck cancer
  and periodontal disease

## Environment

R at `/home/software/anaconda3/envs/renv/bin/Rscript` with ggplot2, dplyr,
patchwork, aplot, cowplot, ggrepel, MaAsLin3. Python at the corresponding
`python` with pandas, numpy, openpyxl and biom. `scripts/15_alt_layout.R` also
needs ImageMagick (`convert`) and `pdflatex`. `01` needs TreeCluster
(`/home/huangxc/.local/bin/TreeCluster.py`) and the organelle step needs
vsearch.
