# Global oral 16S catalogue — the analysis behind the four figures

This pipeline answers two questions about a global collection of oral 16S
samples:

1. how much of the observed ASV catalogue eHOMD represents at species level;
2. which lineages fall below that boundary, and whether they form coherent,
   widely shared sequence clusters.

It does this by dereplicating every ASV exactly, comparing the unique sequences
with eHOMD, removing chimeras, organelles and contaminants, inserting the
remaining candidate sequences into a reference tree, clustering them, and then
testing whether the resulting clusters are supported by independent projects,
subjects, geography and disease status.

The directory is ordered by the steps of that analysis, so the numbered
directories can be read as a methods section.

## Deliverables

All four PDFs are in `19_figures/final/`:

| File | Content |
|---|---|
| `fig_v3_master_alt.pdf` | main figure: A data-set flow, B eHOMD coverage, C per-genus identity curves, D multi-ASV profile |
| `fig_v3_supported_singleton_top40_profile.pdf` | supplement: the 40 most prevalent supported singleton clusters |
| `fig_head_and_neck_cancer_ctrl_only.pdf` | supplement: head-and-neck cancer, control cohorts only |
| `fig_periodontal_disease_ctrl_only.pdf` | supplement: periodontal disease, control cohorts only |

The intermediate panels that are composited into the main figure are in
`19_figures/panels/`; they are regenerable and can be ignored when reading the
results.

### Figure → script → table

| Figure | Script | Inputs |
|---|---|---|
| main panel A | `19_figures/19e_schematic_flow.R` | none (counts are transcribed from the run) |
| main panel B | `19_figures/19a_plot_homd_phylum_coverage_panel.R` (sources `19b`) | `08_homd_phylum_coverage.tsv`, `08_homd_phylum_lt_fractions.tsv`, `08_homd_total_lt_fractions.tsv`, `08_homd_phylum_band_composition.tsv` |
| main panel C | `19_figures/19c_plot_homd_genus_identity_curves.R` | `08_homd_genus_identity_curves.tsv` |
| main panel D | `19_figures/19d_plot_supported_asv_profiles.R` | `13_anchors_supported.tsv`, `12_master_asv_table.tsv`, `13_headline_profile.tsv` |
| main figure | `19_figures/19f_master_alt_layout.R` | the four panels above |
| singleton supplement | `19_figures/19d_plot_supported_asv_profiles.R` | same as panel D |
| ctrl-only supplements | `19_figures/19g_plot_ctrl_only.R` | `14_cluster_relabund.tsv`, `14_cluster_sample_metadata.tsv`, `17_ctrl_only_results.tsv`, `17_ctrl_only_samples.tsv`, `14_supported_cluster_annotation.tsv` |

## The steps

Every step directory holds its scripts and a `results/` subdirectory. Steps 01,
06, 09 and 10 are SLURM jobs on the cluster; the rest are short Python or R
scripts. Counts in the last column are the ones the published figures carry.

| Step | Scripts | Purpose | Key result |
|---|---|---|---|
| 01 data assembly | `01a_process_data.sh`, `01b_process_forwards.R`, `01c_collect_asvs_and_abundance.py` | DADA2 per project, pool the ASV tables, dereplicate exactly | `01_derep_summary.tsv`, `01_derep_representatives.fa`, `01_asv_abundance_metadata.tsv` |
| 02 eHOMD mapping | `02a_homd_mapping.sbatch` | vsearch against eHOMD 16.03 per amplicon region; captured = ≥97 % identity and ≥99 % query coverage | `02_derep_homd_annotated.tsv` |
| 03 chimera removal | `03a_chimera_removal.sbatch` | `uchime_ref` against eHOMD + GG2 + SILVA, calls merged | `03_derep_chimera_annotated.tsv`, 2,099,459 non-chimeric |
| 04 candidate selection | `04a_chimera_filter.py` | keep ASVs that are both eHOMD-uncaptured and non-chimeric | `04_candidates.tsv`, 1,136,142 |
| 05 organelle removal | `05a_extract_fasta.py`, `05b_organelle_map.sbatch`, `05c_organelle_annotate.py` | vsearch against GG2 organelle references; flags are written standalone | `05_organelle_hits.tsv`, `05_organelle_flagged.tsv`, `05_organelle_excluded.tsv` |
| 06 control catalog | `06a`–`06i` | find control samples, run DADA2 on them, keep genuine controls, build the contaminant ASV catalogue | `inputs/contaminant_asv_catalog.tsv` |
| 07 contamination firewall | `07a_contamination_firewall.py`, `07b_contamination_firewall.sbatch` | Stage 1: 97 % clustering of candidates with every control ASV. Stage 2: GG2 genus against a curated contaminant list. Organelle flags from step 05 are merged afterwards | `07_contamination_flagged.tsv` |
| 08 coverage analysis | `08a`–`08d` | species and phylum coverage, identity fractions, identity-band composition, per-genus identity curves | `08_homd_*.tsv` |
| 09 support floor + SEPP | `09a_support_floor.py`, `09b_sepp_place.sbatch`, `09c_sepp_missing.sh`, `09d_subject_prevalence.py` | ≥2 BioProjects and ≥50 reads; SEPP insertion; per-subject prevalence | `09_support_floor_passed.tsv`, 54,290; `09_subject_prevalence.tsv` |
| 10 clustering | `10a_treecluster_cluster.py`, `10b_treecluster_cluster.sbatch`, `10c_cluster_aggregate.py` | TreeCluster `max_clade` at t = 0.03 on the SEPP tree; aggregate cluster members | `10_asv_clusters.tsv`, `10_clusters.tsv`, `10_clusters_aggregated.tsv`, 54,277 placed |
| 11 geography | `11a_cluster_geography.py` | site and continent distribution of each cluster | `11_candidate_novel_with_phylum.tsv` |
| 12 master table | `12a_master_table.py`, `12b_add_organelle_flags.py` | one row per ASV; attach the cluster-gate organelle flag | `12_master_asv_table.tsv` |
| 13 supported clusters | `13a_supported_clusters.py`, `13b_build_anchors.py`, `13c_build_headline_profile.py` | subject prevalence > 50 %, ≥2 BioProjects, ≥2 continents, no contaminant and no organelle member | `13_supported_clusters.tsv`, 486 |
| 14 feature profile | `14a_cluster_feature_profile.py` | cluster counts, relative abundance, sample metadata, cluster annotation | `14_cluster_counts.tsv`, `14_cluster_relabund.tsv`, `14_cluster_sample_metadata.tsv`, `14_supported_cluster_annotation.tsv` |
| 15 disease cohorts | `15a_select_cohorts.R` | disease cohorts with more than 100 samples and HC samples in the same projects; the 15 largest are kept | `15_cohorts.tsv` |
| 16 MaAsLin3 | `16a_maaslin3_disease.R` | disease models for each cohort selected in step 15 | `16_maaslin3_all_results.tsv`, `16_maaslin3_significant.tsv` |
| 17 ctrl-only | `17a_ctrl_only_fit.R` | refit the significant features in cohorts that have controls | `17_ctrl_only_results.tsv`, `17_ctrl_only_samples.tsv` |
| 18 cluster ids | `18a_rename_cluster_ids.py` | rename `21907 → c21907`, `-28239 → s28239` in every table; textual and idempotent | `18_cluster_id_rename_map.tsv` |
| 19 figures | `19a`–`19g` | all drawing | the four PDFs |

The data-set flow drawn in panel A is: 2,844,174 non-redundant ASVs →
2,099,459 non-chimeric → 2,083,093 organelle-free → 1,871,022 contaminant-free
→ 965,488 at ≥97 % identity and 905,534 below it. The candidate chain that
feeds the insertion analysis is separate: 1,136,142 eHOMD-uncaptured,
non-chimeric ASVs → 54,290 above the support floor → 54,277 placed in the tree
→ 486 supported clusters.

## Two organelle screens, kept apart

There are two organelle screens in this analysis and they must not be
substituted for one another.

- **Coverage screen (step 05).** Every candidate is searched against the GG2
  organelle references; any hit at the run threshold flags the ASV. This
  produces 16,366 flagged ASVs and is used only to remove them from the eHOMD
  coverage statistics in step 08.
- **Cluster gate (step 12).** `inputs/organelle_cluster_gate_mapping.tsv` is a
  frozen mapping of 1,562 ASVs. It sets `any_organelle` in the cluster and
  master tables, and step 13 uses it to exclude clusters. It comes from an
  earlier, separate screen.

All 1,562 gate ASVs are inside the support floor, but only 1,533 of them are in
the 16,366 step-05 flags, and no threshold on the step-05 hits reproduces
1,562. Replacing the frozen mapping with the step-05 table would change
`any_organelle` and therefore the 486 supported clusters. The mapping is
therefore kept as an input and is not recomputed.

Step 05 runs before step 07 in the paper's order, but this does not change any
number: step 05 flags are computed on all 1,136,142 candidates and written to
their own tables; step 07 runs its 97 % clustering on the same complete
candidate set and merges the organelle flags only afterwards. `vsearch
--cluster_fast` is sensitive to the order and membership of its input, so the
clustering input must stay complete — filtering organelle rows out beforehand
would change cluster membership.

## Running it

`run_order.sh` lists the commands in order. It is meant to be read and run step
by step, not executed end to end: steps 01, 06, 09 and 10 are SLURM jobs with
multi-hour runtimes, and steps 02 and 03 run vsearch against multi-gigabyte
references.

Environment:

| Tool | Path |
|---|---|
| Python (pandas) | `/home/software/anaconda3/envs/renv/bin/python` |
| R (ggplot2, patchwork, maaslin3) | `/home/software/anaconda3/envs/renv/bin/Rscript` |
| QIIME2 (DADA2) | `conda activate qiime2-2022.2` |
| vsearch | `/home/software/anaconda3/envs/qiime2-2022.2/bin/vsearch` |
| TreeCluster | `/home/huangxc/.local/bin/TreeCluster.py` |

Each script carries its own absolute paths; there is no config file to edit.

## Inputs, outputs and what is not regenerated here

`inputs/` holds the few files several steps share and that are expensive or
impossible to rebuild from this repository alone:

| File | Why it is frozen |
|---|---|
| `sepp_placements_parsed.tsv`, `tree.nwk` | SEPP placements and the reference tree used for insertion and clustering |
| `gg2_2022.10_tree_taxonomy.tsv.gz` | reference taxonomy for tree labels and phylum assignment |
| `contaminant_asv_catalog.tsv`, `contaminant_genera.tsv` | control-derived contaminant catalogue and the curated genus list |
| `gg2_organelle_ref.fasta` | organelle reference for the coverage screen |
| `organelle_cluster_gate_mapping.tsv` | the frozen 1,562-ASV cluster gate (see above) |

Raw reads and the per-project outputs of step 01 live outside this repository.
Steps 01–18 are not re-run for the paper: they need the raw reads, multi-hour
vsearch runs and SEPP, and their outputs are moved here rather than
regenerated. The large upstream tables in `results/` are regenerable from the
raw reads but are too big for version control and are listed in `../.gitignore`;
everything else under `results/` is tracked.

The figure chain (step 19) is re-run from those tables and reproduces the
committed PDFs. `19d` and `19g` fix a seed before their jittered point layers,
so the scatter positions are stable between runs. The PDFs differ from the
versions generated at submission time only in that jitter realisation.

`unrenamed_backup/` holds the same scripts under their original filenames, with
`MAPPING.tsv` giving the old path, the new path and the purpose of each.
