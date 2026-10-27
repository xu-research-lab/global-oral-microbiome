#!/bin/bash
# Command order for the global oral 16S analysis.
#
# This file is a reading guide, not a script to run end to end.  Steps 01, 06,
# 09 and 10 are SLURM jobs with multi-hour runtimes; steps 02 and 03 run
# vsearch against multi-gigabyte references.  Run the stage you need, check its
# results/ directory, then continue.
#
# Every command below is written relative to this directory:
#   cd /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

PY=/home/software/anaconda3/envs/renv/bin/python
R=/home/software/anaconda3/envs/renv/bin/Rscript

## ---- 01 data assembly -----------------------------------------------------
# 01a is a batch template: one job per project, taking the project directory.
# sbatch 01_data_assembly/01a_process_data.sh <project_dir>
$PY 01_data_assembly/01c_collect_asvs_and_abundance.py

## ---- 02 eHOMD mapping -----------------------------------------------------
sbatch 02_homd_mapping/02a_homd_mapping.sbatch

## ---- 03 chimera removal ---------------------------------------------------
sbatch 03_chimera_removal/03a_chimera_removal.sbatch

## ---- 04 candidate selection -----------------------------------------------
$PY 04_candidate_selection/04a_chimera_filter.py

## ---- 05 organelle removal (coverage screen) -------------------------------
$PY 05_organelle_removal/05a_extract_fasta.py
sbatch 05_organelle_removal/05b_organelle_map.sbatch   # runs 05c when vsearch finishes

## ---- 06 control catalog ---------------------------------------------------
$PY 06_control_catalog/06a_find_control_samples.py
$PY 06_control_catalog/06b_extract_control_sample_ids.py
sbatch 06_control_catalog/06c_convert_sra_to_fastq.sbatch
sbatch 06_control_catalog/06d_dada2_controls.sbatch
sbatch 06_control_catalog/06e_dada2_PRJEB11419.sbatch
sbatch 06_control_catalog/06f_dada2_extra_controls.sbatch
sbatch 06_control_catalog/06g_control_quality_viz.sbatch
sbatch 06_control_catalog/06h_filter_genuine_controls.sbatch
sbatch 06_control_catalog/06i_build_contaminant_catalog.sbatch

## ---- 07 contamination firewall --------------------------------------------
sbatch 07_contamination_firewall/07b_contamination_firewall.sbatch

## ---- 08 coverage analysis -------------------------------------------------
$PY 08_coverage_analysis/08a_homd_coverage_clean.py
$PY 08_coverage_analysis/08b_homd_lt_fractions.py
$PY 08_coverage_analysis/08c_homd_phylum_band_composition.py
$PY 08_coverage_analysis/08d_homd_genus_identity_curves.py

## ---- 09 support floor and SEPP --------------------------------------------
$PY 09_support_floor/09a_support_floor.py
sbatch 09_support_floor/09b_sepp_place.sbatch
bash 09_support_floor/09c_sepp_missing.sh              # direct, not sbatch
$PY 09_support_floor/09d_subject_prevalence.py

## ---- 10 clustering --------------------------------------------------------
sbatch 10_clustering/10b_treecluster_cluster.sbatch    # runs 10a
$PY 10_clustering/10c_cluster_aggregate.py

## ---- 11 geography ---------------------------------------------------------
$PY 11_geography/11a_cluster_geography.py

## ---- 12 master table ------------------------------------------------------
$PY 12_master_table/12a_master_table.py
$PY 12_master_table/12b_add_organelle_flags.py

## ---- 13 supported clusters ------------------------------------------------
$PY 13_supported_clusters/13a_supported_clusters.py
$PY 13_supported_clusters/13b_build_anchors.py
$PY 13_supported_clusters/13c_build_headline_profile.py

## ---- 14 feature profile ---------------------------------------------------
$PY 14_feature_profile/14a_cluster_feature_profile.py

## ---- 15 disease cohorts ---------------------------------------------------
$R 15_cohorts/15a_select_cohorts.R

## ---- 16 MaAsLin3 ----------------------------------------------------------
$R 16_maaslin3/16a_maaslin3_disease.R

## ---- 17 control-only refit ------------------------------------------------
$R 17_ctrl_only/17a_ctrl_only_fit.R

## ---- 18 cluster ids -------------------------------------------------------
$PY 18_cluster_ids/18a_rename_cluster_ids.py

## ---- 19 figures -----------------------------------------------------------
$R 19_figures/19a_plot_homd_phylum_coverage_panel.R    # panel B, sources 19b
$R 19_figures/19c_plot_homd_genus_identity_curves.R    # panel C
$R 19_figures/19d_plot_supported_asv_profiles.R        # panel D and the singleton supplement
$R 19_figures/19e_schematic_flow.R                     # panel A
$R 19_figures/19f_master_alt_layout.R                  # composites A-D into the main figure
$R 19_figures/19g_plot_ctrl_only.R                     # the two ctrl-only supplements
