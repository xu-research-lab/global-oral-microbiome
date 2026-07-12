"""Shared paths and conventions for the novel-oral-taxa (not-in-HOMD) pipeline.

Cluster defaults follow ANALYSIS_METHODS_AND_FILTERS.md. Override any of these on the
command line (every script exposes the relevant --flags) so the scripts remain runnable
outside the cluster / on test data.
"""

from __future__ import annotations

# --- reference databases -----------------------------------------------------
HOMD_FASTA = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta"
GG2_FASTA = "/beegfs/dongbiao/greengene2/2024.09.seqs.fna"
GG2_TAXONOMY = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"
# SILVA is used only for the multi-DB novelty cross-check (Phase 4), not for placement.
SILVA_FASTA = ""  # set to a SILVA 138 SSU NR99 fasta if available

# --- per-project inputs ------------------------------------------------------
ORAL_METADATA = "/beegfs/huangxiaochang/global_oral/metadata_v2_merge/metadata-v2-merge-oral.csv"
PROJECT_ASV_FASTA = "/beegfs/dongbiao/oral/{bioproject}/results/ASVs.fa"
PROJECT_ASV_COUNTS = "/beegfs/dongbiao/oral/{bioproject}/results/ASVs_counts.tsv"

# --- upstream pipeline outputs (scripts 14-15) -------------------------------
PIPELINE_DIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
# Annotated low-HOMD table (script 15). NOTE: this is already prevalence-filtered
# (subject_prevalence >= 50). To keep the rare tail (Phase 1) prefer the pre-prevalence
# table `low_identity_lt97_subject_prevalence_filtered.tsv` or the full derep+HOMD table.
LOW_HOMD_ANNOTATED = PIPELINE_DIR + "/low_identity_lt97_subject_prevalence_filtered_gg2_annotated.tsv"

# --- tools -------------------------------------------------------------------
VSEARCH = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"

# --- pipeline output root ----------------------------------------------------
OUT_ROOT = PIPELINE_DIR + "/novel_oral_taxa"

# --- HOMD capture / low-HOMD thresholds --------------------------------------
# "captured" = HOMD identity >= HOMD_ID_CAPTURE and query coverage >= QCOV_MIN.
# Identity is GLOBAL (over the full query), so identity <= qcov and the qcov clause is
# largely redundant; confirm vsearch --iddef is query-global (e.g. --iddef 1).
HOMD_ID_CAPTURE = 97.0   # relaxed capture threshold (percent)
QCOV_MIN = 99.0          # query coverage (percent)

# --- relaxed support floor (KEEP the rare tail; do NOT require >=50 subjects) --
MIN_BIOPROJECTS = 2
MIN_TOTAL_READS = 50

# --- support tiers (label, do not delete) ------------------------------------
TIER_HIGH_SUBJECTS = 50
TIER_MED_SUBJECTS = 10

CONTINENTS = ["North America", "Europe", "Oceania", "Asia", "Africa", "South America"]
