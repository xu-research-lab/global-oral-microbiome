#!/bin/bash
#SBATCH --job-name=da_sex_cont
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --output=results/logs/da_sex_continent_%j.log
#SBATCH --error=results/logs/da_sex_continent_%j.err
#SBATCH --exclude=cu01

# ==============================================================================
# run_all_da.sbatch  (sex-within-continent da_v2 — tsv based, da_v4 aligned)
# Full SERIAL run of the three-method DA consensus pipeline for SEX WITHIN
# CONTINENT (saliva site only): per-continent F vs M comparisons (main effect =
# sex_cleaned; covariates age_group, platform, seq_region, bmi_group —
# continent_group is the stratum and is excluded; conditional
# (1|bioproject_id) random effect). HC scope only (status_cleaned == "HC").
#
# Stages:
#   01-run-three-method-da.R   (LinDA + Maaslin2 + ANCOMBC2, up to 5 continents)
#   02-merge-and-consensus.R   (merge + n2 consensus + summary + enhanced + fulltable)
#   03-plot-consensus-bubble.R (figures A/B/C)
#
# For the PARALLEL version (ANCOMBC2 array + afterok merge) use instead:
#   bash differential_abundance/sex_differential_continent/submit_sex_da_continent.sh
#
# Submit from project root:
#   sbatch differential_abundance/sex_differential_continent/run_all_da.sbatch
# Monitor:
#   squeue -j <jobid> ; tail -f results/logs/da_sex_continent_<jobid>.log
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

SEX=differential_abundance/sex_differential_continent

echo "============================================================"
echo "Sex DA run (da_v2: tsv based; per-continent F vs M, saliva; HC only)"
echo "Robust = >=2 methods (n2)"
echo "Results root: ${DA_RESULTS_ROOT}/sex_continent_hc"
echo "Start: $(date)"
echo "Job ${SLURM_JOB_ID} on ${SLURM_NODELIST}"
echo "============================================================"

# Sanity: tsv must exist before any R driver runs.
if [ ! -f data/metadata-and-biom-v2-qc/genus_table_counts.tsv ]; then
  echo "[FATAL] genus_table_counts.tsv missing in data/metadata-and-biom-v2-qc/."
  echo "        Run: sbatch script/prepare_genus_tables.sbatch"
  exit 1
fi

SCOPE=hc
echo ""
echo "############## SCOPE: ${SCOPE} ##############"
echo "==== [1] three-method per-continent DA (LinDA + Maaslin2 + ANCOMBC2; SLOW) ===="
$RSCRIPT $SEX/01-run-three-method-da.R --scope=${SCOPE}
echo "==== [2] merge + consensus (robust = >=2 methods) ===="
$RSCRIPT $SEX/02-merge-and-consensus.R --scope=${SCOPE}
echo "==== [3] consensus figures ===="
$RSCRIPT $SEX/03-plot-consensus-bubble.R --scope=${SCOPE}
echo "Scope ${SCOPE} done: $(date)"

echo ""
echo "============================================================"
echo "ALL DONE: $(date)"
echo "HC-scope consensus: ${DA_RESULTS_ROOT}/sex_continent_hc/sex_continent_da_consensus/"
echo "============================================================"
