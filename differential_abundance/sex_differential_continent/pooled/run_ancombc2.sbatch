#!/bin/bash
#SBATCH --job-name=da_pool_ancom
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --output=results/logs/da_sex_pooled_ancom_%j.log
#SBATCH --error=results/logs/da_sex_pooled_ancom_%j.err
#SBATCH --exclude=cu01

# ==============================================================================
# run_ancombc2.sbatch  (sex_differential_continent/pooled — da_v2 tsv based)
# ANCOMBC2 for the SINGLE pooled comparison saliva_pooled_F_vs_M (all
# continents merged, continent_group covariate). No SLURM array: there is
# exactly ONE comparison (the parent per-continent version arrayed over 5).
# ~10.5k HC saliva samples in one job => 3-day limit (parent: 1 day x 5).
#
# Results: differential_abundance/result/sex_pooled_hc/sex_ancombc2_pairwise/
# Submit from project root:
#   sbatch differential_abundance/sex_differential_continent/pooled/run_ancombc2.sbatch
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

POOLED=differential_abundance/sex_differential_continent/pooled

# Sanity: counts tsv must exist.
if [ ! -f data/metadata-and-biom-v2-qc/genus_table_counts.tsv ]; then
  echo "[FATAL] genus_table_counts.tsv missing in data/metadata-and-biom-v2-qc/."
  echo "        Run: sbatch script/prepare_genus_tables.sbatch"
  exit 1
fi

echo "============================================================"
echo "ANCOMBC2 pooled sex DA (scope: hc, saliva F vs M)"
echo "  label:  saliva_pooled_F_vs_M"
echo "  driver: $POOLED/01-run-three-method-da.R"
echo "  args:   --method=ancombc2 --scope=hc"
echo "  root:   $DA_RESULTS_ROOT"
echo "  node:   $SLURM_NODELIST"
echo "  start:  $(date)"
echo "============================================================"

$RSCRIPT "$POOLED/01-run-three-method-da.R" --method=ancombc2 --scope=hc
RC=$?

echo "============================================================"
echo "pooled ancombc2 exit=$RC  end=$(date)"
echo "============================================================"
exit $RC