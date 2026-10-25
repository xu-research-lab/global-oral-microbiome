#!/bin/bash
#SBATCH --job-name=da_contova_ancom
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --output=results/logs/da_contova_ancom_%A_%a.log
#SBATCH --error=results/logs/da_contova_ancom_%A_%a.err
#SBATCH --exclude=cu01
#SBATCH --array=1-6

# ==============================================================================
# run_ancombc2_array.sbatch  (continent OVA da_v2 — tsv based)
# Parallelize every ANCOMBC2 continent OVA comparison across SLURM array
# tasks. Each array task = ONE continent OVA comparison (valid continent vs
# the pooled "Other" of the other valid continents; ANCOMBC2 only).
#
# SCOPE is selected by the SCOPE env var (all | hc, default hc); pass it at
# submit time via --export=SCOPE=all. The --only index resolves WITHIN the
# chosen scope. valid continents <= 6 in both scopes (min 1 sample each); the
# array cap is 6 — one task per comparison, no oversubscription in either scope.
#
# Already-completed comparisons are auto-skipped (analysis.done sentinel).
#
# Submit from project root:
#   sbatch differential_abundance/continent_differential_ova/run_ancombc2_array.sbatch                        # hc (default)
#   sbatch --export=SCOPE=all differential_abundance/continent_differential_ova/run_ancombc2_array.sbatch     # all
# Monitor:
#   squeue -j <jobid> ; tail -f results/logs/da_contova_ancom_<jobid>_<taskid>.log
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

SCOPE=${SCOPE:-hc}
case "$SCOPE" in all|hc) ;; *) echo "[FATAL] SCOPE must be all|hc, got: $SCOPE"; exit 1 ;; esac

CONTOVA=differential_abundance/continent_differential_ova
TASKID=${SLURM_ARRAY_TASK_ID}

# Sanity: counts tsv must exist.
if [ ! -f data/metadata-and-biom-v2-qc/genus_table_counts.tsv ]; then
  echo "[FATAL] genus_table_counts.tsv missing in data/metadata-and-biom-v2-qc/."
  echo "        Run: sbatch script/prepare_genus_tables.sbatch"
  exit 1
fi

echo "============================================================"
echo "ANCOMBC2 continent-OVA (scope=$SCOPE) array task $TASKID/$SLURM_ARRAY_TASK_MAX  (da_v2 tsv)"
echo "  driver: $CONTOVA/01-run-three-method-da.R"
echo "  args:   --scope=$SCOPE --method=ancombc2 --only=$TASKID"
echo "  root:   $DA_RESULTS_ROOT"
echo "  node:   $SLURM_NODELIST"
echo "  start:  $(date)"
echo "============================================================"

$RSCRIPT "$CONTOVA/01-run-three-method-da.R" --scope=$SCOPE --method=ancombc2 --only=${TASKID}
RC=$?

echo "============================================================"
echo "task $TASKID exit=$RC  end=$(date)"
echo "============================================================"
exit $RC
