#!/bin/bash
#SBATCH --job-name=da_ageova_ancom
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --output=results/logs/da_ageova_ancom_%A_%a.log
#SBATCH --error=results/logs/da_ageova_ancom_%A_%a.err
#SBATCH --exclude=cu01
#SBATCH --array=1-6

# ==============================================================================
# run_ancombc2_array.sbatch  (age OVA da_v2 — tsv based, HC scope ONLY)
# Parallelize every ANCOMBC2 age OVA comparison across SLURM array tasks.
# Each array task = ONE saliva OVA comparison (valid age group vs the pooled
# "Other" of the other valid age groups; ANCOMBC2 only).
#
# HC scope only: every task passes --scope=hc --only=$SLURM_ARRAY_TASK_ID
# (the --only index resolves WITHIN the hc scope). saliva has <=6 valid age
# groups, so the array cap is 6 — one task per comparison, no oversubscription.
#
# Already-completed comparisons are auto-skipped (analysis.done sentinel).
#
# Submit from project root:
#   sbatch differential_abundance/age_differential_ova/run_ancombc2_array.sbatch
# Monitor:
#   squeue -j <jobid> ; tail -f results/logs/da_ageova_ancom_<jobid>_<taskid>.log
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

AGE=differential_abundance/age_differential_ova
TASKID=${SLURM_ARRAY_TASK_ID}

# Sanity: counts tsv must exist.
if [ ! -f data/metadata-and-biom-v2-qc/genus_table_counts.tsv ]; then
  echo "[FATAL] genus_table_counts.tsv missing in data/metadata-and-biom-v2-qc/."
  echo "        Run: sbatch script/prepare_genus_tables.sbatch"
  exit 1
fi

echo "============================================================"
echo "ANCOMBC2 age-OVA (HC) array task $TASKID/$SLURM_ARRAY_TASK_MAX  (da_v2 tsv)"
echo "  driver: $AGE/01-run-three-method-da.R"
echo "  args:   --scope=hc --method=ancombc2 --only=$TASKID"
echo "  root:   $DA_RESULTS_ROOT"
echo "  node:   $SLURM_NODELIST"
echo "  start:  $(date)"
echo "============================================================"

$RSCRIPT "$AGE/01-run-three-method-da.R" --scope=hc --method=ancombc2 --only=${TASKID}
RC=$?

echo "============================================================"
echo "task $TASKID exit=$RC  end=$(date)"
echo "============================================================"
exit $RC
