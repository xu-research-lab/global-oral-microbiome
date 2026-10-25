#!/bin/bash
#SBATCH --job-name=da_sc_ancom
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --output=results/logs/da_sex_cont_ancom_%A_%a.log
#SBATCH --error=results/logs/da_sex_cont_ancom_%A_%a.err
#SBATCH --exclude=cu01
#SBATCH --array=1-5

# ==============================================================================
# run_ancombc2_array.sbatch  (sex-within-continent da_v2 — tsv based)
# Parallelize every ANCOMBC2 sex-within-continent comparison across SLURM array
# tasks. Each array task = ONE per-continent F_vs_M comparison (ANCOMBC2 only),
# saliva site, HC scope.
#
# Task table (1-5): index -> --only=<continent index>  (order = CONTINENT_ORDER)
#   1 -> Asia,  2 -> Europe,  3 -> North America,
#   4 -> Oceania, 5 -> South America
# (Africa has 0 male saliva-HC samples and is absent from CONTINENT_ORDER.)
#
# Already-completed comparisons are auto-skipped (analysis.done sentinel).
#
# Submit from project root:
#   sbatch differential_abundance/sex_differential_continent/run_ancombc2_array.sbatch
# Monitor:
#   squeue -j <jobid> ; tail -f results/logs/da_sex_cont_ancom_<jobid>_<taskid>.log
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

SEX=differential_abundance/sex_differential_continent

TASKID=${SLURM_ARRAY_TASK_ID}

# Sanity: counts tsv must exist.
if [ ! -f data/metadata-and-biom-v2-qc/genus_table_counts.tsv ]; then
  echo "[FATAL] genus_table_counts.tsv missing in data/metadata-and-biom-v2-qc/."
  echo "        Run: sbatch script/prepare_genus_tables.sbatch"
  exit 1
fi

# Dispatch table: TASKID -> "--only|label"  (order = CONTINENT_ORDER)
dispatch() {
  case "$1" in
    1)  echo "--only=1|hc_saliva_Asia_F_vs_M";;
    2)  echo "--only=2|hc_saliva_Europe_F_vs_M";;
    3)  echo "--only=3|hc_saliva_North_America_F_vs_M";;
    4)  echo "--only=4|hc_saliva_Oceania_F_vs_M";;
    5)  echo "--only=5|hc_saliva_South_America_F_vs_M";;
    *)  echo "INVALID";;
  esac
}

ROW=$(dispatch "$TASKID")
if [ "$ROW" = "INVALID" ]; then
  echo "[task $TASKID] invalid task id"; exit 1
fi

ONLY=$(echo "$ROW" | cut -d'|' -f1)
LABEL=$(echo "$ROW" | cut -d'|' -f2)
ARGS="--scope=hc ${ONLY}"

echo "============================================================"
echo "ANCOMBC2 sex-within-continent array task $TASKID/$SLURM_ARRAY_TASK_MAX  (da_v2 tsv)"
echo "  label:  $LABEL"
echo "  driver: $SEX/01-run-three-method-da.R"
echo "  args:   --method=ancombc2 --scope=hc ${ONLY}"
echo "  root:   $DA_RESULTS_ROOT"
echo "  node:   $SLURM_NODELIST"
echo "  start:  $(date)"
echo "============================================================"

$RSCRIPT "$SEX/01-run-three-method-da.R" --method=ancombc2 $ARGS
RC=$?

echo "============================================================"
echo "task $TASKID ($LABEL) exit=$RC  end=$(date)"
echo "============================================================"
exit $RC
