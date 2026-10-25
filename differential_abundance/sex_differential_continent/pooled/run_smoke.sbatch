#!/bin/bash
#SBATCH --job-name=da_pool_smoke
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --output=results/logs/da_sex_pooled_smoke_%j.log
#SBATCH --error=results/logs/da_sex_pooled_smoke_%j.err
#SBATCH --exclude=cu01

# ==============================================================================
# run_smoke.sbatch  (sex_differential_continent/pooled — END-TO-END SMOKE)
# Subsampled run (a few BioProjects x capped samples per sex, see
# pooled_helpers.R) exercising the WHOLE pooled chain on a compute node:
#   1. 01-run-three-method-da.R --smoke  (all three methods, hc scope)
#   2. 02-merge-and-consensus.R          (over the smoke results)
# Results land under a SEPARATE root (differential_abundance/result/da_v2_smoke_pooled) so real
# results are never touched. NEVER run the driver directly on the login node.
#
# Submit from project root:
#   sbatch differential_abundance/sex_differential_continent/pooled/run_smoke.sbatch
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result_smoke_pooled
mkdir -p results/logs

POOLED=differential_abundance/sex_differential_continent/pooled

# Sanity: counts tsv must exist.
if [ ! -f data/metadata-and-biom-v2-qc/genus_table_counts.tsv ]; then
  echo "[FATAL] genus_table_counts.tsv missing in data/metadata-and-biom-v2-qc/."
  exit 1
fi

echo "============================================================"
echo "POOLED sex DA SMOKE (all 3 methods, hc, subsampled)"
echo "Start: $(date)  Job ${SLURM_JOB_ID} on ${SLURM_NODELIST}"
echo "  root: $DA_RESULTS_ROOT"
echo "============================================================"

echo ""
echo "==== [1/2] 01 driver (smoke) ===="
$RSCRIPT $POOLED/01-run-three-method-da.R --scope=hc --method=linda,maaslin2,ancombc2 --smoke
RC1=$?
echo "01 exit=$RC1  $(date)"
if [ $RC1 -ne 0 ]; then
  echo "[FATAL] smoke 01 failed; aborting before 02."
  exit 1
fi

echo ""
echo "==== [2/2] 02 merge + consensus (smoke) ===="
$RSCRIPT $POOLED/02-merge-and-consensus.R --scope=hc
RC2=$?
echo "02 exit=$RC2  $(date)"

echo ""
echo "============================================================"
echo "SMOKE DONE (01=$RC1, 02=$RC2): $(date)"
echo "Smoke consensus: ${DA_RESULTS_ROOT}/sex_pooled_hc/sex_pooled_da_consensus/"
echo "Clean up with: rm -rf differential_abundance/result/da_v2_smoke_pooled"
echo "============================================================"
exit $(( RC1 + RC2 ))