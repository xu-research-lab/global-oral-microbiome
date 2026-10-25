#!/bin/bash
#SBATCH --job-name=da_pool_merge
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=06:00:00
#SBATCH --output=results/logs/da_sex_pooled_merge_%j.log
#SBATCH --error=results/logs/da_sex_pooled_merge_%j.err
#SBATCH --exclude=cu01

# ==============================================================================
# run_merge_consensus.sbatch  (sex_differential_continent/pooled — da_v2 tsv)
# Merge + n2 consensus (robust = significant in >=2 of LinDA / Maaslin2 /
# ANCOMBC2) over the SINGLE pooled comparison saliva_pooled_F_vs_M, HC scope.
# No bubble-plot stage: the parent 03-plot-consensus-bubble.R is per-continent
# and is not part of the pooled pipeline.
#
# MUST be submitted with --dependency=afterok:<lm_job>:<ancom_job> so it only
# starts once both methods are done (the merge stage scans analysis.done
# sentinels).
#
# Submit from project root:
#   sbatch --dependency=afterok:<lm_jobid>:<ancom_jobid> \
#     differential_abundance/sex_differential_continent/pooled/run_merge_consensus.sbatch
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

POOLED=differential_abundance/sex_differential_continent/pooled

echo "============================================================"
echo "Pooled sex DA merge + consensus (robust = >=2 methods; hc, saliva)"
echo "Start: $(date)  Job ${SLURM_JOB_ID} on ${SLURM_NODELIST}"
echo "============================================================"

echo ""
echo "==== merge + consensus + summary + enhanced + full table (hc) ===="
$RSCRIPT $POOLED/02-merge-and-consensus.R --scope=hc
echo "hc pooled consensus done: $(date)"

echo ""
echo "============================================================"
echo "ALL DONE: $(date)"
echo "Consensus: ${DA_RESULTS_ROOT}/sex_pooled_hc/sex_pooled_da_consensus/"
echo "============================================================"