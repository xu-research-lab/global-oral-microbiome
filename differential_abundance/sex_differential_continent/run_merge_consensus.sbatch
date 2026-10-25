#!/bin/bash
#SBATCH --job-name=da_sc_merge
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=06:00:00
#SBATCH --output=results/logs/da_sex_cont_merge_%j.log
#SBATCH --error=results/logs/da_sex_cont_merge_%j.err
#SBATCH --exclude=cu01

# ==============================================================================
# run_merge_consensus.sbatch  (sex-within-continent da_v2 — tsv based)
# Final stage: merge + n2 consensus (robust = significant in >=2 of
# LinDA / Maaslin2 / ANCOMBC2) + figures, HC scope only (saliva).
#
# MUST be submitted with --dependency=afterok:<ancom_array_job>:<linda_maas_job>
# so it only starts once every comparison is done (the merge stage scans for
# analysis.done sentinels).
#
# Submit from project root:
#   sbatch --dependency=afterok:<ancom_array_jobid>:<linda_maas_jobid> \
#     differential_abundance/sex_differential_continent/run_merge_consensus.sbatch
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

SEX=differential_abundance/sex_differential_continent

echo "============================================================"
echo "Sex DA (continent) merge + consensus (robust = >=2 methods; hc, saliva)"
echo "Start: $(date)  Job ${SLURM_JOB_ID} on ${SLURM_NODELIST}"
echo "============================================================"

SCOPE=hc
echo ""
echo "############## SCOPE: ${SCOPE} ##############"
echo "==== [1] merge + consensus + summary + enhanced + full table ===="
$RSCRIPT $SEX/02-merge-and-consensus.R --scope=${SCOPE}
echo "==== [2] consensus figures ===="
$RSCRIPT $SEX/03-plot-consensus-bubble.R --scope=${SCOPE}
echo "Scope ${SCOPE} consensus done: $(date)"

echo ""
echo "============================================================"
echo "ALL DONE: $(date)"
echo "HC-scope consensus: ${DA_RESULTS_ROOT}/sex_continent_hc/sex_continent_da_consensus/"
echo "============================================================"
