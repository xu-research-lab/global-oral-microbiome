#!/bin/bash
#SBATCH --job-name=da_ageova_merge
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=06:00:00
#SBATCH --output=results/logs/da_ageova_merge_%j.log
#SBATCH --error=results/logs/da_ageova_merge_%j.err
#SBATCH --exclude=cu01

# ==============================================================================
# run_merge_consensus.sbatch  (age OVA da_v2 — tsv based, HC scope only)
# Final stage: merge + n2 consensus (robust = significant in >=2 of
# LinDA / Maaslin2 / ANCOMBC2) for the HC scope.
# Figures use up to 10 genera per unique global abundance peak, selected
# by MaAsLin2 coefficient (03-plot-saliva-ova-bubble.R).
#
# MUST be submitted with --dependency=afterok:<ancom_array_job>:<linda_maas_job>
# so it only starts once every comparison is done (the merge stage scans for
# analysis.done sentinels).
#
# Submit from project root:
#   sbatch --dependency=afterok:<ancom_array_jobid>:<linda_maas_jobid> \
#     differential_abundance/age_differential_ova/run_merge_consensus.sbatch
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

AGE=differential_abundance/age_differential_ova

echo "============================================================"
echo "Age OVA DA merge + consensus (HC scope; robust = >=2 methods)"
echo "Start: $(date)  Job ${SLURM_JOB_ID} on ${SLURM_NODELIST}"
echo "============================================================"

$RSCRIPT $AGE/02-merge-and-consensus.R --scope=hc
RC=$?
if [ $RC -ne 0 ]; then exit $RC; fi
echo "==== [2] consensus figures (saliva OVA bubble) ===="
$RSCRIPT $AGE/03-plot-saliva-ova-bubble.R --scope=hc --select=target_topn --each-n=10
RC=$?

echo ""
echo "============================================================"
echo "ALL DONE (exit=$RC): $(date)"
echo "HC-scope consensus: ${DA_RESULTS_ROOT}/age_hc_ova/age_da_consensus/"
echo "============================================================"
exit $RC
