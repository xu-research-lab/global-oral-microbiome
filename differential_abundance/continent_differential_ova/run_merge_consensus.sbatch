#!/bin/bash
#SBATCH --job-name=da_contova_merge
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=06:00:00
#SBATCH --output=results/logs/da_contova_merge_%j.log
#SBATCH --error=results/logs/da_contova_merge_%j.err
#SBATCH --exclude=cu01

# ==============================================================================
# run_merge_consensus.sbatch  (continent OVA da_v2 — tsv based)
# Final stage: merge + n2 consensus (robust = significant in >=2 of
# LinDA / Maaslin2 / ANCOMBC2) + figures (03-plot-continent-ova-bubble.R).
# SCOPE is selected by the SCOPE env var (all | hc, default hc); pass it at
# submit time via --export=SCOPE=all.
#
# MUST be submitted with --dependency=afterok:<ancom_array_job>:<linda_maas_job>
# (both for the SAME scope) so it only starts once every comparison is done
# (the merge stage scans for analysis.done sentinels).
#
# Submit from project root:
#   sbatch --dependency=afterok:<A>:<B> differential_abundance/continent_differential_ova/run_merge_consensus.sbatch
#   sbatch --export=SCOPE=all --dependency=afterok:<A>:<B> \
#     differential_abundance/continent_differential_ova/run_merge_consensus.sbatch
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

SCOPE=${SCOPE:-hc}
case "$SCOPE" in all|hc) ;; *) echo "[FATAL] SCOPE must be all|hc, got: $SCOPE"; exit 1 ;; esac

CONTOVA=differential_abundance/continent_differential_ova
RES_TAG=$( [ "$SCOPE" = hc ] && echo continent_hc_ova || echo continent_ova )

echo "============================================================"
echo "Continent OVA DA merge + consensus (scope=$SCOPE; robust = >=2 methods)"
echo "Start: $(date)  Job ${SLURM_JOB_ID} on ${SLURM_NODELIST}"
echo "============================================================"

$RSCRIPT $CONTOVA/02-merge-and-consensus.R --scope=$SCOPE
RC=$?
if [ $RC -ne 0 ]; then exit $RC; fi
echo "==== [2] consensus figures (continent OVA bubble) ===="
$RSCRIPT $CONTOVA/03-plot-continent-ova-bubble.R --scope=$SCOPE --select=target_topn --each-n=10
RC=$?

echo ""
echo "============================================================"
echo "ALL DONE (exit=$RC): $(date)"
echo "${SCOPE}-scope consensus: ${DA_RESULTS_ROOT}/${RES_TAG}/continent_da_consensus/"
echo "============================================================"
exit $RC
