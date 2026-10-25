#!/bin/bash
#SBATCH --job-name=da_sc_lm
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --output=results/logs/da_sex_cont_lm_%j.log
#SBATCH --error=results/logs/da_sex_cont_lm_%j.err
#SBATCH --exclude=cu01

# ==============================================================================
# run_linda_maaslin2.sbatch  (sex-within-continent da_v2 — tsv based)
# LinDA + Maaslin2 for the per-continent sex DA (saliva site), HC scope only.
# ANCOMBC2 runs separately via run_ancombc2_array.sbatch (SLURM array, 5 tasks).
#
# Submit from project root:
#   sbatch differential_abundance/sex_differential_continent/run_linda_maaslin2.sbatch
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

SEX=differential_abundance/sex_differential_continent

echo "============================================================"
echo "Sex DA (continent): LinDA + Maaslin2 (scope: hc, saliva)"
echo "Start: $(date)  Job ${SLURM_JOB_ID} on ${SLURM_NODELIST}"
echo "============================================================"

if [ ! -f data/metadata-and-biom-v2-qc/genus_table_counts.tsv ]; then
  echo "[FATAL] genus_table_counts.tsv missing in data/metadata-and-biom-v2-qc/."
  exit 1
fi

SCOPE=hc
echo ""
echo "############## SCOPE: ${SCOPE} — LinDA + Maaslin2 ##############"
$RSCRIPT $SEX/01-run-three-method-da.R --scope=${SCOPE} --method=linda,maaslin2
echo "Scope ${SCOPE} linda+maaslin2 done: $(date)"

echo ""
echo "============================================================"
echo "ALL DONE: $(date)"
echo "Results: ${DA_RESULTS_ROOT}/sex_continent_hc/{sex_linda,sex_maaslin2}_pairwise/"
echo "============================================================"
