#!/bin/bash
#SBATCH --job-name=da_ageova_lm
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --output=results/logs/da_ageova_lm_%j.log
#SBATCH --error=results/logs/da_ageova_lm_%j.err
#SBATCH --exclude=cu01

# ==============================================================================
# run_linda_maaslin2.sbatch  (age OVA da_v2 — tsv based, HC scope only)
# LinDA + Maaslin2 for the saliva age DA (OVA: each valid age group vs the
# pooled "Other"), HC scope only. ANCOMBC2 runs separately via
# run_ancombc2_array.sbatch (SLURM array, one comparison per task).
#
# Submit from project root:
#   sbatch differential_abundance/age_differential_ova/run_linda_maaslin2.sbatch
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

AGE=differential_abundance/age_differential_ova

echo "============================================================"
echo "Age OVA DA (saliva, HC): LinDA + Maaslin2"
echo "Start: $(date)  Job ${SLURM_JOB_ID} on ${SLURM_NODELIST}"
echo "============================================================"

if [ ! -f data/metadata-and-biom-v2-qc/genus_table_counts.tsv ]; then
  echo "[FATAL] genus_table_counts.tsv missing in data/metadata-and-biom-v2-qc/."
  exit 1
fi

$RSCRIPT $AGE/01-run-three-method-da.R --scope=hc --method=linda,maaslin2
RC=$?

echo ""
echo "============================================================"
echo "DONE (exit=$RC): $(date)"
echo "Results: ${DA_RESULTS_ROOT}/age_hc_ova/{age_linda,age_maaslin2}_ova/"
echo "============================================================"
exit $RC
