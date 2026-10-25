#!/bin/bash
#SBATCH --job-name=da_siteova_lm
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --output=results/logs/da_siteova_lm_%j.log
#SBATCH --error=results/logs/da_siteova_lm_%j.err
#SBATCH --exclude=cu01

# ==============================================================================
# run_linda_maaslin2.sbatch  (site OVA da_v2 — tsv based)
# LinDA + Maaslin2 for the site DA (OVA: each valid site vs the pooled
# "Other"). SCOPE is selected by the SCOPE env var (all | hc, default hc);
# pass it at submit time via --export=SCOPE=all. ANCOMBC2 runs separately via
# run_ancombc2_array.sbatch (SLURM array, one comparison per task).
#
# Submit from project root:
#   sbatch differential_abundance/site_differential_ova/run_linda_maaslin2.sbatch                        # hc (default)
#   sbatch --export=SCOPE=all differential_abundance/site_differential_ova/run_linda_maaslin2.sbatch      # all
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

SCOPE=${SCOPE:-hc}
case "$SCOPE" in all|hc) ;; *) echo "[FATAL] SCOPE must be all|hc, got: $SCOPE"; exit 1 ;; esac

SITEOVA=differential_abundance/site_differential_ova
RES_TAG=$( [ "$SCOPE" = hc ] && echo site_hc_ova || echo site_ova )

echo "============================================================"
echo "Site OVA DA (scope=$SCOPE): LinDA + Maaslin2"
echo "Start: $(date)  Job ${SLURM_JOB_ID} on ${SLURM_NODELIST}"
echo "============================================================"

if [ ! -f data/metadata-and-biom-v2-qc/genus_table_counts.tsv ]; then
  echo "[FATAL] genus_table_counts.tsv missing in data/metadata-and-biom-v2-qc/."
  exit 1
fi

$RSCRIPT $SITEOVA/01-run-three-method-da.R --scope=$SCOPE --method=linda,maaslin2
RC=$?

echo ""
echo "============================================================"
echo "DONE (exit=$RC): $(date)"
echo "Results: ${DA_RESULTS_ROOT}/${RES_TAG}/{site_linda,site_maaslin2}_ova/"
echo "============================================================"
exit $RC
