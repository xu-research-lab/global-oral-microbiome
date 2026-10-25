#!/bin/bash
#SBATCH --job-name=da_pool_lm
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --output=results/logs/da_sex_pooled_lm_%j.log
#SBATCH --error=results/logs/da_sex_pooled_lm_%j.err
#SBATCH --exclude=cu01

# ==============================================================================
# run_linda_maaslin2.sbatch  (sex_differential_continent/pooled — da_v2 tsv)
# LinDA + Maaslin2 for the POOLED sex DA: ONE saliva F_vs_M comparison with
# ALL continents merged; continent_group enters as a COVARIATE. HC scope.
#
# Results: differential_abundance/result/sex_pooled_hc/{sex_linda,sex_maaslin2}_pairwise/
# Submit from project root:
#   sbatch differential_abundance/sex_differential_continent/pooled/run_linda_maaslin2.sbatch
# ==============================================================================

set -uo pipefail

cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../../.."
export RSCRIPT="${RSCRIPT:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
export DA_RESULTS_ROOT=differential_abundance/result
mkdir -p results/logs

POOLED=differential_abundance/sex_differential_continent/pooled

echo "============================================================"
echo "Sex DA (pooled continents): LinDA + Maaslin2 (saliva F vs M, hc)"
echo "Start: $(date)  Job ${SLURM_JOB_ID} on ${SLURM_NODELIST}"
echo "============================================================"

# Sanity: counts tsv must exist.
if [ ! -f data/metadata-and-biom-v2-qc/genus_table_counts.tsv ]; then
  echo "[FATAL] genus_table_counts.tsv missing in data/metadata-and-biom-v2-qc/."
  exit 1
fi

echo ""
echo "############ SCOPE: hc — LinDA + Maaslin2 (pooled) ############"
$RSCRIPT $POOLED/01-run-three-method-da.R --scope=hc --method=linda,maaslin2
echo "hc linda+maaslin2 (pooled) done: $(date)"

echo ""
echo "============================================================"
echo "ALL DONE: $(date)"
echo "Results: ${DA_RESULTS_ROOT}/sex_pooled_hc/{sex_linda,sex_maaslin2}_pairwise/"
echo "============================================================"