#!/bin/bash
#SBATCH -J age_rand_split
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH -t 2-00:00:00
#SBATCH -o age_prediction/logs/age_rand_split_%j.out
#SBATCH -e age_prediction/logs/age_rand_split_%j.err

# Random (sample-level) split age RF CV (4 oral sites) + 1x4 scatter
# figure, the companion of the leave-studies-out version:
#   1) run_age_predict_random_split.py — per-fold prediction CSVs to
#      age_prediction/result/random_split_rf/ (KFold(5) over SAMPLES,
#      studies mixed; same pool / RF params as the leave-studies run)
#   2) plot_random_split_scatter.py — 1x4 figure to
#      age_prediction/result/figures/
#      age_random_split_scatter_by_site.{png,pdf}

set -euo pipefail

SCRIPT_DIR="$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
LOG_DIR=${ROOT}/logs
mkdir -p "${ROOT}/age_prediction/logs"

source /home/xuhuan/miniconda3/etc/profile.d/conda.sh
conda activate global_oral

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-local} cpus=${SLURM_CPUS_PER_TASK:-16}"
python -c "import biom, sklearn, pandas; print('sklearn', sklearn.__version__)"

cd "${ROOT}"
python "${SCRIPT_DIR}/run_age_predict_random_split.py"
python "${SCRIPT_DIR}/plot_random_split_scatter.py"
echo "[$(date)] all done"