#!/bin/bash
#SBATCH -J age_accel_cv
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH -t 2-00:00:00
#SBATCH -o age_prediction/logs/age_accel_cv_%j.out
#SBATCH -e age_prediction/logs/age_accel_cv_%j.err

# Leave-studies-out age RF CV (4 oral sites) + one-row replot of the
# disease microbiome-age-acceleration figure.
#   1) run_age_acceleration_cv.py  — per-fold prediction CSVs to
#      age_prediction/result/leave_studies_out_rf/
#   2) plot_age_acceleration_1row.py — stats + 1x4 figure to
#      age_prediction/result/figures/
#      age_disease_microbiome_age_acceleration_1row.{png,pdf}
#      (original figure is left untouched)

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
python "${SCRIPT_DIR}/run_age_acceleration_cv.py"
python "${SCRIPT_DIR}/plot_age_acceleration_1row.py"
echo "[$(date)] all done"