#!/bin/bash
#SBATCH -J rf_oral_gut_cv
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 7-00:00:00
#SBATCH --exclude=cu01
#SBATCH -o continent_rf/logs/rf_cv_repeat5_%j.out
#SBATCH -e continent_rf/logs/rf_cv_repeat5_%j.err

# Run oral+gut project-level RF OVA/OVO (5-fold x 5 repeats).
# The script hardcodes datasets=['oral','gut'] and uses N_JOBS=8.

set -euo pipefail

SCRIPT_DIR="$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
LOG_DIR=${ROOT}/continent_rf/logs
mkdir -p "${LOG_DIR}"

source /home/xuhuan/miniconda3/etc/profile.d/conda.sh
conda activate global_oral

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-local} cpus=${SLURM_CPUS_PER_TASK:-8}"
echo "[$(date)] python=$(which python)"
python -c "import biom, sklearn, pandas, numpy; print('sklearn', sklearn.__version__)"

cd "${ROOT}"

echo "[$(date)] START project-level CV (oral + gut)"
python -u "${SCRIPT_DIR}/run_oral_gut_project_cv_repeat5.py"
echo "[$(date)] DONE project-level CV"

echo "[$(date)] START plots"
python -u "${SCRIPT_DIR}/plot_oral_gut_project_cv.py"
echo "[$(date)] DONE plots"

echo "[$(date)] ALL DONE"
