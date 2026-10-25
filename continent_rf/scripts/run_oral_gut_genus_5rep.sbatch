#!/bin/bash
#SBATCH -J rf_genus_5rep
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 7-00:00:00
#SBATCH --exclude=cu01
#SBATCH -o continent_rf/logs/rf_genus_5rep_%j.out
#SBATCH -e continent_rf/logs/rf_genus_5rep_%j.err

# Genus-level oral+gut project-level RF OVA/OVO (5-fold x 5 repeats).
# Writes ONLY to continent_rf/result_genus_5rep/ (original result/ untouched).
# Comparison table vs the feature-level 5-rep results is generated at the end.

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

echo "[$(date)] START genus-level CV (oral + gut, 5 repeats)"
python -u "${SCRIPT_DIR}/run_oral_gut_genus_5rep.py"
echo "[$(date)] DONE genus-level CV"

echo "[$(date)] ALL DONE"