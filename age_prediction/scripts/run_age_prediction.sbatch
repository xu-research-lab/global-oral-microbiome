#!/bin/bash
#SBATCH -J age_prediction
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH -t 2-00:00:00
#SBATCH -o age_prediction/logs/age_prediction_%j.out
#SBATCH -e age_prediction/logs/age_prediction_%j.err

# Submit from the repository root; arguments are forwarded to the runner.
set -euo pipefail
ROOT="${GLOBAL_ORAL_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
source /home/xuhuan/miniconda3/etc/profile.d/conda.sh
conda activate global_oral
cd "${ROOT}"
python -u age_prediction/scripts/run_age_prediction.py --jobs "${SLURM_CPUS_PER_TASK:-16}" "$@"
