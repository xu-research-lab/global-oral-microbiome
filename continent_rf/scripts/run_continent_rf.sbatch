#!/bin/bash
#SBATCH -J continent_rf
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 7-00:00:00
#SBATCH --exclude=cu01
#SBATCH -o continent_rf/logs/continent_rf_%j.out
#SBATCH -e continent_rf/logs/continent_rf_%j.err

# Submit from the repository root; CLI arguments are forwarded to the runner.
set -euo pipefail

# Slurm copies this file into its spool; use the submission directory there.
ROOT="${GLOBAL_ORAL_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
source /home/xuhuan/miniconda3/etc/profile.d/conda.sh
conda activate global_oral
cd "${ROOT}"
python -u continent_rf/scripts/run_continent_rf.py --jobs "${SLURM_CPUS_PER_TASK:-8}" "$@"
