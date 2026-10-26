#!/usr/bin/env bash
#SBATCH --job-name=oral_da
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --exclude=cu01
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
set -euo pipefail
project_root="${GLOBAL_ORAL_ROOT:-${SLURM_SUBMIT_DIR:-}}"
: "${project_root:?Set GLOBAL_ORAL_ROOT or submit from the global_oral project root}"
cd "$project_root"
export LC_ALL=en_US.UTF-8
rscript="${RSCRIPT_BIN:-/home/xuhuan/miniconda3/envs/R4/bin/Rscript}"
exec "$rscript" differential_abundance/scripts/run_da.R "$@"
