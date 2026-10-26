#!/usr/bin/env bash
# Submit from network/ with --array=1-N; N is the study-list line count.
#SBATCH --job-name=sparcc_paired
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cu
#SBATCH --mem=8G
#SBATCH --output=logs/sparcc_%A_%a.out
#SBATCH --error=logs/sparcc_%A_%a.err
set -euo pipefail

NETWORK_DIR=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
PYTHON_BIN=${PYTHON_BIN:-python}
cd "$NETWORK_DIR"
"$PYTHON_BIN" scripts/run_network.py sparcc --unit "${SLURM_ARRAY_TASK_ID:?Set a Slurm array index}"
