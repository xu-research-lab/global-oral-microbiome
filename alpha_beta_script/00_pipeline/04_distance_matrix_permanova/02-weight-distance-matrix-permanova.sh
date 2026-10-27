#!/bin/bash
#SBATCH --job-name=permanova-unifrac-weight-distance-matrix-0909
#SBATCH --partition=cu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=250G
#SBATCH --time=72:00:00
#SBATCH --output=/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/output/02-weight-distance-matrix-permanova_%j.out
#SBATCH --error=/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/output/02-weight-distance-matrix-permanova_%j.err
#SBATCH --exclude=cu01

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

source /home/hexincheng/miniconda3/bin/activate /home/hexincheng/miniconda3/envs/py39

echo "Running on node: $(hostname)"
echo "Start time: $(date)"
echo "Python: $(which python)"
python -c "import sys, pandas, skbio; print('Python', sys.version.split()[0]); print('pandas', pandas.__version__); print('scikit-bio', skbio.__version__)"

/usr/bin/time -v python /home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/02-weight-distance-matrix-permanova.py

echo "End time: $(date)"