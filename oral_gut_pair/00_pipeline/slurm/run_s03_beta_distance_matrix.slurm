#!/bin/bash
#SBATCH -J ts_beta_s03
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=120G
#SBATCH --exclude=cu01,cu02
#SBATCH -t 48:00:00
#SBATCH -o logs/s03_beta_diversity_distance_matrix/s03_beta_%j.out
#SBATCH -e logs/s03_beta_diversity_distance_matrix/s03_beta_%j.err

set -euo pipefail

echo "Job started at: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Node: $(hostname)"
echo "Workdir before cd: $(pwd)"

PROJECT_ROOT="/home/zhangcaiyi/project/oral_gut_pair_TimeSeries_11bp"
cd "$PROJECT_ROOT"

echo "Workdir after cd: $(pwd)"

module load miniconda/4.9.2 || true
source activate oralgut_rf || source activate base

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "Python: $(which python)"
python --version

echo "Checking required input files..."
test -s scripts/s03_beta_diversity_distance_matrix.py
test -s 01.filter_data/filter_metadata.csv
test -s 01.filter_data/filter.biom
test -s 01.filter_data/filter_ra.biom
test -s 00.data/phylogeny_time_series_11bp.nwk

echo "Running beta distance matrix calculation..."

python scripts/s03_beta_diversity_distance_matrix.py \
  --filter-dir 01.filter_data \
  --metadata 01.filter_data/filter_metadata.csv \
  --count-biom 01.filter_data/filter.biom \
  --ra-biom 01.filter_data/filter_ra.biom \
  --tree 00.data/phylogeny_time_series_11bp.nwk \
  --out-dir 03.beta_diversity/distance_matrix \
  --sample-col sample_id \
  --bioproject-col bioproject \
  --site-col site_group_cleaned \
  --pair-col pair_id \
  --metrics braycurtis,weighted_unifrac,unweighted_unifrac \
  --midpoint-root-tree \
  --overwrite

echo "Checking output files..."
ls -lh 03.beta_diversity/distance_matrix || true

echo "Job finished at: $(date)"
