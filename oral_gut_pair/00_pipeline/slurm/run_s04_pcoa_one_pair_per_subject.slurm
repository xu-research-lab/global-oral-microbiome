#!/bin/bash
#SBATCH -J ts_s04_pcoa_1pair
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=80G
#SBATCH --exclude=cu01,cu02
#SBATCH -o logs/s04_pcoa_paired_adonis2/s04_one_pair_per_subject_%j.out
#SBATCH -e logs/s04_pcoa_paired_adonis2/s04_one_pair_per_subject_%j.err

set -euo pipefail
set -x

PROJECT_ROOT="/home/zhangcaiyi/project/oral_gut_pair_TimeSeries_11bp"
cd "$PROJECT_ROOT"

echo "===== START s04 one-pair-per-subject sensitivity ====="
date
hostname
pwd

module load miniconda/4.9.2
source activate oralgut_rf

module load R/4.1.3
export R_LIBS_USER=~/Rlibs

echo "Python:"
which python
python --version

echo "Rscript:"
which Rscript
Rscript --version

echo "Check R packages:"
Rscript -e 'library(vegan); library(permute); cat("R packages OK\n")'

echo "Input check:"
ls -lh 01.filter_data/filter_metadata_one_pair_per_subject.csv
ls -lh 03.beta_diversity/distance_matrix/braycurtis_filter_ra_distance_matrix.tsv
ls -lh 03.beta_diversity/distance_matrix/weighted_unifrac_filter_distance_matrix.tsv
ls -lh 03.beta_diversity/distance_matrix/unweighted_unifrac_filter_distance_matrix.tsv

python -u scripts/s04_beta_diversity_pcoa_permanova_paired_adonis2.py \
  --distance-dir 03.beta_diversity/distance_matrix \
  --metadata 01.filter_data/filter_metadata_one_pair_per_subject.csv \
  --out-dir 03.beta_diversity/pcoa_paired_adonis2_one_pair_per_subject \
  --sample-col sample_id \
  --bioproject-col bioproject \
  --site-col site_group_cleaned \
  --pair-col pair_id \
  --status-col status_cleaned \
  --metrics braycurtis,weighted_unifrac,unweighted_unifrac \
  --permutations 999 \
  --seed 2026 \
  --adonis2-by terms \
  --paired true \
  --strict-one-to-one true \
  --grid-n-cols 4 \
  --overwrite

echo "Output check:"
ls -lh 03.beta_diversity/pcoa_paired_adonis2_one_pair_per_subject

echo "===== DONE s04 one-pair-per-subject sensitivity ====="
date
