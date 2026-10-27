#!/bin/bash
#SBATCH --job-name=unweight_pcoa_csv
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p cu
#SBATCH --mem=50G
#SBATCH -o /home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/output/01-export-unweight-pcoa-hdf5-to-csv.out
#SBATCH -e /home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/output/01-export-unweight-pcoa-hdf5-to-csv.err
#SBATCH --exclude=cu01

source /home/hexincheng/miniconda3/bin/activate /home/hexincheng/miniconda3/envs/py39

python /home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/01-export-unweight-pcoa-hdf5-to-csv.py