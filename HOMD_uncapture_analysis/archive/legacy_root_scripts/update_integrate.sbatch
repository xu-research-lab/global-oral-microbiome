#!/bin/bash
#SBATCH --job-name=upd_integrate
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/upd_integrate_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/upd_integrate_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=2:00:00

set -e
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv

python3 /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/update_integrated_annotations.py

echo "Done at $(date)"
