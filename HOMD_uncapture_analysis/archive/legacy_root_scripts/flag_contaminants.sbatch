#!/bin/bash
#SBATCH --job-name=flag_contam
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/flag_contam_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/flag_contam_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=4:00:00

set -e
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv

python3 /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/flag_contaminants.py

echo "Done at $(date)"
