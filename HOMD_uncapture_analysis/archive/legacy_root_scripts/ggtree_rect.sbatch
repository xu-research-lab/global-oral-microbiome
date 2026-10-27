#!/bin/bash
#SBATCH --job-name=ggtree_phy
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/ggtree_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/ggtree_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=4:00:00

set -e
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv

Rscript /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/ggtree_phylotypes_rect.R

echo "Done at $(date)"
