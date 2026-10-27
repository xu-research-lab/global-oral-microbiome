#!/bin/bash
#SBATCH --job-name=build_subtree
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/build_subtree_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/build_subtree_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=6:00:00

set -e
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate mpa4.1

python3 /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/05_sepp/build_phylotype_subtree.py

echo "Done at $(date)"
