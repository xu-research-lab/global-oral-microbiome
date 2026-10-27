#!/bin/bash
#SBATCH --job-name=pv2_self97
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/pv2_self97_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/pv2_self97_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=120G
#SBATCH --time=12:00:00

set -e
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate mpa4.1

cd /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis
python3 pipeline_v2/03_self_cluster.py

echo "Done at $(date)"
