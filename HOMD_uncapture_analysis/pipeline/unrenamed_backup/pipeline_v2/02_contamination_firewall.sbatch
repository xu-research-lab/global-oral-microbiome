#!/bin/bash
#SBATCH --job-name=contamination_firewall
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/contamination_firewall_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/contamination_firewall_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=180G
#SBATCH --time=24:00:00

set -e
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate mpa4.1

PIPE=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline

cd "$PIPE/07_contamination_firewall"
python3 07a_contamination_firewall.py

echo "Done at $(date)"
