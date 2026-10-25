#!/bin/bash
#SBATCH --job-name=pv3_treecluster
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/pv3_treecluster_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/pv3_treecluster_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=12:00:00

# ============================================================
# pipeline_v3 Step 01: TreeCluster OTU-like clustering + HOMD coverage
# ============================================================
set -e
# uses base anaconda python3 (has pandas + treeswift via pip --user)
PYTHON=/home/software/anaconda3/bin/python3

PIPE_V3="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"

"$PYTHON" "$PIPE_V3/01_treecluster_cluster.py" \
  --tree /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/05_sepp/tree_export/tree.nwk \
  --threshold 0.03 \
  --method max_clade \
  --outdir "$PIPE_V3"

echo "Done at $(date)"
