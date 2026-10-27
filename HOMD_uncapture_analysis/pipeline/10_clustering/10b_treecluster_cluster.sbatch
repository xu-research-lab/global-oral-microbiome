#!/bin/bash
#SBATCH --job-name=treecluster
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/pv3_treecluster_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/pv3_treecluster_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=12:00:00

# ============================================================
# Step 10b: TreeCluster clustering on the SEPP tree
# ============================================================
set -e
# uses base anaconda python3 (has pandas + treeswift via pip --user)
PYTHON=/home/software/anaconda3/bin/python3

PIPE="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline"
PIPE_V3="$PIPE"

"$PYTHON" "$PIPE_V3/10_clustering/10a_treecluster_cluster.py" \
  --tree "$PIPE_V3/inputs/tree.nwk" \
  --threshold 0.03 \
  --method max_clade \
  --outdir "$PIPE_V3/10_clustering/results"

echo "Done at $(date)"
