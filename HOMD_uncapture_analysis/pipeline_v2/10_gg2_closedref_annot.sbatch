#!/bin/bash
#SBATCH --job-name=pv2_gg2_annot
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/pv2_gg2_annot_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/pv2_gg2_annot_%j.err
#SBATCH --nodelist=cu03
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=120G
#SBATCH --time=48:00:00

# ============================================================
# pipeline_v2 Step 10: GG2 2024.09 closed-reference annotation +
#                     edge-based phylotypes (edge x species label)
# ============================================================
set -e

PIPE_V2="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2"
PYTHON="/home/huangxc/.conda/envs/rachis-qiime2-2026.4/bin/python3"

echo "=== GG2 closed-reference annotation started at $(date) ==="
echo ""

"$PYTHON" "$PIPE_V2/10_gg2_closedref_annot.py"

echo ""
echo "=== Done at $(date) ==="
