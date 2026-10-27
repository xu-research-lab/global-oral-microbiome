#!/bin/bash
#SBATCH --job-name=sepp_place
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/sepp_place_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/sepp_place_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=120G
#SBATCH --time=24:00:00

# ============================================================
# Phase 4: True SEPP placement onto GG2 backbone
# Uses the official GG2 SEPP reference database
# ============================================================

set -e

source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate rachis-qiime2-2026.4

CANDIDATES_FA=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/02_chimera/candidates_nonchimeric.fasta
SEPP_REF=/beegfs/huangxiaochang/global_oral/2022.10.backbone.sepp-reference.qza
OUTDIR=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/05_sepp
THREADS=$SLURM_CPUS_PER_TASK

mkdir -p "$OUTDIR"

echo "=== SEPP placement started at $(date) ==="
echo "Candidates: $CANDIDATES_FA"
echo "SEPP reference: $SEPP_REF"
echo "Threads: $THREADS"

# Step 1: Import candidate sequences
echo ""
echo "[Step 1] Importing candidate sequences..."
qiime tools import \
  --input-path "$CANDIDATES_FA" \
  --output-path "$OUTDIR/candidates.qza" \
  --type 'FeatureData[Sequence]' 2>&1

# Step 2: Run SEPP
echo ""
echo "[Step 2] Running SEPP placement (this is the heavy step)..."
qiime fragment-insertion sepp \
  --i-representative-sequences "$OUTDIR/candidates.qza" \
  --i-reference-database "$SEPP_REF" \
  --p-threads $THREADS \
  --o-tree "$OUTDIR/sepp_insertion_tree.qza" \
  --o-placements "$OUTDIR/sepp_placements.qza" 2>&1

# Step 3: Export tree and placements
echo ""
echo "[Step 3] Exporting tree and placements..."
qiime tools export \
  --input-path "$OUTDIR/sepp_insertion_tree.qza" \
  --output-path "$OUTDIR/tree_export" 2>&1

qiime tools export \
  --input-path "$OUTDIR/sepp_placements.qza" \
  --output-path "$OUTDIR/placements_export" 2>&1

echo ""
echo "=== SEPP placement finished at $(date) ==="
echo ""
echo "Outputs:"
ls -lh "$OUTDIR"/sepp_*.qza
echo ""
echo "Tree: $OUTDIR/tree_export/tree.nwk"
echo "Placements: $OUTDIR/placements_export/placements.json"
ls -lh "$OUTDIR/tree_export/" 2>/dev/null
ls -lh "$OUTDIR/placements_export/" 2>/dev/null
