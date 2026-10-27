#!/bin/bash
# pipeline_v2 Step 04c: SEPP placement for the candidates that step 04b could
# not reuse.  Run directly (not sbatch); the query is a few hundred sequences.

set -e

source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate rachis-qiime2-2026.4

PIPE_V2="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2"
MISSING_FA="$PIPE_V2/candidates_missing_sepp.fasta"
SEPP_REF="/beegfs/huangxiaochang/global_oral/2022.10.backbone.sepp-reference.qza"
OUTDIR="$PIPE_V2/sepp_missing"

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

N_QUERY=$(grep -c "^>" "$MISSING_FA")
echo "=== SEPP for $N_QUERY missing ASVs started at $(date) ==="

# Step 1: Import
echo "[1] Importing sequences..."
qiime tools import \
  --input-path "$MISSING_FA" \
  --output-path "$OUTDIR/candidates.qza" \
  --type 'FeatureData[Sequence]' 2>&1

# Step 2: Run SEPP
echo "[2] Running SEPP..."
qiime fragment-insertion sepp \
  --i-representative-sequences "$OUTDIR/candidates.qza" \
  --i-reference-database "$SEPP_REF" \
  --p-threads 28 \
  --o-tree "$OUTDIR/sepp_insertion_tree.qza" \
  --o-placements "$OUTDIR/sepp_placements.qza" 2>&1

# Step 3: Export
echo "[3] Exporting..."
qiime tools export \
  --input-path "$OUTDIR/sepp_placements.qza" \
  --output-path "$OUTDIR/placements_export" 2>&1

echo "=== SEPP done at $(date) ==="
echo "Placements: $OUTDIR/placements_export/placements.json"
ls -lh "$OUTDIR/placements_export/"
