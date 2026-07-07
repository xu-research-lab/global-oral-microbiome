#!/bin/bash
# db-RDA 分析运行脚本
# 1. 从 QIIME2 PCoA qza 文件提取 PCoA 坐标
# 2. 运行 R db-RDA 分析

set -e

QZA_FILE="/beegfs/huangxiaochang/global_oral/bray_curtis_pcoa_results.qza"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Activate R environment
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv

echo "=== db-RDA Analysis Pipeline ==="
echo "Working directory: $SCRIPT_DIR"
echo "PCoA file: $QZA_FILE"

# Step 1: Extract PCoA coordinates from qza (skip if already done)
if [ -f "$SCRIPT_DIR/pcoa_coordinates.tsv" ] && [ -f "$SCRIPT_DIR/pcoa_eigenvalues.tsv" ]; then
    echo ""
    echo "[1/2] PCoA coordinates already extracted, skipping..."
else
    echo ""
    echo "[1/2] Extracting PCoA coordinates from QIIME2 qza file..."
    python3 "$SCRIPT_DIR/extract_pcoa_from_qza.py" "$QZA_FILE" "$SCRIPT_DIR" 1000
fi

# Step 2: Run R db-RDA analysis
echo ""
echo "[2/2] Running db-RDA analysis in R..."
cd "$SCRIPT_DIR"
Rscript dbRDA_analysis.R

echo ""
echo "=== Pipeline Complete ==="
echo "Results saved in: $SCRIPT_DIR"
