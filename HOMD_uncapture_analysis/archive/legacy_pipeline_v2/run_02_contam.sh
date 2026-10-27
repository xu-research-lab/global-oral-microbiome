#!/bin/bash
# Run contamination firewall on cu05 (foreground)
# Usage: bash pipeline_v2/run_02_contam.sh

set -e

source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate mpa4.1

cd /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis

echo "=== pipeline_v2 Step 02: Contamination Firewall ==="
echo "Started at $(date)"
echo "Host: $(hostname)"
echo ""

python3 pipeline_v2/02_contamination_firewall.py

echo ""
echo "=== Done at $(date) ==="
