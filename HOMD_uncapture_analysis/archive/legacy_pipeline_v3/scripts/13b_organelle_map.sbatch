#!/bin/bash
#SBATCH --job-name=pv3_organelle
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/pv3_organelle_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/pv3_organelle_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=250G
#SBATCH --time=24:00:00

set -e
VSEARCH=/home/software/anaconda3/envs/qiime2-2022.2/bin/vsearch
PYTHON=/home/software/anaconda3/envs/renv/bin/python
PIPE_V3="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
QUERY="$PIPE_V3/13b_uncaptured_asvs.fasta"
REF="$PIPE_V3/gg2_organelle_ref.fasta"
OUT="$PIPE_V3/13b_organelle.hits.tsv"

echo "=== vsearch organelle mapping (28 threads) ==="
"$VSEARCH" --usearch_global "$QUERY" \
  --db "$REF" \
  --id 0.97 \
  --threads 28 \
  --userout "$OUT" \
  --userfields query+target+id+alnlen+mism+qcov+tcov \
  --wordlength 12 \
  --notrunclabels \
  --quiet

echo "=== annotate 02_contamination_flagged.tsv ==="
"$PYTHON" "$PIPE_V3/scripts/13c_organelle_annotate.py"

echo "Done at $(date)"
