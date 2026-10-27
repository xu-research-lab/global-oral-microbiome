#!/bin/bash
#SBATCH --job-name=organelle_map
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/organelle_map_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/organelle_map_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=250G
#SBATCH --time=24:00:00

set -e
VSEARCH=/home/software/anaconda3/envs/qiime2-2022.2/bin/vsearch
PYTHON=/home/software/anaconda3/envs/renv/bin/python
PIPE="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline"
QUERY="$PIPE/05_organelle_removal/results/05_candidates.fasta"
REF="$PIPE/inputs/gg2_organelle_ref.fasta"
HITS="$PIPE/05_organelle_removal/results/05_organelle_hits.tsv"

echo "=== vsearch organelle mapping (28 threads) ==="
"$VSEARCH" --usearch_global "$QUERY" \
  --db "$REF" \
  --id 0.97 \
  --threads 28 \
  --userout "$HITS" \
  --userfields query+target+id+alnlen+mism+qcov+tcov \
  --wordlength 12 \
  --notrunclabels \
  --quiet

echo "=== annotate organelle hits ==="
"$PYTHON" "$PIPE/05_organelle_removal/05c_organelle_annotate.py"

echo "Done at $(date)"
