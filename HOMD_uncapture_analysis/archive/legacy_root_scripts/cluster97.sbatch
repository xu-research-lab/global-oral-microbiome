#!/bin/bash
#SBATCH --job-name=cluster97_filter
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/cluster97_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/cluster97_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=6:00:00

set -e
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv

BASE="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"
CATALOG="${BASE}/contaminant_catalog/contaminant_catalog.fasta"
CANDIDATES="${BASE}/pipeline/02_chimera/candidates_nonchimeric.fasta"
TMPDIR="${BASE}/pipeline/03_decontam_genus/tmp_97cluster"
OUTDIR="${BASE}/pipeline/03_decontam_genus"
mkdir -p "$TMPDIR"

echo "=== 97% Clustering: Control ASVs + Candidate ASVs ==="
echo "Started: $(date)"

# 1. Merge control + candidate FASTA
cat "$CATALOG" "$CANDIDATES" > "$TMPDIR/merged.fa"
grep -c "^>" "$TMPDIR/merged.fa"

# 2. vsearch 97% cluster
vsearch --cluster_fast "$TMPDIR/merged.fa" \
  --id 0.97 \
  --uc "$TMPDIR/clusters.uc" \
  --threads 16 \
  --sizeout \
  --fasta_width 0

# 3. Check results
python3 << 'PYINNER'
# Count how many clusters contain both control and candidate ASVs
# And how many candidates are removed
ctrls = set()
cands = set()
current_clust = None
n_total_clusters = 0
n_mixed = 0
n_cand_removed = 0

with open("'"$TMPDIR/clusters.uc"'") as f:
    for line in f:
        parts = line.strip().split("\t")
        if parts[0] == "C":
            n_total_clusters += 1
            if current_clust and ctrls and cands:
                n_mixed += 1
                n_cand_removed += len(cands)
            current_clust = parts[1]
            ctrls = set()
            cands = set()
        elif parts[0] == "H":
            seq_id = parts[8]
            if seq_id.startswith("ASV"):
                ctrls.add(seq_id)
            else:
                cands.add(seq_id)

# Last cluster
if current_clust and ctrls and cands:
    n_mixed += 1
    n_cand_removed += len(cands)

print(f"\n=== 97% Clustering Results ===")
print(f"Total clusters: {n_total_clusters}")
print(f"Mixed clusters (ctrl+cand): {n_mixed}")
print(f"Candidates removed: {n_cand_removed}")
PYINNER

rm -rf "$TMPDIR"
echo "Done at $(date)"
