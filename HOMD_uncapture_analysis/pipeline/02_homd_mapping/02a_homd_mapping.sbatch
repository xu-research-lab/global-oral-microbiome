#!/bin/bash
#SBATCH --job-name=homd_map_uncaptured
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/homd_mapping_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/homd_mapping_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=240G
#SBATCH --time=48:00:00

# ============================================================
# Step 02: map exact-derep ASVs to eHOMD per amplicon region
#
# HOMD-uncaptured definition (per region, per query):
#   - No HOMD hit, OR
#   - Global identity < 97%, OR
#   - Query coverage (qcov) < 99%
#
# Uses --iddef 1 (query-global identity) per the pipeline plan.
# --maxrejects 32 ensures more thorough search than the default.
# ============================================================

set -e

source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv  # has python with pandas

VSEARCH="/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"
HOMD_DB="/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta"
DEREP_FA="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/01_data_assembly/results/01_derep_representatives.fa"
OUTDIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/02_homd_mapping/results"
DEREP_TABLE="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/01_data_assembly/results/01_derep_summary.tsv"

mkdir -p "$OUTDIR/logs" "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs"

echo "=== HOMD mapping started at $(date) ==="
echo "Query FASTA: $DEREP_FA"
echo "HOMD DB: $HOMD_DB"
echo "Threads: $SLURM_CPUS_PER_TASK"
echo ""

# ---- Step 1: Run vsearch against HOMD ----
HITS_OUT="$OUTDIR/02_derep_vs_homd.hits.tsv"

if [ -s "$HITS_OUT" ]; then
    echo "[skip] HOMD hits file already exists: $HITS_OUT"
else
    echo "[vsearch] running usearch_global (--maxrejects 32, --iddef 1)..."
    $VSEARCH \
        --usearch_global "$DEREP_FA" \
        --db "$HOMD_DB" \
        --id 0.0 \
        --iddef 1 \
        --maxaccepts 1 \
        --maxrejects 32 \
        --strand both \
        --userout "$HITS_OUT" \
        --userfields "query+target+id+alnlen+mism+gaps+qlo+qhi+tlo+thi+ql+tl+qcov+tcov" \
        --threads $SLURM_CPUS_PER_TASK \
        --quiet

    echo "[vsearch] done: $(wc -l < "$HITS_OUT") hits"
fi

# ---- Step 2: Parse hits, flag HOMD-uncaptured ----
echo ""
echo "[python] flagging HOMD-uncaptured ASVs..."

export OUTDIR DEREP_TABLE HITS_OUT

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os
import json

OUTDIR = os.environ['OUTDIR']
DEREP_TABLE = os.environ['DEREP_TABLE']
HITS_OUT = os.environ['HITS_OUT']

print(f"Loading derep table: {DEREP_TABLE}")
derep = pd.read_csv(DEREP_TABLE, sep='\t', dtype=str)
derep['n_bioprojects'] = pd.to_numeric(derep['n_bioprojects'], errors='coerce').fillna(0).astype(int)
derep['n_asvs'] = pd.to_numeric(derep['n_asvs'], errors='coerce').fillna(0).astype(int)
derep['total_reads'] = pd.to_numeric(derep['total_reads'], errors='coerce').fillna(0.0)
derep['prevalence_samples'] = pd.to_numeric(derep['prevalence_samples'], errors='coerce').fillna(0).astype(int)
print(f"  {len(derep):,} derep sequences")

# Load HOMD hits
cols = ['derep_id','homd_target','identity','alnlen','mismatches','gaps',
        'qlo','qhi','tlo','thi','query_length','target_length','query_cov','target_cov']
hits = pd.read_csv(HITS_OUT, sep='\t', names=cols)
for c in ['identity','query_cov','alnlen']:
    hits[c] = pd.to_numeric(hits[c], errors='coerce').fillna(0)
print(f"  {len(hits):,} HOMD hit rows")

# Best hit per query (should be one since --maxaccepts 1)
hits['rank'] = hits.groupby('derep_id')['identity'].rank(ascending=False, method='first')
best = hits[hits['rank'] == 1].drop(columns=['rank']).copy()

# Merge with derep table
merged = derep.merge(
    best[['derep_id','homd_target','identity','alnlen','query_cov','target_cov']],
    on='derep_id', how='left'
)
merged['identity'] = merged['identity'].fillna(0)
merged['query_cov'] = merged['query_cov'].fillna(0)
merged['has_homd_hit'] = merged['homd_target'].notna()

# ---- HOMD-uncaptured definition ----
# Per plan: no HOMD hit, OR identity < 97%, OR qcov < 99%
merged['homd_uncaptured'] = (
    (~merged['has_homd_hit']) |
    (merged['identity'] < 97.0) |
    (merged['query_cov'] < 99.0)
)

# Sub-categories for reporting
merged['uncaptured_reason'] = 'captured'  # default: well-represented in HOMD
merged.loc[~merged['has_homd_hit'], 'uncaptured_reason'] = 'no_hit'
merged.loc[(merged['has_homd_hit']) & (merged['identity'] < 97.0), 'uncaptured_reason'] = 'low_identity'
merged.loc[(merged['has_homd_hit']) & (merged['identity'] >= 97.0) & (merged['query_cov'] < 99.0), 'uncaptured_reason'] = 'low_qcov'
merged.loc[(merged['has_homd_hit']) & (merged['identity'] >= 97.0) & (merged['query_cov'] >= 99.0), 'uncaptured_reason'] = 'captured'

# Save
outpath = os.path.join(OUTDIR, '02_derep_homd_annotated.tsv')
merged.to_csv(outpath, sep='\t', index=False)
print(f"  Wrote {outpath} ({len(merged):,} rows)")

# Summary
total = len(merged)
uncaptured = merged['homd_uncaptured'].sum()
captured = total - uncaptured
print(f"\n=== HOMD MAPPING SUMMARY ===")
print(f"Total derep sequences: {total:,}")
print(f"HOMD-captured (id>=97% AND qcov>=99%): {captured:,} ({captured/total*100:.1f}%)")
print(f"HOMD-uncaptured: {uncaptured:,} ({uncaptured/total*100:.1f}%)")

reason_counts = merged['uncaptured_reason'].value_counts()
for reason, count in reason_counts.items():
    print(f"  {reason}: {count:,} ({count/total*100:.1f}%)")

# Per-region summary
if 'amplicon_regions' in merged.columns:
    print(f"\n=== PER-REGION SUMMARY ===")
    for region, group in merged.groupby('amplicon_regions'):
        n_uncaptured = group['homd_uncaptured'].sum()
        print(f"  {region}: {n_uncaptured:,} / {len(group):,} uncaptured ({n_uncaptured/len(group)*100:.1f}%)")

# Save summary stats
import json
stats = {
    'total_derep_sequences': int(total),
    'homd_captured': int(captured),
    'homd_uncaptured': int(uncaptured),
    'by_reason': {k: int(v) for k, v in reason_counts.items()},
}
with open(os.path.join(OUTDIR, '02_homd_mapping_stats.json'), 'w') as f:
    json.dump(stats, f, indent=2)
print(f"Stats saved to " + os.path.join(OUTDIR, '02_homd_mapping_stats.json'))
PYEOF

echo ""
echo "=== HOMD mapping finished at $(date) ==="
echo "Outputs:"
ls -lh "$OUTDIR"/*.tsv "$OUTDIR"/*.json
