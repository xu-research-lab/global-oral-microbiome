#!/bin/bash
#SBATCH --job-name=chimera_ref
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/chimera_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/chimera_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=240G
#SBATCH --time=72:00:00

# ============================================================
# Step 03: reference-based chimera removal
#
# Use vsearch --uchime_ref against combined HOMD + GG2 + SILVA
# to flag cross-study/cross-taxon chimeras that DADA2 misses.
#
# Reference databases:
#   HOMD:  /beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta  (6,600 seqs)
#   GG2:   /beegfs/dongbiao/greengene2/backbone_fasta/dna-sequences.fasta              (337,506 seqs)
#   SILVA: /beegfs/dongbiao/silva_138.2/SILVA_138.2_SSURef_NR99.fasta                  (510,495 seqs)
#   Total combined: ~854,601 sequences
# ============================================================

set -e

source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2

VSEARCH="$(which vsearch)"
DEREP_FA="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/01_data_assembly/results/01_derep_representatives.fa"
OUTDIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/03_chimera_removal/results"
HOMD_ANNOTATED="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/02_homd_mapping/results/02_derep_homd_annotated.tsv"

# Reference databases
HOMD_DB="/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta"
GG2_DB="/beegfs/dongbiao/greengene2/backbone_fasta/dna-sequences.fasta"
SILVA_DB="/beegfs/dongbiao/silva_138.2/SILVA_138.2_SSURef_NR99.fasta"
COMBINED_DB="$OUTDIR/03_combined_reference.fasta"

mkdir -p "$OUTDIR" "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs"

echo "=== Chimera removal started at $(date) ==="
echo "VSEARCH: $VSEARCH"
echo "Query FASTA: $DEREP_FA"
echo "Combined DB: $COMBINED_DB"
echo "Threads: $SLURM_CPUS_PER_TASK"
echo ""

# ---- Step 0: Build combined reference database ----
if [ -s "$COMBINED_DB" ]; then
    echo "[skip] Combined reference database already exists: $COMBINED_DB"
else
    echo "[build] Combining reference databases (adding prefixes to avoid ID collisions)..."

    # Use sed to add prefix: >HMT_xxx -> >HOMD|HMT_xxx
    sed 's/^>/&HOMD|/' "$HOMD_DB" > "$COMBINED_DB"
    echo "  HOMD: done ($(python3 -c "print(sum(1 for l in open('$COMBINED_DB') if l.startswith('>')))") sequences so far)"

    sed 's/^>/&GG2|/' "$GG2_DB" >> "$COMBINED_DB"
    echo "  GG2: appended ($(python3 -c "print(sum(1 for l in open('$COMBINED_DB') if l.startswith('>')))") sequences so far)"

    sed 's/^>/&SILVA|/' "$SILVA_DB" >> "$COMBINED_DB"
    echo "  SILVA: appended"

    N_COMBINED=$(python3 -c "print(sum(1 for l in open('$COMBINED_DB') if l.startswith('>')))")
    echo "[build] Combined database: $N_COMBINED sequences"
    echo "[build] Size: $(du -h "$COMBINED_DB" | cut -f1)"
fi

# ---- Step 1: Run uchime_ref ----
CHIMERAS_OUT="$OUTDIR/03_chimeras.fa"
NONCHIMERAS_OUT="$OUTDIR/03_nonchimeras.fa"
UCHIME_OUT="$OUTDIR/03_uchime_detailed.tsv"

if [ -s "$UCHIME_OUT" ]; then
    echo "[skip] uchime_ref output already exists: $UCHIME_OUT"
else
    N_QUERY=$(python3 -c "print(sum(1 for l in open('$DEREP_FA') if l.startswith('>')))")
    N_REF=$(python3 -c "print(sum(1 for l in open('$COMBINED_DB') if l.startswith('>')))")
    echo "[uchime_ref] Starting at $(date)..."
    echo "  Query: $N_QUERY derep sequences"
    echo "  Reference: $N_REF sequences"

    $VSEARCH \
        --uchime_ref "$DEREP_FA" \
        --db "$COMBINED_DB" \
        --chimeras "$CHIMERAS_OUT" \
        --nonchimeras "$NONCHIMERAS_OUT" \
        --uchimeout "$UCHIME_OUT" \
        --threads $SLURM_CPUS_PER_TASK

    echo "[uchime_ref] Done at $(date)"
    N_CHIM=$(python3 -c "print(sum(1 for l in open('$CHIMERAS_OUT') if l.startswith('>')))")
    N_NONCHIM=$(python3 -c "print(sum(1 for l in open('$NONCHIMERAS_OUT') if l.startswith('>')))")
    echo "  Chimeras: $N_CHIM sequences"
    echo "  Non-chimeras: $N_NONCHIM sequences"
fi

# ---- Step 2: Parse results and update annotation ----
echo ""
echo "[python] Parsing chimera results and updating annotations..."

export OUTDIR HOMD_ANNOTATED UCHIME_OUT

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os
import json

OUTDIR = os.environ['OUTDIR']
HOMD_ANNOTATED = os.environ['HOMD_ANNOTATED']
UCHIME_OUT = os.environ['UCHIME_OUT']

# Load uchime_ref output
# Format: query_id<TAB>score<TAB>Y/N<TAB>parent_A<TAB>parent_B<TAB>...
uchime = pd.read_csv(UCHIME_OUT, sep='\t', header=None)
uchime.columns = [
    'derep_id', 'score', 'is_chimera', 'parent_A', 'parent_B',
    'id_QM', 'id_QA', 'id_QB', 'id_QC',
    'id_AB', 'id_AC', 'id_QT', 'id_AT', 'id_BT', 'id_CT',
    'divergence'
] + [f'extra_{i}' for i in range(max(0, uchime.shape[1] - 16))]

for c in ['score', 'divergence']:
    if c in uchime.columns:
        uchime[c] = pd.to_numeric(uchime[c], errors='coerce').fillna(0)

n_chimeric = (uchime['is_chimera'] == 'Y').sum()
n_nonchim = (uchime['is_chimera'] == 'N').sum()
n_other = len(uchime) - n_chimeric - n_nonchim
print(f"Uchime results:")
print(f"  Chimeric (Y): {n_chimeric:,}")
print(f"  Non-chimeric (N): {n_nonchim:,}")
print(f"  Other: {n_other:,}")
print(f"  Total: {len(uchime):,}")

# Score distribution for chimeras
chimeras_only = uchime[uchime['is_chimera'] == 'Y']
if len(chimeras_only) > 0:
    print(f"  Chimera score range: {chimeras_only['score'].min():.3f} - {chimeras_only['score'].max():.3f}")
    print(f"  Chimera score median: {chimeras_only['score'].median():.3f}")

# Load HOMD annotated table
print(f"\nLoading HOMD annotated table...")
homd = pd.read_csv(HOMD_ANNOTATED, sep='\t', dtype=str)
homd['n_bioprojects'] = pd.to_numeric(homd['n_bioprojects'], errors='coerce').fillna(0).astype(int)
homd['n_asvs'] = pd.to_numeric(homd['n_asvs'], errors='coerce').fillna(0).astype(int)
homd['total_reads'] = pd.to_numeric(homd['total_reads'], errors='coerce').fillna(0.0)
homd['prevalence_samples'] = pd.to_numeric(homd['prevalence_samples'], errors='coerce').fillna(0).astype(int)
homd['identity'] = pd.to_numeric(homd['identity'], errors='coerce').fillna(0)
homd['query_cov'] = pd.to_numeric(homd['query_cov'], errors='coerce').fillna(0)
print(f"  {len(homd):,} sequences")

# Merge chimera flag
chimera_flag = uchime[['derep_id', 'is_chimera', 'score', 'parent_A', 'parent_B', 'divergence']].copy()
chimera_flag.columns = ['derep_id', 'uchime_flag', 'uchime_score', 'uchime_parent_A', 'uchime_parent_B', 'uchime_divergence']
chimera_flag['is_chimera'] = chimera_flag['uchime_flag'] == 'Y'

merged = homd.merge(chimera_flag, on='derep_id', how='left')
merged['is_chimera'] = merged['is_chimera'].fillna(False)
merged['uchime_flag'] = merged['uchime_flag'].fillna('NO_RESULT')
merged['uchime_score'] = merged['uchime_score'].fillna(-1)
merged['uchime_divergence'] = merged['uchime_divergence'].fillna(-1)

n_chim = merged['is_chimera'].sum()
print(f"\n  Chimeric sequences: {n_chim:,} ({n_chim/len(merged)*100:.2f}%)")

# Among HOMD-uncaptured, how many are chimeric?
is_uncaptured = (merged['homd_uncaptured'].astype(str).str.lower() == 'true')
uncaptured = merged[is_uncaptured]
n_uncaptured_chim = uncaptured['is_chimera'].sum()
print(f"  HOMD-uncaptured: {len(uncaptured):,}")
print(f"  HOMD-uncaptured AND chimeric: {n_uncaptured_chim:,} ({n_uncaptured_chim/len(uncaptured)*100:.2f}% of uncaptured)")

# New column: candidate for novelty (uncaptured + non-chimera)
merged['passed_chimera_filter'] = ~merged['is_chimera']
# Post-chimera: uncaptured AND non-chimeric
merged['is_candidate'] = is_uncaptured & (~merged['is_chimera'])
n_candidate = merged['is_candidate'].sum()
print(f"  Candidates (uncaptured + non-chimera): {n_candidate:,} ({n_candidate/len(merged)*100:.2f}%)")

# Save
outpath = os.path.join(OUTDIR, '03_derep_chimera_annotated.tsv')
merged.to_csv(outpath, sep='\t', index=False)
print(f"\n  Wrote {outpath} ({len(merged):,} rows)")

# ---- Per-region summary ----
if 'amplicon_regions' in merged.columns:
    print(f"\n=== PER-REGION CHIMERA & CANDIDATE SUMMARY ===")
    region_stats = {}
    main_regions = ['V3-V4', 'V4', 'V1-V3', 'V4-V5', 'V3-V5', 'V6-V8', 'V5-V7', 'V1-V2', 'V2-V3', 'V1-V3,V7-V9', 'V2,V3,V4,V6-V7,V8,V9']
    for region in main_regions:
        if region in merged['amplicon_regions'].values:
            group = merged[merged['amplicon_regions'] == region]
            n_chim_r = group['is_chimera'].sum()
            n_uncapt = (group['homd_uncaptured'].astype(str).str.lower() == 'true').sum()
            n_cand = group['is_candidate'].sum()
            print(f"  {region}: {len(group):,} total | {n_chim_r:,} chimeras ({n_chim_r/len(group)*100:.1f}%) | {n_uncapt:,} uncaptured | {n_cand:,} candidates")
            region_stats[region] = {
                'total': len(group),
                'chimeras': int(n_chim_r),
                'chimera_pct': round(n_chim_r/len(group)*100, 1),
                'homd_uncaptured': int(n_uncapt),
                'candidates': int(n_cand),
            }

# ---- Save stats ----
total = len(merged)
summary = {
    'total_derep_sequences': int(total),
    'chimera_flagged': int(n_chim),
    'chimera_pct': round(n_chim/total*100, 2),
    'non_chimera': int(total - n_chim),
    'homd_uncaptured_total': int(len(uncaptured)),
    'homd_uncaptured_chimera': int(n_uncaptured_chim),
    'candidates_post_chimera': int(n_candidate),
    'candidate_pct': round(n_candidate/total*100, 2),
}
if 'amplicon_regions' in merged.columns:
    summary['per_region'] = region_stats

summary_path = os.path.join(OUTDIR, '03_chimera_stats.json')
with open(summary_path, 'w') as f:
    json.dump(summary, f, indent=2)
print(f"\nStats saved to {summary_path}")
print("Done!")
PYEOF

echo ""
echo "=== Chimera removal finished at $(date) ==="
echo "Outputs:"
ls -lh "$OUTDIR"/*.tsv "$OUTDIR"/*.json "$OUTDIR"/*.fa 2>/dev/null
