#!/bin/bash
#SBATCH --job-name=filter_ctrl
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/filter_ctrl_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/filter_ctrl_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=4:00:00

# ============================================================
# Filter genuine control samples:
#   Keep: n_asvs < 30 AND oral_frac < 20%
#   Exclude: everything else
# Visualize with scatter plot + per-BP bar chart
# ============================================================

set -e
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv

BASE_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"
TEMP_DIR="${BASE_DIR}/control_quality_viz/temp"
OUT_DIR="${BASE_DIR}/control_quality_viz"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({'font.size': 11, 'axes.titlesize': 13, 'axes.labelsize': 11,
                     'figure.dpi': 150, 'savefig.dpi': 150, 'savefig.bbox': 'tight'})

BASE_DIR = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"
TEMP_DIR = f"{BASE_DIR}/control_quality_viz/temp"
OUT_DIR = f"{BASE_DIR}/control_quality_viz"

CATEGORY_ORDER = ['extraction_blank', 'blank', 'ntc_no_template', 'negative_control', 'reagent_water']
CATEGORY_LABELS = {
    'extraction_blank': 'Extraction Blank', 'blank': 'Blank',
    'ntc_no_template': 'NTC', 'negative_control': 'Negative Control',
    'reagent_water': 'Reagent Water'
}
CAT_PALETTE = {
    'extraction_blank': '#d62728', 'blank': '#1f77b4',
    'ntc_no_template': '#2ca02c', 'negative_control': '#ff7f0e',
    'reagent_water': '#9467bd',
}

# ============================================================
# Core oral genera (high specificity for oral cavity)
# ============================================================
ORAL_GENERA = {
    'Streptococcus', 'Prevotella', 'Veillonella', 'Haemophilus',
    'Neisseria', 'Rothia', 'Fusobacterium', 'Leptotrichia',
    'Porphyromonas', 'Actinomyces', 'Oribacterium', 'Capnocytophaga',
    'Parvimonas', 'Alloprevotella', 'Gemella', 'Granulicatella',
    'Selenomonas', 'Campylobacter', 'Tannerella',
    'Lachnoanaerobaculum', 'Peptostreptococcus', 'Schaalia',
    'Kingella', 'Eikenella', 'Aggregatibacter', 'Filifactor',
    'Cardiobacterium', 'Abiotrophia', 'Treponema',
}

# ============================================================
# Load data
# ============================================================
print("Loading data...")

# Parse SILVA taxonomy
cols = ['query', 'target', 'id', 'alnlen', 'mism', 'gaps',
        'qlo', 'qhi', 'tlo', 'thi', 'ql', 'tl', 'qcov', 'tcov']
hits = pd.read_csv(f"{TEMP_DIR}/silva_hits.txt", sep="\t", names=cols)
hits['true_id'] = (hits['alnlen'] - hits['mism'] - hits['gaps']) / hits['ql'] * 100

silva_tax = {}
with open("/beegfs/dongbiao/silva_138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta") as f:
    for line in f:
        if line.startswith(">"):
            full_line = line.strip().lstrip(">")
            space_idx = full_line.find(" ")
            if space_idx > 0:
                silva_tax[full_line[:space_idx]] = full_line[space_idx+1:]

def parse_taxonomy(target_id):
    tax_str = silva_tax.get(target_id, "Unknown")
    ranks = [r.strip() for r in tax_str.split(";")]
    result = {}
    for i, rn in enumerate(['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
        result[rn] = ranks[i] if i < len(ranks) else 'Unknown'
    return result

# ASV -> genus
asv_genus = {}
for _, row in hits.iterrows():
    if row['true_id'] >= 80 and row['query'] not in asv_genus:
        asv_genus[row['query']] = parse_taxonomy(row['target'])['genus']

# Load sample data
with open(f"{TEMP_DIR}/sample_data.json") as f:
    sample_data = json.load(f)

# ============================================================
# Compute per-sample oral fraction
# ============================================================
print("Computing per-sample oral fractions...")
results = []

for s in sample_data:
    total = s['total_reads']
    n_asvs = s['n_asvs']
    category = s['category']
    if category not in CATEGORY_ORDER:
        continue

    oral_reads = 0
    genus_reads = defaultdict(int)

    for asv_id, count in s['asv_counts'].items():
        genus = asv_genus.get(asv_id, 'Unclassified')
        genus_reads[genus] += count
        if genus in ORAL_GENERA:
            oral_reads += count

    oral_frac = oral_reads / total * 100 if total > 0 else 0
    top3 = sorted(genus_reads.items(), key=lambda x: -x[1])[:3]

    results.append({
        'sample': s['sample'], 'sample_base': s['sample_base'],
        'bioproject': s['bioproject'], 'category': category,
        'total_reads': total, 'n_asvs': n_asvs,
        'oral_reads': oral_reads, 'oral_frac': oral_frac,
        'top3_genera': top3,
    })

df = pd.DataFrame(results)
df['category'] = pd.Categorical(df['category'], categories=CATEGORY_ORDER, ordered=True)

# ============================================================
# Apply threshold: n_asvs < 30 AND oral_frac < 20%
# ============================================================
THRESHOLD_ASV = 1000   # hard upper bound: >1000 ASVs is never a genuine control
THRESHOLD_ORAL = 20

df['keep'] = (df['n_asvs'] <= THRESHOLD_ASV) & (df['oral_frac'] < THRESHOLD_ORAL)

# ============================================================
# Summary
# ============================================================
print("\n" + "="*70)
print(f"THRESHOLD: oral genera < {THRESHOLD_ORAL}% AND ASVs ≤ {THRESHOLD_ASV}")
print("="*70)
for cat in CATEGORY_ORDER:
    sub = df[df['category'] == cat]
    if len(sub) == 0:
        continue
    keep_n = sub['keep'].sum()
    excl_n = (~sub['keep']).sum()
    print(f"  {cat:20s}: keep={keep_n:4d}, exclude={excl_n:4d} ({(excl_n/len(sub)*100):4.1f}% excluded)")

total_keep = df['keep'].sum()
total_excl = (~df['keep']).sum()
print(f"\n  OVERALL: keep={total_keep}, exclude={total_excl} ({total_excl/len(df)*100:.1f}% excluded)")

# Show excluded samples with high oral signal
print("\n" + "="*70)
print("EXCLUDED SAMPLES (sorted by oral fraction)")
print("="*70)
excluded = df[~df['keep']].sort_values('oral_frac', ascending=False)
for _, row in excluded.iterrows():
    top3_str = ', '.join([f"{g}({c})" for g, c in row['top3_genera']])
    print(f"  {row['sample_base'][:28]:28s} | {row['category'][:18]:18s} | {row['bioproject']:15s} | "
          f"ASVs={row['n_asvs']:5d} | oral={row['oral_frac']:5.1f}% | reads={row['total_reads']:>8,} | "
          f"top: {top3_str[:100]}")

# ============================================================
# PLOT 1: Scatter plot with threshold lines
# ============================================================
print("\n[Plotting] Scatter: ASV count vs oral fraction (with density margins)...")

# Filter data to visible range
plot_df = df[df['n_asvs'] <= 200].copy()
n_offscale = (df['n_asvs'] > 200).sum()

# Create figure with marginal axes
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4],
                      left=0.08, right=0.95, bottom=0.08, top=0.93,
                      hspace=0.05, wspace=0.05)

ax = fig.add_subplot(gs[1, 0])        # main scatter
ax_histx = fig.add_subplot(gs[0, 0])  # top histogram (x-axis: ASV count)
ax_histy = fig.add_subplot(gs[1, 1])  # right histogram (y-axis: oral frac)

# ---- Main scatter ----
for cat in CATEGORY_ORDER:
    sub = plot_df[plot_df['category'] == cat]
    if len(sub) == 0:
        continue
    kept = sub[sub['keep']]
    excl = sub[~sub['keep']]
    if len(kept) > 0:
        ax.scatter(kept['n_asvs'], kept['oral_frac'],
                  c=CAT_PALETTE[cat], alpha=0.25, s=12,
                  edgecolors='none')
    if len(excl) > 0:
        ax.scatter(excl['n_asvs'], excl['oral_frac'],
                  c=CAT_PALETTE[cat], alpha=0.85, s=40,
                  edgecolors='black', linewidth=0.5, marker='D',
                  label=f"{CATEGORY_LABELS.get(cat, cat)} excluded (n={len(excl)})")

# Threshold line for oral fraction only (ASV threshold at 1000 is off-scale)
ax.axhline(y=THRESHOLD_ORAL, color='#7f7f7f', linestyle='-', linewidth=2, alpha=0.7)

# Shade keep region (oral < 20%, entire x range visible)
ax.fill_between([0, 200], 0, THRESHOLD_ORAL, alpha=0.1, color='#2ca02c')
ax.text(100, THRESHOLD_ORAL/2, 'KEEP\nn=%d' % total_keep,
       ha='center', va='center', fontsize=11, color='#2ca02c', fontweight='bold')

ax.text(100, 65,
       f'EXCLUDE\n(oral ≥ {THRESHOLD_ORAL}% or ASV > {THRESHOLD_ASV})',
       ha='center', va='center', fontsize=9, color='#d62728', fontweight='bold')

# Note samples outside visible range
n_outside = (df['n_asvs'] > 200).sum()
n_excl_offscale = ((df['n_asvs'] > 200) & (~df['keep'])).sum()
if n_outside > 0:
    ax.text(195, 98, f'{n_outside} samples\nASVs > 200\n({n_excl_offscale} excluded)',
           ha='right', va='top', fontsize=8, color='gray')

ax.set_xlim(-2, 200)
ax.set_ylim(-2, 105)
ax.set_xlabel('Number of ASVs')
ax.set_ylabel('Oral Genera Relative Abundance (%)')

handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), fontsize=8, loc='upper right',
         framealpha=0.9, ncol=1)

# ---- Top histogram: ASV count distribution ----
# Use log-spaced bins for ASV count to see the low end clearly
x_bins = np.concatenate([np.linspace(0, 50, 51), np.linspace(55, 200, 30)])

ax_histx.hist(plot_df[plot_df['keep']]['n_asvs'], bins=x_bins,
             color='#2ca02c', alpha=0.5, density=True, label='Kept')
ax_histx.hist(plot_df[~plot_df['keep']]['n_asvs'], bins=x_bins,
             color='#d62728', alpha=0.5, density=True, label='Excluded')

ax_histx.set_xlim(-2, 200)
ax_histx.set_ylabel('Density', fontsize=9)
ax_histx.set_title(f'Control Sample Filtering: oral < {THRESHOLD_ORAL}% & ASV ≤ {THRESHOLD_ASV}\n'
                   f'(n={total_excl}/{len(df)} excluded, ASV threshold at {THRESHOLD_ASV} off-scale →)',
                   fontsize=11)
ax_histx.legend(fontsize=7, loc='upper right')
ax_histx.tick_params(labelbottom=False)

# ---- Right histogram: oral fraction distribution ----
y_bins = np.linspace(0, 100, 51)

ax_histy.hist(df[df['keep']]['oral_frac'].clip(0, 100), bins=y_bins,
             color='#2ca02c', alpha=0.5, density=True, orientation='horizontal',
             label='Kept')
ax_histy.hist(df[~df['keep']]['oral_frac'].clip(0, 100), bins=y_bins,
             color='#d62728', alpha=0.5, density=True, orientation='horizontal',
             label='Excluded')

ax_histy.axhline(y=THRESHOLD_ORAL, color='#7f7f7f', linestyle='-', linewidth=2, alpha=0.7)
ax_histy.set_ylim(-2, 105)
ax_histy.set_xlabel('Density', fontsize=9)
ax_histy.tick_params(labelleft=False)

plt.tight_layout()
fig.savefig(f"{OUT_DIR}/03_oral_vs_asv_scatter.png", dpi=150, bbox_inches='tight')
plt.close()
print("  -> 03_oral_vs_asv_scatter.png")

# ============================================================
# PLOT 2: Per-BioProject keep/exclude bar chart
# ============================================================
print("[Plotting] Per-BioProject keep/exclude bar chart...")

bp_order = sorted(df['bioproject'].unique(),
                  key=lambda b: (~df[df['bioproject']==b]['keep']).sum())

fig, ax = plt.subplots(figsize=(14, 8))
y_pos = np.arange(len(bp_order))

keep_counts = []
excl_counts = []
for bp in bp_order:
    sub = df[df['bioproject'] == bp]
    keep_counts.append(sub['keep'].sum())
    excl_counts.append((~sub['keep']).sum())

# Stacked bars
bar_keep = ax.barh(y_pos, keep_counts, 0.7, color='#2ca02c', alpha=0.85, label='Keep')
bar_excl = ax.barh(y_pos, excl_counts, 0.7, left=keep_counts, color='#d62728', alpha=0.85, label='Exclude')

# Labels
for i, (bp, k, e) in enumerate(zip(bp_order, keep_counts, excl_counts)):
    total = k + e
    if e > 0:
        ax.text(total + total*0.05, i, f'{k}/{total}', va='center', fontsize=8, fontweight='bold',
               color='#d62728' if e > k else '#333333')

ax.set_yticks(y_pos)
ax.set_yticklabels(bp_order, fontsize=9)
ax.set_xlabel('Number of Samples')
ax.set_title(f'Control Sample Filtering by BioProject\n(oral < {THRESHOLD_ORAL}% & ASV ≤ {THRESHOLD_ASV})')
ax.legend(fontsize=10, loc='lower right')

# Add total annotation
ax.text(0.98, 0.02, f'Total: {total_keep} keep / {total_excl} exclude',
       transform=ax.transAxes, fontsize=9, color='gray', ha='right', va='bottom')

plt.tight_layout()
fig.savefig(f"{OUT_DIR}/04_bioproject_keep_exclude.png", dpi=150, bbox_inches='tight')
plt.close()
print("  -> 04_bioproject_keep_exclude.png")

# ============================================================
# Export
# ============================================================
print("\n[Export] Filtered sample lists...")

# Full results
df_export = df.copy()
df_export['top3_str'] = df_export['top3_genera'].apply(
    lambda x: '; '.join([f"{g}:{c}" for g, c in x]))
df_export = df_export.drop(columns=['top3_genera'])
df_export = df_export.sort_values(['keep', 'oral_frac'], ascending=[False, False])
df_export.to_csv(f"{OUT_DIR}/control_classification.csv", index=False, float_format='%.2f')

# Keep list
keep_list = df[df['keep']][['sample', 'sample_base', 'bioproject', 'category',
                              'total_reads', 'n_asvs', 'oral_frac']]
keep_list = keep_list.sort_values(['bioproject', 'n_asvs'])
keep_list.to_csv(f"{OUT_DIR}/genuine_controls_keep.csv", index=False, float_format='%.2f')

# Exclude list
excl_list = df[~df['keep']][['sample', 'sample_base', 'bioproject', 'category',
                               'total_reads', 'n_asvs', 'oral_frac']]
excl_list = excl_list.sort_values('oral_frac', ascending=False)
excl_list.to_csv(f"{OUT_DIR}/controls_excluded.csv", index=False, float_format='%.2f')

print(f"  -> control_classification.csv ({len(df_export)} samples)")
print(f"  -> genuine_controls_keep.csv ({len(keep_list)} samples)")
print(f"  -> controls_excluded.csv ({len(excl_list)} samples)")

print("\n=== Done! ===")
PYEOF

echo ""
echo "=== Filtering done at $(date) ==="
ls -lh "${OUT_DIR}/"control_classification.csv \
      "${OUT_DIR}/"genuine_controls_keep.csv \
      "${OUT_DIR}/"controls_excluded.csv \
      "${OUT_DIR}/"03_oral_vs_asv_scatter.png \
      "${OUT_DIR}/"04_bioproject_keep_exclude.png 2>/dev/null
