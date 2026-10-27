#!/bin/bash
#SBATCH --job-name=ctrl_tax_viz
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/ctrl_tax_viz_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/ctrl_tax_viz_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=120G
#SBATCH --time=12:00:00

# ============================================================
# Control Sample Visualization
# Plot 1: Boxplot of ASV richness per BioProject (with jitter)
# Plot 2: Per-BP stacked bar chart at genus level
# ============================================================

set -e

source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv

VSEARCH="/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"
SILVA_DB="/beegfs/dongbiao/silva_138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta"
BASE_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"
VIZ_DIR="${BASE_DIR}/control_quality_viz"
TEMP_DIR="${VIZ_DIR}/temp"
mkdir -p "${VIZ_DIR}" "${TEMP_DIR}"

echo "=== Control Sample Visualization ==="
echo "Started: $(date)"

# ============================================================
# Step 1: Collect all control ASVs & sample metadata
# ============================================================
echo "[Step 1] Collecting ASVs and sample data..."

python3 << 'PYEOF'
import pandas as pd
import os, json
from collections import defaultdict, Counter

BASE_DIR = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"
TEMP_DIR = f"{BASE_DIR}/control_quality_viz/temp"

# Load control sample metadata
ctrl = pd.read_csv(f"{BASE_DIR}/control_sample_ids.tsv", sep="\t")
ctrl = ctrl[ctrl['category'] != 'mock_community']
sample_cat = dict(zip(ctrl['sample_id'], ctrl['category']))
sample_bp = dict(zip(ctrl['sample_id'], ctrl['bioproject']))

# Collect ASVs
all_asvs = {}
asv_id_to_seq = {}
sample_data = []
asv_id_counter = 0

result_dirs = [f"{BASE_DIR}/control_dada2", f"{BASE_DIR}/control_dada2_extra"]

for result_dir in result_dirs:
    if not os.path.isdir(result_dir):
        continue
    for bioproject in sorted(os.listdir(result_dir)):
        proj_path = os.path.join(result_dir, bioproject)
        if not os.path.isdir(proj_path):
            continue
        counts_file = os.path.join(proj_path, "results", "ASVs_counts.tsv")
        fasta_file = os.path.join(proj_path, "results", "ASVs.fa")
        if not os.path.exists(counts_file) or not os.path.exists(fasta_file):
            continue

        # Parse FASTA
        bioproject_asv_map = {}
        with open(fasta_file) as f:
            header = None; seq_parts = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header is not None:
                        bioproject_asv_map[header.lstrip(">")] = "".join(seq_parts).upper()
                    header = line; seq_parts = []
                else:
                    seq_parts.append(line)
            if header is not None:
                bioproject_asv_map[header.lstrip(">")] = "".join(seq_parts).upper()

        # Register globally
        local_to_global = {}
        for local_id, seq in bioproject_asv_map.items():
            if seq not in all_asvs:
                all_asvs[seq] = f"ASV_{asv_id_counter:06d}"
                asv_id_to_seq[f"ASV_{asv_id_counter:06d}"] = seq
                asv_id_counter += 1
            local_to_global[local_id] = all_asvs[seq]

        # Parse counts
        counts_df = pd.read_csv(counts_file, sep="\t", index_col=0)
        for sample_name in counts_df.columns:
            read_counts = counts_df[sample_name]
            nonzero = read_counts[read_counts > 0]
            total_reads = int(nonzero.sum())
            n_asvs = len(nonzero)

            global_counts = {}
            for local_asv, count in nonzero.items():
                if local_asv in local_to_global:
                    global_counts[local_to_global[local_asv]] = int(count)

            sample_base = sample_name.replace(".fastq.gz", "").replace(".R1.filtered.fastq.gz", "")
            for suffix in ["_1", "_2", "_R1", "_R2"]:
                if sample_base.endswith(suffix):
                    sample_base = sample_base[:-len(suffix)]

            category = sample_cat.get(sample_base, "unknown")
            bp = sample_bp.get(sample_base, bioproject)

            sample_data.append({
                'sample': sample_name,
                'sample_base': sample_base,
                'bioproject': bp,
                'category': category,
                'total_reads': total_reads,
                'n_asvs': n_asvs,
                'global_counts': global_counts,
            })

print(f"Total unique ASVs: {len(all_asvs)}")
print(f"Total samples: {len(sample_data)}")

# Print per-BP summary
bp_summary = defaultdict(lambda: {'n_samples': 0, 'total_asvs': 0, 'categories': set()})
for s in sample_data:
    bp = s['bioproject']
    bp_summary[bp]['n_samples'] += 1
    bp_summary[bp]['total_asvs'] += s['n_asvs']
    bp_summary[bp]['categories'].add(s['category'])

print("\nBioProject summary:")
for bp in sorted(bp_summary.keys()):
    info = bp_summary[bp]
    cats = ','.join(sorted(info['categories']))
    print(f"  {bp}: {info['n_samples']} samples, {info['total_asvs']} total ASVs, categories: {cats}")

# Save ASV fasta
with open(f"{TEMP_DIR}/all_control_asvs.fa", "w") as f:
    for asv_id, seq in sorted(asv_id_to_seq.items()):
        f.write(f">{asv_id}\n{seq}\n")
print(f"\nWrote {len(asv_id_to_seq)} ASVs")

# Save sample data
serializable = []
for s in sample_data:
    sd = {k: v for k, v in s.items() if k != 'global_counts'}
    sd['asv_ids'] = list(s['global_counts'].keys())
    sd['asv_counts'] = s['global_counts']
    serializable.append(sd)
with open(f"{TEMP_DIR}/sample_data.json", "w") as f:
    json.dump(serializable, f)

print("Step 1 done.")
PYEOF

# ============================================================
# Step 2: Taxonomic assignment via vsearch against SILVA
# ============================================================
echo ""
echo "[Step 2] Running taxonomic assignment against SILVA..."

HITS_FILE="${TEMP_DIR}/silva_hits.txt"
if [ -s "$HITS_FILE" ] && [ "$(wc -l < "$HITS_FILE")" -gt 0 ]; then
    echo "[skip] SILVA hits already exist: $(wc -l < "$HITS_FILE") hits"
else
    "$VSEARCH" \
        --usearch_global "${TEMP_DIR}/all_control_asvs.fa" \
        --db "$SILVA_DB" \
        --id 0.0 \
        --iddef 1 \
        --maxaccepts 1 \
        --maxrejects 32 \
        --userout "$HITS_FILE" \
        --userfields "query+target+id+alnlen+mism+gaps+qlo+qhi+tlo+thi+ql+tl+qcov+tcov" \
        --threads 28 \
        --maxhits 1 \
        --quiet
    echo "SILVA hits: $(wc -l < "$HITS_FILE")"
fi

# ============================================================
# Step 3: Parse taxonomy & generate plots
# ============================================================
echo ""
echo "[Step 3] Parsing taxonomy and generating plots..."

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import json, os, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"
TEMP_DIR = f"{BASE_DIR}/control_quality_viz/temp"
VIZ_DIR = f"{BASE_DIR}/control_quality_viz"

sns.set_style("whitegrid")
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 13, 'axes.labelsize': 11,
                     'figure.dpi': 150, 'savefig.dpi': 150, 'savefig.bbox': 'tight'})

CATEGORY_ORDER = ['extraction_blank', 'blank', 'ntc_no_template', 'negative_control', 'reagent_water']
CATEGORY_LABELS = {
    'extraction_blank': 'Extraction Blank', 'blank': 'Blank',
    'ntc_no_template': 'NTC', 'negative_control': 'Negative Control',
    'reagent_water': 'Reagent Water'
}
CAT_PALETTE = {
    'extraction_blank': '#d62728', 'blank': '#1f77b4',
    'ntc_no_template': '#2ca02c', 'negative_control': '#ff7f0e',
    'reagent_water': '#9467bd', 'unknown': '#7f7f7f'
}

# ============================================================
# Load & parse SILVA taxonomy
# ============================================================
print("Parsing SILVA taxonomy headers...")
silva_tax = {}
with open("/beegfs/dongbiao/silva_138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta") as f:
    for line in f:
        if line.startswith(">"):
            full_line = line.strip().lstrip(">")
            space_idx = full_line.find(" ")
            if space_idx > 0:
                tax_str = full_line[space_idx+1:]
                seq_id = full_line[:space_idx]
                silva_tax[seq_id] = tax_str
print(f"  Parsed {len(silva_tax):,} SILVA entries")

# Load SILVA hits
cols = ['query', 'target', 'id', 'alnlen', 'mism', 'gaps',
        'qlo', 'qhi', 'tlo', 'thi', 'ql', 'tl', 'qcov', 'tcov']
hits = pd.read_csv(f"{TEMP_DIR}/silva_hits.txt", sep="\t", names=cols)
hits['true_id'] = (hits['alnlen'] - hits['mism'] - hits['gaps']) / hits['ql'] * 100

# Parse taxonomy for each ASV (best hit ≥80% id)
def parse_taxonomy(target_id, silva_tax_map):
    tax_str = silva_tax_map.get(target_id, "Unknown")
    ranks = [r.strip() for r in tax_str.split(";")]
    result = {}
    for i, rank_name in enumerate(['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
        result[rank_name] = ranks[i] if i < len(ranks) else 'Unknown'
    return result

asv_taxonomy = {}
for _, row in hits.iterrows():
    asv_id = row['query']
    if row['true_id'] < 80:
        continue
    if asv_id not in asv_taxonomy:
        asv_taxonomy[asv_id] = parse_taxonomy(row['target'], silva_tax)

print(f"ASVs with taxonomy (≥80% id): {len(asv_taxonomy)}")

# Load sample data
with open(f"{TEMP_DIR}/sample_data.json") as f:
    sample_data = json.load(f)

# ============================================================
# Build sample-level data with genus composition
# ============================================================
df_rows = []
for s in sample_data:
    asv_counts = s['asv_counts']
    genus_reads = defaultdict(int)
    unclassified = 0
    for asv_id, count in asv_counts.items():
        if asv_id in asv_taxonomy:
            genus_reads[asv_taxonomy[asv_id]['genus']] += count
        else:
            unclassified += count

    df_rows.append({
        'sample': s['sample'], 'sample_base': s['sample_base'],
        'bioproject': s['bioproject'], 'category': s['category'],
        'total_reads': s['total_reads'], 'n_asvs': s['n_asvs'],
        'genus_reads': dict(genus_reads), 'unclassified_reads': unclassified,
    })

df = pd.DataFrame(df_rows)
df['category'] = pd.Categorical(df['category'], categories=CATEGORY_ORDER, ordered=True)

# ============================================================
# Global top genera
# ============================================================
all_genus_reads = defaultdict(int)
for _, row in df.iterrows():
    for g, c in row['genus_reads'].items():
        all_genus_reads[g] += c
top_genera_global = [g for g, _ in sorted(all_genus_reads.items(), key=lambda x: -x[1])[:15]]

# Distinct colors for top genera
GENUS_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#bcbd22', '#17becf', '#aec7e8',
    '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94',
]

# ============================================================
# PLOT 1: Boxplot of ASV richness per BioProject
# ============================================================
print("\n[Plot 1] ASV richness per BioProject...")

# Order BioProjects by median ASV richness
bp_order = df.groupby('bioproject')['n_asvs'].median().sort_values().index.tolist()
present_cats = [c for c in CATEGORY_ORDER if len(df[df['category']==c]) > 0]

fig, ax = plt.subplots(figsize=(20, 8))

data_by_bp = [df[df['bioproject']==bp]['n_asvs'].values for bp in bp_order]
# Color each box by the dominant category
bp_dominant_cat = []
for bp in bp_order:
    sub = df[df['bioproject']==bp]
    dom_cat = sub['category'].value_counts().index[0]
    bp_dominant_cat.append(dom_cat)

bp = ax.boxplot(data_by_bp, patch_artist=True, widths=0.6,
                medianprops={'color': 'black', 'linewidth': 1.2},
                flierprops={'marker': 'o', 'markersize': 3, 'alpha': 0.5})

for i, (box, dom_cat) in enumerate(zip(bp['boxes'], bp_dominant_cat)):
    box.set_facecolor(CAT_PALETTE.get(dom_cat, '#7f7f7f'))
    box.set_alpha(0.75)

# Add jittered points
for i, (bp_name, data) in enumerate(zip(bp_order, data_by_bp)):
    jitter = np.random.normal(0, 0.08, size=len(data))
    sub = df[df['bioproject']==bp_name]
    colors = [CAT_PALETTE.get(c, '#7f7f7f') for c in sub['category']]
    ax.scatter(np.full(len(data), i+1) + jitter, data, alpha=0.35, s=6,
              c=colors, edgecolors='none', zorder=10)

# Add n labels
for i, (bp_name, data) in enumerate(zip(bp_order, data_by_bp)):
    ax.annotate(f'n={len(data)}', xy=(i+1, ax.get_ylim()[1]*0.97),
               ha='center', fontsize=6, color='gray')

ax.set_xticks(np.arange(1, len(bp_order)+1))
ax.set_xticklabels(bp_order, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Number of ASVs')
ax.set_title(f'ASV Richness per BioProject (Control Samples, n={len(df)} samples)')
ax.set_yscale('log')

# Category legend
cat_patches = [mpatches.Patch(color=CAT_PALETTE[c], label=CATEGORY_LABELS.get(c, c))
               for c in present_cats]
ax.legend(handles=cat_patches, fontsize=8, loc='upper left', title='Control Type')

plt.tight_layout()
fig.savefig(f"{VIZ_DIR}/01_asv_richness_per_bioproject.png", dpi=150, bbox_inches='tight')
plt.close()
print("  -> 01_asv_richness_per_bioproject.png saved")

# ============================================================
# PLOT 2: Per-BioProject stacked bar chart (genus level)
# ============================================================
print("[Plot 2] Per-BioProject genus composition stacked bars...")

# For each BP, determine which top genera to show (top 10 within that BP)
for bp_name in bp_order:
    sub = df[df['bioproject']==bp_name]
    if len(sub) == 0:
        continue

    # Aggregate genus reads for this BP
    bp_genus = defaultdict(int)
    bp_unclassified = 0
    for _, row in sub.iterrows():
        for g, c in row['genus_reads'].items():
            bp_genus[g] += c
        bp_unclassified += row['unclassified_reads']

    # Top 10 genera in this BP + Unclassified + Other
    top_bp_genera = [g for g, _ in sorted(bp_genus.items(), key=lambda x: -x[1])[:10]]
    show_genera = top_bp_genera + ['Unclassified']

    # Build color map
    genus_color_map = {}
    for i, g in enumerate(top_bp_genera):
        if g in top_genera_global:
            idx = top_genera_global.index(g)
            genus_color_map[g] = GENUS_COLORS[idx % len(GENUS_COLORS)]
        else:
            genus_color_map[g] = GENUS_COLORS[i % len(GENUS_COLORS)]
    genus_color_map['Unclassified'] = '#d3d3d3'

    # Sort samples by total reads (descending)
    sub_sorted = sub.sort_values('total_reads', ascending=False)
    n_samples = len(sub_sorted)

    # Build matrix: samples x genera
    sample_labels = []
    sample_colors_cat = []
    matrix = np.zeros((n_samples, len(show_genera)))

    for i, (_, row) in enumerate(sub_sorted.iterrows()):
        # Truncate long sample names
        label = row['sample_base']
        if len(label) > 18:
            # Try to get meaningful part
            label = label[:15] + '...'
        sample_labels.append(label)
        sample_colors_cat.append(CAT_PALETTE.get(row['category'], '#7f7f7f'))
        total = row['total_reads']
        if total == 0:
            continue
        for j, genus in enumerate(show_genera):
            if genus == 'Unclassified':
                matrix[i, j] = row['unclassified_reads'] / total * 100
            else:
                matrix[i, j] = row['genus_reads'].get(genus, 0) / total * 100

    # Calculate figure height based on number of samples
    fig_height = max(6, n_samples * 0.18)
    fig, ax = plt.subplots(figsize=(16, fig_height))

    # Stacked bars
    bottoms = np.zeros(n_samples)
    bar_handles = []
    for j, genus in enumerate(show_genera):
        bars = ax.barh(np.arange(n_samples), matrix[:, j], 0.85,
                      left=bottoms, color=genus_color_map[genus],
                      label=genus, alpha=0.85, edgecolor='none')
        bar_handles.append(bars)
        bottoms += matrix[:, j]

    # Sample labels colored by category
    for i, (label, cat_color) in enumerate(zip(sample_labels, sample_colors_cat)):
        ax.text(101, i, label, va='center', fontsize=6, color=cat_color)

    ax.set_yticks([])
    ax.set_xlabel('Relative Abundance (%)')
    ax.set_title(f'{bp_name} — Genus-Level Composition ({n_samples} control samples)')
    ax.set_xlim(0, 120)  # Extra space for labels

    # Legend
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=8, title='Genus',
             title_fontsize=9, ncol=1)

    # Add a small category legend
    cat_text = ', '.join([f"{CATEGORY_LABELS.get(c, c)}={len(sub[sub['category']==c])}"
                         for c in present_cats if len(sub[sub['category']==c]) > 0])
    ax.text(0.5, -0.01, f'Control types: {cat_text}', transform=ax.transAxes,
           fontsize=7, color='gray', ha='center')

    plt.tight_layout()
    fig.savefig(f"{VIZ_DIR}/02_genus_bars_{bp_name}.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> 02_genus_bars_{bp_name}.png ({n_samples} samples)")

# ============================================================
# Summary statistics table
# ============================================================
print("\n[Export] Summary statistics...")

summary_rows = []
for bp_name in bp_order:
    sub = df[df['bioproject']==bp_name]
    summary_rows.append({
        'bioproject': bp_name,
        'n_samples': len(sub),
        'categories': ','.join(sorted(str(c) for c in sub['category'].unique())),
        'median_asvs': sub['n_asvs'].median(),
        'mean_asvs': sub['n_asvs'].mean(),
        'min_asvs': sub['n_asvs'].min(),
        'max_asvs': sub['n_asvs'].max(),
        'median_reads': sub['total_reads'].median(),
        'mean_reads': sub['total_reads'].mean(),
        'mean_pct_unclassified': sub.apply(
            lambda r: r['unclassified_reads']/r['total_reads']*100 if r['total_reads']>0 else 0, axis=1
        ).mean(),
    })
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(f"{VIZ_DIR}/bioproject_summary.csv", index=False, float_format='%.1f')
print(f"  -> bioproject_summary.csv")

# Per-sample export
df_export = df[['sample', 'sample_base', 'bioproject', 'category', 'total_reads', 'n_asvs']].copy()
df_export['pct_unclassified'] = df.apply(
    lambda r: r['unclassified_reads']/r['total_reads']*100 if r['total_reads']>0 else 0, axis=1
)
df_export = df_export.sort_values(['bioproject', 'total_reads'], ascending=[True, False])
df_export.to_csv(f"{VIZ_DIR}/sample_details.csv", index=False, float_format='%.2f')
print(f"  -> sample_details.csv")

print("\n=== All done! ===")
print(f"Figures saved in: {VIZ_DIR}")
import glob
for f in sorted(glob.glob(f"{VIZ_DIR}/*.png")):
    size_mb = os.path.getsize(f) / 1024 / 1024
    print(f"  {os.path.basename(f)} ({size_mb:.1f} MB)")
PYEOF

echo ""
echo "=== Finished at $(date) ==="
ls -lh "${VIZ_DIR}/"*.png 2>/dev/null | head -30
ls -lh "${VIZ_DIR}/"*.csv 2>/dev/null
