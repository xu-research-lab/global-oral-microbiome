#!/usr/bin/env python3
"""
Control Sample QC v2 — Stricter multi-dimensional filtering for genuine blanks.

Improvements over v1 (C08_filter_genuine_controls.sbatch):
- Adds total_reads threshold (genuine blanks must have low biomass)
- Tightens n_asvs threshold (≤100 instead of ≤1000)
- Excludes entire BioProjects where ALL controls are high-read (PRJEB5726)
- Excludes BioProjects lacking oral samples (PRJNA649980, PRJNA662100)
- Records explicit per-sample exclusion reasons
- Multi-panel publication-quality figure

Thresholds:
  total_reads  : 1 – 5,000  (too low = empty file; too high = real sample contamination)
  n_asvs       : ≤ 100       (genuine blanks have few taxa)
  oral_frac    : < 20%       (high oral signal → sample swap / mislabel)
  Per-study    : ≥ 3 reliable controls + ≥ 5 oral samples in metadata

Inputs (from C07_control_quality_viz):
  - control_quality_viz/temp/sample_data.json
  - control_quality_viz/temp/silva_hits.txt
  - control_quality_viz/temp/asv_mapping.json (optional, for FASTA export)

Outputs:
  - control_quality_viz/genuine_controls_keep_v2.csv
  - control_quality_viz/controls_excluded_v2.csv
  - control_quality_viz/control_classification_v2.csv
  - control_quality_viz/05_control_qc_v2.png  (multi-panel figure)
"""

import pandas as pd
import numpy as np
import json
import os
import sys
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Config
# ============================================================
BASE_DIR = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
TEMP_DIR = BASE_DIR / "control_quality_viz" / "temp"
OUT_DIR = BASE_DIR / "control_quality_viz"
ORAL_METADATA_CSV = "/beegfs/huangxiaochang/global_oral/metadata_v2_merge/metadata-v2-merge-oral.csv"
SILVA_FASTA = "/beegfs/dongbiao/silva_138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta"

# ---- Thresholds ----
MIN_READS = 1          # exclude empty files
MAX_READS = 5_000      # genuine blanks have low biomass
MAX_ASVS = 100          # genuine blanks have few taxa
MAX_ORAL_FRAC = 20.0    # % oral genera

# ---- Study-level filters ----
MIN_CONTROLS_PER_STUDY = 3   # need enough blanks for decontam
MIN_ORAL_PER_STUDY = 5       # need enough oral samples for comparison

# ---- Excluded BioProjects (all controls are high-read or no oral samples) ----
EXCLUDED_BIOPROJECTS = {
    "PRJEB5726": "All 29 controls >10K reads (extraction_blank/NTC/reagent_water median ~45K) — clearly not genuine blanks",
    "PRJNA649980": "No oral samples in metadata (5 blanks only)",
    "PRJNA662100": "No oral samples in metadata (2 blanks only)",
}

# ---- Oral genera for oral_frac calculation ----
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

# ---- Plotting ----
CATEGORY_ORDER = ['extraction_blank', 'blank', 'ntc_no_template', 'negative_control', 'reagent_water']
CATEGORY_LABELS = {
    'extraction_blank': 'Extraction Blank', 'blank': 'Blank',
    'ntc_no_template': 'NTC', 'negative_control': 'Negative Control',
    'reagent_water': 'Reagent Water',
}
CAT_PALETTE = {
    'extraction_blank': '#d62728', 'blank': '#1f77b4',
    'ntc_no_template': '#2ca02c', 'negative_control': '#ff7f0e',
    'reagent_water': '#9467bd', 'unknown': '#7f7f7f',
}
CAT_MARKERS = {
    'extraction_blank': 's', 'blank': 'o',
    'ntc_no_template': '^', 'negative_control': 'D',
    'reagent_water': 'v',
}

# Visual settings
plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 10,
    'figure.dpi': 200, 'savefig.dpi': 200, 'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})


def log(msg: str) -> None:
    print(f"[QC-v2] {msg}", flush=True)


# ============================================================
# 1. Load SILVA taxonomy
# ============================================================
def load_silva_taxonomy(hits_file: Path, silva_fasta: str) -> dict:
    """Build ASV ID -> {genus, phylum, ...} from SILVA hits."""
    log("Loading SILVA taxonomy...")

    cols = ['query', 'target', 'id', 'alnlen', 'mism', 'gaps',
            'qlo', 'qhi', 'tlo', 'thi', 'ql', 'tl', 'qcov', 'tcov']
    hits = pd.read_csv(hits_file, sep="\t", names=cols)
    hits['true_id'] = (hits['alnlen'] - hits['mism'] - hits['gaps']) / hits['ql'] * 100

    # Parse SILVA headers: >accession taxonomy_string
    silva_tax = {}
    log(f"  Parsing SILVA headers from {silva_fasta}...")
    with open(silva_fasta) as f:
        for line in f:
            if line.startswith(">"):
                full_line = line.strip().lstrip(">")
                space_idx = full_line.find(" ")
                if space_idx > 0:
                    silva_tax[full_line[:space_idx]] = full_line[space_idx + 1:]
    log(f"  Parsed {len(silva_tax):,} SILVA entries")

    def parse_taxonomy(target_id):
        tax_str = silva_tax.get(target_id, "Unknown")
        ranks = [r.strip() for r in tax_str.split(";")]
        result = {}
        for i, rn in enumerate(['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
            result[rn] = ranks[i] if i < len(ranks) else 'Unknown'
        return result

    asv_taxonomy = {}
    for _, row in hits.iterrows():
        if row['true_id'] >= 80 and row['query'] not in asv_taxonomy:
            asv_taxonomy[row['query']] = parse_taxonomy(row['target'])

    log(f"  ASVs with taxonomy (≥80% id): {len(asv_taxonomy):,}")
    return asv_taxonomy


# ============================================================
# 2. Load & classify samples
# ============================================================
def classify_controls(sample_data: list, asv_taxonomy: dict) -> pd.DataFrame:
    """Compute oral_frac and apply filtering thresholds."""
    log("Classifying control samples...")

    rows = []
    for s in sample_data:
        total = s['total_reads']
        n_asvs = s['n_asvs']
        category = s.get('category', 'unknown')

        if category not in CATEGORY_ORDER and category != 'mock_community':
            continue
        if category == 'mock_community':
            continue

        oral_reads = 0
        genus_reads = defaultdict(int)

        for asv_id, count in s['asv_counts'].items():
            genus = asv_taxonomy.get(asv_id, {}).get('genus', 'Unclassified')
            genus_reads[genus] += count
            if genus in ORAL_GENERA:
                oral_reads += count

        oral_frac = oral_reads / total * 100 if total > 0 else 0.0
        top3 = sorted(genus_reads.items(), key=lambda x: -x[1])[:3]

        rows.append({
            'sample': s['sample'],
            'sample_base': s['sample_base'],
            'bioproject': s['bioproject'],
            'category': category,
            'total_reads': total,
            'n_asvs': n_asvs,
            'oral_reads': oral_reads,
            'oral_frac': round(oral_frac, 2),
            'top3_genera': top3,
        })

    df = pd.DataFrame(rows)
    df['category'] = pd.Categorical(df['category'], categories=CATEGORY_ORDER, ordered=True)
    log(f"  Total control samples (non-mock): {len(df)}")
    return df


# ============================================================
# 3. Apply multi-dimensional filters
# ============================================================
def apply_filters(df: pd.DataFrame) -> pd.DataFrame:
    """Apply filtering and record per-sample exclusion reasons."""
    log("Applying filters...")

    reasons = []  # list of reasons (empty = keep)
    study_excluded = set()

    for _, row in df.iterrows():
        bp = row['bioproject']
        sample_reasons = []

        # Study-level exclusion (override everything)
        if bp in EXCLUDED_BIOPROJECTS:
            reasons.append(f"STUDY_EXCLUDED: {EXCLUDED_BIOPROJECTS[bp]}")
            study_excluded.add(bp)
            continue

        # Total reads
        if row['total_reads'] < MIN_READS:
            sample_reasons.append(f"reads={row['total_reads']}<{MIN_READS}")
        elif row['total_reads'] > MAX_READS:
            sample_reasons.append(f"reads={row['total_reads']}>{MAX_READS}")

        # ASV count
        if row['n_asvs'] > MAX_ASVS:
            sample_reasons.append(f"ASVs={row['n_asvs']}>{MAX_ASVS}")

        # Oral fraction
        if row['oral_frac'] >= MAX_ORAL_FRAC:
            sample_reasons.append(f"oral_frac={row['oral_frac']:.1f}%>={MAX_ORAL_FRAC}%")

        if sample_reasons:
            reasons.append("; ".join(sample_reasons))
        else:
            reasons.append("")  # empty = keep

    df['exclusion_reason'] = reasons
    df['keep'] = df['exclusion_reason'] == ""

    # Summary
    n_keep = df['keep'].sum()
    n_excl = (~df['keep']).sum()
    n_study_excl = len([r for r in reasons if r.startswith("STUDY_EXCLUDED")])

    log(f"  KEPT: {n_keep}  |  EXCLUDED: {n_excl}")
    log(f"    of which study-level exclusions: {n_study_excl}")
    log(f"    remaining individual exclusions: {n_excl - n_study_excl}")
    log(f"  Excluded BioProjects: {sorted(study_excluded)}")

    return df


# ============================================================
# 4. Post-filter: per-study viability
# ============================================================
def check_study_viability(df: pd.DataFrame) -> pd.DataFrame:
    """Check if kept controls meet per-study minimums, load oral sample counts."""
    log("Checking per-study viability...")

    oral_md = pd.read_csv(ORAL_METADATA_CSV, low_memory=False)
    # Count oral samples per BioProject in metadata
    oral_counts = oral_md.groupby('bioproject')['sample_id'].nunique().to_dict()

    kept = df[df['keep']]
    study_summary = []

    for bp in sorted(kept['bioproject'].unique()):
        sub = kept[kept['bioproject'] == bp]
        n_ctrl = len(sub)
        n_oral = oral_counts.get(bp, 0)
        viable = n_ctrl >= MIN_CONTROLS_PER_STUDY and n_oral >= MIN_ORAL_PER_STUDY
        study_summary.append({
            'bioproject': bp,
            'n_reliable_controls': n_ctrl,
            'n_oral_in_metadata': n_oral,
            'viable_for_decontam': viable,
        })

    study_df = pd.DataFrame(study_summary)
    viable_bps = set(study_df[study_df['viable_for_decontam']]['bioproject'])
    nonviable_bps = set(study_df[~study_df['viable_for_decontam']]['bioproject'])

    log(f"  Studies with ≥{MIN_CONTROLS_PER_STUDY} controls & ≥{MIN_ORAL_PER_STUDY} oral: {len(viable_bps)}")
    if nonviable_bps:
        log(f"  Non-viable studies: {sorted(nonviable_bps)}")

    # Add viability column to main df
    bp_viable_map = dict(zip(study_df['bioproject'], study_df['viable_for_decontam']))
    df['study_viable'] = df['bioproject'].map(bp_viable_map).fillna(False)
    df['n_oral_in_metadata'] = df['bioproject'].map(oral_counts).fillna(0).astype(int)

    return df, study_df


# ============================================================
# 5. Multi-panel figure
# ============================================================
def make_figure(df: pd.DataFrame, study_df: pd.DataFrame):
    """Generate a 3-panel summary figure."""
    log("Generating figure...")

    kept = df[df['keep']]
    excluded = df[~df['keep']]

    # Prepare plot data (cap for display)
    plot_df = df.copy()
    plot_df['display_reads'] = plot_df['total_reads'].clip(1, MAX_READS * 2)
    plot_df['display_asvs'] = plot_df['n_asvs'].clip(0.5, MAX_ASVS * 2 + 1)

    fig = plt.figure(figsize=(20, 16))
    gs = GridSpec(2, 3, figure=fig, width_ratios=[3, 1.2, 1.8],
                  height_ratios=[1, 3], hspace=0.28, wspace=0.25,
                  left=0.06, right=0.96, bottom=0.06, top=0.94)

    # ============================
    # Panel A: Main scatter (large, center-left)
    # ============================
    ax_main = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_main)

    # --- Main scatter ---
    # Plot excluded first (background), kept on top
    for cat in CATEGORY_ORDER:
        sub_excl = plot_df[(~plot_df['keep']) & (plot_df['category'] == cat)]
        sub_kept = plot_df[plot_df['keep'] & (plot_df['category'] == cat)]

        color = CAT_PALETTE.get(cat, '#7f7f7f')
        marker = CAT_MARKERS.get(cat, 'o')

        if len(sub_excl) > 0:
            ax_main.scatter(sub_excl['n_asvs'], sub_excl['oral_frac'],
                          c=color, alpha=0.15, s=25, marker=marker,
                          edgecolors='none', zorder=2)
        if len(sub_kept) > 0:
            ax_main.scatter(sub_kept['n_asvs'], sub_kept['oral_frac'],
                          c=color, alpha=0.7, s=18, marker=marker,
                          edgecolors='none', zorder=5)

    # Threshold lines
    ax_main.axvline(x=MAX_ASVS, color='#e74c3c', linestyle='--', linewidth=2, alpha=0.8,
                    label=f'n_ASVs = {MAX_ASVS}')
    ax_main.axhline(y=MAX_ORAL_FRAC, color='#e74c3c', linestyle='--', linewidth=2, alpha=0.8,
                    label=f'Oral ≤ {MAX_ORAL_FRAC}%')

    # Shade the "keep" region
    ax_main.fill_between([0, MAX_ASVS], 0, MAX_ORAL_FRAC,
                        alpha=0.08, color='#27ae60', zorder=0)
    ax_main.text(MAX_ASVS / 2, MAX_ORAL_FRAC / 2,
                f'KEEP ZONE\nn={df["keep"].sum():,}',
                ha='center', va='center', fontsize=13, fontweight='bold',
                color='#27ae60', alpha=0.7)

    # Exclude zones labels
    ax_main.text(MAX_ASVS * 1.5, MAX_ORAL_FRAC / 2,
                f'High ASVs\n(n={(df["n_asvs"] > MAX_ASVS).sum():,})',
                ha='center', va='center', fontsize=9, color='#c0392b', alpha=0.6)
    ax_main.text(MAX_ASVS / 2, MAX_ORAL_FRAC + (100 - MAX_ORAL_FRAC) / 2,
                f'High Oral\n(n={(df["oral_frac"] >= MAX_ORAL_FRAC).sum():,} excl.)',
                ha='center', va='center', fontsize=9, color='#c0392b', alpha=0.6)

    # Read-size annotation: mark points above the read threshold
    high_read_mask = (plot_df['total_reads'] > MAX_READS) & (plot_df['keep'])
    if high_read_mask.sum() > 0:
        # These are already excluded by reads, just note it
        pass

    ax_main.set_xscale('log')
    ax_main.set_xlim(0.8, max(1000, plot_df['n_asvs'].max() * 1.2))
    ax_main.set_ylim(-3, 103)
    ax_main.set_xlabel('Number of ASVs (log scale)', fontsize=11)
    ax_main.set_ylabel('Oral Genera Relative Abundance (%)', fontsize=11)

    # Legend: category markers + threshold lines
    cat_handles = []
    for cat in CATEGORY_ORDER:
        if len(plot_df[plot_df['category'] == cat]) > 0:
            cat_handles.append(
                mpatches.Patch(color=CAT_PALETTE[cat], alpha=0.7,
                              label=f"{CATEGORY_LABELS.get(cat, cat)} (n={len(plot_df[plot_df['category']==cat]):,})")
            )
    threshold_handles = [
        plt.Line2D([0], [0], color='#e74c3c', linestyle='--', linewidth=2, alpha=0.8,
                   label=f'Thresholds: ASV≤{MAX_ASVS}, Oral<{MAX_ORAL_FRAC}%'),
        plt.Line2D([0], [0], color='none', label=f'Reads: {MIN_READS}–{MAX_READS:,}'),
    ]
    legend1 = ax_main.legend(handles=cat_handles, fontsize=8, loc='upper left',
                           title='Control Type', title_fontsize=9, framealpha=0.9)
    ax_main.add_artist(legend1)
    ax_main.legend(handles=threshold_handles, fontsize=8, loc='upper right',
                  framealpha=0.9)

    # --- Top histogram (ASV count distribution) ---
    x_bins = np.logspace(0, np.log10(max(1000, plot_df['n_asvs'].max())), 60)
    ax_histx.hist(plot_df[plot_df['keep']]['n_asvs'], bins=x_bins,
                 color='#27ae60', alpha=0.6, density=False, label=f'Kept (n={df["keep"].sum():,})')
    ax_histx.hist(plot_df[~plot_df['keep']]['n_asvs'], bins=x_bins,
                 color='#e74c3c', alpha=0.5, density=False, label=f'Excluded (n={(~df["keep"]).sum():,})')
    ax_histx.axvline(x=MAX_ASVS, color='#e74c3c', linestyle='--', linewidth=1.5, alpha=0.7)
    ax_histx.set_xscale('log')
    ax_histx.set_ylabel('Count', fontsize=9)
    ax_histx.set_title(
        f'Control Sample QC v2 — Multidimensional Filtering\n'
        f'(reads∈[{MIN_READS},{MAX_READS:,}], ASVs≤{MAX_ASVS}, oral<{MAX_ORAL_FRAC}%, '
        f'≥{MIN_CONTROLS_PER_STUDY} controls & ≥{MIN_ORAL_PER_STUDY} oral per study)',
        fontsize=12, fontweight='bold')
    ax_histx.legend(fontsize=7, loc='upper right')
    ax_histx.tick_params(labelbottom=False)

    # ============================
    # Panel B: Per-BioProject stacked bar (right column, full height)
    # ============================
    ax_bp = fig.add_subplot(gs[:, 2])

    # Build per-BP keep/exclude counts
    bp_order_data = []
    for bp in sorted(df['bioproject'].unique()):
        sub = df[df['bioproject'] == bp]
        n_keep = sub['keep'].sum()
        n_excl_reads = (sub['exclusion_reason'].str.contains('reads=', na=False) & ~sub['exclusion_reason'].str.contains('STUDY_EXCLUDED', na=False)).sum()
        n_excl_asvs = (sub['exclusion_reason'].str.contains('ASVs=', na=False) & ~sub['exclusion_reason'].str.contains('STUDY_EXCLUDED', na=False)).sum()
        n_excl_oral = (sub['exclusion_reason'].str.contains('oral_frac=', na=False) & ~sub['exclusion_reason'].str.contains('STUDY_EXCLUDED', na=False)).sum()
        n_excl_study = sub['exclusion_reason'].str.contains('STUDY_EXCLUDED', na=False).sum()
        # Count those excluded by multiple reasons (counted once in other categories)
        n_excl_multi = max(0, (~sub['keep']).sum() - n_excl_reads - n_excl_asvs - n_excl_oral - n_excl_study)
        bp_order_data.append({
            'bioproject': bp,
            'keep': n_keep,
            'excl_reads': n_excl_reads,
            'excl_asvs': n_excl_asvs,
            'excl_oral': n_excl_oral,
            'excl_study': n_excl_study,
            'excl_multi': n_excl_multi,
            'total': len(sub),
        })

    bp_df = pd.DataFrame(bp_order_data)
    # Sort by total excluded (descending)
    bp_df['n_excl_total'] = bp_df['total'] - bp_df['keep']
    bp_df = bp_df.sort_values('n_excl_total', ascending=True)

    y_positions = np.arange(len(bp_df))
    bar_height = 0.7

    # Stack: keep (green), excl_study (dark red), excl_reads (orange), excl_asvs (purple),
    #         excl_oral (pink), excl_multi (gray)
    left = np.zeros(len(bp_df))
    bars_data = [
        ('keep', '#27ae60', 'Kept'),
        ('excl_study', '#8B0000', 'Study Excluded'),
        ('excl_reads', '#e67e22', f'Reads >{MAX_READS}'),
        ('excl_asvs', '#8e44ad', f'ASVs >{MAX_ASVS}'),
        ('excl_oral', '#e74c3c', f'Oral ≥{MAX_ORAL_FRAC}%'),
        ('excl_multi', '#95a5a6', 'Multiple Reasons'),
    ]

    for col, color, label in bars_data:
        vals = bp_df[col].values
        ax_bp.barh(y_positions, vals, bar_height, left=left, color=color,
                  alpha=0.85, label=label, edgecolor='none')
        left += vals

    # Annotate with keep/total
    for i, (_, row) in enumerate(bp_df.iterrows()):
        total_val = int(row['total'])
        keep_val = int(row['keep'])
        bp_name = row['bioproject']

        # Color the label red if study-excluded
        label_color = '#8B0000' if row['excl_study'] > 0 else '#333333'
        label_weight = 'bold' if row['excl_study'] > 0 else 'normal'

        ax_bp.text(total_val + total_val * 0.02, i,
                  f'{bp_name}  ({keep_val}/{total_val})',
                  va='center', fontsize=7.5, fontweight=label_weight,
                  color=label_color)

    ax_bp.set_yticks(y_positions)
    ax_bp.set_yticklabels([])
    ax_bp.set_xlabel('Number of Samples', fontsize=11)
    ax_bp.set_title('Per-BioProject: Kept vs Excluded', fontsize=11, fontweight='bold')

    # Legend
    handles, labels = ax_bp.get_legend_handles_labels()
    # Reverse to match stacking order
    ax_bp.legend(reversed(handles), reversed(labels), fontsize=7, loc='lower right',
                framealpha=0.9, ncol=1)

    # Remove top/right spines
    ax_bp.spines['top'].set_visible(False)
    ax_bp.spines['right'].set_visible(False)

    # ============================
    # Panel C: Total reads distribution (bottom center)
    # ============================
    ax_reads = fig.add_subplot(gs[1, 1])

    # Log-spaced bins for reads
    max_r = max(MAX_READS * 2, plot_df['total_reads'].max())
    r_bins = np.logspace(0, np.log10(max_r + 1), 60)

    ax_reads.hist(plot_df[plot_df['keep']]['total_reads'].clip(1, max_r), bins=r_bins,
                 color='#27ae60', alpha=0.6, density=False, label='Kept')
    ax_reads.hist(plot_df[~plot_df['keep']]['total_reads'].clip(1, max_r), bins=r_bins,
                 color='#e74c3c', alpha=0.5, density=False, label='Excluded')

    # Threshold line
    ax_reads.axvline(x=MAX_READS, color='#e74c3c', linestyle='--', linewidth=2, alpha=0.8,
                    label=f'Max reads = {MAX_READS:,}')

    # Annotation for samples above threshold
    n_high_reads = (df['total_reads'] > MAX_READS).sum()
    ax_reads.annotate(f'{n_high_reads} samples\nabove {MAX_READS:,} reads',
                     xy=(MAX_READS * 1.5, ax_reads.get_ylim()[1] * 0.7 if ax_reads.get_ylim()[1] > 0 else 10),
                     fontsize=8, color='#c0392b', fontweight='bold',
                     arrowprops=dict(arrowstyle='->', color='#c0392b', alpha=0.5))

    ax_reads.set_xscale('log')
    ax_reads.set_xlabel('Total Reads (log scale)', fontsize=11)
    ax_reads.set_ylabel('Count', fontsize=11)
    ax_reads.set_title(f'Total Reads Distribution\n(threshold: {MIN_READS}–{MAX_READS:,})',
                      fontsize=10)
    ax_reads.legend(fontsize=7, loc='upper right')
    ax_reads.spines['top'].set_visible(False)
    ax_reads.spines['right'].set_visible(False)

    # ============================
    # Final layout & save
    # ============================
    out_path = OUT_DIR / "05_control_qc_v2.png"
    fig.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    log(f"  Figure saved: {out_path}")
    return out_path


# ============================================================
# 6. Export tables
# ============================================================
def export_tables(df: pd.DataFrame, study_df: pd.DataFrame):
    """Export filtered sample lists."""
    log("Exporting tables...")

    # Add top3_str to main df for export
    df['top3_str'] = df['top3_genera'].apply(
        lambda x: '; '.join([f"{g}:{c}" for g, c in x]) if x else '')

    # Full classification
    df_export = df.drop(columns=['top3_genera'])
    df_export = df_export.sort_values(['keep', 'oral_frac'], ascending=[False, False])

    classify_path = OUT_DIR / "control_classification_v2.csv"
    df_export.to_csv(classify_path, index=False, float_format='%.2f')
    log(f"  -> {classify_path} ({len(df_export)} samples)")

    # Kept list
    keep_cols = ['sample', 'sample_base', 'bioproject', 'category',
                 'total_reads', 'n_asvs', 'oral_frac', 'study_viable',
                 'n_oral_in_metadata', 'top3_str']
    keep_df = df[df['keep']][keep_cols].sort_values(['bioproject', 'n_asvs'])
    keep_path = OUT_DIR / "genuine_controls_keep_v2.csv"
    keep_df.to_csv(keep_path, index=False, float_format='%.2f')
    log(f"  -> {keep_path} ({len(keep_df)} samples)")

    # Excluded list
    excl_cols = ['sample', 'sample_base', 'bioproject', 'category',
                 'total_reads', 'n_asvs', 'oral_frac', 'exclusion_reason', 'top3_str']
    excl_df = df[~df['keep']][excl_cols].sort_values(['bioproject', 'oral_frac'], ascending=[True, False])
    excl_path = OUT_DIR / "controls_excluded_v2.csv"
    excl_df.to_csv(excl_path, index=False, float_format='%.2f')
    log(f"  -> {excl_path} ({len(excl_df)} samples)")

    # Study viability summary
    study_path = OUT_DIR / "study_decontam_viability_v2.csv"
    study_df.to_csv(study_path, index=False)
    log(f"  -> {study_path} ({len(study_df)} BioProjects)")

    # Summary JSON
    summary = {
        'total_input_samples': len(df),
        'kept_samples': int(df['keep'].sum()),
        'excluded_samples': int((~df['keep']).sum()),
        'exclusion_reasons': df.loc[~df['keep'], 'exclusion_reason'].value_counts().to_dict(),
        'viable_studies_for_decontam': int(study_df['viable_for_decontam'].sum()),
        'thresholds': {
            'min_reads': MIN_READS,
            'max_reads': MAX_READS,
            'max_asvs': MAX_ASVS,
            'max_oral_frac_pct': MAX_ORAL_FRAC,
            'min_controls_per_study': MIN_CONTROLS_PER_STUDY,
            'min_oral_per_study': MIN_ORAL_PER_STUDY,
        },
        'excluded_bioprojects': EXCLUDED_BIOPROJECTS,
    }
    with open(OUT_DIR / "control_qc_v2_summary.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    log(f"  -> control_qc_v2_summary.json")


# ============================================================
# Main
# ============================================================
def main():
    log("=" * 60)
    log("Control Sample QC v2 — Multidimensional Filtering")
    log(f"  Reads: {MIN_READS} – {MAX_READS:,}")
    log(f"  ASVs: ≤ {MAX_ASVS}")
    log(f"  Oral fraction: < {MAX_ORAL_FRAC}%")
    log(f"  Per-study: ≥ {MIN_CONTROLS_PER_STUDY} controls & ≥ {MIN_ORAL_PER_STUDY} oral")
    log(f"  Excluded BioProjects: {list(EXCLUDED_BIOPROJECTS.keys())}")
    log("=" * 60)

    # 1. Load SILVA taxonomy
    asv_taxonomy = load_silva_taxonomy(
        TEMP_DIR / "silva_hits.txt", SILVA_FASTA)

    # 2. Load sample data
    log("Loading sample data...")
    with open(TEMP_DIR / "sample_data.json") as f:
        sample_data = json.load(f)
    log(f"  Loaded {len(sample_data):,} samples")

    # 3. Classify
    df = classify_controls(sample_data, asv_taxonomy)

    # 4. Apply filters
    df = apply_filters(df)

    # 5. Check per-study viability
    df, study_df = check_study_viability(df)

    # 6. Print detailed summary
    _print_summary(df, study_df)

    # 7. Generate figure
    make_figure(df, study_df)

    # 8. Export tables
    export_tables(df, study_df)

    log("Done!")
    return df, study_df


def _print_summary(df, study_df):
    """Print a detailed filtering summary."""
    kept = df[df['keep']]
    excl = df[~df['keep']]

    print()
    print("=" * 70)
    print("FILTERING SUMMARY")
    print("=" * 70)
    print(f"  Input samples:          {len(df):,}")
    print(f"  KEPT (reliable):        {len(kept):,} ({len(kept)/len(df)*100:.1f}%)")
    print(f"  EXCLUDED:               {len(excl):,} ({len(excl)/len(df)*100:.1f}%)")
    print()

    print("  Exclusion reasons breakdown:")
    reason_counts = excl['exclusion_reason'].apply(
        lambda x: x.split(':')[0] if ':' in str(x) else str(x)).value_counts()
    for reason, count in reason_counts.items():
        print(f"    {reason}: {count:,}")
    print()

    print("  Per-category:")
    for cat in CATEGORY_ORDER:
        sub = df[df['category'] == cat]
        sub_k = sub[sub['keep']]
        if len(sub) == 0:
            continue
        print(f"    {cat:20s}: {len(sub_k):5d} kept / {len(sub):5d} total "
              f"({len(sub_k)/len(sub)*100:5.1f}%)")

    print()
    viable = study_df[study_df['viable_for_decontam']]
    print(f"  Viable studies for decontam (≥{MIN_CONTROLS_PER_STUDY} ctrl, ≥{MIN_ORAL_PER_STUDY} oral): {len(viable)}")
    for _, row in viable.iterrows():
        print(f"    {row['bioproject']:15s}: {row['n_reliable_controls']:3d} controls, {row['n_oral_in_metadata']:4d} oral")

    nonviable = study_df[~study_df['viable_for_decontam']]
    if len(nonviable) > 0:
        print(f"  Non-viable studies: {len(nonviable)}")
        for _, row in nonviable.iterrows():
            print(f"    {row['bioproject']:15s}: {row['n_reliable_controls']:3d} controls, {row['n_oral_in_metadata']:4d} oral")

    print()
    print("=" * 70)


if __name__ == "__main__":
    main()
