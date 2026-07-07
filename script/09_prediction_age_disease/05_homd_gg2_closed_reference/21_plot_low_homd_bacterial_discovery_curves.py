#!/usr/bin/env python3
"""Cohort-level discovery curves for robust low-HOMD GG2-bacterial ASVs."""

from __future__ import annotations

import argparse
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"

sns.set_theme(style="whitegrid", context="paper", font_scale=1.12)
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 150, "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False})

PALETTE = {
    "BP>=3, Subj>=50": "#1B5E20",
    "BP>=2, Subj>=20": "#2F7D73",
    "BP>=2, Subj>=10": "#4477AA",
    "BP>=2, Subj>=5": "#CC6677",
    "BP>=2, Subj>=3": "#AA3377",
}


def split_projects(x):
    if pd.isna(x):
        return []
    return [v for v in str(x).split(';') if v and v not in {'NA', 'nan'}]


def display_label(row):
    return f"BP>={int(row.min_bioprojects)}, Subj>={int(row.min_subject_prevalence)}"


def build_discovery(df, label, bin_size=20):
    project_to_asvs = defaultdict(set)
    for r in df[['derep_id', 'bioprojects']].itertuples(index=False):
        for bp in split_projects(r.bioprojects):
            project_to_asvs[bp].add(r.derep_id)
    # Deterministic and interpretable order: cohorts contributing more candidate ASVs enter first.
    order = sorted(project_to_asvs, key=lambda bp: (-len(project_to_asvs[bp]), bp))
    seen = set()
    rows = []
    bin_rows = []
    previous_bin_seen = 0
    for i, bp in enumerate(order, start=1):
        before = len(seen)
        seen.update(project_to_asvs[bp])
        new_count = len(seen) - before
        rows.append({
            'threshold': label,
            'step': i,
            'bioproject': bp,
            'cumulative_bioprojects': i,
            'cumulative_asvs': len(seen),
            'new_asvs_this_bioproject': new_count,
            'project_asv_count': len(project_to_asvs[bp]),
            'total_asvs': len(df),
            'fraction_discovered': len(seen) / len(df) if len(df) else 0,
        })
        if i % bin_size == 0 or i == len(order):
            bin_rows.append({
                'threshold': label,
                'bin_end_bioprojects': i,
                'bin_start_bioprojects': max(1, i - bin_size + 1),
                'new_asvs_in_bin': len(seen) - previous_bin_seen,
                'cumulative_asvs': len(seen),
                'total_asvs': len(df),
                'fraction_discovered': len(seen) / len(df) if len(df) else 0,
            })
            previous_bin_seen = len(seen)
    return pd.DataFrame(rows), pd.DataFrame(bin_rows)


def savefig(fig, path_base):
    fig.tight_layout()
    fig.savefig(path_base + '.pdf', bbox_inches='tight')
    fig.savefig(path_base + '.png', bbox_inches='tight', dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', default=DEFAULT_OUTDIR)
    parser.add_argument('--bin-size', type=int, default=20)
    args = parser.parse_args()

    sens_dir = os.path.join(args.outdir, 'threshold_sensitivity')
    figdir = os.path.join(args.outdir, 'figures')
    os.makedirs(figdir, exist_ok=True)
    summary = pd.read_csv(os.path.join(sens_dir, 'low_homd_threshold_sensitivity_summary.tsv'), sep='\t')
    summary['threshold_display'] = summary.apply(display_label, axis=1)
    order = summary['threshold_display'].tolist()

    curve_parts = []
    bin_parts = []
    for row in summary.itertuples(index=False):
        df = pd.read_csv(row.bacterial_table, sep='\t', usecols=['derep_id', 'bioprojects', 'subject_prevalence', 'total_reads', 'subject_geo_regions'])
        curve, bins = build_discovery(df, row.threshold_display, bin_size=args.bin_size)
        curve_parts.append(curve)
        bin_parts.append(bins)
    curves = pd.concat(curve_parts, ignore_index=True)
    bins = pd.concat(bin_parts, ignore_index=True)
    curves['threshold'] = pd.Categorical(curves['threshold'], categories=order, ordered=True)
    bins['threshold'] = pd.Categorical(bins['threshold'], categories=order, ordered=True)
    curves.to_csv(os.path.join(sens_dir, 'low_homd_bacterial_cohort_discovery_curves.tsv'), sep='\t', index=False)
    bins.to_csv(os.path.join(sens_dir, 'low_homd_bacterial_cohort_discovery_bins.tsv'), sep='\t', index=False)

    fig, axes = plt.subplots(2, 2, figsize=(14.2, 9.0))
    ax = axes[0, 0]
    sns.lineplot(data=curves, x='cumulative_bioprojects', y='cumulative_asvs', hue='threshold', linewidth=2.2, ax=ax, palette=PALETTE)
    ax.set_title('A. Cumulative discovery of low-HOMD bacterial ASVs', loc='left', fontweight='bold')
    ax.set_xlabel('Accumulated BioProjects')
    ax.set_ylabel('Cumulative unique ASVs')
    ax.legend(title='Threshold', frameon=False)

    ax = axes[0, 1]
    sns.lineplot(data=curves, x='cumulative_bioprojects', y='fraction_discovered', hue='threshold', linewidth=2.2, ax=ax, palette=PALETTE, legend=False)
    ax.axhline(0.9, color='#555555', linestyle='--', linewidth=1)
    ax.set_title('B. Fraction of final ASV set discovered', loc='left', fontweight='bold')
    ax.set_xlabel('Accumulated BioProjects')
    ax.set_ylabel('Fraction of final threshold-specific set')
    ax.set_ylim(0, 1.02)
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))

    ax = axes[1, 0]
    plot_bins = bins.copy()
    sns.lineplot(data=plot_bins, x='bin_end_bioprojects', y='new_asvs_in_bin', hue='threshold', marker='o', linewidth=2.0, ax=ax, palette=PALETTE, legend=False)
    ax.set_title(f'C. Newly discovered ASVs per {args.bin_size} BioProjects', loc='left', fontweight='bold')
    ax.set_xlabel('Accumulated BioProjects')
    ax.set_ylabel(f'New ASVs in previous {args.bin_size} BioProjects')

    ax = axes[1, 1]
    final = curves.sort_values('cumulative_bioprojects').groupby('threshold', as_index=False).tail(1).copy()
    first20 = bins[bins['bin_end_bioprojects'].le(args.bin_size)].groupby('threshold', as_index=False)['new_asvs_in_bin'].sum().rename(columns={'new_asvs_in_bin': 'first_bin_asvs'})
    lastbin = bins.sort_values('bin_end_bioprojects').groupby('threshold', as_index=False).tail(1)[['threshold', 'new_asvs_in_bin']].rename(columns={'new_asvs_in_bin': 'last_bin_asvs'})
    bar = final[['threshold', 'total_asvs']].merge(first20, on='threshold').merge(lastbin, on='threshold')
    bar_long = bar.melt(id_vars=['threshold'], value_vars=['first_bin_asvs', 'last_bin_asvs'], var_name='bin', value_name='new_asvs')
    bar_long['bin'] = bar_long['bin'].map({'first_bin_asvs': f'First {args.bin_size} BioProjects', 'last_bin_asvs': 'Last bin'})
    sns.barplot(data=bar_long, x='threshold', y='new_asvs', hue='bin', ax=ax, palette=['#555555', '#C44E52'])
    ax.set_title('D. Early vs late discovery rate', loc='left', fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('New ASVs')
    ax.tick_params(axis='x', rotation=35)
    ax.legend(title='', frameon=False)

    fig.suptitle('Cohort-level discovery curves for GG2-bacterial low-HOMD ASVs\nDeterministic order: BioProjects sorted by candidate-ASV contribution', y=1.02, fontweight='bold')
    savefig(fig, os.path.join(figdir, 'fig14_low_homd_bacterial_cohort_discovery_curves'))

    # A simplified main figure using the recommended intermediate thresholds.
    keep = ['BP>=3, Subj>=50', 'BP>=2, Subj>=20', 'BP>=2, Subj>=10']
    sub_curves = curves[curves['threshold'].isin(keep)].copy()
    sub_bins = bins[bins['threshold'].isin(keep)].copy()
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.6))
    sns.lineplot(data=sub_curves, x='cumulative_bioprojects', y='cumulative_asvs', hue='threshold', linewidth=2.4, ax=axes[0], palette=PALETTE)
    axes[0].set_title('A. Cumulative discovery', loc='left', fontweight='bold')
    axes[0].set_xlabel('Accumulated BioProjects')
    axes[0].set_ylabel('Cumulative GG2-bacterial low-HOMD ASVs')
    axes[0].legend(title='Threshold', frameon=False)
    sns.lineplot(data=sub_bins, x='bin_end_bioprojects', y='new_asvs_in_bin', hue='threshold', marker='o', linewidth=2.2, ax=axes[1], palette=PALETTE, legend=False)
    axes[1].set_title(f'B. Discovery rate per {args.bin_size} BioProjects', loc='left', fontweight='bold')
    axes[1].set_xlabel('Accumulated BioProjects')
    axes[1].set_ylabel('New ASVs')
    fig.suptitle('Low-HOMD bacterial ASVs continue to accrue across cohorts', y=1.02, fontweight='bold')
    savefig(fig, os.path.join(figdir, 'fig14b_low_homd_bacterial_cohort_discovery_selected_thresholds'))

    print(os.path.join(figdir, 'fig14_low_homd_bacterial_cohort_discovery_curves.pdf'))
    print(os.path.join(figdir, 'fig14b_low_homd_bacterial_cohort_discovery_selected_thresholds.pdf'))


if __name__ == '__main__':
    main()
