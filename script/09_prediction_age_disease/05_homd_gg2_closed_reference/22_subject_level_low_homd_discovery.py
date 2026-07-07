#!/usr/bin/env python3
"""Subject-level discovery curves for GG2-bacterial low-HOMD ASVs across threshold sets."""

from __future__ import annotations

import argparse
import json
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
DEFAULT_METADATA = "/beegfs/huangxiaochang/global_oral/metadata_v2_merge/metadata-v2-merge-oral.csv"
DEFAULT_ASV_ROOT = "/beegfs/dongbiao/oral"

sns.set_theme(style="whitegrid", context="paper", font_scale=1.12)
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 150, "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False})

PALETTE = {
    "BP>=3, Subj>=50": "#1B5E20",
    "BP>=2, Subj>=20": "#2F7D73",
    "BP>=2, Subj>=10": "#4477AA",
    "BP>=2, Subj>=5": "#CC6677",
    "BP>=2, Subj>=3": "#AA3377",
}


def log(msg):
    print(f"[subject-discovery] {msg}", flush=True)


def country_to_region(country, continent):
    country = str(country)
    continent = str(continent)
    if country in {'China', 'Japan', 'Republic of Korea', 'India', 'Thailand', 'Bangladesh', 'Israel', 'Qatar', 'Turkey'}:
        return 'Asia'
    if country in {'United States', 'Canada', 'Mexico'}:
        return 'North America'
    if country in {'Italy', 'Spain', 'Germany', 'United Kingdom', 'Sweden', 'Finland', 'Austria', 'Russia', 'Netherlands', 'France', 'Denmark', 'Norway', 'Poland', 'Switzerland'}:
        return 'Europe'
    if country in {'Brazil', 'Chile', 'Argentina', 'Colombia', 'Peru'}:
        return 'South America'
    if country in {'Australia', 'New Zealand'}:
        return 'Oceania'
    if continent in {'Asia', 'North America', 'Europe', 'South America', 'Oceania', 'Africa'}:
        return continent
    return 'Other/NA'


def load_metadata(path):
    cols = ['bioproject', 'sample_id', 'subject_id', 'country_cleaned', 'continent_group']
    md = pd.read_csv(path, usecols=cols, low_memory=False)
    for col in cols:
        md[col] = md[col].fillna('NA').astype(str)
    md['subject_id_clean'] = md['subject_id']
    missing = md['subject_id_clean'].isin(['', 'nan', 'NA'])
    md.loc[missing, 'subject_id_clean'] = md.loc[missing, 'sample_id']
    md['geo_region'] = [country_to_region(c, cont) for c, cont in zip(md['country_cleaned'], md['continent_group'])]
    return md


def display_label(row):
    return f"BP>={int(row.min_bioprojects)}, Subj>={int(row.min_subject_prevalence)}"


def build_presence_table(args, sens_dir, summary):
    out_path = os.path.join(sens_dir, 'low_homd_bacterial_bp2_subj3_reads100_subject_presence_long.tsv')
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0 and not args.force_presence:
        log(f'using existing subject presence table: {out_path}')
        return out_path

    relaxed = summary[summary['threshold_label'].eq('bp2_subj3_reads100')].iloc[0]
    target = pd.read_csv(relaxed['bacterial_table'], sep='\t', usecols=['derep_id', 'sequence'])
    target_seq = dict(zip(target['sequence'], target['derep_id']))
    target_dereps = set(target['derep_id'])
    log(f'building subject presence for max relaxed bacterial set: {len(target_dereps):,} ASVs')

    log('loading metadata')
    md = load_metadata(args.metadata)
    md_by_bp = {bp: sub.copy() for bp, sub in md.groupby('bioproject')}

    log('recovering original project ASV IDs from asv_abundance_metadata.tsv')
    asv_meta_path = os.path.join(args.outdir, 'asv_abundance_metadata.tsv')
    bp_to_asvs = defaultdict(dict)
    for i, chunk in enumerate(pd.read_csv(asv_meta_path, sep='\t', usecols=['bioproject', 'original_asv_id', 'sequence'], chunksize=200000), start=1):
        hit = chunk[chunk['sequence'].isin(target_seq)]
        for r in hit.itertuples(index=False):
            derep_id = target_seq[r.sequence]
            bp_to_asvs[r.bioproject].setdefault(r.original_asv_id, []).append(derep_id)
        if i % 10 == 0:
            log(f'processed {i * 200000:,} ASV metadata rows')

    log(f'scanning ASV count tables across {len(bp_to_asvs)} BioProjects')
    first = True
    total_edges = 0
    for i, (bp, asv_map) in enumerate(sorted(bp_to_asvs.items()), start=1):
        md_bp = md_by_bp.get(bp)
        counts_path = os.path.join(args.asv_root, bp, 'results', 'ASVs_counts.tsv')
        if md_bp is None or not os.path.exists(counts_path):
            continue
        header = pd.read_csv(counts_path, sep='\t', nrows=0).columns.tolist()
        sample_cols = header[1:]
        sample_to_subject = dict(zip(md_bp['sample_id'], md_bp['subject_id_clean']))
        sample_to_region = dict(zip(md_bp['sample_id'], md_bp['geo_region']))
        use_samples = [s for s in sample_cols if s in sample_to_subject]
        if not use_samples:
            continue
        wanted = set(asv_map)
        rows = []
        for chunk in pd.read_csv(counts_path, sep='\t', index_col=0, chunksize=5000):
            idx = chunk.index.intersection(wanted)
            if len(idx) == 0:
                continue
            sub = chunk.loc[idx, use_samples]
            for asv_id, row in sub.iterrows():
                vals = pd.to_numeric(row, errors='coerce').fillna(0).to_numpy(dtype=float)
                present_samples = row.index[vals > 0]
                if len(present_samples) == 0:
                    continue
                # Collapse sample-level detections to subject-level detections within the BioProject.
                subj_regions = {}
                for sample in present_samples:
                    subj = sample_to_subject[sample]
                    region = sample_to_region[sample]
                    subj_regions[subj] = region
                for derep_id in asv_map.get(asv_id, []):
                    for subj, region in subj_regions.items():
                        rows.append({
                            'subject_key': f'{bp}|{subj}',
                            'bioproject': bp,
                            'subject_id_clean': subj,
                            'geo_region': region,
                            'derep_id': derep_id,
                        })
        if rows:
            out = pd.DataFrame(rows).drop_duplicates(['subject_key', 'derep_id'])
            out.to_csv(out_path, sep='\t', index=False, mode='w' if first else 'a', header=first)
            first = False
            total_edges += len(out)
        if i % 25 == 0:
            log(f'processed {i}/{len(bp_to_asvs)} BioProjects; subject-ASV edges={total_edges:,}')
    log(f'wrote subject presence edges: {total_edges:,}')
    return out_path


def summarize_curves(curves):
    return curves.groupby(['threshold', 'accumulated_subjects'], as_index=False)['cumulative_asvs'].agg(
        median='median',
        q025=lambda x: np.quantile(x, 0.025),
        q975=lambda x: np.quantile(x, 0.975),
    )


def run_subject_permutations(presence, selected_dereps, threshold, n_perm=20, seed=17, step_size=500):
    sub = presence[presence['derep_id'].isin(selected_dereps)].copy()
    subject_groups = sub.groupby('subject_key')['derep_id'].apply(lambda x: set(x)).to_dict()
    subjects = np.array(sorted(subject_groups), dtype=object)
    if len(subjects) == 0:
        return pd.DataFrame(), pd.DataFrame()
    checkpoints = set(range(step_size, len(subjects) + 1, step_size))
    checkpoints.add(len(subjects))
    rng = np.random.default_rng(seed)
    curve_rows = []
    bin_rows = []
    for perm in range(n_perm):
        order = subjects.copy()
        rng.shuffle(order)
        seen = set()
        prev_bin = 0
        for i, subj in enumerate(order, start=1):
            seen.update(subject_groups[subj])
            if i in checkpoints:
                cur = len(seen)
                curve_rows.append({'threshold': threshold, 'perm': perm, 'accumulated_subjects': i, 'cumulative_asvs': cur, 'total_subjects': len(subjects), 'total_asvs': len(selected_dereps), 'fraction_discovered': cur / len(selected_dereps) if selected_dereps else 0})
                bin_rows.append({'threshold': threshold, 'perm': perm, 'bin_end_subjects': i, 'bin_start_subjects': max(1, i - step_size + 1), 'new_asvs_in_bin': cur - prev_bin, 'cumulative_asvs': cur, 'total_subjects': len(subjects), 'total_asvs': len(selected_dereps)})
                prev_bin = cur
    return pd.DataFrame(curve_rows), pd.DataFrame(bin_rows)


def savefig(fig, path_base):
    fig.tight_layout()
    fig.savefig(path_base + '.pdf', bbox_inches='tight')
    fig.savefig(path_base + '.png', bbox_inches='tight', dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', default=DEFAULT_OUTDIR)
    parser.add_argument('--metadata', default=DEFAULT_METADATA)
    parser.add_argument('--asv-root', default=DEFAULT_ASV_ROOT)
    parser.add_argument('--n-perm', type=int, default=20)
    parser.add_argument('--step-size', type=int, default=500)
    parser.add_argument('--seed', type=int, default=17)
    parser.add_argument('--force-presence', action='store_true')
    args = parser.parse_args()

    sens_dir = os.path.join(args.outdir, 'threshold_sensitivity')
    figdir = os.path.join(args.outdir, 'figures')
    os.makedirs(figdir, exist_ok=True)
    summary = pd.read_csv(os.path.join(sens_dir, 'low_homd_threshold_sensitivity_summary.tsv'), sep='\t')
    summary['threshold_display'] = summary.apply(display_label, axis=1)
    order = summary['threshold_display'].tolist()

    presence_path = build_presence_table(args, sens_dir, summary)
    log('loading subject presence table')
    presence = pd.read_csv(presence_path, sep='\t')
    log(f'presence edges={len(presence):,}; subjects={presence.subject_key.nunique():,}; ASVs={presence.derep_id.nunique():,}')

    curve_parts = []
    bin_parts = []
    stats = []
    for row in summary.itertuples(index=False):
        btab = pd.read_csv(row.bacterial_table, sep='\t', usecols=['derep_id'])
        selected = set(btab['derep_id'])
        log(f'running subject permutations for {row.threshold_display}: ASVs={len(selected):,}')
        curves, bins = run_subject_permutations(presence, selected, row.threshold_display, n_perm=args.n_perm, seed=args.seed, step_size=args.step_size)
        curve_parts.append(curves)
        bin_parts.append(bins)
        stats.append({'threshold': row.threshold_display, 'asvs': len(selected), 'subjects_with_any_asv': int(presence[presence['derep_id'].isin(selected)].subject_key.nunique())})
    curves = pd.concat(curve_parts, ignore_index=True)
    bins = pd.concat(bin_parts, ignore_index=True)
    curves['threshold'] = pd.Categorical(curves['threshold'], categories=order, ordered=True)
    bins['threshold'] = pd.Categorical(bins['threshold'], categories=order, ordered=True)
    curve_path = os.path.join(sens_dir, 'low_homd_bacterial_subject_discovery_curves.tsv')
    bin_path = os.path.join(sens_dir, 'low_homd_bacterial_subject_discovery_bins.tsv')
    summary_path = os.path.join(sens_dir, 'low_homd_bacterial_subject_discovery_summary.tsv')
    curves.to_csv(curve_path, sep='\t', index=False)
    bins.to_csv(bin_path, sep='\t', index=False)
    pd.DataFrame(stats).to_csv(summary_path, sep='\t', index=False)

    curve_summary = summarize_curves(curves)
    curve_summary['threshold'] = pd.Categorical(curve_summary['threshold'], categories=order, ordered=True)
    final_bins = bins.sort_values('bin_end_subjects').groupby(['threshold', 'perm'], as_index=False).tail(1)
    first_bins = bins.groupby(['threshold', 'perm'], as_index=False).head(1)

    fig, axes = plt.subplots(2, 2, figsize=(14.2, 9.0))
    ax = axes[0, 0]
    for threshold in order:
        sub = curve_summary[curve_summary['threshold'].eq(threshold)].sort_values('accumulated_subjects')
        if sub.empty:
            continue
        x = sub['accumulated_subjects'].to_numpy(float)
        y = sub['median'].to_numpy(float)
        lo = sub['q025'].to_numpy(float)
        hi = sub['q975'].to_numpy(float)
        ax.plot(x, y, linewidth=2.2, color=PALETTE.get(threshold), label=threshold)
        ax.fill_between(x, lo, hi, color=PALETTE.get(threshold), alpha=0.15, linewidth=0)
    ax.set_title('A. Subject-level cumulative discovery', loc='left', fontweight='bold')
    ax.set_xlabel('Accumulated subjects')
    ax.set_ylabel('Cumulative unique ASVs')
    ax.legend(title='Threshold', frameon=False)

    ax = axes[0, 1]
    curves['fraction_discovered'] = curves['cumulative_asvs'] / curves['total_asvs']
    frac_summary = curves.groupby(['threshold', 'accumulated_subjects'], as_index=False)['fraction_discovered'].agg(median='median', q025=lambda x: np.quantile(x, 0.025), q975=lambda x: np.quantile(x, 0.975))
    frac_summary['threshold'] = pd.Categorical(frac_summary['threshold'], categories=order, ordered=True)
    for threshold in order:
        sub = frac_summary[frac_summary['threshold'].eq(threshold)].sort_values('accumulated_subjects')
        if sub.empty:
            continue
        x = sub['accumulated_subjects'].to_numpy(float)
        y = sub['median'].to_numpy(float)
        lo = sub['q025'].to_numpy(float)
        hi = sub['q975'].to_numpy(float)
        ax.plot(x, y, linewidth=2.2, color=PALETTE.get(threshold), label=threshold)
        ax.fill_between(x, lo, hi, color=PALETTE.get(threshold), alpha=0.15, linewidth=0)
    ax.axhline(0.9, color='#555555', linestyle='--', linewidth=1)
    ax.set_title('B. Fraction of final ASV set discovered', loc='left', fontweight='bold')
    ax.set_xlabel('Accumulated subjects')
    ax.set_ylabel('Fraction discovered')
    ax.set_ylim(0, 1.02)
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))

    ax = axes[1, 0]
    bin_summary = bins.groupby(['threshold', 'bin_end_subjects'], as_index=False)['new_asvs_in_bin'].agg(median='median', q025=lambda x: np.quantile(x, 0.025), q975=lambda x: np.quantile(x, 0.975))
    bin_summary['threshold'] = pd.Categorical(bin_summary['threshold'], categories=order, ordered=True)
    for threshold in order:
        sub = bin_summary[bin_summary['threshold'].eq(threshold)].sort_values('bin_end_subjects')
        if sub.empty:
            continue
        ax.plot(sub['bin_end_subjects'].to_numpy(float), sub['median'].to_numpy(float), linewidth=2.0, marker='o', markersize=3, color=PALETTE.get(threshold), label=threshold)
    ax.set_title(f'C. Newly discovered ASVs per {args.step_size} subjects', loc='left', fontweight='bold')
    ax.set_xlabel('Accumulated subjects')
    ax.set_ylabel(f'New ASVs in previous {args.step_size} subjects')

    ax = axes[1, 1]
    first = first_bins.groupby('threshold', as_index=False)['new_asvs_in_bin'].median().rename(columns={'new_asvs_in_bin': 'first_bin_new_asvs'})
    last = final_bins.groupby('threshold', as_index=False)['new_asvs_in_bin'].median().rename(columns={'new_asvs_in_bin': 'last_bin_new_asvs'})
    comp = first.merge(last, on='threshold')
    comp['threshold'] = pd.Categorical(comp['threshold'], categories=order, ordered=True)
    comp_long = comp.melt(id_vars='threshold', value_vars=['first_bin_new_asvs', 'last_bin_new_asvs'], var_name='bin', value_name='new_asvs')
    comp_long['bin'] = comp_long['bin'].map({'first_bin_new_asvs': f'First {args.step_size} subjects', 'last_bin_new_asvs': 'Last bin'})
    sns.barplot(data=comp_long, x='threshold', y='new_asvs', hue='bin', ax=ax, palette=['#555555', '#C44E52'])
    ax.set_title('D. Early vs late subject-level discovery', loc='left', fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Median new ASVs')
    ax.tick_params(axis='x', rotation=35)
    ax.legend(title='', frameon=False)

    fig.suptitle(f'Subject-level discovery curves for GG2-bacterial low-HOMD ASVs ({args.n_perm} random subject orders)', y=1.02, fontweight='bold')
    savefig(fig, os.path.join(figdir, 'fig15_low_homd_bacterial_subject_discovery_curves'))

    keep = ['BP>=3, Subj>=50', 'BP>=2, Subj>=20', 'BP>=2, Subj>=10']
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.6))
    for threshold in keep:
        sub = curve_summary[curve_summary['threshold'].eq(threshold)].sort_values('accumulated_subjects')
        if sub.empty:
            continue
        axes[0].plot(sub['accumulated_subjects'].to_numpy(float), sub['median'].to_numpy(float), linewidth=2.4, color=PALETTE.get(threshold), label=threshold)
        axes[0].fill_between(sub['accumulated_subjects'].to_numpy(float), sub['q025'].to_numpy(float), sub['q975'].to_numpy(float), color=PALETTE.get(threshold), alpha=0.15, linewidth=0)
        bsub = bin_summary[bin_summary['threshold'].eq(threshold)].sort_values('bin_end_subjects')
        axes[1].plot(bsub['bin_end_subjects'].to_numpy(float), bsub['median'].to_numpy(float), linewidth=2.2, marker='o', markersize=3, color=PALETTE.get(threshold), label=threshold)
    axes[0].set_title('A. Cumulative discovery', loc='left', fontweight='bold')
    axes[0].set_xlabel('Accumulated subjects')
    axes[0].set_ylabel('Cumulative unique ASVs')
    axes[0].legend(title='Threshold', frameon=False)
    axes[1].set_title(f'B. Discovery rate per {args.step_size} subjects', loc='left', fontweight='bold')
    axes[1].set_xlabel('Accumulated subjects')
    axes[1].set_ylabel('New ASVs')
    fig.suptitle('Subject-level low-HOMD bacterial ASV discovery', y=1.02, fontweight='bold')
    savefig(fig, os.path.join(figdir, 'fig15b_low_homd_bacterial_subject_discovery_selected_thresholds'))

    stats_path = os.path.join(sens_dir, 'low_homd_bacterial_subject_discovery_stats.json')
    stats = {'presence_table': presence_path, 'curve_table': curve_path, 'bin_table': bin_path, 'summary_table': summary_path, 'n_permutations': args.n_perm, 'step_size_subjects': args.step_size, 'thresholds': stats}
    with open(stats_path, 'w') as handle:
        json.dump(stats, handle, indent=2)
    log(json.dumps(stats, indent=2))


if __name__ == '__main__':
    main()
