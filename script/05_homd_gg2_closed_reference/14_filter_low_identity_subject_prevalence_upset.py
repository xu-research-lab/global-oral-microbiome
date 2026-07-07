#!/usr/bin/env python3
"""Filter low-HOMD-identity ASVs by subject prevalence and draw geographic UpSet plots."""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
DEFAULT_METADATA = "/beegfs/huangxiaochang/global_oral/metadata_v2_merge/metadata-v2-merge-oral.csv"
DEFAULT_ASV_ROOT = "/beegfs/dongbiao/oral"


def log(msg):
    print(f"[subject-low97] {msg}", flush=True)


def split_values(x):
    if pd.isna(x):
        return []
    return [v for v in str(x).split(';') if v and v != 'nan' and v != 'NA']


def load_metadata(path):
    cols = ['bioproject','sample_id','subject_id','country_cleaned','continent_group']
    md = pd.read_csv(path, usecols=cols, low_memory=False)
    for col in cols:
        md[col] = md[col].fillna('NA').astype(str)
    md['subject_id_clean'] = md['subject_id']
    # If subject_id is missing, use sample_id so prevalence remains conservative and computable.
    missing = md['subject_id_clean'].isin(['', 'nan', 'NA'])
    md.loc[missing, 'subject_id_clean'] = md.loc[missing, 'sample_id']
    return md


def country_to_region(country, continent):
    country = str(country)
    continent = str(continent)
    if country in {'China','Japan','Republic of Korea','India','Thailand','Bangladesh','Israel','Qatar','Turkey'}:
        return 'Asia'
    if country in {'United States','Canada','Mexico'}:
        return 'North America'
    if country in {'Italy','Spain','Germany','United Kingdom','Sweden','Finland','Austria','Russia','Netherlands','France','Denmark','Norway','Poland','Switzerland'}:
        return 'Europe'
    if country in {'Brazil','Chile','Argentina','Colombia','Peru'}:
        return 'South America'
    if country in {'Australia','New Zealand'}:
        return 'Oceania'
    if continent in {'Asia','North America','Europe','South America','Oceania','Africa'}:
        return continent
    return 'Other/NA'


def compute_project_subject_presence(bp, asv_ids, md_bp, asv_root):
    counts_path = os.path.join(asv_root, bp, 'results', 'ASVs_counts.tsv')
    if not os.path.exists(counts_path) or not asv_ids:
        return {}
    header = pd.read_csv(counts_path, sep='\t', nrows=0).columns.tolist()
    sample_cols = header[1:]
    sample_to_subject = dict(zip(md_bp['sample_id'], md_bp['subject_id_clean']))
    sample_to_region = dict(zip(md_bp['sample_id'], md_bp['geo_region']))
    use_samples = [s for s in sample_cols if s in sample_to_subject]
    if not use_samples:
        return {}
    wanted = set(asv_ids)
    out = {}
    # Read only needed rows by chunks; ASV tables vary in size.
    for chunk in pd.read_csv(counts_path, sep='\t', index_col=0, chunksize=5000):
        sub = chunk.loc[chunk.index.intersection(wanted), use_samples]
        if sub.empty:
            continue
        for asv_id, row in sub.iterrows():
            present_samples = row.index[row.to_numpy(dtype=float) > 0]
            subjects = {sample_to_subject[s] for s in present_samples}
            regions = {sample_to_region[s] for s in present_samples if sample_to_region[s] != 'Other/NA'}
            out[asv_id] = (subjects, regions)
    return out


def plot_upset(set_counter, output_base, title, top_n=40):
    items = set_counter.most_common(top_n)
    if not items:
        return
    combos = [k for k, v in items]
    counts = [v for k, v in items]
    regions = sorted(set().union(*[set(c) for c in combos]))

    fig = plt.figure(figsize=(max(10, top_n * 0.32), 6.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.4, 1.4], hspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_mat = fig.add_subplot(gs[1], sharex=ax_bar)

    x = np.arange(len(combos))
    ax_bar.bar(x, counts, color='#34495E', width=0.82)
    ax_bar.set_ylabel('Unique ASV count')
    ax_bar.set_title(title)
    ax_bar.tick_params(axis='x', labelbottom=False)
    ax_bar.spines[['top','right']].set_visible(False)

    for i, combo in enumerate(combos):
        present = set(combo)
        ys = []
        for j, region in enumerate(regions):
            if region in present:
                ax_mat.scatter(i, j, s=34, color='#222222', zorder=3)
                ys.append(j)
            else:
                ax_mat.scatter(i, j, s=12, color='#D0D0D0', zorder=2)
        if len(ys) > 1:
            ax_mat.plot([i, i], [min(ys), max(ys)], color='#222222', linewidth=0.8, zorder=1)
    ax_mat.set_yticks(np.arange(len(regions)))
    ax_mat.set_yticklabels(regions)
    ax_mat.set_xlabel('Region intersection')
    ax_mat.set_xlim(-0.6, len(combos) - 0.4)
    ax_mat.spines[['top','right']].set_visible(False)
    ax_mat.set_xticks([])
    fig.tight_layout()
    fig.savefig(output_base + '.pdf', bbox_inches='tight')
    fig.savefig(output_base + '.png', bbox_inches='tight', dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', default=DEFAULT_OUTDIR)
    parser.add_argument('--metadata', default=DEFAULT_METADATA)
    parser.add_argument('--asv-root', default=DEFAULT_ASV_ROOT)
    parser.add_argument('--min-bioprojects', type=int, default=3)
    parser.add_argument('--min-subject-prevalence', type=int, default=50)
    parser.add_argument('--min-reads', type=float, default=100)
    args = parser.parse_args()

    low_path = os.path.join(args.outdir, 'low_identity_lt97_or_nohit_asvs.tsv')
    out_prefix = os.path.join(args.outdir, 'low_identity_lt97_subject_prevalence_filtered')
    out_table = out_prefix + '.tsv'
    stats_path = out_prefix + '.stats.json'
    figdir = os.path.join(args.outdir, 'figures')
    os.makedirs(figdir, exist_ok=True)

    log('loading low-identity ASVs')
    usecols = ['derep_id','sequence','n_asvs','n_bioprojects','total_reads','prevalence_samples','representative_asv_id','bioprojects','identity','query_cov','has_homd_hit','countries','continents','sequencing_regions','sites','statuses']
    low = pd.read_csv(low_path, sep='\t', usecols=usecols)
    for c in ['n_bioprojects','total_reads','prevalence_samples','identity','query_cov']:
        low[c] = pd.to_numeric(low[c], errors='coerce').fillna(0)
    pre = low[(low['n_bioprojects'] >= args.min_bioprojects) & (low['total_reads'] >= args.min_reads)].copy()
    log(f'pre-filtered by n_bioprojects/read count: {len(pre):,}')

    log('loading metadata')
    md = load_metadata(args.metadata)
    md['geo_region'] = [country_to_region(c, cont) for c, cont in zip(md['country_cleaned'], md['continent_group'])]
    md_by_bp = {bp: sub.copy() for bp, sub in md.groupby('bioproject')}

    # Map candidate ASVs to original project ASV IDs.
    project_to_asvs = defaultdict(dict)
    for row in pre[['derep_id','bioprojects','representative_asv_id']].itertuples(index=False):
        # Use exact-sequence membership from bioprojects, but only representative ASV ID is guaranteed here.
        # The original exact sequence can have multiple per-project ASV IDs; recover them from asv_abundance_metadata if needed.
        pass

    # More accurate: recover all source ASV IDs by joining sequence to asv_abundance_metadata in chunks.
    target_seq = dict(zip(pre['sequence'], pre['derep_id'])) if 'sequence' in pre.columns else {}
    if not target_seq:
        raise RuntimeError('sequence column is required in low identity table')
    asv_meta_path = os.path.join(args.outdir, 'asv_abundance_metadata.tsv')
    derep_to_source = defaultdict(list)
    log('recovering source ASV IDs from asv_abundance_metadata.tsv')
    for chunk in pd.read_csv(asv_meta_path, sep='\t', usecols=['bioproject','original_asv_id','sequence'], chunksize=200000):
        hit = chunk[chunk['sequence'].isin(target_seq)]
        for r in hit.itertuples(index=False):
            derep_to_source[target_seq[r.sequence]].append((r.bioproject, r.original_asv_id))

    bp_to_asvs = defaultdict(dict)
    for derep_id, sources in derep_to_source.items():
        for bp, asv_id in sources:
            bp_to_asvs[bp].setdefault(asv_id, []).append(derep_id)

    log(f'computing subject prevalence across {len(bp_to_asvs)} bioprojects')
    derep_subjects = defaultdict(set)
    derep_regions = defaultdict(set)
    for i, (bp, asv_map) in enumerate(sorted(bp_to_asvs.items()), start=1):
        md_bp = md_by_bp.get(bp)
        if md_bp is None:
            continue
        presence = compute_project_subject_presence(bp, list(asv_map.keys()), md_bp, args.asv_root)
        for asv_id, (subjects, regions) in presence.items():
            for derep_id in asv_map.get(asv_id, []):
                derep_subjects[derep_id].update(f'{bp}|{s}' for s in subjects)
                derep_regions[derep_id].update(regions)
        if i % 25 == 0:
            log(f'processed {i}/{len(bp_to_asvs)} bioprojects')

    pre['subject_prevalence'] = pre['derep_id'].map(lambda x: len(derep_subjects.get(x, set())))
    pre['subject_geo_regions'] = pre['derep_id'].map(lambda x: ';'.join(sorted(derep_regions.get(x, set()))) or 'NA')
    pre['n_subject_geo_regions'] = pre['subject_geo_regions'].map(lambda x: 0 if x == 'NA' else len(x.split(';')))
    filt = pre[pre['subject_prevalence'] >= args.min_subject_prevalence].copy()
    filt = filt.sort_values(['subject_prevalence','n_bioprojects','total_reads'], ascending=False)
    filt.to_csv(out_table, sep='\t', index=False)

    region_counter = Counter(tuple(split_values(x)) for x in filt['subject_geo_regions'])
    region_counter.pop(tuple(), None)
    pd.DataFrame([{'regions': ';'.join(k), 'unique_asv_count': v} for k, v in region_counter.most_common()]).to_csv(out_prefix + '.region_intersections.tsv', sep='\t', index=False)
    plot_upset(region_counter, os.path.join(figdir, 'fig8_low_identity_subject_prevalence_region_upset'), 'Low-HOMD-identity ASVs shared across geographic regions')

    stats = {
        'low_identity_input': int(len(low)),
        'after_n_bioprojects_and_reads_filter': int(len(pre)),
        'after_subject_prevalence_filter': int(len(filt)),
        'filters': {'n_bioprojects_min': args.min_bioprojects, 'subject_prevalence_min': args.min_subject_prevalence, 'total_reads_min': args.min_reads},
        'region_intersections_top20': [{'regions': ';'.join(k), 'unique_asv_count': v} for k, v in region_counter.most_common(20)],
        'outputs': {'filtered_table': out_table, 'region_intersections': out_prefix + '.region_intersections.tsv', 'upset_pdf': os.path.join(figdir, 'fig8_low_identity_subject_prevalence_region_upset.pdf')}
    }
    with open(stats_path, 'w') as handle:
        json.dump(stats, handle, indent=2)
    log(json.dumps(stats, indent=2))

if __name__ == '__main__':
    main()
