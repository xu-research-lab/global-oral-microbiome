#!/usr/bin/env python3
"""Threshold sensitivity for recurrent low-HOMD ASVs using subject prevalence and GG2 bacterial QC."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from collections import Counter, defaultdict

import pandas as pd

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
DEFAULT_METADATA = "/beegfs/huangxiaochang/global_oral/metadata_v2_merge/metadata-v2-merge-oral.csv"
DEFAULT_ASV_ROOT = "/beegfs/dongbiao/oral"
DEFAULT_VSEARCH = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"
DEFAULT_GG2_REF = "/beegfs/dongbiao/greengene2/2024.09.seqs.fna"
DEFAULT_GG2_TAX = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"
THRESHOLDS = [
    (3, 50, 100),
    (2, 20, 100),
    (2, 10, 100),
    (2, 5, 100),
    (2, 3, 100),
]


def log(msg):
    print(f"[threshold-sensitivity] {msg}", flush=True)


def label_threshold(nb, sp, reads):
    return f"bp{nb}_subj{sp}_reads{int(reads)}"


def split_values(x):
    if pd.isna(x):
        return []
    return [v for v in str(x).split(';') if v and v not in {'nan', 'NA'}]


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
    for chunk in pd.read_csv(counts_path, sep='\t', index_col=0, chunksize=5000):
        idx = chunk.index.intersection(wanted)
        if len(idx) == 0:
            continue
        sub = chunk.loc[idx, use_samples]
        for asv_id, row in sub.iterrows():
            vals = pd.to_numeric(row, errors='coerce').fillna(0).to_numpy(dtype=float)
            present_samples = row.index[vals > 0]
            subjects = {sample_to_subject[s] for s in present_samples}
            regions = {sample_to_region[s] for s in present_samples if sample_to_region[s] != 'Other/NA'}
            out[asv_id] = (subjects, regions)
    return out


def load_taxonomy(path):
    tax = {}
    with open(path) as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                tax[parts[0]] = parts[1]
    return tax


def parse_rank(taxonomy, prefix):
    if pd.isna(taxonomy):
        return 'Unassigned'
    for part in str(taxonomy).split(';'):
        part = part.strip()
        if part.startswith(prefix):
            val = part.split('__', 1)[1].strip()
            return val if val else 'Unassigned'
    return 'Unassigned'


def classify_qc(taxonomy):
    t = str(taxonomy)
    levels = [x.strip() for x in t.split(';')]
    if not t or t == 'nan' or t == 'Unassigned':
        return 'unassigned'
    if any('Chloroplast' in x for x in levels):
        return 'chloroplast'
    if any('Mitochondria' in x for x in levels):
        return 'mitochondria'
    if levels and levels[0].startswith('d__Bacteria'):
        return 'bacterial'
    if levels and levels[0].startswith('d__Archaea'):
        return 'archaea'
    if levels and levels[0].startswith('d__Eukaryota'):
        return 'eukaryota'
    return 'other_non_bacterial'


def write_fasta(df, path):
    with open(path, 'w') as handle:
        for row in df.itertuples(index=False):
            handle.write(f">{row.derep_id}\n")
            seq = row.sequence
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + '\n')


def run_gg2(df, outdir, args):
    fasta_path = os.path.join(outdir, 'threshold_sensitivity_max_relaxed.fa')
    hit_path = os.path.join(outdir, 'threshold_sensitivity_max_relaxed_vs_gg2_best_hits.tsv')
    ann_path = os.path.join(outdir, 'threshold_sensitivity_max_relaxed_gg2_annotated.tsv')
    write_fasta(df, fasta_path)
    if args.force_gg2 or not (os.path.exists(hit_path) and os.path.getsize(hit_path) > 0):
        cmd = [args.vsearch, '--usearch_global', fasta_path, '--db', args.gg2_ref, '--id', '0.7', '--maxaccepts', '1', '--maxrejects', '0', '--strand', 'both', '--userout', hit_path, '--userfields', 'query+target+id+alnlen+qcov+tcov+ql+tl', '--threads', str(args.threads), '--quiet']
        log('running GG2 alignment for max relaxed threshold set')
        log(' '.join(cmd))
        subprocess.run(cmd, check=True)
    else:
        log(f'using existing GG2 hits: {hit_path}')

    hit_cols = ['derep_id', 'gg2_target', 'gg2_identity', 'gg2_alnlen', 'gg2_qcov', 'gg2_tcov', 'query_length', 'target_length']
    hits = pd.read_csv(hit_path, sep='\t', names=hit_cols) if os.path.getsize(hit_path) else pd.DataFrame(columns=hit_cols)
    tax = load_taxonomy(args.gg2_tax)
    if len(hits):
        hits['gg2_taxonomy'] = hits['gg2_target'].map(tax).fillna('Unassigned')
        hits['has_taxonomy'] = hits['gg2_taxonomy'].ne('Unassigned')
        sort_cols = ['derep_id', 'gg2_identity', 'gg2_qcov', 'gg2_alnlen']
        best_overall = hits.sort_values(sort_cols, ascending=[True, False, False, False]).drop_duplicates('derep_id').rename(columns={
            'gg2_target': 'gg2_best_target', 'gg2_identity': 'gg2_best_identity', 'gg2_alnlen': 'gg2_best_alnlen', 'gg2_qcov': 'gg2_best_qcov', 'gg2_tcov': 'gg2_best_tcov', 'query_length': 'gg2_best_query_length', 'target_length': 'gg2_best_target_length', 'gg2_taxonomy': 'gg2_best_taxonomy', 'has_taxonomy': 'gg2_best_has_taxonomy'})
        best_tax = hits[hits['has_taxonomy']].sort_values(sort_cols, ascending=[True, False, False, False]).drop_duplicates('derep_id')
        hits = best_overall.merge(best_tax[['derep_id', 'gg2_target', 'gg2_identity', 'gg2_alnlen', 'gg2_qcov', 'gg2_tcov', 'query_length', 'target_length', 'gg2_taxonomy']], on='derep_id', how='left')
    else:
        hits['gg2_taxonomy'] = []
    ann = df.merge(hits, on='derep_id', how='left')
    ann['gg2_taxonomy'] = ann['gg2_taxonomy'].fillna('Unassigned')
    ann['gg2_qc_class'] = ann['gg2_taxonomy'].map(classify_qc)
    for col, prefix in [('gg2_domain', 'd__'), ('gg2_phylum', 'p__'), ('gg2_class', 'c__'), ('gg2_order', 'o__'), ('gg2_family', 'f__'), ('gg2_genus', 'g__'), ('gg2_species', 's__')]:
        ann[col] = ann['gg2_taxonomy'].map(lambda x, p=prefix: parse_rank(x, p))
    ann.to_csv(ann_path, sep='\t', index=False)
    return ann, hit_path, ann_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', default=DEFAULT_OUTDIR)
    parser.add_argument('--metadata', default=DEFAULT_METADATA)
    parser.add_argument('--asv-root', default=DEFAULT_ASV_ROOT)
    parser.add_argument('--vsearch', default=DEFAULT_VSEARCH)
    parser.add_argument('--gg2-ref', default=DEFAULT_GG2_REF)
    parser.add_argument('--gg2-tax', default=DEFAULT_GG2_TAX)
    parser.add_argument('--threads', type=int, default=24)
    parser.add_argument('--force-gg2', action='store_true')
    args = parser.parse_args()

    sens_dir = os.path.join(args.outdir, 'threshold_sensitivity')
    os.makedirs(sens_dir, exist_ok=True)
    low_path = os.path.join(args.outdir, 'low_identity_lt97_or_nohit_asvs.tsv')
    max_cache = os.path.join(sens_dir, 'low_identity_lt97_bp2_reads100_subject_prevalence.tsv')

    if os.path.exists(max_cache) and os.path.getsize(max_cache) > 0:
        log(f'using existing subject-prevalence max cache: {max_cache}')
        max_df = pd.read_csv(max_cache, sep='\t')
    else:
        log('loading max relaxed low-HOMD set: n_bioprojects>=2, total_reads>=100')
        usecols = ['derep_id', 'sequence', 'n_asvs', 'n_bioprojects', 'total_reads', 'prevalence_samples', 'representative_asv_id', 'bioprojects', 'identity', 'query_cov', 'has_homd_hit', 'countries', 'continents', 'sequencing_regions', 'sites', 'statuses']
        parts = []
        for chunk in pd.read_csv(low_path, sep='\t', usecols=usecols, chunksize=200000):
            chunk['n_bioprojects'] = pd.to_numeric(chunk['n_bioprojects'], errors='coerce').fillna(0)
            chunk['total_reads'] = pd.to_numeric(chunk['total_reads'], errors='coerce').fillna(0)
            hit = chunk[(chunk['n_bioprojects'] >= 2) & (chunk['total_reads'] >= 100)].copy()
            if len(hit):
                parts.append(hit)
        max_df = pd.concat(parts, ignore_index=True)
        log(f'max relaxed set rows: {len(max_df):,}')

        log('loading metadata')
        md = load_metadata(args.metadata)
        md_by_bp = {bp: sub.copy() for bp, sub in md.groupby('bioproject')}

        target_seq = dict(zip(max_df['sequence'], max_df['derep_id']))
        asv_meta_path = os.path.join(args.outdir, 'asv_abundance_metadata.tsv')
        derep_to_source = defaultdict(list)
        log('recovering source ASV IDs from asv_abundance_metadata.tsv')
        for i, chunk in enumerate(pd.read_csv(asv_meta_path, sep='\t', usecols=['bioproject', 'original_asv_id', 'sequence'], chunksize=200000), start=1):
            hit = chunk[chunk['sequence'].isin(target_seq)]
            for r in hit.itertuples(index=False):
                derep_to_source[target_seq[r.sequence]].append((r.bioproject, r.original_asv_id))
            if i % 10 == 0:
                log(f'processed {i * 200000:,} ASV metadata rows')

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

        max_df['subject_prevalence'] = max_df['derep_id'].map(lambda x: len(derep_subjects.get(x, set())))
        max_df['subject_geo_regions'] = max_df['derep_id'].map(lambda x: ';'.join(sorted(derep_regions.get(x, set()))) or 'NA')
        max_df['n_subject_geo_regions'] = max_df['subject_geo_regions'].map(lambda x: 0 if x == 'NA' else len(x.split(';')))
        max_df = max_df.sort_values(['subject_prevalence', 'n_bioprojects', 'total_reads'], ascending=False)
        max_df.to_csv(max_cache, sep='\t', index=False)

    log('running or loading GG2 annotation for max relaxed set')
    ann, hit_path, ann_path = run_gg2(max_df, sens_dir, args)

    rows = []
    for nb, sp, reads in THRESHOLDS:
        label = label_threshold(nb, sp, reads)
        sub = ann[(ann['n_bioprojects'] >= nb) & (ann['total_reads'] >= reads) & (ann['subject_prevalence'] >= sp)].copy()
        out_prefix = os.path.join(sens_dir, f'low_identity_lt97_{label}')
        sub.to_csv(out_prefix + '.gg2_annotated.tsv', sep='\t', index=False)
        bacterial = sub[sub['gg2_qc_class'].eq('bacterial')].copy()
        bacterial.to_csv(out_prefix + '.gg2_bacterial_only.tsv', sep='\t', index=False)
        qc = sub.groupby('gg2_qc_class', dropna=False).agg(
            derep_sequences=('derep_id', 'count'),
            source_asvs=('n_asvs', 'sum'),
            total_reads=('total_reads', 'sum'),
            subject_prevalence=('subject_prevalence', 'sum'),
            max_subject_prevalence=('subject_prevalence', 'max'),
            max_bioprojects=('n_bioprojects', 'max'),
        ).reset_index().sort_values('subject_prevalence', ascending=False)
        qc.to_csv(out_prefix + '.gg2_qc_summary.tsv', sep='\t', index=False)
        region_counter = Counter(tuple(split_values(x)) for x in bacterial['subject_geo_regions'])
        region_counter.pop(tuple(), None)
        pd.DataFrame([{'regions': ';'.join(k), 'bacterial_asv_count': v} for k, v in region_counter.most_common()]).to_csv(out_prefix + '.bacterial_region_intersections.tsv', sep='\t', index=False)
        counts = qc.set_index('gg2_qc_class')['derep_sequences'].to_dict() if len(qc) else {}
        rows.append({
            'threshold_label': label,
            'min_bioprojects': nb,
            'min_subject_prevalence': sp,
            'min_total_reads': reads,
            'filtered_low_homd_asvs': int(len(sub)),
            'gg2_bacterial_asvs': int(len(bacterial)),
            'gg2_unassigned_asvs': int(counts.get('unassigned', 0)),
            'gg2_chloroplast_asvs': int(counts.get('chloroplast', 0)),
            'gg2_mitochondria_asvs': int(counts.get('mitochondria', 0)),
            'gg2_archaea_asvs': int(counts.get('archaea', 0)),
            'bacterial_source_asvs': int(bacterial['n_asvs'].sum()) if len(bacterial) else 0,
            'bacterial_total_reads': float(bacterial['total_reads'].sum()) if len(bacterial) else 0,
            'bacterial_subject_prevalence_burden': int(bacterial['subject_prevalence'].sum()) if len(bacterial) else 0,
            'bacterial_max_subject_prevalence': int(bacterial['subject_prevalence'].max()) if len(bacterial) else 0,
            'bacterial_max_bioprojects': int(bacterial['n_bioprojects'].max()) if len(bacterial) else 0,
            'top_region_intersection': ';'.join(region_counter.most_common(1)[0][0]) if region_counter else 'NA',
            'top_region_intersection_count': int(region_counter.most_common(1)[0][1]) if region_counter else 0,
            'annotated_table': out_prefix + '.gg2_annotated.tsv',
            'bacterial_table': out_prefix + '.gg2_bacterial_only.tsv',
            'qc_summary': out_prefix + '.gg2_qc_summary.tsv',
            'region_intersections': out_prefix + '.bacterial_region_intersections.tsv',
        })
        log(f'{label}: total={len(sub):,}, bacterial={len(bacterial):,}')

    summary = pd.DataFrame(rows)
    summary_path = os.path.join(sens_dir, 'low_homd_threshold_sensitivity_summary.tsv')
    summary.to_csv(summary_path, sep='\t', index=False)
    stats = {'thresholds': rows, 'max_cache': max_cache, 'gg2_hit_path': hit_path, 'gg2_annotated_max_table': ann_path, 'summary': summary_path}
    stats_path = os.path.join(sens_dir, 'low_homd_threshold_sensitivity_stats.json')
    with open(stats_path, 'w') as handle:
        json.dump(stats, handle, indent=2)
    log(json.dumps(stats, indent=2))


if __name__ == '__main__':
    main()
