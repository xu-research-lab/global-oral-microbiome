#!/usr/bin/env python3
"""Analyze observed ASVs with low HOMD identity (<97% or no hit), prioritizing recurrent/prevalent signals."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
DEFAULT_VSEARCH = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"
DEFAULT_GG2_REF = "/beegfs/dongbiao/greengene2/2024.09.seqs.fna"
DEFAULT_GG2_TAX = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"


def log(msg):
    print(f"[low97] {msg}", flush=True)


def split_values(x):
    if pd.isna(x):
        return []
    return [v for v in str(x).split(';') if v and v != 'nan']


def write_fasta(df, path):
    with open(path, 'w') as handle:
        for row in df.itertuples(index=False):
            handle.write(f">{row.derep_id}\n")
            seq = row.sequence
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + '\n')


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


def add_geo_counts(df):
    df = df.copy()
    df['n_countries'] = df['countries'].map(lambda x: len(set(split_values(x))))
    df['n_continents'] = df['continents'].map(lambda x: len(set(split_values(x))))
    df['n_regions'] = df['sequencing_regions'].map(lambda x: len(set(split_values(x))))
    return df


def expand_group_tax(df, group_col, tax_col):
    rows = []
    for row in df.itertuples(index=False):
        for group in split_values(getattr(row, group_col)):
            rows.append({
                'grouping': group_col,
                'level': group,
                'taxon_rank': tax_col.replace('gg2_', ''),
                'taxon': getattr(row, tax_col),
                'derep_sequences': 1,
                'source_asvs': row.n_asvs,
                'total_reads': row.total_reads,
                'prevalence_samples': row.prevalence_samples,
                'n_bioprojects_sum': row.n_bioprojects,
            })
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).groupby(['grouping','level','taxon_rank','taxon'], dropna=False).agg(
        derep_sequences=('derep_sequences','sum'),
        source_asvs=('source_asvs','sum'),
        total_reads=('total_reads','sum'),
        prevalence_samples=('prevalence_samples','sum'),
        n_bioprojects_sum=('n_bioprojects_sum','sum'),
    ).reset_index().sort_values(['grouping','level','total_reads'], ascending=[True,True,False])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', default=DEFAULT_OUTDIR)
    parser.add_argument('--vsearch', default=DEFAULT_VSEARCH)
    parser.add_argument('--gg2-ref', default=DEFAULT_GG2_REF)
    parser.add_argument('--gg2-tax', default=DEFAULT_GG2_TAX)
    parser.add_argument('--threads', type=int, default=24)
    parser.add_argument('--top-n-annotate', type=int, default=50000)
    parser.add_argument('--min-bioprojects-priority', type=int, default=3)
    parser.add_argument('--min-prevalence-priority', type=int, default=10)
    parser.add_argument('--min-reads-priority', type=float, default=1000)
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()

    cov_path = os.path.join(args.outdir, 'qcov_sensitivity_by_derep_sequence.tsv')
    low_path = os.path.join(args.outdir, 'low_identity_lt97_or_nohit_asvs.tsv')
    priority_path = os.path.join(args.outdir, 'low_identity_lt97_priority_asvs.tsv')
    fasta_path = os.path.join(args.outdir, 'low_identity_lt97_priority_asvs.fa')
    gg2_hit_path = os.path.join(args.outdir, 'low_identity_lt97_priority_vs_gg2_best_hits.tsv')
    annot_path = os.path.join(args.outdir, 'low_identity_lt97_priority_gg2_annotated.tsv')
    qc_path = os.path.join(args.outdir, 'low_identity_lt97_priority_gg2_qc_summary.tsv')
    tax_path = os.path.join(args.outdir, 'low_identity_lt97_priority_gg2_taxonomy_summary.tsv')
    group_tax_path = os.path.join(args.outdir, 'low_identity_lt97_priority_group_taxonomy_summary.tsv')
    stats_path = os.path.join(args.outdir, 'low_identity_lt97_analysis_stats.json')

    usecols = ['derep_id','sequence','length','n_asvs','n_bioprojects','total_reads','prevalence_samples','representative_asv_id','bioprojects','identity','query_cov','has_homd_hit','sequencing_regions','platforms','continents','countries','sites','disease_studies','statuses']
    log('loading coverage table')
    df = pd.read_csv(cov_path, sep='\t', usecols=usecols)
    for c in ['identity','query_cov','n_asvs','n_bioprojects','total_reads','prevalence_samples']:
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0)
    df['has_homd_hit'] = df['has_homd_hit'].astype(bool)
    low = df[(~df['has_homd_hit']) | (df['identity'] < 97) | (df['query_cov'] < 99)].copy()
    low = add_geo_counts(low)
    low['identity_class'] = np.select(
        [~low['has_homd_hit'], low['identity'] < 80, low['identity'] < 90, low['identity'] < 95, low['identity'] < 97, low['query_cov'] < 99],
        ['no_hit','below_80','80_90','90_95','95_97','qcov_lt99'],
        default='other_low'
    )
    low['priority_score'] = (
        np.log10(low['prevalence_samples'] + 1)
        * np.log10(low['n_bioprojects'] + 1)
        * np.log10(low['total_reads'] + 1)
        * (1 + np.log10(low['n_countries'] + 1))
    )
    low.to_csv(low_path, sep='\t', index=False)

    priority = low[(low['n_bioprojects'] >= args.min_bioprojects_priority) | (low['prevalence_samples'] >= args.min_prevalence_priority) | (low['total_reads'] >= args.min_reads_priority)].copy()
    priority = priority.sort_values(['priority_score','n_bioprojects','prevalence_samples','total_reads'], ascending=False).head(args.top_n_annotate)
    priority.to_csv(priority_path, sep='\t', index=False)
    write_fasta(priority, fasta_path)

    log(f'low identity/no-hit ASVs: {len(low):,}; priority annotated set: {len(priority):,}')
    if args.force or not (os.path.exists(gg2_hit_path) and os.path.getsize(gg2_hit_path) > 0):
        cmd = [args.vsearch, '--usearch_global', fasta_path, '--db', args.gg2_ref, '--id', '0.7', '--maxaccepts', '1', '--maxrejects', '0', '--strand', 'both', '--userout', gg2_hit_path, '--userfields', 'query+target+id+alnlen+qcov+tcov+ql+tl', '--threads', str(args.threads), '--quiet']
        log('running GG2 annotation for priority low-identity ASVs')
        log(' '.join(cmd))
        subprocess.run(cmd, check=True)
    else:
        log(f'using existing GG2 hits: {gg2_hit_path}')

    hit_cols = ['derep_id','gg2_target','gg2_identity','gg2_alnlen','gg2_qcov','gg2_tcov','query_length','target_length']
    hits = pd.read_csv(gg2_hit_path, sep='\t', names=hit_cols) if os.path.getsize(gg2_hit_path) else pd.DataFrame(columns=hit_cols)
    if len(hits):
        hits = hits.sort_values(['derep_id','gg2_identity','gg2_qcov','gg2_alnlen'], ascending=[True,False,False,False]).drop_duplicates('derep_id')
    tax = load_taxonomy(args.gg2_tax)
    hits['gg2_taxonomy'] = hits['gg2_target'].map(tax).fillna('Unassigned') if len(hits) else []
    ann = priority.merge(hits, on='derep_id', how='left')
    ann['gg2_taxonomy'] = ann['gg2_taxonomy'].fillna('Unassigned')
    ann['gg2_qc_class'] = ann['gg2_taxonomy'].map(classify_qc)
    for col, prefix in [('gg2_domain','d__'),('gg2_phylum','p__'),('gg2_class','c__'),('gg2_order','o__'),('gg2_family','f__'),('gg2_genus','g__'),('gg2_species','s__')]:
        ann[col] = ann['gg2_taxonomy'].map(lambda x, p=prefix: parse_rank(x, p))
    ann.to_csv(annot_path, sep='\t', index=False)

    qc = ann.groupby('gg2_qc_class', dropna=False).agg(derep_sequences=('derep_id','count'), source_asvs=('n_asvs','sum'), total_reads=('total_reads','sum'), prevalence_samples=('prevalence_samples','sum'), max_bioprojects=('n_bioprojects','max')).reset_index().sort_values('total_reads', ascending=False)
    qc.to_csv(qc_path, sep='\t', index=False)

    tax_rows = []
    for rank in ['gg2_phylum','gg2_class','gg2_order','gg2_family','gg2_genus','gg2_species']:
        tmp = ann.groupby(rank, dropna=False).agg(derep_sequences=('derep_id','count'), source_asvs=('n_asvs','sum'), total_reads=('total_reads','sum'), prevalence_samples=('prevalence_samples','sum'), max_bioprojects=('n_bioprojects','max'), max_countries=('n_countries','max')).reset_index().rename(columns={rank:'taxon'})
        tmp.insert(0,'rank',rank.replace('gg2_',''))
        tax_rows.append(tmp)
    pd.concat(tax_rows, ignore_index=True).sort_values(['rank','total_reads'], ascending=[True,False]).to_csv(tax_path, sep='\t', index=False)

    bacteria = ann[ann['gg2_qc_class'].eq('bacterial')].copy()
    parts = []
    for group in ['continents','countries','statuses','disease_studies','sequencing_regions','sites']:
        part = expand_group_tax(bacteria, group, 'gg2_genus')
        if len(part): parts.append(part)
    pd.concat(parts, ignore_index=True).to_csv(group_tax_path, sep='\t', index=False) if parts else pd.DataFrame().to_csv(group_tax_path, sep='\t', index=False)

    stats = {
        'total_derep_sequences': int(len(df)),
        'low_identity_or_nohit_sequences': int(len(low)),
        'priority_annotated_sequences': int(len(priority)),
        'priority_filters': {'min_bioprojects_or': args.min_bioprojects_priority, 'min_prevalence_or': args.min_prevalence_priority, 'min_reads_or': args.min_reads_priority},
        'identity_class_counts': low['identity_class'].value_counts().to_dict(),
        'qc_summary': qc.to_dict(orient='records'),
        'outputs': {'low_all': low_path, 'priority': priority_path, 'annotated': annot_path, 'qc': qc_path, 'taxonomy': tax_path, 'group_taxonomy': group_tax_path}
    }
    with open(stats_path, 'w') as handle:
        json.dump(stats, handle, indent=2)
    log(f'wrote {stats_path}')

if __name__ == '__main__':
    main()
