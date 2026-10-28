#!/usr/bin/env python3
"""
Generate QC metadata tables and BIOM/QZA subsets.

Filtering — dataset-level only, evaluated on oral samples:
  0. CR=0 removal: samples with closed-reference mapping rate == 0 (absent from BIOM)
  1. chimera: remove studies where >30% oral samples have chimera >25%
  2. closed-reference: remove studies where oral aggregate CR mapping <70%

CR=0 filtering ensures metadata matches the closed-reference BIOM table. ASV-based filtering removed (v4).

Outputs (all in data/metadata-and-biom-v2-qc/):
  - metadata_all_qc.tsv, metadata_oral_qc.tsv, metadata_oral_qc_dedup.tsv
  - biom_all_qc.biom, biom_oral_qc.biom, biom_oral_qc_dedup.biom
  - biom_all_qc.qza, genus_table_all_qc.qza, genus_table_relabund_all_qc.qza
  - QC/filter_funnel_summary.tsv, QC/removed_studies_by_step.tsv
"""
import csv
import os
import subprocess
from collections import defaultdict, Counter
from pathlib import Path

import biom
from biom.util import biom_open

BASE = '/beegfs/huangxiaochang/global_oral/global-oral-microbiome'
DATA_DIR = f'{BASE}/data'
OUT_DIR = f'{DATA_DIR}/metadata-and-biom-v2-qc'
RAW = f'{OUT_DIR}/metadata_all_for_qc_raw.tsv'
ALL_QC = f'{OUT_DIR}/metadata_all_qc.tsv'
ORAL_QC = f'{OUT_DIR}/metadata_oral_qc.tsv'
ORAL_DEDUP = f'{OUT_DIR}/metadata_oral_qc_dedup.tsv'
SOURCE_BIOM = '/beegfs/dongbiao/oral/table_gg2_id99_61235_115338_v2.biom'
BIOM_ALL_QC = f'{OUT_DIR}/biom_all_qc.biom'
BIOM_ORAL_QC = f'{OUT_DIR}/biom_oral_qc.biom'
BIOM_ORAL_DEDUP = f'{OUT_DIR}/biom_oral_qc_dedup.biom'
QZA_ALL_QC = f'{OUT_DIR}/biom_all_qc.qza'
QZA_GENUS = f'{OUT_DIR}/genus_table_all_qc.qza'
QZA_GENUS_RELABUND = f'{OUT_DIR}/genus_table_relabund_all_qc.qza'
QC_DIR = f'{OUT_DIR}/QC'
FUNNEL_SUMMARY = f'{QC_DIR}/filter_funnel_summary.tsv'
REMOVED_STUDIES = f'{QC_DIR}/removed_studies_by_step.tsv'


def count_subjects(rows):
    """Number of individuals: distinct subject_id per study, summed across studies.

    subject_id is only unique within a study (the same value recurs in other
    studies), so a global unique count would undercount individuals.
    """
    return len({(r['bioproject_id'], r['subject_id']) for r in rows
                if r.get('subject_id') not in (None, '', 'NA')})


def dedupe_fieldnames(fieldnames):
    seen = set()
    deduped = []
    dropped = []
    for name in fieldnames:
        if name in seen:
            dropped.append(name)
            continue
        seen.add(name)
        deduped.append(name)
    return deduped, dropped


def load_rows():
    rows = []
    with open(RAW, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        raw_fieldnames = reader.fieldnames or []
        fieldnames, dropped_columns = dedupe_fieldnames(raw_fieldnames)
        for row in reader:
            rows.append({name: row.get(name, '') for name in fieldnames})
    return rows, fieldnames, dropped_columns


def count_state(rows):
    return {
        'datasets': len({r['bioproject_id'] for r in rows if r.get('bioproject_id')}),
        'samples_total': len(rows),
        'samples_oral': sum(1 for r in rows if r.get('site_group') == 'oral'),
    }


def add_funnel_row(out, step_order, phase, step_name, criterion, before, after, removed_studies=None):
    b = count_state(before)
    a = count_state(after)
    removed_studies = removed_studies or set()
    out.append({
        'step_order': step_order,
        'phase': phase,
        'step_name': step_name,
        'criterion': criterion,
        'datasets_before': b['datasets'],
        'datasets_after': a['datasets'],
        'datasets_removed': b['datasets'] - a['datasets'],
        'samples_before': b['samples_total'],
        'samples_after': a['samples_total'],
        'samples_removed': b['samples_total'] - a['samples_total'],
        'oral_samples_before': b['samples_oral'],
        'oral_samples_after': a['samples_oral'],
        'oral_samples_removed': b['samples_oral'] - a['samples_oral'],
        'removed_study_count': len(removed_studies),
        'removed_studies': ','.join(sorted(removed_studies)),
    })


def write_tsv(path, rows, fieldnames):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def write_sample_ids(rows, output_path):
    sample_ids = [r['sample_id'] for r in rows]
    Path(output_path).write_text('\n'.join(sample_ids) + '\n')


def subset_biom(sample_rows, id_path, output_biom, full_table):
    biom_ids = set(full_table.ids(axis='sample'))
    sample_keys = []
    valid_rows = []
    for r in sample_rows:
        sid = r['sample_id']
        if sid in biom_ids:
            sample_keys.append(sid)
            valid_rows.append(r)
    n_missing = len(sample_rows) - len(valid_rows)
    if n_missing > 0:
        print(f'    Warning: {n_missing}/{len(sample_rows)} sample IDs not in BIOM')
    if not sample_keys:
        print(f'    ERROR: No valid sample IDs for {output_biom}, skipping')
        return
    Path(id_path).write_text('\n'.join(sample_keys) + '\n')
    subset_table = full_table.filter(sample_keys, axis='sample', inplace=False)
    sorted_table = subset_table.sort_order(sample_keys, axis='sample')
    with biom_open(output_biom, 'w') as handle:
        sorted_table.to_hdf5(handle, 'metadata_qc_pipeline_v4')
    print(f'    Wrote {output_biom}: {len(sample_keys)} samples')


def regenerate_biom_tables(all_rows, oral_rows, dedup_rows):
    Path(QC_DIR).mkdir(parents=True, exist_ok=True)
    print(f'  Loading source BIOM: {SOURCE_BIOM}')
    full_table = biom.load_table(SOURCE_BIOM)
    outputs = [
        (all_rows, f'{QC_DIR}/sample_ids_all_qc.txt', BIOM_ALL_QC),
        (oral_rows, f'{QC_DIR}/sample_ids_oral_qc.txt', BIOM_ORAL_QC),
        (dedup_rows, f'{QC_DIR}/sample_ids_oral_qc_dedup.txt', BIOM_ORAL_DEDUP),
    ]
    for sample_rows, id_path, output_biom in outputs:
        print(f'  Regenerating {output_biom} from {len(sample_rows)} metadata samples')
        subset_biom(sample_rows, id_path, output_biom, full_table)


def remove_studies(rows, flag_col):
    remove = {r['bioproject_id'] for r in rows if r.get(flag_col) == 'yes'}
    after = [r for r in rows if r['bioproject_id'] not in remove]
    return after, remove


def generate_qza(input_biom, output_qza, qza_type='FeatureTable[Frequency]'):
    """Convert BIOM to QZA using qiime2."""
    if os.path.exists(output_qza):
        print(f'    Removing existing {output_qza}')
        os.remove(output_qza)
    cmd = ['qiime', 'tools', 'import',
           '--input-path', input_biom,
           '--output-path', output_qza,
           '--type', qza_type]
    print(f'    Running: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)
    print(f'    Wrote {output_qza}')


def main():
    Path(QC_DIR).mkdir(parents=True, exist_ok=True)
    rows, fieldnames, dropped_columns = load_rows()
    print(f'Read {len(rows)} rows from {RAW}')
    if dropped_columns:
        print(f'Dropped {len(dropped_columns)} duplicated columns from input header')

    funnel = []
    removed_study_rows = []
    current = rows
    add_funnel_row(funnel, 0, 'start', 'Raw merged metadata (DADA2-processed, CR>0)',
                   'Samples with valid CR mapping rate > 0', current, current)

    # ── Dataset-level filters (2 steps) ──
    dataset_steps = [
        (1, 'Dataset chimera filter', 'study_chim_gt25_high',
         'Remove studies where >30% oral samples have chimera >25%'),
        (2, 'Dataset close-reference filter', 'study_cr_mapping_low',
         'Remove studies where oral aggregate CR mapping rate <70%'),
    ]
    for step_order, name, flag, criterion in dataset_steps:
        before = current
        current, removed = remove_studies(current, flag)
        add_funnel_row(funnel, step_order, 'dataset', name, criterion, before, current, removed)
        for study in sorted(removed):
            samples = [r for r in before if r['bioproject_id'] == study]
            removed_study_rows.append({
                'step_order': step_order,
                'step_name': name,
                'bioproject_id': study,
                'samples_total': len(samples),
                'samples_oral': sum(1 for r in samples if r.get('site_group') == 'oral'),
                'study_chim_gt25_high': samples[0].get('study_chim_gt25_high', 'NA') if samples else 'NA',
                'study_cr_mapping_rate': samples[0].get('study_cr_mapping_rate', 'NA') if samples else 'NA',
                'study_cr_mapping_low': samples[0].get('study_cr_mapping_low', 'NA') if samples else 'NA',
            })
        print(f'  {name}: {len(removed)} datasets removed; {len(current)} samples remain')

    # ── No further sample-level filters (CR=0 already removed upstream) ──
    all_qc_rows = current
    write_tsv(ALL_QC, all_qc_rows, fieldnames)

    oral_rows = [r for r in all_qc_rows if r.get('site_group') == 'oral']
    write_tsv(ORAL_QC, oral_rows, fieldnames)

    groups = defaultdict(list)
    for r in oral_rows:
        groups[(r['bioproject_id'], r['subject_id'], r.get('site_cleaned', r.get('site', '')))].append(r)
    dedup_rows = [v[0] for v in groups.values()]
    write_tsv(ORAL_DEDUP, dedup_rows, fieldnames)

    # ── Funnel & removed-studies tables ──
    funnel_fields = list(funnel[0].keys())
    write_tsv(FUNNEL_SUMMARY, funnel, funnel_fields)
    removed_fields = ['step_order', 'step_name', 'bioproject_id', 'samples_total',
                      'samples_oral', 'study_chim_gt25_high', 'study_cr_mapping_rate', 'study_cr_mapping_low']
    write_tsv(REMOVED_STUDIES, removed_study_rows, removed_fields)

    n_all = len(all_qc_rows)
    n_oral = len(oral_rows)
    n_dedup = len(dedup_rows)
    n_subjects = count_subjects(dedup_rows)
    n_studies_all = len(set(r['bioproject_id'] for r in all_qc_rows))
    n_studies_oral = len(set(r['bioproject_id'] for r in oral_rows))
    oral_sites = Counter(r.get('site_cleaned', r.get('site', '')) for r in oral_rows)
    d_sites = Counter(r.get('site_cleaned', r.get('site', '')) for r in dedup_rows)

    print('\n=== Generated ===')
    print(f'  metadata_all_qc.tsv:        {n_all:>7} samples, {n_studies_all} studies')
    print(f'  metadata_oral_qc.tsv:       {n_oral:>7} samples, {n_studies_oral} studies')
    print(f'    site_cleaned distribution: {dict(oral_sites.most_common(10))}')
    print(f'  metadata_oral_qc_dedup.tsv: {n_dedup:>7} samples, {n_subjects} subjects (removed {n_oral - n_dedup} duplicates)')
    print(f'    site_cleaned distribution: {dict(d_sites.most_common(10))}')
    print(f'  Funnel summary: {FUNNEL_SUMMARY}')
    print(f'  Removed studies: {REMOVED_STUDIES}')

    # ── BIOM tables ──
    print(f'\n=== Regenerating BIOM tables ===')
    regenerate_biom_tables(all_qc_rows, oral_rows, dedup_rows)

    # ── QZA files (run separately under qiime2-2022.2 environment) ──
    print(f'\n=== QZA generation ===')
    print(f'  SKIPPED in this run. Activate qiime2-2022.2 and run generate_qza.py separately.')
    print(f'  (biom_all_qc.qza, genus_table_all_qc.qza, genus_table_relabund_all_qc.qza)')


if __name__ == '__main__':
    main()
