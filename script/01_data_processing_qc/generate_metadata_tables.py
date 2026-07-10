#!/usr/bin/env python3
"""
Generate QC metadata tables and BIOM subsets from metadata_all_for_qc_raw.tsv.

Filtering is intentionally funnel-shaped:
  Dataset-level filters, evaluated on oral samples only:
    1. chimera: remove studies where >30% oral samples have chimera >25%
    2. ASV: remove studies where >50% oral samples have ASV <= 30
    3. closed-reference: remove studies where aggregate GG2 CR mapping <30%

  Sample-level filters, applied after dataset filters:
    4. chimera: keep only samples with chimera flag == no
    5. ASV: keep only samples with ASV flag == no
    6. closed-reference: keep only samples with GG2 CR mapping flag == no

Outputs:
  - metadata_all_qc.tsv
  - metadata_oral_qc.tsv
  - metadata_oral_qc_dedup.tsv
  - biom_all_qc.biom, biom_oral_qc.biom, biom_oral_qc_dedup.biom
  - QC/filter_funnel_summary.tsv
  - QC/removed_studies_by_step.tsv
"""
import csv
from collections import defaultdict, Counter
from pathlib import Path

import biom
from biom.util import biom_open

BASE = '/beegfs/huangxiaochang/global_oral/new_analysis'
RAW = f'{BASE}/metadata_all_for_qc_raw.tsv'
ALL_QC = f'{BASE}/metadata_all_qc.tsv'
ORAL_QC = f'{BASE}/metadata_oral_qc.tsv'
ORAL_DEDUP = f'{BASE}/metadata_oral_qc_dedup.tsv'
SOURCE_BIOM = '/beegfs/dongbiao/oral/table_gg2_id99_61235_115338_v2.biom'
BIOM_ALL_QC = f'{BASE}/biom_all_qc.biom'
BIOM_ORAL_QC = f'{BASE}/biom_oral_qc.biom'
BIOM_ORAL_DEDUP = f'{BASE}/biom_oral_qc_dedup.biom'
QC_DIR = f'{BASE}/QC'
FUNNEL_SUMMARY = f'{QC_DIR}/filter_funnel_summary.tsv'
REMOVED_STUDIES = f'{QC_DIR}/removed_studies_by_step.tsv'
CR_THRESHOLD = 0.30


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
    sample_ids = [r['sample_id'] for r in sample_rows]
    biom_ids = set(full_table.ids(axis='sample'))
    valid_ids = [sid for sid in sample_ids if sid in biom_ids]
    n_missing = len(sample_ids) - len(valid_ids)
    if n_missing > 0:
        print(f'    Warning: {n_missing} sample IDs not found in BIOM, will be excluded from {output_biom}')
    if not valid_ids:
        print(f'    ERROR: No valid sample IDs for {output_biom}, skipping')
        return
    Path(id_path).write_text('\n'.join(valid_ids) + '\n')
    subset_table = full_table.filter(valid_ids, axis='sample', inplace=False)
    sorted_table = subset_table.sort_order(valid_ids, axis='sample')
    with biom_open(output_biom, 'w') as handle:
        sorted_table.to_hdf5(handle, 'metadata_qc_pipeline_gg2_v2_cr30')
    print(f'    Wrote {output_biom}: {len(valid_ids)} samples')

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


def main():
    Path(QC_DIR).mkdir(parents=True, exist_ok=True)
    rows, fieldnames, dropped_columns = load_rows()
    print(f'Read {len(rows)} rows from {RAW}')
    if dropped_columns:
        print(f'Dropped {len(dropped_columns)} duplicated columns from input header')

    funnel = []
    removed_study_rows = []
    current = rows
    add_funnel_row(funnel, 0, 'start', 'Raw merged metadata', 'No filtering', current, current)

    dataset_steps = [
        (1, 'Dataset chimera filter', 'study_chim_gt25_high', 'Remove datasets where >30% oral samples have chimera >25%'),
        (2, 'Dataset ASV filter', 'study_asv_lt30_high', 'Remove datasets where >50% oral samples have ASV <=30'),
        (3, 'Dataset close-reference filter', 'study_cr_mapping_low', f'Remove datasets with aggregate oral GG2 CR mapping <{CR_THRESHOLD:.0%}'),
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
                'study_asv_lt30_high': samples[0].get('study_asv_lt30_high', 'NA') if samples else 'NA',
                'study_cr_mapping_rate': samples[0].get('study_cr_mapping_rate', 'NA') if samples else 'NA',
                'study_cr_mapping_low': samples[0].get('study_cr_mapping_low', 'NA') if samples else 'NA',
            })
        print(f'  {name}: {len(removed)} datasets removed; {len(current)} samples remain')

    sample_steps = [
        (4, 'Sample chimera filter', lambda r: r.get('chim_gt_25') == 'no', 'Keep only samples with chimera <=25% and non-missing chimera data'),
        (5, 'Sample ASV filter', lambda r: r.get('asv_le_30') == 'no', 'Keep only samples with ASV >30 and non-missing ASV count'),
        (6, 'Sample close-reference filter', lambda r: r.get('sample_cr_mapping_low') == 'no', f'Keep only samples with GG2 CR mapping >={CR_THRESHOLD:.0%} and non-missing CR data'),
    ]
    for step_order, name, keep_fn, criterion in sample_steps:
        before = current
        current = [r for r in current if keep_fn(r)]
        add_funnel_row(funnel, step_order, 'sample', name, criterion, before, current)
        print(f'  {name}: {len(before) - len(current)} samples removed; {len(current)} samples remain')

    all_qc_rows = current
    write_tsv(ALL_QC, all_qc_rows, fieldnames)

    oral_rows = [r for r in all_qc_rows if r.get('site_group') == 'oral']
    write_tsv(ORAL_QC, oral_rows, fieldnames)

    groups = defaultdict(list)
    for r in oral_rows:
        groups[(r['subject_id'], r['site_cleaned'])].append(r)
    dedup_rows = [v[0] for v in groups.values()]
    write_tsv(ORAL_DEDUP, dedup_rows, fieldnames)

    funnel_fields = list(funnel[0].keys())
    write_tsv(FUNNEL_SUMMARY, funnel, funnel_fields)
    removed_fields = ['step_order', 'step_name', 'bioproject_id', 'samples_total', 'samples_oral', 'study_chim_gt25_high', 'study_asv_lt30_high', 'study_cr_mapping_rate', 'study_cr_mapping_low']
    write_tsv(REMOVED_STUDIES, removed_study_rows, removed_fields)

    n_all = len(all_qc_rows)
    n_oral = len(oral_rows)
    n_dedup = len(dedup_rows)
    n_subjects = len(set(r['subject_id'] for r in dedup_rows))
    n_studies_all = len(set(r['bioproject_id'] for r in all_qc_rows))
    n_studies_oral = len(set(r['bioproject_id'] for r in oral_rows))
    oral_sites = Counter(r['site_cleaned'] for r in oral_rows)
    d_sites = Counter(r['site_cleaned'] for r in dedup_rows)

    print('\n=== Generated ===')
    print(f'  metadata_all_qc.tsv:        {n_all:>7} samples, {n_studies_all} studies')
    print(f'  metadata_oral_qc.tsv:       {n_oral:>7} samples, {n_studies_oral} studies')
    print(f'    site_cleaned distribution: {dict(oral_sites.most_common(10))}')
    print(f'  metadata_oral_qc_dedup.tsv: {n_dedup:>7} samples, {n_subjects} subjects (removed {n_oral - n_dedup} duplicates)')
    print(f'    site_cleaned distribution: {dict(d_sites.most_common(10))}')
    print(f'  Funnel summary: {FUNNEL_SUMMARY}')
    print(f'  Removed studies: {REMOVED_STUDIES}')

    print(f'\n=== Regenerating BIOM tables from exact-derep close-reference BIOM ===')
    regenerate_biom_tables(all_qc_rows, oral_rows, dedup_rows)


if __name__ == '__main__':
    main()
