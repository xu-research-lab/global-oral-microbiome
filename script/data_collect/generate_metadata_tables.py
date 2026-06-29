#!/usr/bin/env python3
"""
Generate three derived metadata tables from metadata_all_for_qc_raw.tsv:
  1. metadata_all_qc.tsv        - quality-filtered, all body sites
  2. metadata_oral_qc.tsv       - oral subset (by site_group='oral'), not deduplicated
  3. metadata_oral_qc_dedup.tsv - oral subset, deduplicated by (subject_id, site_cleaned)

Filters applied:
  - Exclude abnormal studies where ASV <= 30 in >50% samples OR
    chimera >25% in >30% samples OR study-level closed-reference mapping rate <30%
  - Exclude remaining samples with asv_count <= 30
  - Exclude remaining samples with closed-reference sample mapping rate <30%
  - Exclude remaining samples with NA Observed (alpha diversity) or NA subject_id
  - Oral tables: include all samples where site_group = 'oral'
"""
import csv
import os
import subprocess
from collections import defaultdict, Counter
from pathlib import Path

import biom
from biom.util import biom_open

BASE = '/beegfs/huangxiaochang/global_oral/new_analysis'
RAW = f'{BASE}/metadata_all_for_qc_raw.tsv'
ALL_QC = f'{BASE}/metadata_all_qc.tsv'
ORAL_QC = f'{BASE}/metadata_oral_qc.tsv'
ORAL_DEDUP = f'{BASE}/metadata_oral_qc_dedup.tsv'
SOURCE_BIOM = '/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile/biom_exports/gg2_reference_id_profile.biom'
BIOM_ALL_QC = f'{BASE}/biom_all_qc.biom'
BIOM_ORAL_QC = f'{BASE}/biom_oral_qc.biom'
BIOM_ORAL_DEDUP = f'{BASE}/biom_oral_qc_dedup.biom'
QC_DIR = f'{BASE}/QC'


def dedupe_fieldnames(fieldnames):
    """Keep the first occurrence of duplicated columns."""
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


def get_abnormal_studies(rows):
    """Identify studies to exclude under the current rule:
    remove any study where (based on ORAL samples only):
    >50% oral samples have ASV <= 30, OR
    >30% oral samples have chimera >25%, OR
    study-level closed-reference mapping rate (oral samples) <30%.
    """
    study_data = defaultdict(list)
    for r in rows:
        study_data[r['bioproject_id']].append(r)

    abnormal = set()
    for study, samples in study_data.items():
        if len(samples) == 0:
            continue
        # Use only oral samples for threshold computation
        oral_samples = [r for r in samples if r.get('site_group', '') == 'oral']
        n_oral = len(oral_samples)
        if n_oral == 0:
            continue
        low_asv = sum(1 for r in oral_samples if r.get('asv_le_30') == 'yes') / n_oral
        high_chim = sum(1 for r in oral_samples if r.get('chim_gt_25') == 'yes') / n_oral
        low_cr_mapping = any(r.get('study_cr_mapping_low') == 'yes' for r in samples)
        if low_asv > 0.5 or high_chim > 0.3 or low_cr_mapping:
            abnormal.add(study)
    return abnormal



def write_sample_ids(rows, output_path):
    # BIOM uses "bioproject|sample_id" format; construct keys for subsetting
    biom_keys = [f"{r['bioproject_id']}|{r['sample_id']}" for r in rows]
    Path(output_path).write_text('\n'.join(biom_keys) + '\n')


def subset_biom(sample_rows, id_path, output_biom):
    # BIOM sample IDs are in "bioproject|sample_id" format
    # Build BIOM-compatible keys for each sample
    biom_keys = [f"{r['bioproject_id']}|{r['sample_id']}" for r in sample_rows]

    # Load BIOM first to get valid sample IDs, then intersect with our list
    full_table = biom.load_table(SOURCE_BIOM)
    biom_ids = set(full_table.ids(axis='sample'))

    # Filter to only samples present in the BIOM table
    valid_keys = [k for k in biom_keys if k in biom_ids]
    n_missing = len(biom_keys) - len(valid_keys)
    if n_missing > 0:
        print(f'    Warning: {n_missing} sample IDs not found in BIOM, will be excluded from {output_biom}')

    if not valid_keys:
        print(f'    ERROR: No valid sample IDs for {output_biom}, skipping')
        return

    # Write the BIOM-format keys for subsetting
    Path(id_path).write_text('\n'.join(valid_keys) + '\n')
    tmp_biom = f'{output_biom}.tmp'
    sorted_tmp_biom = f'{output_biom}.sorted.tmp'
    subprocess.run([
        'biom', 'subset-table',
        '-i', SOURCE_BIOM,
        '-a', 'sample',
        '-s', id_path,
        '-o', tmp_biom,
    ], check=True)

    table = biom.load_table(tmp_biom)
    sorted_table = table.sort_order(valid_keys, axis='sample')
    with biom_open(sorted_tmp_biom, 'w') as handle:
        sorted_table.to_hdf5(handle, 'metadata_qc_pipeline')
    os.replace(sorted_tmp_biom, output_biom)
    os.remove(tmp_biom)


def regenerate_biom_tables(step3, oral_rows, dedup_rows):
    Path(QC_DIR).mkdir(parents=True, exist_ok=True)
    outputs = [
        (step3, f'{QC_DIR}/sample_ids_all_qc.txt', BIOM_ALL_QC),
        (oral_rows, f'{QC_DIR}/sample_ids_oral_qc.txt', BIOM_ORAL_QC),
        (dedup_rows, f'{QC_DIR}/sample_ids_oral_qc_dedup.txt', BIOM_ORAL_DEDUP),
    ]
    for sample_rows, id_path, output_biom in outputs:
        print(f'  Regenerating {output_biom} from {len(sample_rows)} samples')
        subset_biom(sample_rows, id_path, output_biom)

def main():
    rows = []
    with open(RAW, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        raw_fieldnames = reader.fieldnames or []
        fieldnames, dropped_columns = dedupe_fieldnames(raw_fieldnames)
        for row in reader:
            rows.append({name: row.get(name, '') for name in fieldnames})
    print(f"Read {len(rows)} rows from {RAW}")
    if dropped_columns:
        print(f"Dropped {len(dropped_columns)} duplicated columns from input header")

    # Identify abnormal studies
    abnormal = get_abnormal_studies(rows)
    print(f"Abnormal studies to exclude: {abnormal}")

    # Filter 1: exclude abnormal studies
    step1 = [r for r in rows if r['bioproject_id'] not in abnormal]
    n_removed_study = len(rows) - len(step1)
    print(f"  After removing abnormal studies: {len(step1)} samples (-{n_removed_study})")

    # Filter 2: exclude asv <= 30
    step2 = [r for r in step1 if r['asv_le_30'] != 'yes']
    n_removed_asv = len(step1) - len(step2)
    print(f"  After removing asv<=30: {len(step2)} samples (-{n_removed_asv})")

    # Filter 3: exclude samples with closed-reference mapping rate < 30%
    step3 = [r for r in step2 if r['sample_cr_mapping_low'] != 'yes']
    n_removed_cr = len(step2) - len(step3)
    print(f"  After removing sample_cr_mapping<30%: {len(step3)} samples (-{n_removed_cr})")

    # Filter 4: exclude NA alpha or NA subject_id
    step4 = [r for r in step3 if r['Observed'] != 'NA' and r['subject_id'] != 'NA']
    n_removed_na = len(step3) - len(step4)
    print(f"  After removing NA alpha/subject_id: {len(step4)} samples (-{n_removed_na})")

    # Table 1: all sites, QC
    with open(ALL_QC, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(step4)

    # Table 2: oral subset (by site_group='oral'), not deduplicated
    oral_rows = [r for r in step4 if r['site_group'] == 'oral']
    with open(ORAL_QC, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(oral_rows)

    # Table 3: oral, deduplicated by (subject_id, site_cleaned)
    groups = defaultdict(list)
    for r in oral_rows:
        key = (r['subject_id'], r['site_cleaned'])
        groups[key].append(r)
    dedup_rows = [v[0] for v in groups.values()]
    with open(ORAL_DEDUP, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(dedup_rows)

    # Summary
    n_all = len(step4)
    n_oral = len(oral_rows)
    n_dedup = len(dedup_rows)
    n_subjects = len(set(r['subject_id'] for r in dedup_rows))
    n_studies_all = len(set(r['bioproject_id'] for r in step4))
    n_studies_oral = len(set(r['bioproject_id'] for r in oral_rows))

    oral_sites = Counter(r['site_cleaned'] for r in oral_rows)
    d_sites = Counter(r['site_cleaned'] for r in dedup_rows)

    print(f"\n=== Generated ===")
    print(f"  metadata_all_qc.tsv:        {n_all:>7} samples, {n_studies_all} studies")
    print(f"  metadata_oral_qc.tsv:       {n_oral:>7} samples, {n_studies_oral} studies")
    print(f"    site_cleaned distribution: {dict(oral_sites.most_common(10))}")
    print(f"  metadata_oral_qc_dedup.tsv: {n_dedup:>7} samples, {n_subjects} subjects (removed {n_oral - n_dedup} duplicates)")
    print(f"    site_cleaned distribution: {dict(d_sites.most_common(10))}")

    print(f"\n=== Regenerating BIOM tables from {SOURCE_BIOM} ===")
    regenerate_biom_tables(step4, oral_rows, dedup_rows)


if __name__ == '__main__':
    main()
