#!/usr/bin/env python3
"""
Merge all input sources into metadata_all_for_qc_raw.tsv.

Authoritative base table:
  - metadata-v2-merge-cleaned.csv

All other sources are left-joined onto the base metadata by
(bioproject, sample_id). The output must therefore preserve the same
row count, unique sample count, unique subject count, and unique study
count as the base metadata table.
"""
import csv
import re
from collections import defaultdict

BASE = '/beegfs/huangxiaochang/global_oral/new_analysis'
META = '/beegfs/huangxiaochang/global_oral/metadata_v2_merge/metadata-v2-merge-cleaned.csv'
CR_TABLE = '/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile/all_samples_gg2_closed_reference_mapping_fraction.tsv'
OUT = f'{BASE}/metadata_all_for_qc_raw.tsv'


def norm(value):
    if value is None:
        return ''
    return str(value).strip()


def strip_fastq_suffix(sample_id):
    sample_id = norm(sample_id)
    sample_id = re.sub(r'_1\.fastq$', '', sample_id)
    sample_id = re.sub(r'\.fastq$', '', sample_id)
    sample_id = re.sub(r'_R1\.fastq\.gz$', '', sample_id)
    return sample_id


def load_tsv(path):
    with open(path, 'r', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def load_csv(path):
    with open(path, 'r', newline='') as f:
        return list(csv.DictReader(f))


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


def unique_nonempty(values):
    return len({norm(v) for v in values if norm(v) and norm(v) != 'NA'})


def main():
    print('=' * 60)
    print('Step 1: Load authoritative base metadata')
    meta_rows = load_csv(META)
    meta_fieldnames = [f for f in meta_rows[0].keys() if f not in ('bioproject', 'sample_id', 'subject_id')]
    base_studies = unique_nonempty(r.get('bioproject') for r in meta_rows)
    base_samples = unique_nonempty(r.get('sample_id') for r in meta_rows)
    base_subjects = unique_nonempty(r.get('subject_id') for r in meta_rows)
    print(f'  base rows: {len(meta_rows)}')
    print(f'  base studies: {base_studies}')
    print(f'  base samples: {base_samples}')
    print(f'  base subjects: {base_subjects}')

    print('\nStep 2: Load ASV counts')
    asv_rows = load_tsv(f'{BASE}/sample_asv_counts_new.tsv')
    asv_map = {(norm(r['bioproject_id']), norm(r['sample_id'])): norm(r['asv_count']) for r in asv_rows}
    print(f'  {len(asv_rows)} entries, {len(asv_map)} unique')

    print('\nStep 3: Load chimera data')
    chim_rows = load_csv(f'{BASE}/oral_chim_perc_YL_260519.csv')
    chim_map = {}
    for r in chim_rows:
        key = (norm(r['Project']), strip_fastq_suffix(r['SampleID']))
        chim_map[key] = norm(r['chim_perc'])
    print(f'  {len(chim_rows)} entries, {len(chim_map)} unique after stripping suffix')

    print('\nStep 4: Load raw reads')
    reads_rows = load_tsv(f'{BASE}/raw_reads_stats.tsv')
    reads_map = {(norm(r['Project']), norm(r['Sample'])): norm(r['Total_Reads']) for r in reads_rows}
    print(f'  {len(reads_rows)} entries, {len(reads_map)} unique')

    print('\nStep 5: Load alpha diversity data')
    alpha_rows = load_csv(f'{BASE}/alpha_diversity_data.csv')
    meta_column_set = set(meta_rows[0].keys())
    alpha_fieldnames = [f for f in alpha_rows[0].keys() if f not in meta_column_set and f not in ('bioproject', 'sample_id')]
    alpha_map = {}
    for r in alpha_rows:
        key = (norm(r['bioproject']), norm(r['sample_id']))
        alpha_map[key] = {col: norm(r[col]) for col in alpha_fieldnames}
    print(f'  {len(alpha_rows)} entries, {len(alpha_map)} unique')
    print(f'  alpha-only fields kept: {alpha_fieldnames}')

    print('\nStep 6: Load sample-level closed-reference mapping data')
    cr_rows = load_tsv(CR_TABLE)
    cr_sample_map = {}  # (bioproject, sample_id) -> dict of all CR columns
    for r in cr_rows:
        key = (norm(r['bioproject']), norm(r['sample_id']))
        cr_sample_map[key] = r
    print(f'  {len(cr_rows)} entries, {len(cr_sample_map)} unique sample keys')

    print('\nStep 7: Left join all inputs onto base metadata')
    merged = []
    matched_asv = matched_chim = matched_reads = matched_alpha = matched_cr = 0

    for meta in meta_rows:
        bioproject = norm(meta.get('bioproject'))
        sample_id = norm(meta.get('sample_id'))
        key = (bioproject, sample_id)

        row = {'bioproject_id': bioproject, 'sample_id': sample_id, 'subject_id': norm(meta.get('subject_id')) or 'NA'}

        asv_val = asv_map.get(key, 'NA')
        if asv_val != 'NA':
            matched_asv += 1
        row['asv_count'] = asv_val

        chim_val = chim_map.get(key, 'NA')
        if chim_val != 'NA':
            matched_chim += 1
        row['chim_perc'] = chim_val

        reads_val = reads_map.get(key, 'NA')
        if reads_val != 'NA':
            matched_reads += 1
        row['total_reads'] = reads_val

        cr_data = cr_sample_map.get(key)
        if cr_data:
            matched_cr += 1
            row['sample_cr_mapping_rate'] = norm(cr_data['mapped_fraction'])
        else:
            row['sample_cr_mapping_rate'] = 'NA'

        alpha_data = alpha_map.get(key, {})
        if alpha_data:
            matched_alpha += 1
        for col in alpha_fieldnames:
            row[col] = alpha_data.get(col, 'NA')

        for col in meta_fieldnames:
            meta_value = norm(meta.get(col, 'NA'))
            row[col] = meta_value if meta_value else 'NA'

        merged.append(row)

    print(f'  matched ASV:              {matched_asv}/{len(merged)}')
    print(f'  matched chimera:          {matched_chim}/{len(merged)}')
    print(f'  matched reads:            {matched_reads}/{len(merged)}')
    print(f'  matched alpha:            {matched_alpha}/{len(merged)}')
    print(f'  matched CR mapping:       {matched_cr}/{len(merged)}')

    print('\nStep 8: Compute QC flags (study-level thresholds based on ORAL samples only)')
    study_samples = defaultdict(list)
    for row in merged:
        study_samples[row['bioproject_id']].append(row)

    # Compute study-level aggregate CR mapping rate from ORAL sample-level data only
    study_cr_mapping_rate = {}
    for study, samples in study_samples.items():
        oral_samples = [s for s in samples if s.get('site_group', '') == 'oral']
        total_reads_sum = 0.0
        cr_reads_sum = 0.0
        for s in oral_samples:
            key = (study, s['sample_id'])
            cr_data = cr_sample_map.get(key)
            if cr_data:
                try:
                    total_reads_sum += float(cr_data['total_reads'])
                    cr_reads_sum += float(cr_data['closed_reference_reads'])
                except (ValueError, TypeError):
                    pass
        if total_reads_sum > 0:
            study_cr_mapping_rate[study] = cr_reads_sum / total_reads_sum
        else:
            study_cr_mapping_rate[study] = None

    study_asv_flag = {}
    study_chim_flag = {}
    study_cr_flag = {}

    for study, samples in study_samples.items():
        # Use only oral samples for study-level threshold computation
        oral_samples = [s for s in samples if s.get('site_group', '') == 'oral']
        n_oral = len(oral_samples)

        n_asv_le_30 = sum(1 for s in oral_samples if s['asv_count'] != 'NA' and float(s['asv_count']) <= 30)
        n_chim_gt_25 = sum(1 for s in oral_samples if s['chim_perc'] != 'NA' and float(s['chim_perc']) > 25)

        study_asv_flag[study] = 'yes' if n_oral and (n_asv_le_30 / n_oral > 0.5) else 'no'
        study_chim_flag[study] = 'yes' if n_oral and (n_chim_gt_25 / n_oral > 0.3) else 'no'

        cr_rate = study_cr_mapping_rate.get(study)
        if cr_rate is not None:
            study_cr_flag[study] = 'yes' if cr_rate < 0.3 else 'no'
        else:
            study_cr_flag[study] = 'NA'

    for row in merged:
        study = row['bioproject_id']
        try:
            row['asv_le_30'] = 'yes' if float(row['asv_count']) <= 30 else 'no'
        except (ValueError, TypeError):
            row['asv_le_30'] = 'NA'
        row['study_asv_lt30_high'] = study_asv_flag.get(study, 'NA')
        try:
            row['chim_gt_25'] = 'yes' if float(row['chim_perc']) > 25 else 'no'
        except (ValueError, TypeError):
            row['chim_gt_25'] = 'NA'
        row['study_chim_gt25_high'] = study_chim_flag.get(study, 'NA')

        # Sample-level CR mapping flag (new)
        try:
            row['sample_cr_mapping_low'] = 'yes' if float(row['sample_cr_mapping_rate']) < 0.3 else 'no'
        except (ValueError, TypeError):
            row['sample_cr_mapping_low'] = 'NA'

        # Study-level aggregate CR mapping rate
        try:
            row['study_cr_mapping_rate'] = '{:.6f}'.format(study_cr_mapping_rate[study])
        except (KeyError, TypeError):
            row['study_cr_mapping_rate'] = 'NA'
        row['study_cr_mapping_low'] = study_cr_flag.get(study, 'NA')

    base_fields = ['bioproject_id', 'sample_id', 'subject_id', 'chim_perc', 'asv_count', 'total_reads', 'sample_cr_mapping_rate', 'study_cr_mapping_rate']
    qc_fields = ['asv_le_30', 'study_asv_lt30_high', 'chim_gt_25', 'study_chim_gt25_high', 'sample_cr_mapping_low', 'study_cr_mapping_low']
    out_fields_raw = base_fields + alpha_fieldnames + meta_fieldnames + qc_fields
    out_fields, dropped_out_fields = dedupe_fieldnames(out_fields_raw)

    with open(OUT, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(merged)

    out_rows = len(merged)
    out_studies = unique_nonempty(r.get('bioproject_id') for r in merged)
    out_samples = unique_nonempty(r.get('sample_id') for r in merged)
    out_subjects = unique_nonempty(r.get('subject_id') for r in merged)

    print(f'\nStep 9: Written {out_rows} rows to {OUT}')
    print(f'  Columns: {len(out_fields)}')
    if dropped_out_fields:
        print(f'  Dropped duplicate output columns: {len(dropped_out_fields)}')

    print('\n=== Output/Base Consistency Check ===')
    print(f'  base rows / output rows: {len(meta_rows)} / {out_rows}')
    print(f'  base studies / output studies: {base_studies} / {out_studies}')
    print(f'  base samples / output samples: {base_samples} / {out_samples}')
    print(f'  base subjects / output subjects: {base_subjects} / {out_subjects}')

    if (len(meta_rows), base_studies, base_samples, base_subjects) != (out_rows, out_studies, out_samples, out_subjects):
        raise SystemExit('ERROR: output counts do not match authoritative metadata base table')

    n_asv_le_30 = sum(1 for r in merged if r['asv_le_30'] == 'yes')
    n_chim_gt_25 = sum(1 for r in merged if r['chim_gt_25'] == 'yes')
    n_sample_cr_low = sum(1 for r in merged if r['sample_cr_mapping_low'] == 'yes')
    n_oral_raw = sum(1 for r in merged if r.get('site_group', '') == 'oral')

    print('\n=== QC Summary ===')
    print(f'  Total samples:       {out_rows}')
    print(f'  Total studies:       {out_studies}')
    print(f'  asv <= 30:           {n_asv_le_30} ({100 * n_asv_le_30 / out_rows:.1f}%)')
    print(f'  chim > 25%:          {n_chim_gt_25} ({100 * n_chim_gt_25 / out_rows:.1f}%)')
    print(f'  sample cr_map<30%:   {n_sample_cr_low} ({100 * n_sample_cr_low / out_rows:.1f}%)')
    print(f'  studies asv_lt30>50%:{sum(1 for v in study_asv_flag.values() if v == "yes")}')
    print(f'  studies chim>25>30%: {sum(1 for v in study_chim_flag.values() if v == "yes")}')
    print(f'  studies cr_map<30%:  {sum(1 for v in study_cr_flag.values() if v == "yes")}')
    print(f'  Oral raw (site_group=oral): {n_oral_raw}')


if __name__ == '__main__':
    main()
