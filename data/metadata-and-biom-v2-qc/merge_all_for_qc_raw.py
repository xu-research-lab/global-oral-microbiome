#!/usr/bin/env python3
"""
Merge all input sources into metadata_all_for_qc_raw.tsv.

Authoritative base table:
  - metadata-v2-merge-cleaned.csv (repo)

All other sources are left-joined onto the base metadata by
(bioproject, sample_id). After merging, only samples with valid
sample_cr_mapping_rate (DADA2-processed) are retained.

Study-level QC flags (chimera and CR mapping) are computed from
oral samples only (site_group='oral').
"""
import csv
import re
from collections import defaultdict

BASE = '/beegfs/huangxiaochang/global_oral/global-oral-microbiome'
DATA_DIR = f'{BASE}/data'
OUT_DIR = f'{DATA_DIR}/metadata-and-biom-v2-qc'
META = f'{DATA_DIR}/metadata-v2-merge-cleaned.csv'
CR_TABLE = '/beegfs/huangxiaochang/global_oral/new_analysis/sample_closed_reference_mapping_rate_v2.tsv'
STUDY_CR_TABLE = '/beegfs/huangxiaochang/global_oral/new_analysis/study_closed_reference_mapping_rate_v2.tsv'
ASV_TABLE = '/beegfs/huangxiaochang/global_oral/new_analysis/sample_asv_counts_new.tsv'
CHIM_TABLE = '/beegfs/huangxiaochang/global_oral/new_analysis/oral_chim_perc_YL_260519.csv'
READS_TABLE = '/beegfs/huangxiaochang/global_oral/new_analysis/raw_reads_stats.tsv'
ALPHA_TABLE = '/beegfs/huangxiaochang/global_oral/new_analysis/alpha_diversity_data.csv'
OUT = f'{OUT_DIR}/metadata_all_for_qc_raw.tsv'

# ── QC thresholds ───────────────────────────────────────────────────
CR_THRESHOLD = 0.70          # study-level oral CR mapping rate
CHIM_STUDY_RATIO = 0.30      # max fraction of oral samples with chim >25%


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


def to_float(value):
    try:
        if norm(value) in {'', 'NA'}:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    import os
    os.makedirs(OUT_DIR, exist_ok=True)

    print('=' * 60)
    print('Step 1: Load authoritative base metadata')
    meta_rows = load_csv(META)
    meta_fieldnames = [f for f in meta_rows[0].keys()
                       if f not in ('bioproject', 'sample_id', 'subject_id')]
    base_studies = unique_nonempty(r.get('bioproject') for r in meta_rows)
    base_samples = unique_nonempty(r.get('sample_id') for r in meta_rows)
    base_subjects = unique_nonempty(r.get('subject_id') for r in meta_rows)
    print(f'  base rows: {len(meta_rows)}')
    print(f'  base studies: {base_studies}')
    print(f'  base samples: {base_samples}')
    print(f'  base subjects: {base_subjects}')

    print('\nStep 2: Load ASV counts')
    asv_rows = load_tsv(ASV_TABLE)
    asv_map = {(norm(r['bioproject_id']), norm(r['sample_id'])): norm(r['asv_count'])
               for r in asv_rows}
    print(f'  {len(asv_rows)} entries, {len(asv_map)} unique')

    print('\nStep 3: Load chimera data')
    chim_rows = load_csv(CHIM_TABLE)
    chim_map = {}
    for r in chim_rows:
        key = (norm(r['Project']), strip_fastq_suffix(r['SampleID']))
        chim_map[key] = norm(r['chim_perc'])
    print(f'  {len(chim_rows)} entries, {len(chim_map)} unique after stripping suffix')

    print('\nStep 4: Load raw reads')
    reads_rows = load_tsv(READS_TABLE)
    reads_map = {(norm(r['Project']), norm(r['Sample'])): norm(r['Total_Reads'])
                 for r in reads_rows}
    print(f'  {len(reads_rows)} entries, {len(reads_map)} unique')

    print('\nStep 5: Load alpha diversity data')
    alpha_rows = load_csv(ALPHA_TABLE)
    meta_column_set = set(meta_rows[0].keys())
    alpha_fieldnames = [f for f in alpha_rows[0].keys()
                        if f not in meta_column_set
                        and f not in ('bioproject', 'sample_id')]
    alpha_map = {}
    for r in alpha_rows:
        key = (norm(r['bioproject']), norm(r['sample_id']))
        alpha_map[key] = {col: norm(r[col]) for col in alpha_fieldnames}
    print(f'  {len(alpha_rows)} entries, {len(alpha_map)} unique')
    print(f'  alpha-only fields kept: {alpha_fieldnames}')

    print('\nStep 6: Load GG2 v2 closed-reference mapping data')
    cr_rows = load_tsv(CR_TABLE)
    cr_sample_map = {}
    for r in cr_rows:
        key = (norm(r['bioproject']), norm(r['sample_id']))
        cr_sample_map[key] = r
    # Study-level CR table (pre-computed from all samples, for reference only)
    study_cr_rows = load_tsv(STUDY_CR_TABLE)
    study_cr_all = {norm(r['bioproject']): r for r in study_cr_rows}
    print(f'  {len(cr_rows)} sample entries, {len(cr_sample_map)} unique sample keys')
    print(f'  {len(study_cr_rows)} study entries (all-sample aggregate, reference only)')

    print('\nStep 7: Left join all inputs onto base metadata')
    merged = []
    matched_asv = matched_chim = matched_reads = matched_alpha = matched_cr = 0

    for meta in meta_rows:
        bioproject = norm(meta.get('bioproject'))
        sample_id = norm(meta.get('sample_id'))
        key = (bioproject, sample_id)
        row = {'bioproject_id': bioproject, 'sample_id': sample_id,
               'subject_id': norm(meta.get('subject_id')) or 'NA'}

        asv_val = asv_map.get(key, 'NA')
        matched_asv += asv_val != 'NA'
        row['asv_count'] = asv_val

        chim_val = chim_map.get(key, 'NA')
        matched_chim += chim_val != 'NA'
        row['chim_perc'] = chim_val

        reads_val = reads_map.get(key, 'NA')
        matched_reads += reads_val != 'NA'
        row['total_reads'] = reads_val

        cr_data = cr_sample_map.get(key)
        if cr_data:
            matched_cr += 1
            row['sample_cr_mapping_rate'] = norm(cr_data.get('mapped_fraction'))
        else:
            row['sample_cr_mapping_rate'] = 'NA'

        alpha_data = alpha_map.get(key, {})
        matched_alpha += bool(alpha_data)
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

    # ── Filter to DADA2-processed samples (have CR mapping data) ──
    print('\nStep 8: Retain only samples with valid CR mapping data (DADA2-processed)')
    merged_cr = [r for r in merged if r['sample_cr_mapping_rate'] != 'NA']
    n_excluded = len(merged) - len(merged_cr)
    print(f'  retained: {len(merged_cr)} samples')
    print(f'  excluded (no CR data): {n_excluded}')
    merged = merged_cr

    # ── Compute study-level QC flags (oral samples only) ──
    print('\nStep 9: Compute study-level QC flags (oral samples only)')
    study_samples = defaultdict(list)
    for row in merged:
        study_samples[row['bioproject_id']].append(row)

    # Compute study-level oral CR mapping rate from sample-level data
    study_cr_mapping_rate = {}
    for study, samples in study_samples.items():
        oral_samples = [s for s in samples if s.get('site_group', '') == 'oral']
        total_reads_sum = 0.0
        cr_reads_sum = 0.0
        for s in oral_samples:
            cr_data = cr_sample_map.get((study, s['sample_id']))
            if cr_data:
                tr = to_float(cr_data.get('total_reads_dada2'))
                cr = to_float(cr_data.get('closed_reference_reads'))
                if tr is not None and cr is not None:
                    total_reads_sum += tr
                    cr_reads_sum += cr
        if total_reads_sum > 0:
            study_cr_mapping_rate[study] = cr_reads_sum / total_reads_sum
        else:
            study_cr_mapping_rate[study] = None

    study_chim_flag = {}
    study_cr_flag = {}
    for study, samples in study_samples.items():
        oral_samples = [s for s in samples if s.get('site_group', '') == 'oral']
        n_oral = len(oral_samples)
        n_chim_gt_25 = sum(1 for s in oral_samples
                          if to_float(s['chim_perc']) is not None
                          and to_float(s['chim_perc']) > 25)
        study_chim_flag[study] = 'yes' if n_oral and (n_chim_gt_25 / n_oral > CHIM_STUDY_RATIO) else 'no'
        cr_rate = study_cr_mapping_rate.get(study)
        study_cr_flag[study] = ('yes' if cr_rate < CR_THRESHOLD else 'no') if cr_rate is not None else 'NA'

    # ── Per-row flags ──
    for row in merged:
        study = row['bioproject_id']
        chim = to_float(row['chim_perc'])
        cr = to_float(row['sample_cr_mapping_rate'])
        row['chim_gt_25'] = 'yes' if chim is not None and chim > 25 else ('NA' if chim is None else 'no')
        row['study_chim_gt25_high'] = study_chim_flag.get(study, 'NA')
        row['sample_cr_mapping_low'] = 'yes' if cr is not None and cr < CR_THRESHOLD else ('NA' if cr is None else 'no')
        rate = study_cr_mapping_rate.get(study)
        row['study_cr_mapping_rate'] = f'{rate:.6f}' if rate is not None else 'NA'
        row['study_cr_mapping_low'] = study_cr_flag.get(study, 'NA')

    # ── Write ──
    base_fields = ['bioproject_id', 'sample_id', 'subject_id',
                   'chim_perc', 'asv_count', 'total_reads',
                   'sample_cr_mapping_rate', 'study_cr_mapping_rate']
    qc_fields = ['chim_gt_25', 'study_chim_gt25_high',
                 'sample_cr_mapping_low', 'study_cr_mapping_low']
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

    print(f'\nStep 10: Written {out_rows} rows to {OUT}')
    print(f'  Columns: {len(out_fields)}')
    if dropped_out_fields:
        print(f'  Dropped duplicate output columns: {len(dropped_out_fields)}')

    n_chim_gt_25 = sum(1 for r in merged if r['chim_gt_25'] == 'yes')
    n_sample_cr_low = sum(1 for r in merged if r['sample_cr_mapping_low'] == 'yes')
    n_oral_raw = sum(1 for r in merged if r.get('site_group', '') == 'oral')
    print('\n=== QC Summary ===')
    print(f'  Total samples:              {out_rows}')
    print(f'  Total studies:              {out_studies}')
    print(f'  Total subjects:             {out_subjects}')
    print(f'  chim > 25%:                 {n_chim_gt_25} ({100*n_chim_gt_25/out_rows:.1f}%)')
    print(f'  sample cr_map<{CR_THRESHOLD*100:.0f}%:          {n_sample_cr_low} ({100*n_sample_cr_low/out_rows:.1f}%)')
    print(f'  studies chim>25>30% (oral): {sum(1 for v in study_chim_flag.values() if v == "yes")}')
    print(f'  studies cr_map<{CR_THRESHOLD*100:.0f}% (oral):   {sum(1 for v in study_cr_flag.values() if v == "yes")}')
    print(f'  Oral raw (site_group=oral): {n_oral_raw}')


if __name__ == '__main__':
    main()
