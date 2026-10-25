#!/usr/bin/env python3
"""
Extract actual Run/sample IDs for control samples from raw SraRunTable.csv files.

Uses the same regex categories as 25_negative_control_audit.py but outputs
per-sample IDs instead of per-project counts.

Output:
- control_sample_ids.tsv: sample_id, bioproject, category, matched_value, column, file_source
- control_samples_by_category.tsv: summary counts
"""
import pandas as pd
import numpy as np
from pathlib import Path
import csv
import re
import json
import glob

PROJECT_DIR = Path('/beegfs/huangxiaochang/global_oral/global-oral-microbiome/data/projects')
METADATA_CSV = '/beegfs/huangxiaochang/global_oral/global-oral-microbiome/data/metadata-v2-merge-oral.csv'
OUTDIR = Path('/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis')

# Same categories as 25_negative_control_audit.py
CATEGORIES = [
    ("extraction_blank", re.compile(r"\bext(raction)?[\s_\-\.]*(blank|b\d|ctrl)|extblank|extractionblank", re.I)),
    ("blank",            re.compile(r"\bblank\b|\bblnk\b", re.I)),
    ("ntc_no_template",  re.compile(r"\bntc\b|no[\s\-_]*template", re.I)),
    ("negative_control", re.compile(r"negative[\s_\-]*control|neg[\s_\-\.]*ctrl|neg[\s_\-\.]*control", re.I)),
    ("mock_community",   re.compile(r"mock[\s_\-]*(community|comm|control)", re.I)),
    ("reagent_water",    re.compile(r"reagent[\s_\-]*(control|blank)|sterile[\s\-_]*water|water[\s_\-]*(control|blank)|no[\s\-_]*dna", re.I)),
]

STRONG_CATEGORIES = {"extraction_blank", "blank", "ntc_no_template"}

# Columns most likely to contain the Run/sample ID
ID_COLUMNS = ['Run', 'run', 'Sample Name', 'Sample_Name', 'Sample_name',
              'sample_name', 'sample_id', 'SampleId', 'library_id']

# Columns that give additional context about the sample type
CONTEXT_COLUMNS = ['Sample_Name', 'Sample_name', 'Sample Name',
                   'sample_type', 'SampleType', 'SampleDescription',
                   'isolation_source', 'IsolationSource', 'Isolation_source',
                   'LibrarySource', 'LibrarySelection',
                   'Library Name', 'library_name',
                   'host_sex', 'Host', 'HOST',
                   'Submitter_Id', 'submitted_subject_id',
                   'material_type', 'material_type (exp)',
                   'source_material_id', 'BioSample',
                   'description', 'Description', 'title', 'Title',
                   'miscellaneous_parameter', 'study_design']

OUTDIR.mkdir(parents=True, exist_ok=True)

# Load metadata for cross-reference
print("Loading metadata...")
md = pd.read_csv(METADATA_CSV, dtype=str, low_memory=False)
md_samples = set(zip(md['bioproject'].astype(str), md['sample_id'].astype(str)))
print(f"  {len(md):,} metadata samples")

results = []
summary_by_category = {cat: 0 for cat, _ in CATEGORIES}
projects_with_controls = set()

for proj_dir in sorted(PROJECT_DIR.iterdir()):
    if not proj_dir.is_dir():
        continue
    bp = proj_dir.name

    # Find all relevant CSV files
    files_to_scan = []
    for name in ['SraRunTable.csv', 'metadata.csv']:
        f = proj_dir / name
        if f.exists():
            files_to_scan.append(f)
    for f in glob.glob(str(proj_dir / '*[Ss]ra*un*able*.csv')):
        if Path(f) not in files_to_scan:
            files_to_scan.append(Path(f))

    for fpath in files_to_scan:
        try:
            df = pd.read_csv(fpath, dtype=str, low_memory=False)
        except Exception:
            continue

        # Find the Run/sample ID column
        id_col = None
        for col in ID_COLUMNS:
            if col in df.columns:
                id_col = col
                break
        if id_col is None:
            # Use first column as fallback
            id_col = df.columns[0]

        for idx, row in df.iterrows():
            sample_id = str(row.get(id_col, f'row_{idx}'))
            line = ' '.join(str(v) for v in row.values if pd.notna(v))

            for cat, pat in CATEGORIES:
                m = pat.search(line)
                if m:
                    # Find which column(s) matched
                    matched_cols = []
                    matched_vals = []
                    for col in df.columns:
                        val = str(row.get(col, ''))
                        if pat.search(val):
                            matched_cols.append(col)
                            matched_vals.append(val[:200])  # truncate

                    # Get context
                    context = {}
                    for col in CONTEXT_COLUMNS:
                        if col in df.columns:
                            val = str(row[col])
                            if val and val.lower() != 'nan':
                                context[col] = val[:300]

                    in_md = (bp, sample_id) in md_samples

                    results.append({
                        'bioproject': bp,
                        'sample_id': sample_id,
                        'category': cat,
                        'matched_value': m.group(0),
                        'matched_columns': '|'.join(matched_cols[:5]),
                        'matched_full_values': '|'.join(matched_vals[:3]),
                        'file_source': fpath.name,
                        'in_metadata': in_md,
                        'context': json.dumps(context, ensure_ascii=False),
                    })

                    summary_by_category[cat] += 1
                    projects_with_controls.add(bp)
                    break  # one category per row

# Save results
detail_df = pd.DataFrame(results)
# Sort: strong categories first, then by bioproject
strong_order = {'extraction_blank': 0, 'blank': 1, 'ntc_no_template': 2,
                'negative_control': 3, 'mock_community': 4, 'reagent_water': 5}
detail_df['cat_order'] = detail_df['category'].map(strong_order)
detail_df = detail_df.sort_values(['cat_order', 'bioproject', 'sample_id']).drop(columns=['cat_order'])

detail_path = OUTDIR / 'control_sample_ids.tsv'
detail_df.to_csv(detail_path, sep='\t', index=False)
print(f'\nWrote {detail_path}: {len(detail_df)} control samples from {len(projects_with_controls)} BioProjects')

# Summary by category and project
summary_rows = []
for bp in sorted(projects_with_controls):
    bp_data = detail_df[detail_df['bioproject'] == bp]
    row = {'bioproject': bp, 'total': len(bp_data),
           'in_metadata': bp_data['in_metadata'].sum()}
    for cat, _ in CATEGORIES:
        row[cat] = (bp_data['category'] == cat).sum()
    summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
summary_path = OUTDIR / 'control_samples_by_category.tsv'
summary_df.to_csv(summary_path, sep='\t', index=False)

# Print overview
print(f"\n=== CONTROL SAMPLE AUDIT ===")
print(f"BioProjects with controls: {len(projects_with_controls)}")
print(f"Total control samples: {len(detail_df)}")
print(f"In metadata (should be 0): {detail_df['in_metadata'].sum()}")
print(f"\nBy category:")
for cat, _ in CATEGORIES:
    cnt = summary_by_category[cat]
    n_bp = summary_df[summary_df[cat] > 0].shape[0]
    strong_mark = ' **' if cat in STRONG_CATEGORIES else ''
    print(f"  {cat}: {cnt} samples across {n_bp} BioProjects{strong_mark}")

print(f"\nHigh-confidence controls (extraction_blank + blank + NTC):")
strong = detail_df[detail_df['category'].isin(STRONG_CATEGORIES)]
strong_bps = set(strong['bioproject'])
for bp in sorted(strong_bps):
    bp_strong = strong[strong['bioproject'] == bp]
    cats = bp_strong['category'].value_counts().to_dict()
    cat_str = ', '.join(f'{k}={v}' for k, v in cats.items())
    # Show a few example IDs
    examples = bp_strong['sample_id'].head(5).tolist()
    print(f"  {bp}: {len(bp_strong)} controls ({cat_str})")
    print(f"    Examples: {', '.join(examples)}")

# Check if control BIOM data exists
print(f"\n=== BIOM AVAILABILITY ===")
biom_root = Path('/beegfs/dongbiao/oral')
for bp in sorted(strong_bps):
    biom_path = biom_root / bp / 'results' / 'DADA2.biom'
    if biom_path.exists():
        print(f"  {bp}: BIOM exists ✓")
    else:
        print(f"  {bp}: BIOM NOT FOUND ✗")
