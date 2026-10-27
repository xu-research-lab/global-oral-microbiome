#!/usr/bin/env python3
"""
Scan all project SraRunTable.csv files to identify blank/negative control/NTC/mock samples.

Output:
- control_samples_detailed.tsv: per-sample annotation of all potential controls
- control_samples_summary.tsv: per-BioProject summary
"""
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import json
import re

PROJECT_DIR = Path('/beegfs/huangxiaochang/global_oral/global-oral-microbiome/data/projects')
METADATA_CSV = '/beegfs/huangxiaochang/global_oral/global-oral-microbiome/data/metadata-v2-merge-oral.csv'
OUTDIR = Path('/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/06_control_catalog/results')

def log(msg):
    print(f'[find-controls] {msg}', flush=True)

# ---- Keywords for identifying control samples ----
# These are in lowercase for case-insensitive matching
BLANK_KEYWORDS = [
    'blank', 'ntc', 'no template', 'no-template', 'template-free',
    'negative control', 'neg control', 'negctrl', 'neg-ctrl',
    'extraction control', 'extraction blank', 'extraction negative',
    'pcr control', 'pcr negative', 'pcr blank',
    'mock', 'mock community', 'mock-community',
    'buffer', 'buffer control', 'elution buffer',
    'water', 'water control', 'molecular grade water', 'pcr water',
    'kit control', 'kit blank', 'reagent control', 'reagent blank',
    'empty', 'unused', 'unused swab', 'empty tube',
    'control blank', 'process control', 'processing control',
    'field blank', 'collection blank', 'sampling blank',
    'dna extraction', 'extraction kit',
]

# For sample name matching - more specific patterns
CONTROL_SAMPLE_NAME_PATTERNS = [
    re.compile(p, re.IGNORECASE) for p in [
        r'\bblank\b', r'\bNTC\b', r'\bntc\b',
        r'\bmock\b', r'\bmock[-\s]?community\b',
        r'\bnegative[-\s]?control\b', r'\bneg[-\s]?ctrl\b',
        r'\bextraction[-\s]?(blank|control|negative)\b',
        r'\bPCR[-\s]?(blank|control|negative|water)\b',
        r'\bbuffer\b', r'\bwater[-\s]?(control|blank)?\b',
        r'\bkit[-\s]?(control|blank|extract)\b',
        r'\bempty\b', r'\bunused\b',
        r'\btemplate[-\s]?free\b',
        r'\bcontrol[-\s]?(blank|sample)\b',
    ]
]

# Fields to search in SraRunTable (varies by project)
SAMPLE_NAME_FIELDS = [
    'Sample Name', 'Sample_Name', 'Sample_name', 'sample_name',
    'sample_id', 'SampleId', 'SampleName',
    'Submitter_Id', 'submitter_id',
    'source_material_id', 'Source_material_id',
    'Library Name', 'library_name',
    'isolate', 'Isolate',
    'title', 'Title',
]

DESCRIPTION_FIELDS = [
    'isolation_source', 'Isolation_source', 'IsolationSource',
    'material_type', 'material_type (exp)',
    'miscellaneous_parameter', 'study_design',
    'description', 'Description', 'SampleDescription',
    'sample_type', 'SampleType',
]

HOST_FIELDS = [
    'host_sex', 'Host_sex', 'HOST', 'host', 'Host',
    'submitted_subject_id', 'subject_id',
]

LIBRARY_FIELDS = [
    'LibrarySource', 'LibrarySelection', 'library_strategy (exp)',
]


def scan_sra_table(project_path: Path) -> dict:
    """Scan one project's SraRunTable.csv for control samples."""
    sra_file = project_path / 'SraRunTable.csv'
    if not sra_file.exists():
        return {'status': 'no_sra_table'}

    try:
        df = pd.read_csv(sra_file, dtype=str, low_memory=False)
    except Exception as e:
        return {'status': f'read_error: {e}'}

    results = []

    for idx, row in df.iterrows():
        sample_id = row.get('Run', f'row_{idx}')
        if pd.isna(sample_id):
            sample_id = f'row_{idx}'
        sample_id = str(sample_id)

        evidence = []
        evidence_details = {}

        # 1. Check Sample Name fields for blank/NTC keywords
        for field in SAMPLE_NAME_FIELDS:
            val = str(row.get(field, '')).strip()
            if val and val.lower() != 'nan':
                for pattern in CONTROL_SAMPLE_NAME_PATTERNS:
                    if pattern.search(val):
                        evidence.append(f'sample_name:{field}')
                        evidence_details[f'sample_name_{field}'] = val
                        break

        # 2. Check description/isolation fields
        for field in DESCRIPTION_FIELDS:
            val = str(row.get(field, '')).strip()
            if val and val.lower() != 'nan':
                val_lower = val.lower()
                for kw in BLANK_KEYWORDS:
                    if kw in val_lower:
                        evidence.append(f'description:{field}:{kw}')
                        evidence_details[f'desc_{field}'] = val
                        break

        # 3. Check host fields for missing/non-standard values
        # (weak signal - blanks have no host)
        host_values = []
        for field in HOST_FIELDS:
            val = str(row.get(field, '')).strip()
            if val and val.lower() not in ('nan', ''):
                host_values.append(val)

        # 4. Check LibrarySource for non-standard values
        for field in LIBRARY_FIELDS:
            val = str(row.get(field, '')).strip()
            if val and val.lower() != 'nan':
                # Genomic DNA / Metagenomic are normal
                if val.lower() not in ('genomic dna', 'metagenomic', 'genomic', 'pcr', 'other'):
                    if any(kw in val.lower() for kw in BLANK_KEYWORDS):
                        evidence.append(f'library:{field}:{val}')
                        evidence_details[f'lib_{field}'] = val

        if evidence:
            results.append({
                'bioproject': project_path.name,
                'sample_id': sample_id,
                'evidence': '|'.join(evidence),
                'evidence_details': json.dumps(evidence_details),
                'host_values': ';'.join(host_values) if host_values else 'NA',
                'in_metadata': 'unknown',  # will be filled later
            })

    bp = project_path.name
    return {
        'status': 'ok',
        'bioproject': bp,
        'total_sra_samples': len(df),
        'control_candidates': len(results),
        'samples': results,
    }


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    log('loading metadata for cross-reference')
    md = pd.read_csv(METADATA_CSV, dtype=str, low_memory=False)
    md_samples = set(zip(md['bioproject'].astype(str), md['sample_id'].astype(str)))
    md_bps = set(md['bioproject'].astype(str))
    log(f'{len(md):,} metadata samples across {len(md_bps)} BioProjects')

    log(f'scanning {len(list(PROJECT_DIR.iterdir()))} project directories')

    all_controls = []
    project_summaries = []

    for proj_dir in sorted(PROJECT_DIR.iterdir()):
        if not proj_dir.is_dir():
            continue
        bp = proj_dir.name
        result = scan_sra_table(proj_dir)

        if result['status'] == 'no_sra_table':
            continue
        elif result['status'].startswith('read_error'):
            log(f'  {bp}: {result["status"]}')
            continue

        if result['control_candidates'] > 0:
            log(f'  {bp}: {result["control_candidates"]} control candidates / {result["total_sra_samples"]} total')
            for s in result['samples']:
                s['in_metadata'] = (s['bioproject'], s['sample_id']) in md_samples
                all_controls.append(s)

            project_summaries.append({
                'bioproject': bp,
                'total_sra_samples': result['total_sra_samples'],
                'control_candidates': result['control_candidates'],
            })
        else:
            project_summaries.append({
                'bioproject': bp,
                'total_sra_samples': result['total_sra_samples'],
                'control_candidates': 0,
            })

    # Save detailed output
    detail_df = pd.DataFrame(all_controls)
    detail_path = OUTDIR / 'control_samples_detailed.tsv'
    detail_df.to_csv(detail_path, sep='\t', index=False)
    log(f'wrote {detail_path} ({len(detail_df)} control candidates)')

    # Summary
    summary_df = pd.DataFrame(project_summaries)
    summary_path = OUTDIR / 'control_samples_summary.tsv'
    summary_df.to_csv(summary_path, sep='\t', index=False)

    # Print overview
    projects_with_controls = summary_df[summary_df['control_candidates'] > 0]
    total_controls = summary_df['control_candidates'].sum()
    log(f'\n=== SUMMARY ===')
    log(f'Projects with control candidates: {len(projects_with_controls)} / {len(summary_df)}')
    log(f'Total control candidate samples: {total_controls}')
    log(f'In metadata: {detail_df["in_metadata"].sum() if len(detail_df) > 0 else 0}')

    for _, row in projects_with_controls.iterrows():
        bp_controls = [s for s in all_controls if s['bioproject'] == row['bioproject']]
        in_md = sum(1 for s in bp_controls if s['in_metadata'])
        log(f'  {row["bioproject"]}: {row["control_candidates"]} found, {in_md} in metadata')
        for s in bp_controls[:5]:
            log(f'    sample={s["sample_id"]}, evidence={s["evidence"][:120]}')


if __name__ == '__main__':
    main()
