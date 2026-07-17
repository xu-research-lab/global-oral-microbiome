#!/usr/bin/env python3
"""
Generate QC pipeline report as Markdown with embedded figures.
"""
import csv
import os
from collections import defaultdict, Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

OUT_DIR = '/beegfs/huangxiaochang/global_oral/global-oral-microbiome/data/metadata-and-biom-v2-qc'
RAW = f'{OUT_DIR}/metadata_all_for_qc_raw.tsv'
ALL_QC = f'{OUT_DIR}/metadata_all_qc.tsv'
ORAL_QC = f'{OUT_DIR}/metadata_oral_qc.tsv'
ORAL_DEDUP = f'{OUT_DIR}/metadata_oral_qc_dedup.tsv'
FUNNEL = f'{OUT_DIR}/QC/filter_funnel_summary.tsv'
REMOVED = f'{OUT_DIR}/QC/removed_studies_by_step.tsv'
OUT_MD = f'{OUT_DIR}/data_qc_pipeline.md'
FIGS_DIR = f'{OUT_DIR}/figures'
CR_THRESHOLD = 0.70  # must match merge_all_for_qc_raw.py

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 11,
    'axes.titlesize': 14, 'axes.labelsize': 12,
    'figure.dpi': 150, 'savefig.dpi': 150,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.1,
})


def save_fig(fig, name):
    os.makedirs(FIGS_DIR, exist_ok=True)
    path = f'{FIGS_DIR}/{name}'
    fig.savefig(path)
    plt.close(fig)
    return f'figures/{name}'


def load_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


print("Loading data...")
raw_rows = load_tsv(RAW)
qc_rows = load_tsv(ALL_QC)
oral_rows = load_tsv(ORAL_QC)
dedup_rows = load_tsv(ORAL_DEDUP)
funnel = load_tsv(FUNNEL)
removed_studies = load_tsv(REMOVED)

n_raw = len(raw_rows)
n_qc = len(qc_rows)
n_oral = len(oral_rows)
n_dedup = len(dedup_rows)
n_studies_raw = len(set(r['bioproject_id'] for r in raw_rows))
n_studies_qc = len(set(r['bioproject_id'] for r in qc_rows))
n_subjects = len(set(r['subject_id'] for r in dedup_rows))
n_oral_raw = sum(1 for r in raw_rows if r.get('site_group') == 'oral')

n_chim_studies = sum(1 for r in removed_studies if 'chimera' in r['step_name'].lower())
n_cr_studies = sum(1 for r in removed_studies if 'close' in r['step_name'].lower())
removed_chim_samples = sum(int(r['samples_total']) for r in removed_studies if 'chimera' in r['step_name'].lower())
removed_cr_samples = sum(int(r['samples_total']) for r in removed_studies if 'close' in r['step_name'].lower())

# ── Plot 1: Sample Flow Waterfall ──
print("Plot 1: Sample flow...")
fig, ax = plt.subplots(figsize=(12, 6))
labels = ['Raw\n(DADA2)', 'All\nQC', 'Oral\nQC', 'Oral\nDedup']
vals = [n_raw, n_qc, n_oral, n_dedup]
colors = ['#2196F3', '#4CAF50', '#00BCD4', '#795548']
bars = ax.bar(labels, vals, color=colors, width=0.5, edgecolor='white')
for bar, v in zip(bars, vals):
    ax.text(bar.get_x() + bar.get_width()/2, v + 800, f'{v:,}', ha='center', fontsize=12, fontweight='bold')
ax.set_ylabel('Number of Samples')
ax.set_title(f'Sample Count Progression (CR threshold = {CR_THRESHOLD*100:.0f}%)', fontweight='bold')
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x/1000:.0f}k'))
ax.grid(axis='y', alpha=0.3)
plot1 = save_fig(fig, 'sample_flow.png')

# ── Plot 2: Filter Impact ──
print("Plot 2: Filter impact...")
fig, ax = plt.subplots(figsize=(10, 6))
step_labels = ['Chimera\n(>30% oral chim>25%)', f'CR Mapping\n(oral <{CR_THRESHOLD*100:.0f}%)']
step_samples = [removed_chim_samples, removed_cr_samples]
step_counts = [n_chim_studies, n_cr_studies]
step_colors = ['#FF9800', '#9C27B0']
bars = ax.bar(step_labels, step_samples, color=step_colors, width=0.4, edgecolor='white')
for bar, samples, studies in zip(bars, step_samples, step_counts):
    ax.text(bar.get_x() + bar.get_width()/2, samples + 200,
            f'{studies} studies\n{samples:,} samples', ha='center', fontsize=12, fontweight='bold')
ax.set_ylabel('Total Samples Removed')
ax.set_title('Dataset-Level Filter Impact', fontweight='bold')
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x/1000:.0f}k'))
ax.grid(axis='y', alpha=0.3)
plot2 = save_fig(fig, 'filter_impact.png')

# ── Plot 3: CR Mapping Rate Distributions ──
print("Plot 3: CR mapping rates...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

sample_rates = []
for r in raw_rows:
    try: sample_rates.append(float(r.get('sample_cr_mapping_rate', 'NA')))
    except (ValueError, TypeError): pass

study_rates = {}
for r in raw_rows:
    s = r['bioproject_id']
    if s not in study_rates:
        try: study_rates[s] = float(r.get('study_cr_mapping_rate', 'NA'))
        except (ValueError, TypeError): pass

ax1.hist(sample_rates, bins=80, color='#2196F3', edgecolor='white', alpha=0.85)
ax1.axvline(x=CR_THRESHOLD, color='red', linestyle='--', linewidth=1.5, label=f'{CR_THRESHOLD*100:.0f}% (study filter)')
ax1.set_xlabel('Sample CR Mapping Rate')
ax1.set_ylabel('Samples')
ax1.set_title(f'Sample-Level CR Mapping (n={len(sample_rates):,})', fontweight='bold')
ax1.legend(fontsize=9)
ax1.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x/1000:.0f}k'))
ax1.grid(axis='y', alpha=0.3)

study_vals = list(study_rates.values())
ax2.hist(study_vals, bins=40, color='#4CAF50', edgecolor='white', alpha=0.85)
ax2.axvline(x=CR_THRESHOLD, color='red', linestyle='--', linewidth=1.5, label=f'{CR_THRESHOLD*100:.0f}% (filter)')
n_below = sum(1 for v in study_vals if v < CR_THRESHOLD)
ax2.set_xlabel('Study CR Mapping Rate (oral aggregate)')
ax2.set_ylabel('Studies')
ax2.set_title(f'Study-Level CR Mapping (n={len(study_vals)}, {n_below} < {CR_THRESHOLD*100:.0f}%)', fontweight='bold')
ax2.legend(fontsize=9)
ax2.grid(axis='y', alpha=0.3)

plt.tight_layout()
plot3 = save_fig(fig, 'cr_mapping_rates.png')

# ── Plot 4: Oral site distribution ──
print("Plot 4: Oral sites...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

oral_sites = Counter(r.get('site_cleaned', r.get('site', '')) for r in oral_rows)
d_sites = Counter(r.get('site_cleaned', r.get('site', '')) for r in dedup_rows)

for ax, sites, title in [(ax1, oral_sites, f'Oral QC (n={n_oral:,})'),
                            (ax2, d_sites, f'Oral Dedup (n={n_dedup:,}, {n_subjects:,} subjects)')]:
    top = sites.most_common(8)
    labels_v, vals_v = zip(*top) if top else ([], [])
    ax.barh(range(len(labels_v)), vals_v, color='#00BCD4', edgecolor='white')
    for i, v in enumerate(vals_v):
        ax.text(v + 30, i, f'{v:,}', va='center', fontsize=9, fontweight='bold')
    ax.set_yticks(range(len(labels_v)))
    ax.set_yticklabels(labels_v, fontsize=10)
    ax.set_title(title, fontweight='bold')
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plot4 = save_fig(fig, 'oral_sites.png')

# ── Removed studies table ──
removed_md = ''
for r in removed_studies:
    icon = '🟠' if 'chimera' in r['step_name'].lower() else '🟣'
    removed_md += f"| {icon} {r['step_name']} | {r['bioproject_id']} | {r['samples_total']} | {r['samples_oral']} | {r.get('study_cr_mapping_rate', 'NA')} |\n"

# ── Output file sizes ──
def size_mb(path):
    try: return os.path.getsize(path) / 1024**2
    except OSError: return 0

# ── Assemble Markdown ──
print("Assembling Markdown report...")
md = f"""# Global Oral Microbiome — Metadata QC Pipeline v4

**Date:** 2026-07-17
**Directory:** `data/metadata-and-biom-v2-qc/`
**Scripts:** `merge_all_for_qc_raw.py` → `generate_metadata_tables.py` → `generate_qza.py`

---

## 1. Key Statistics

| Metric | Value |
|--------|------:|
| Raw samples (DADA2-processed) | {n_raw:,} |
| Raw studies | {n_studies_raw} |
| **QC-passed samples** | **{n_qc:,}** |
| QC-passed studies | {n_studies_qc} |
| Studies removed | {n_studies_raw - n_studies_qc} |
| Oral QC samples | {n_oral:,} |
| Oral dedup samples | {n_dedup:,} |
| Unique subjects (dedup) | {n_subjects:,} |
| Oral samples in raw | {n_oral_raw:,} |

## 2. Filter Pipeline

Filtering is **dataset-level only**. All thresholds computed on **oral samples** (`site_group='oral'`).

| Step | Criterion | Datasets Removed | Samples Removed | Remaining Samples |
|------|-----------|:---:|:---:|---:|
| 0 — Raw merge | DADA2-processed samples only (`sample_cr_mapping_rate != NA`) | — | — | {n_raw:,} |
| 1 — Chimera filter | >30% oral samples with chimera >25% | {n_chim_studies} | {removed_chim_samples:,} | {int(funnel[1]['samples_after']):,} |
| 2 — **CR mapping filter** | **Oral aggregate GG2 CR mapping < {CR_THRESHOLD*100:.0f}%** | {n_cr_studies} | {removed_cr_samples:,} | **{n_qc:,}** |
| 3 — Oral subset | `site_group='oral'` | — | — | {n_oral:,} |
| 4 — Dedup | 1 sample per (subject_id × site_cleaned) | — | — | **{n_dedup:,}** |

> **v4 Changes:** Dataset-level only. ASV filter removed. CR threshold = {CR_THRESHOLD*100:.0f}%. No sample-level filters.

## 3. Sample Flow

![Sample Flow]({plot1})

## 4. Filter Impact

![Filter Impact]({plot2})

## 5. CR Mapping Rate Distributions

![CR Mapping Rates]({plot3})

*Left: sample-level CR mapping rates. Right: study-level oral aggregate CR mapping rates. Red dashed line = {CR_THRESHOLD*100:.0f}% filter threshold.*

## 6. Oral Site Distribution

![Oral Sites]({plot4})

## 7. Removed Studies Detail

| Step | Study | Total Samples | Oral Samples | CR Rate |
|------|-------|:---:|:---:|---:|
{removed_md}

## 8. Output Files

| File | Samples | Features | Studies | Subjects | Size |
|------|:---:|:---:|:---:|:---:|:---:|
| `metadata_all_for_qc_raw.tsv` | {n_raw:,} | — | {n_studies_raw} | — | {size_mb(RAW):.0f} MB |
| `metadata_all_qc.tsv` | {n_qc:,} | — | {n_studies_qc} | — | {size_mb(ALL_QC):.0f} MB |
| `metadata_oral_qc.tsv` | {n_oral:,} | — | {n_studies_qc} | — | {size_mb(ORAL_QC):.0f} MB |
| `metadata_oral_qc_dedup.tsv` | {n_dedup:,} | — | {n_studies_qc} | {n_subjects:,} | {size_mb(ORAL_DEDUP):.0f} MB |
| `biom_all_qc.biom` | 108,691 | 61,235 | — | — | {size_mb(f'{OUT_DIR}/biom_all_qc.biom'):.0f} MB |
| `biom_oral_qc.biom` | 80,788 | 61,235 | — | — | {size_mb(f'{OUT_DIR}/biom_oral_qc.biom'):.0f} MB |
| `biom_oral_qc_dedup.biom` | 46,223 | 61,235 | — | — | {size_mb(f'{OUT_DIR}/biom_oral_qc_dedup.biom'):.0f} MB |
| `biom_all_qc.qza` | — | — | — | — | {size_mb(f'{OUT_DIR}/biom_all_qc.qza'):.0f} MB |
| `genus_table_all_qc.qza` | — | — | — | — | {size_mb(f'{OUT_DIR}/genus_table_all_qc.qza'):.0f} MB |
| `genus_table_relabund_all_qc.qza` | — | — | — | — | {size_mb(f'{OUT_DIR}/genus_table_relabund_all_qc.qza'):.0f} MB |

### Auxiliary Files

| File | Description |
|------|-------------|
| `QC/filter_funnel_summary.tsv` | Per-step sample/study counts |
| `QC/removed_studies_by_step.tsv` | Studies removed at each step |
| `figures/` | Report figures (PNG) |
| `merge_all_for_qc_raw.py` | Merge metadata sources → compute QC flags |
| `generate_metadata_tables.py` | Dataset filters → metadata + BIOM tables |
| `generate_qc_report.py` | This report |
| `generate_qza.py` | BIOM → QZA + genus collapse (requires qiime2-2022.2) |
| `README.md` | Directory documentation |

---

*Generated 2026-07-17 by `generate_qc_report.py` | Filter: dataset-only, Chimera >25% in >30% oral, CR <{CR_THRESHOLD*100:.0f}% (oral aggregate)*
"""

with open(OUT_MD, 'w') as f:
    f.write(md)
print(f'Report: {OUT_MD} ({os.path.getsize(OUT_MD)/1024:.1f} KB)')
print(f'Figures: {FIGS_DIR}/')
print('Done!')
