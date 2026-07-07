#!/usr/bin/env python3
"""
Generate a comprehensive QC report HTML with embedded plots.

Reads metadata_all_for_qc_raw.tsv, computes all filter step details,
generates plots as base64-embedded PNGs, and writes data_qc_pipeline.html.
"""
import csv
import io
import base64
import os
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

BASE = '/beegfs/huangxiaochang/global_oral/new_analysis'
RAW = f'{BASE}/metadata_all_for_qc_raw.tsv'
OUT_HTML = f'{BASE}/data_qc_pipeline.html'

# ── matplotlib style ──────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

# ── helpers ────────────────────────────────────────────────────────
def fig_to_b64(fig):
    """Convert a matplotlib figure to a base64 PNG <img> tag."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return f'<img src="data:image/png;base64,{b64}" style="max-width:100%;margin:8px 0;">'

def load_rows(path):
    with open(path, 'r') as f:
        return list(csv.DictReader(f, delimiter='\t'))

# ── Load data ──────────────────────────────────────────────────────
print("Loading raw metadata...")
rows = load_rows(RAW)
print(f"  {len(rows)} samples")

# Group by study
study_data = defaultdict(list)
for r in rows:
    study_data[r['bioproject_id']].append(r)

# ── Step-by-step filtering ─────────────────────────────────────────
print("Computing filter steps...")

# Identify abnormal studies (study-level thresholds based on ORAL samples only)
abnormal_studies = set()
study_reasons = {}  # study -> [reasons]
for study, samples in study_data.items():
    # Use only oral samples for threshold computation
    oral_samples = [s for s in samples if s.get('site_group', '') == 'oral']
    n_oral = len(oral_samples)
    if n_oral == 0:
        continue
    reasons = []
    low_asv_frac = sum(1 for s in oral_samples if s.get('asv_le_30') == 'yes') / n_oral
    high_chim_frac = sum(1 for s in oral_samples if s.get('chim_gt_25') == 'yes') / n_oral
    low_cr = any(s.get('study_cr_mapping_low') == 'yes' for s in samples)
    if low_asv_frac > 0.5:
        reasons.append('asv')
    if high_chim_frac > 0.3:
        reasons.append('chim')
    if low_cr:
        reasons.append('cr')
    if reasons:
        abnormal_studies.add(study)
        study_reasons[study] = reasons

# Step 1: remove abnormal studies
step1 = [r for r in rows if r['bioproject_id'] not in abnormal_studies]
n_step0 = len(rows)
n_step1 = len(step1)
n_removed_study = n_step0 - n_step1
# Count oral in step0
oral_raw = sum(1 for r in rows if r.get('site_group', '') == 'oral')
oral_step1 = sum(1 for r in step1 if r.get('site_group', '') == 'oral')

# Step 2: remove asv <= 30
step2 = [r for r in step1 if r['asv_le_30'] != 'yes']
n_step2 = len(step2)
n_removed_asv = n_step1 - n_step2

# Step 3: remove sample cr mapping < 30%
step3 = [r for r in step2 if r['sample_cr_mapping_low'] != 'yes']
n_step3 = len(step3)
n_removed_cr_sample = n_step2 - n_step3

# Step 4: remove NA alpha / NA subject_id
step4 = [r for r in step3 if r['Observed'] != 'NA' and r['subject_id'] != 'NA']
n_step4 = len(step4)
n_removed_na = n_step3 - n_step4

# Oral subset
oral_rows = [r for r in step4 if r['site_group'] == 'oral']
n_oral = len(oral_rows)

# Dedup
groups = defaultdict(list)
for r in oral_rows:
    key = (r['subject_id'], r['site_cleaned'])
    groups[key].append(r)
dedup_rows = [v[0] for v in groups.values()]
n_dedup = len(dedup_rows)
n_removed_dup = n_oral - n_dedup

n_studies_raw = len(study_data)
n_studies_qc = len(set(r['bioproject_id'] for r in step4))
n_subjects_dedup = len(set(r['subject_id'] for r in dedup_rows))

print(f"  Step0 (raw):          {n_step0} samples, {n_studies_raw} studies")
print(f"  Step1 (-abnormal):    {n_step1} samples (-{n_removed_study}), {n_studies_raw - len(abnormal_studies)} studies")
print(f"  Step2 (-ASV<=30):     {n_step2} samples (-{n_removed_asv})")
print(f"  Step3 (-CR sample):   {n_step3} samples (-{n_removed_cr_sample})")
print(f"  Step4 (-NA/alpha):    {n_step4} samples (-{n_removed_na}), {n_studies_qc} studies")
print(f"  Oral:                 {n_oral} samples")
print(f"  Oral dedup:           {n_dedup} samples, {n_subjects_dedup} subjects")

# ── PLOT 1: Sample Flow Waterfall ──────────────────────────────────
print("Generating Plot 1: Sample flow waterfall...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# All sites waterfall
labels_all = ['Raw\nMerge', '-Abnormal\nStudies', '-ASV\n≤30', '-CR map\n<30%', '-NA alpha\n/subject_id']
values_all = [n_step0, n_step1, n_step2, n_step3, n_step4]
deltas_all = [0, n_step0 - n_step1, n_step1 - n_step2, n_step2 - n_step3, n_step3 - n_step4]
colors_all = ['#2196F3', '#FF9800', '#f44336', '#9C27B0', '#4CAF50']
# waterfall style: show as connected bars
for i in range(len(values_all)):
    ax1.bar(i, values_all[i], color=colors_all[i], width=0.6, edgecolor='white', linewidth=0.5)
    ax1.text(i, values_all[i] + 800, f'{values_all[i]:,}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    if deltas_all[i] > 0:
        ax1.text(i, values_all[i] / 2, f'-{deltas_all[i]:,}', ha='center', va='center', fontsize=8, color='white', fontweight='bold')
ax1.set_xticks(range(len(labels_all)))
ax1.set_xticklabels(labels_all, fontsize=9)
ax1.set_ylabel('Number of Samples', fontsize=12)
ax1.set_title('All Body Sites — Filter Cascade', fontsize=14, fontweight='bold')
ax1.set_ylim(0, max(values_all) * 1.12)
ax1.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x/1000:.0f}k'))
ax1.grid(axis='y', alpha=0.3)

# Oral subset waterfall
oral_after_study = oral_step1
# We don't have per-step oral counts without recomputing. Let's approximate from the all-sites numbers.
# Actually let's recompute properly
oral_step2_val = sum(1 for r in step2 if r.get('site_group', '') == 'oral')
oral_step3_val = sum(1 for r in step3 if r.get('site_group', '') == 'oral')
oral_step4_val = sum(1 for r in step4 if r.get('site_group', '') == 'oral')

oral_vals = [oral_raw, oral_after_study, oral_step2_val, oral_step3_val, oral_step4_val, n_oral, n_dedup]
oral_labels = ['Raw\nOral', '-Study\nFilter', '-ASV\n≤30', '-CR map\n<30%', '-NA\nalpha/subj', 'Oral\nSubset', '-Dedup']
oral_colors = ['#2196F3', '#FF9800', '#f44336', '#9C27B0', '#4CAF50', '#00BCD4', '#795548']

for i in range(len(oral_vals)):
    ax2.bar(i, oral_vals[i], color=oral_colors[i], width=0.6, edgecolor='white', linewidth=0.5)
    ax2.text(i, oral_vals[i] + 400, f'{oral_vals[i]:,}', ha='center', va='bottom', fontsize=9, fontweight='bold')
ax2.set_xticks(range(len(oral_labels)))
ax2.set_xticklabels(oral_labels, fontsize=8)
ax2.set_ylabel('Number of Samples', fontsize=12)
ax2.set_title('Oral Sites — Filter Cascade', fontsize=14, fontweight='bold')
ax2.set_ylim(0, max(oral_vals) * 1.12)
ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x/1000:.0f}k'))
ax2.grid(axis='y', alpha=0.3)

plt.tight_layout()
plot1_html = fig_to_b64(fig)

# ── PLOT 2: Study Filter Venn/Upset-style ─────────────────────────
print("Generating Plot 2: Study filter criteria overview...")
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(16, 6))

# Count studies by filter category
only_asv = set()
only_chim = set()
only_cr = set()
asv_chim = set()
asv_cr = set()
chim_cr = set()
all3 = set()

for study in abnormal_studies:
    reasons = set(study_reasons[study])
    if reasons == {'asv'}:
        only_asv.add(study)
    elif reasons == {'chim'}:
        only_chim.add(study)
    elif reasons == {'cr'}:
        only_cr.add(study)
    elif reasons == {'asv', 'chim'}:
        asv_chim.add(study)
    elif reasons == {'asv', 'cr'}:
        asv_cr.add(study)
    elif reasons == {'chim', 'cr'}:
        chim_cr.add(study)
    elif reasons == {'asv', 'chim', 'cr'}:
        all3.add(study)

categories = ['ASV only', 'Chimera\nonly', 'CR only', 'ASV+\nChimera', 'ASV+\nCR', 'Chimera+\nCR', 'All 3']
counts = [len(only_asv), len(only_chim), len(only_cr), len(asv_chim), len(asv_cr), len(chim_cr), len(all3)]
bar_colors = ['#f44336', '#FF9800', '#9C27B0', '#E91E63', '#673AB7', '#FF5722', '#333333']

bars = ax_a.barh(categories, counts, color=bar_colors, edgecolor='white')
for bar, count in zip(bars, counts):
    ax_a.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2, str(count),
              va='center', fontsize=11, fontweight='bold')
ax_a.set_xlabel('Number of Studies', fontsize=12)
ax_a.set_title(f'Abnormal Studies by Filter Criteria\n(Total: {len(abnormal_studies)} removed)', fontsize=13, fontweight='bold')
ax_a.set_xlim(0, max(counts) * 1.3)
ax_a.grid(axis='x', alpha=0.3)

# Sample count impact by criterion
# How many samples are removed due to each criterion
sample_impact = {}
for study in abnormal_studies:
    n_s = len(study_data[study])
    for reason in study_reasons[study]:
        sample_impact[reason] = sample_impact.get(reason, 0) + n_s

impact_labels = ['ASV >50%\nlow-ASV', 'Chimera >25%\nin >30% samples', 'CR mapping\n<30%']
impact_colors = ['#f44336', '#FF9800', '#9C27B0']
impact_values = [sample_impact.get('asv', 0), sample_impact.get('chim', 0), sample_impact.get('cr', 0)]

bars2 = ax_b.bar(impact_labels, impact_values, color=impact_colors, edgecolor='white', width=0.5)
for bar, val in zip(bars2, impact_values):
    ax_b.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 200, f'{val:,}',
              ha='center', fontsize=12, fontweight='bold')
ax_b.set_ylabel('Total Samples in Affected Studies', fontsize=12)
ax_b.set_title('Sample Impact by Filter Criterion\n(non-mutually-exclusive)', fontsize=13, fontweight='bold')
ax_b.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x/1000:.0f}k'))
ax_b.grid(axis='y', alpha=0.3)

plt.tight_layout()
plot2_html = fig_to_b64(fig)

# ── PLOT 3: CR Mapping Rate Distributions ──────────────────────────
print("Generating Plot 3: CR mapping rate distributions...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Sample-level CR mapping rates
sample_rates = []
for r in rows:
    val = r.get('sample_cr_mapping_rate', 'NA')
    if val != 'NA':
        try:
            sample_rates.append(float(val))
        except ValueError:
            pass

ax1.hist(sample_rates, bins=80, color='#2196F3', edgecolor='white', alpha=0.85, linewidth=0.3)
ax1.axvline(x=0.3, color='red', linestyle='--', linewidth=2, label='30% threshold')
n_below = sum(1 for v in sample_rates if v < 0.3)
ax1.text(0.15, ax1.get_ylim()[1] * 0.85 if ax1.get_ylim()[1] > 0 else 10000,
         f'{n_below:,} samples\n<30% ({100*n_below/len(sample_rates):.1f}%)',
         fontsize=11, color='red', fontweight='bold',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax1.set_xlabel('Sample CR Mapping Rate', fontsize=12)
ax1.set_ylabel('Number of Samples', fontsize=12)
ax1.set_title(f'Sample-Level CR Mapping Rate Distribution\n(n={len(sample_rates):,})', fontsize=13, fontweight='bold')
ax1.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x/1000:.0f}k'))
ax1.legend(fontsize=10)
ax1.grid(axis='y', alpha=0.3)

# Study-level CR mapping rates
study_rates = {}
for study, samples in study_data.items():
    for s in samples:
        val = s.get('study_cr_mapping_rate', 'NA')
        if val != 'NA':
            try:
                study_rates[study] = float(val)
            except ValueError:
                pass
            break

study_rate_vals = list(study_rates.values())
ax2.hist(study_rate_vals, bins=40, color='#4CAF50', edgecolor='white', alpha=0.85, linewidth=0.3)
ax2.axvline(x=0.3, color='red', linestyle='--', linewidth=2, label='30% threshold')
n_study_below = sum(1 for v in study_rate_vals if v < 0.3)
ax2.text(0.15, ax2.get_ylim()[1] * 0.85 if ax2.get_ylim()[1] > 0 else 10,
         f'{n_study_below} studies\n<30% ({100*n_study_below/len(study_rate_vals):.1f}%)',
         fontsize=11, color='red', fontweight='bold',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax2.set_xlabel('Study-Level CR Mapping Rate (aggregate)', fontsize=12)
ax2.set_ylabel('Number of Studies', fontsize=12)
ax2.set_title(f'Study-Level CR Mapping Rate Distribution\n(n={len(study_rate_vals)})', fontsize=13, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(axis='y', alpha=0.3)

plt.tight_layout()
plot3_html = fig_to_b64(fig)

# ── PLOT 4: Study CR mapping rate sorted bar (top/bottom 20) ──────
print("Generating Plot 4: Study CR mapping rate overview...")
fig, ax = plt.subplots(figsize=(18, 10))

sorted_studies = sorted(study_rates.items(), key=lambda x: x[1])
# Show all below 30% + bottom part
n_show = min(80, len(sorted_studies))
show_studies = sorted_studies[:n_show]  # lowest rates
study_names = [s[0] for s in show_studies]
study_vals_list = [s[1] for s in show_studies]

colors = ['#f44336' if v < 0.3 else '#4CAF50' for v in study_vals_list]
bars = ax.barh(range(len(study_names)), study_vals_list, color=colors, edgecolor='white', linewidth=0.3, height=0.75)
ax.axvline(x=0.3, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
ax.set_yticks(range(len(study_names)))
ax.set_yticklabels(study_names, fontsize=7, fontfamily='monospace')
ax.set_xlabel('Study-Level Aggregate CR Mapping Rate', fontsize=12)
ax.set_title(f'Study-Level CR Mapping Rate (Lowest {n_show} Studies)\nRed = below 30% threshold, Green = above',
             fontsize=14, fontweight='bold')
ax.set_xlim(0, max(study_vals_list) * 1.05)
ax.invert_yaxis()
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plot4_html = fig_to_b64(fig)

# ── PLOT 5: Oral site distribution ─────────────────────────────────
print("Generating Plot 5: Oral site_cleaned distribution...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

oral_sites = Counter(r['site_cleaned'] for r in oral_rows)
d_sites = Counter(r['site_cleaned'] for r in dedup_rows)

# Before dedup
top_n = 10
top_oral = oral_sites.most_common(top_n)
labels_o = [t[0] for t in top_oral]
vals_o = [t[1] for t in top_oral]
ax1.barh(range(len(labels_o)), vals_o, color='#00BCD4', edgecolor='white')
for i, v in enumerate(vals_o):
    ax1.text(v + 50, i, f'{v:,}', va='center', fontsize=9, fontweight='bold')
ax1.set_yticks(range(len(labels_o)))
ax1.set_yticklabels(labels_o, fontsize=10)
ax1.set_title(f'Oral QC — site_cleaned Distribution\n(n={n_oral:,} samples)', fontsize=13, fontweight='bold')
ax1.invert_yaxis()
ax1.grid(axis='x', alpha=0.3)

# After dedup
top_dedup = d_sites.most_common(top_n)
labels_d = [t[0] for t in top_dedup]
vals_d = [t[1] for t in top_dedup]
ax2.barh(range(len(labels_d)), vals_d, color='#795548', edgecolor='white')
for i, v in enumerate(vals_d):
    ax2.text(v + 20, i, f'{v:,}', va='center', fontsize=9, fontweight='bold')
ax2.set_yticks(range(len(labels_d)))
ax2.set_yticklabels(labels_d, fontsize=10)
ax2.set_title(f'Oral Dedup — site_cleaned Distribution\n(n={n_dedup:,} samples, {n_subjects_dedup:,} subjects)',
              fontsize=13, fontweight='bold')
ax2.invert_yaxis()
ax2.grid(axis='x', alpha=0.3)

plt.tight_layout()
plot5_html = fig_to_b64(fig)

# ── PLOT 6: Filter step summary donut ──────────────────────────────
print("Generating Plot 6: Summary donut charts...")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Donut 1: Study filter - which criteria flagged studies
criteria_counts = [len(abnormal_studies), n_studies_raw - len(abnormal_studies)]
wedges1, texts1, autotexts1 = axes[0].pie(
    criteria_counts, labels=['Removed', 'Kept'],
    colors=['#f44336', '#4CAF50'], autopct='%1.1f%%',
    startangle=90, pctdistance=0.6,
    wedgeprops=dict(width=0.35, edgecolor='white'))
for at in autotexts1:
    at.set_fontweight('bold')
axes[0].set_title(f'Study Filter\n({len(abnormal_studies)}/{n_studies_raw} removed)', fontsize=12, fontweight='bold')

# Donut 2: Sample removals at each step
sample_removal_labels = ['Kept\n(final)', 'Abnormal\nstudies', 'ASV\n≤30', 'CR map\n<30%', 'NA alpha\n/subject']
sample_removal_vals = [n_step4, n_removed_study, n_removed_asv, n_removed_cr_sample, n_removed_na]
sample_removal_colors = ['#4CAF50', '#FF9800', '#f44336', '#9C27B0', '#607D8B']
wedges2, texts2, autotexts2 = axes[1].pie(
    sample_removal_vals, labels=sample_removal_labels, colors=sample_removal_colors,
    autopct='%1.1f%%', startangle=90, pctdistance=0.75,
    wedgeprops=dict(width=0.35, edgecolor='white'))
for i, at in enumerate(autotexts2):
    at.set_fontsize(9)
    at.set_fontweight('bold')
axes[1].set_title(f'Sample Filter Breakdown\n(Total: {n_step0:,})', fontsize=12, fontweight='bold')

# Donut 3: Oral vs non-oral in final QC
non_oral_qc = n_step4 - n_oral
oral_dedup_vals = [n_dedup, n_oral - n_dedup, non_oral_qc]
oral_dedup_labels = ['Oral\n(dedup)', 'Oral\n(duplicates)', 'Non-oral']
oral_dedup_colors = ['#00BCD4', '#B2EBF2', '#9E9E9E']
wedges3, texts3, autotexts3 = axes[2].pie(
    oral_dedup_vals, labels=oral_dedup_labels, colors=oral_dedup_colors,
    autopct='%1.1f%%', startangle=90, pctdistance=0.75,
    wedgeprops=dict(width=0.35, edgecolor='white'))
for at in autotexts3:
    at.set_fontweight('bold')
axes[2].set_title(f'Final QC Breakdown\n(Total: {n_step4:,})', fontsize=12, fontweight='bold')

plt.tight_layout()
plot6_html = fig_to_b64(fig)

# ── Embed R-generated QC panel plots ─────────────────────────────────
print("Loading R-generated QC panel plots...")

def png_to_b64_html(path, max_width='100%'):
    """Read a PNG file and return a base64 <img> tag."""
    with open(path, 'rb') as f:
        b64 = base64.b64encode(f.read()).decode()
    return f'<img src="data:image/png;base64,{b64}" style="max-width:{max_width};margin:8px 0;">'

r_study_panel = png_to_b64_html(
    f'{BASE}/QC/metadata_raw_qc_panels_v3_study_level.png')
r_sample_panel = png_to_b64_html(
    f'{BASE}/QC/metadata_raw_qc_panels_v3_sample_level.png')
r_sample_post_filter = png_to_b64_html(
    f'{BASE}/QC/sample_level_qc_panels_v3.png')

# ── Tables: Abnormal studies detail ─────────────────────────────────
print("Generating abnormal studies detail table...")
# Build detailed table rows for all abnormal studies
abnormal_table_rows = []
for study in sorted(abnormal_studies):
    samples = study_data[study]
    n = len(samples)
    reasons = study_reasons[study]
    low_asv_frac = sum(1 for s in samples if s.get('asv_le_30') == 'yes') / n
    high_chim_frac = sum(1 for s in samples if s.get('chim_gt_25') == 'yes') / n
    cr_rate = study_rates.get(study, None)
    cr_str = f'{cr_rate*100:.1f}%' if cr_rate is not None else 'NA'

    reason_strs = []
    if 'asv' in reasons:
        reason_strs.append(f'ASV≤30 in {low_asv_frac*100:.0f}%')
    if 'chim' in reasons:
        reason_strs.append(f'Chim>25% in {high_chim_frac*100:.0f}%')
    if 'cr' in reasons:
        reason_strs.append(f'CR map = {cr_str}')

    tr_class = ''
    if len(reasons) >= 2:
        tr_class = ' style="background:#fff3e0;"'

    abnormal_table_rows.append(
        f'<tr{tr_class}><td>{study}</td><td>{n}</td><td>{" + ".join(reason_strs)}</td></tr>'
    )

# ── Assemble HTML ──────────────────────────────────────────────────
print("Assembling HTML report...")

# Count QC flag statistics
n_asv_le_30 = sum(1 for r in rows if r['asv_le_30'] == 'yes')
n_chim_gt_25 = sum(1 for r in rows if r['chim_gt_25'] == 'yes')
n_sample_cr_low = sum(1 for r in rows if r['sample_cr_mapping_low'] == 'yes')
n_study_asv = sum(1 for study in abnormal_studies if 'asv' in study_reasons[study])
n_study_chim = sum(1 for study in abnormal_studies if 'chim' in study_reasons[study])
n_study_cr = sum(1 for study in abnormal_studies if 'cr' in study_reasons[study])

html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Metadata QC Pipeline Documentation (v3)</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; max-width: 1100px; margin: 40px auto; padding: 0 20px; color: #222; line-height: 1.6; background: #fafafa; }}
h1 {{ border-bottom: 3px solid #333; padding-bottom: 10px; color: #111; }}
h2 {{ margin-top: 36px; border-bottom: 2px solid #999; padding-bottom: 6px; color: #333; }}
h3 {{ margin-top: 28px; color: #555; }}
table {{ border-collapse: collapse; width: 100%; margin: 12px 0; background: white; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: left; }}
th {{ background: #f5f5f5; font-weight: 600; }}
tr:nth-child(even) {{ background: #fafafa; }}
code {{ background: #f0f0f0; padding: 1px 5px; border-radius: 3px; font-size: 0.92em; }}
pre {{ background: #f5f5f5; padding: 12px 16px; border-radius: 4px; overflow-x: auto; border: 1px solid #e0e0e0; }}
pre code {{ background: none; padding: 0; }}
.flow {{ background: #e8f5e9; border-left: 4px solid #4caf50; padding: 12px 16px; margin: 16px 0; border-radius: 2px; }}
.note {{ background: #fff3e0; border-left: 4px solid #ff9800; padding: 12px 16px; margin: 16px 0; border-radius: 2px; }}
.stats-box {{ display: flex; flex-wrap: wrap; gap: 12px; margin: 16px 0; }}
.stat-card {{ background: white; border: 1px solid #e0e0e0; border-radius: 8px; padding: 16px 20px; flex: 1; min-width: 140px; text-align: center; box-shadow: 0 2px 4px rgba(0,0,0,0.06); }}
.stat-card .number {{ font-size: 28px; font-weight: 700; color: #1565C0; }}
.stat-card .label {{ font-size: 12px; color: #888; margin-top: 4px; text-transform: uppercase; letter-spacing: 0.5px; }}
.stat-card.red .number {{ color: #c62828; }}
.stat-card.green .number {{ color: #2E7D32; }}
.stat-card.purple .number {{ color: #6A1B9A; }}
.figure-box {{ background: white; border: 1px solid #e0e0e0; border-radius: 8px; padding: 16px; margin: 20px 0; box-shadow: 0 2px 6px rgba(0,0,0,0.06); }}
.figure-box .caption {{ font-size: 12px; color: #999; margin-top: 8px; text-align: center; }}
.filter-step {{ background: white; border-left: 4px solid #2196F3; padding: 10px 16px; margin: 10px 0; border-radius: 2px; box-shadow: 0 1px 2px rgba(0,0,0,0.05); }}
.filter-step.removed {{ border-left-color: #f44336; background: #fff5f5; }}
.filter-step.final {{ border-left-color: #4CAF50; background: #f1f8e9; }}
details {{ margin: 12px 0; }}
details summary {{ cursor: pointer; font-weight: 600; color: #1565C0; font-size: 14px; padding: 4px 0; }}
details summary:hover {{ color: #0D47A1; }}
</style>
</head>
<body>

<h1>🦷 Global Oral Microbiome — Metadata QC Pipeline v3</h1>
<p>Date: 2026-06-26 &nbsp;|&nbsp; Directory: <code>{BASE}</code> &nbsp;|&nbsp;
Scripts: <code>merge_all_for_qc_raw.py</code>, <code>generate_metadata_tables.py</code></p>

<div class="note">
<strong>v3 Updates (2026-06-26):</strong>
<ul>
  <li>CR data source: switched from study-level <code>study_close_reference_keeping_rate_cleaned_all.csv</code>
      to sample-level <code>all_samples_gg2_closed_reference_mapping_fraction.tsv</code></li>
  <li>CR filter threshold: <strong>50% → 30%</strong></li>
  <li>New sample-level CR mapping filter: removes individual samples with mapping rate &lt;30%</li>
  <li>New Source BIOM: <code>gg2_reference_id_profile.biom</code></li>
</ul>
</div>

<!-- ── KEY STATISTICS CARDS ── -->
<h2>1. Key Statistics at a Glance</h2>
<div class="stats-box">
  <div class="stat-card">
    <div class="number">{n_step0:,}</div>
    <div class="label">Raw Samples</div>
  </div>
  <div class="stat-card">
    <div class="number">{n_studies_raw}</div>
    <div class="label">Raw Studies</div>
  </div>
  <div class="stat-card red">
    <div class="number">{len(abnormal_studies)}</div>
    <div class="label">Studies Removed</div>
  </div>
  <div class="stat-card green">
    <div class="number">{n_step4:,}</div>
    <div class="label">Final QC Samples</div>
  </div>
  <div class="stat-card">
    <div class="number">{n_oral:,}</div>
    <div class="label">Oral QC Samples</div>
  </div>
  <div class="stat-card purple">
    <div class="number">{n_dedup:,}</div>
    <div class="label">Oral Dedup Samples</div>
  </div>
  <div class="stat-card">
    <div class="number">{n_subjects_dedup:,}</div>
    <div class="label">Unique Subjects</div>
  </div>
</div>

<!-- ── PLOT 1: SAMPLE FLOW ── -->
<h2>2. Sample Filter Cascade</h2>
<div class="figure-box">
{plot1_html}
<div class="caption">Figure 1: Sample counts at each filter step for all body sites (left) and oral-only sites (right).</div>
</div>

<div style="display: flex; flex-wrap: wrap; gap: 12px;">

<div class="filter-step" style="flex: 1; min-width: 200px;">
  <strong>Step 1 — Merge:</strong><br>
  <code>merge_all_for_qc_raw.py</code> joins 6 input files<br>
  → <strong>{n_step0:,} samples</strong>, {n_studies_raw} studies, {oral_raw:,} oral
</div>

<div class="filter-step removed" style="flex: 1; min-width: 200px;">
  <strong>Step 2 — Study Filter (OR):</strong><br>
  ASV≤30 in &gt;50% OR chim&gt;25% in &gt;30% OR CR map &lt;30%<br>
  → <strong>{len(abnormal_studies)} studies removed</strong> (−{n_removed_study:,} samples)
</div>

<div class="filter-step removed" style="flex: 1; min-width: 200px;">
  <strong>Step 3 — ASV ≤ 30:</strong><br>
  Remove individual samples with asv_count ≤ 30<br>
  → <strong>−{n_removed_asv:,} samples</strong>
</div>

<div class="filter-step removed" style="flex: 1; min-width: 200px;">
  <strong>Step 4 — CR Mapping &lt;30% <em>(NEW v3)</em>:</strong><br>
  Remove individual samples with CR mapping rate &lt;30%<br>
  → <strong>−{n_removed_cr_sample:,} samples</strong>
</div>

<div class="filter-step removed" style="flex: 1; min-width: 200px;">
  <strong>Step 5 — Completeness:</strong><br>
  Remove samples with NA alpha diversity (Observed) or NA subject_id<br>
  → <strong>−{n_removed_na:,} samples</strong>
</div>

<div class="filter-step final" style="flex: 1; min-width: 200px;">
  <strong>Step 6 — Oral Subset:</strong><br>
  Select <code>site_group='oral'</code><br>
  → <strong>{n_oral:,} oral samples</strong>
</div>

<div class="filter-step final" style="flex: 1; min-width: 200px;">
  <strong>Step 7 — Dedup:</strong><br>
  Keep 1 sample per (subject_id × site_cleaned)<br>
  → <strong>{n_dedup:,} samples</strong>, {n_subjects_dedup:,} subjects (−{n_removed_dup:,} duplicates)
</div>

</div>

<!-- ── PLOT 2: STUDY FILTER CRITERIA ── -->
<h2>3. Study-Level Filter Detail</h2>
<div class="figure-box">
{plot2_html}
<div class="caption">Figure 2: Left: Number of abnormal studies by filter criterion combination. Right: Total sample impact per criterion (non-mutually-exclusive, a study may be counted in multiple criteria).</div>
</div>

<h3>3a. All {len(abnormal_studies)} Removed Studies — Detail</h3>
<table style="font-size:13px;">
<tr><th>Study</th><th>Samples</th><th>Reason(s) for Removal</th></tr>
{"".join(abnormal_table_rows)}
</table>

<details>
<summary>Show study reasons breakdown (click to expand)</summary>
<table style="font-size:13px;">
<tr><th>Category</th><th>Count</th><th>Studies</th></tr>
<tr><td>ASV only</td><td>{len(only_asv)}</td><td>{", ".join(sorted(only_asv)) if only_asv else "—"}</td></tr>
<tr><td>Chimera only</td><td>{len(only_chim)}</td><td>{", ".join(sorted(only_chim)) if only_chim else "—"}</td></tr>
<tr><td>CR only</td><td>{len(only_cr)}</td><td>{", ".join(sorted(only_cr)) if only_cr else "—"}</td></tr>
<tr style="background:#fff3e0;"><td>ASV + Chimera</td><td>{len(asv_chim)}</td><td>{", ".join(sorted(asv_chim)) if asv_chim else "—"}</td></tr>
<tr style="background:#fff3e0;"><td>ASV + CR</td><td>{len(asv_cr)}</td><td>{", ".join(sorted(asv_cr)) if asv_cr else "—"}</td></tr>
<tr style="background:#fff3e0;"><td>Chimera + CR</td><td>{len(chim_cr)}</td><td>{", ".join(sorted(chim_cr)) if chim_cr else "—"}</td></tr>
<tr style="background:#ffccbc;"><td>All 3 criteria</td><td>{len(all3)}</td><td>{", ".join(sorted(all3)) if all3 else "—"}</td></tr>
</table>
</details>

<!-- ── R PLOTS: STUDY/SAMPLE QC PANELS ── -->
<h2>4. R-Generated QC Panels — Study-Level Overview</h2>
<div class="figure-box">
{r_study_panel}
<div class="caption">Figure 7 (R): Study-level QC panels. From top: (1) Mean raw reads per study (log10); (2) Ratio of samples with ASV &le; 30 (red = flagged &gt;50%); (3) Closed-reference mapping rate (grey line = 30% threshold); (4) Ratio of samples with chimera &gt; 25% (red = flagged &gt;30%); (5) Sample counts per study (total samples + unique subjects). Studies ordered by sample count descending.</div>
</div>

<h2>5. R-Generated QC Panels — Sample-Level (All Raw Data)</h2>
<div class="figure-box">
{r_sample_panel}
<div class="caption">Figure 8 (R): Sample-level QC panels for all {n_step0:,} samples in the raw merge. From top: (1) Raw reads per sample, colored by study; (2) Chimera percentage (red line = 25%); (3) ASV count per sample (red line = 30). Study boundaries shown as grey vertical lines. Samples ordered by study (descending sample count), then by decreasing raw reads.</div>
</div>

<h2>6. R-Generated QC Panels — Sample-Level (After Study Filter)</h2>
<div class="figure-box">
{r_sample_post_filter}
<div class="caption">Figure 9 (R): Sample-level QC panels for {n_step1:,} samples after removing {len(abnormal_studies)} abnormal studies. From top: (1) Raw reads (bars) + dechimera reads overlay (black curve), red dashed line = 1,000 reads; (2) Chimera percentage (red line = 25%); (3) Observed ASV count (red line = 30). Samples ordered by study rank, then by descending raw reads.</div>
</div>

<!-- ── PLOTS 3&4: CR MAPPING RATE ── -->
<h2>7. Closed-Reference Mapping Rate Distributions</h2>
<div class="figure-box">
{plot3_html}
<div class="caption">Figure 12: Distribution of CR mapping rates. Left: per-sample level (n={len(sample_rates):,}, {n_below:,} below 30%). Right: per-study aggregate (n={len(study_rate_vals)}, {n_study_below} below 30%). Red dashed line = 30% threshold.</div>
</div>

<div class="figure-box">
{plot4_html}
<div class="caption">Figure 13: Study-level aggregate CR mapping rates, sorted low-to-high (bottom {n_show} shown). Red = below 30%, Green = above 30%. The aggregate is computed as sum(closed_reference_reads) / sum(total_reads) across all samples in each study.</div>
</div>

<!-- ── ORAL SITES ── -->
<h2>8. Oral Site Distribution</h2>
<div class="figure-box">
{plot5_html}
<div class="caption">Figure 10: Distribution of <code>site_cleaned</code> categories in the oral QC cohort. Left: before deduplication ({n_oral:,} samples). Right: after deduplication by (subject_id, site_cleaned) ({n_dedup:,} samples, {n_subjects_dedup:,} unique subjects).</div>
</div>

<!-- ── SUMMARY DONUTS ── -->
<h2>9. Summary Overview</h2>
<div class="figure-box">
{plot6_html}
<div class="caption">Figure 11: Donut charts summarizing the pipeline. Left: study filter outcome. Middle: sample removal breakdown by filter step. Right: final QC composition (oral dedup vs oral duplicates vs non-oral).</div>
</div>

<!-- ── FILTER SUMMARY TABLE ── -->
<h2>10. Sample Count Change Summary</h2>
<table>
<tr><th>Stage</th><th>All Sites</th><th>Δ</th><th>Oral</th><th>Δ</th></tr>
<tr><td>Raw merge</td><td style="text-align:right">{n_step0:,}</td><td>—</td><td style="text-align:right">{oral_raw:,}</td><td>—</td></tr>
<tr><td>− Abnormal studies ({len(abnormal_studies)})</td><td style="text-align:right">{n_step1:,}</td><td style="color:#c62828;">−{n_removed_study:,}</td><td style="text-align:right">{oral_step1:,}</td><td style="color:#c62828;">−{oral_raw - oral_step1:,}</td></tr>
<tr><td>− ASV ≤ 30</td><td style="text-align:right">{n_step2:,}</td><td style="color:#c62828;">−{n_removed_asv:,}</td><td style="text-align:right">{oral_step2_val:,}</td><td style="color:#c62828;">−{oral_step1 - oral_step2_val:,}</td></tr>
<tr><td>− Sample CR mapping &lt;30% <em>(NEW v3)</em></td><td style="text-align:right">{n_step3:,}</td><td style="color:#c62828;">−{n_removed_cr_sample:,}</td><td style="text-align:right">{oral_step3_val:,}</td><td style="color:#c62828;">−{oral_step2_val - oral_step3_val:,}</td></tr>
<tr><td>− NA alpha / subject_id</td><td style="text-align:right"><strong>{n_step4:,}</strong></td><td style="color:#c62828;">−{n_removed_na:,}</td><td style="text-align:right">{oral_step4_val:,}</td><td style="color:#c62828;">−{oral_step3_val - oral_step4_val:,}</td></tr>
<tr><td>Oral subset / Dedup</td><td>—</td><td>—</td><td style="text-align:right"><strong>{n_dedup:,}</strong></td><td style="color:#c62828;">−{n_removed_dup:,}</td></tr>
</table>

<!-- ── QC FLAG COLUMNS ── -->
<h2>11. QC Flag Columns in metadata_all_for_qc_raw.tsv</h2>
<table>
<tr><th>Column</th><th>Meaning</th><th>Values</th><th>Samples Flagged</th></tr>
<tr><td><code>asv_le_30</code></td><td>Sample has ≤ 30 ASVs</td><td>yes / no / NA</td><td style="text-align:right">{n_asv_le_30:,} ({100*n_asv_le_30/n_step0:.1f}%)</td></tr>
<tr><td><code>study_asv_lt30_high</code></td><td>Study has &gt;50% samples with ASV ≤ 30</td><td>yes / no / NA</td><td style="text-align:right">{n_study_asv} studies</td></tr>
<tr><td><code>chim_gt_25</code></td><td>Sample has &gt;25% chimera reads</td><td>yes / no / NA</td><td style="text-align:right">{n_chim_gt_25:,} ({100*n_chim_gt_25/n_step0:.1f}%)</td></tr>
<tr><td><code>study_chim_gt25_high</code></td><td>Study has &gt;30% samples with chimera &gt;25%</td><td>yes / no / NA</td><td style="text-align:right">{n_study_chim} studies</td></tr>
<tr><td><code>sample_cr_mapping_rate</code></td><td>Sample-level CR mapping fraction</td><td>0–1 float / NA</td><td>—</td></tr>
<tr><td><code>sample_cr_mapping_low</code> <em>(NEW v3)</em></td><td>Sample CR mapping rate &lt;30%</td><td>yes / no / NA</td><td style="text-align:right">{n_sample_cr_low:,} ({100*n_sample_cr_low/n_step0:.1f}%)</td></tr>
<tr><td><code>study_cr_mapping_rate</code></td><td>Study-level aggregate CR mapping rate</td><td>0–1 float / NA</td><td>—</td></tr>
<tr><td><code>study_cr_mapping_low</code></td><td>Study aggregate CR mapping rate &lt;30%</td><td>yes / no / NA</td><td style="text-align:right">{n_study_cr} studies</td></tr>
</table>

<!-- ── INPUT FILES ── -->
<h2>12. Input Files</h2>
<table>
<tr><th>File</th><th>Source</th><th>Content</th></tr>
<tr><td><code>oral_chim_perc_YL_260519.csv</code></td><td>DADA2 chimera detection</td><td>chimera percentage per sample</td></tr>
<tr><td><code>raw_reads_stats.tsv</code></td><td>count_raw_reads.py</td><td>total raw reads per sample</td></tr>
<tr><td><code>sample_asv_counts_new.tsv</code></td><td>count_sample_asvs.py (DADA2.biom)</td><td>ASV count per sample</td></tr>
<tr><td><code>alpha_diversity_data.csv</code></td><td>closed-reference OTU picking</td><td>alpha diversity + metadata per sample</td></tr>
<tr><td><code>all_samples_gg2_closed_reference_mapping_fraction.tsv</code></td><td>GG2 CR containment derep</td><td>sample-level mapping fraction (119,181 entries)</td></tr>
<tr><td><code>metadata-v2-merge-cleaned.csv</code></td><td>metadata curation</td><td>subject_id / site_group / site_cleaned</td></tr>
<tr><td><code>gg2_reference_id_profile.biom</code></td><td>GG2 closed-reference profile</td><td>OTU table (HDF5, 65 MB)</td></tr>
</table>

<!-- ── OUTPUT FILES ── -->
<h2>13. Output Files</h2>
<table>
<tr><th>File</th><th>Samples</th><th>Studies</th><th>Subjects</th><th>Description</th></tr>
<tr><td><code>metadata_all_for_qc_raw.tsv</code></td><td style="text-align:right">{n_step0:,}</td><td style="text-align:right">{n_studies_raw}</td><td>—</td><td>Raw merge of all sources, no filtering</td></tr>
<tr><td><code>metadata_all_qc.tsv</code></td><td style="text-align:right">{n_step4:,}</td><td style="text-align:right">{n_studies_qc}</td><td>—</td><td>QC-filtered, all body sites</td></tr>
<tr><td><code>metadata_oral_qc.tsv</code></td><td style="text-align:right">{n_oral:,}</td><td style="text-align:right">{n_studies_qc}</td><td>—</td><td>QC-filtered, oral (<code>site_group='oral'</code>)</td></tr>
<tr><td><code>metadata_oral_qc_dedup.tsv</code></td><td style="text-align:right">{n_dedup:,}</td><td style="text-align:right">{n_studies_qc}</td><td style="text-align:right">{n_subjects_dedup:,}</td><td>QC-filtered, oral, deduplicated</td></tr>
</table>

<p style="text-align:center;color:#999;margin-top:48px;font-size:12px;">
Generated 2026-06-26 by <code>generate_qc_report.py</code><br>
CR filter thresholds: study-level aggregate mapping rate &lt;30%, sample-level mapping fraction &lt;30%
</p>

</body>
</html>'''

# ── Write ──────────────────────────────────────────────────────────
with open(OUT_HTML, 'w', encoding='utf-8') as f:
    f.write(html)

out_size = os.path.getsize(OUT_HTML)
print(f"\nReport written to: {OUT_HTML}")
print(f"  Size: {out_size / (1024*1024):.1f} MB")
print("Done!")
