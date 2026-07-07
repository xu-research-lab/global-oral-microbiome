#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon

BASE = Path('/beegfs/huangxiaochang/global_oral/new_analysis/age_prediction/grouped_bioproject_cv')
PRED = BASE / 'age_rf_grouped_cv_predictions.csv'
OUT = BASE / 'plots'
OUT.mkdir(parents=True, exist_ok=True)

sns.set_theme(style='whitegrid', context='notebook')
site_order = ['Tooth', 'Saliva', 'Tongue', 'Oral Mucosa', 'Oral_swab']

pred = pd.read_csv(PRED)
pred['MAD'] = pred['PredictedAge'] - pred['ActualAge']
site_order = [s for s in site_order if s in set(pred['Site'])] + sorted(set(pred['Site']) - set(site_order))
pred['Site'] = pd.Categorical(pred['Site'], categories=site_order, ordered=True)

rows = []
for site, df in pred.groupby('Site', observed=True):
    mad = df['MAD'].dropna().to_numpy(float)
    if len(mad) == 0:
        continue
    try:
        p = wilcoxon(mad, alternative='two-sided').pvalue if len(mad) > 0 and np.any(mad != 0) else np.nan
    except Exception:
        p = np.nan
    rows.append({
        'Site': site,
        'N': len(mad),
        'BioProjects': df['BioProject'].nunique(),
        'MAD_mean': np.mean(mad),
        'MAD_median': np.median(mad),
        'MAD_sd': np.std(mad, ddof=1),
        'MAD_IQR': np.percentile(mad, 75) - np.percentile(mad, 25),
        'MAD_q25': np.percentile(mad, 25),
        'MAD_q75': np.percentile(mad, 75),
        'Wilcoxon_p_vs_0': p,
    })
summary = pd.DataFrame(rows)
summary.to_csv(OUT / 'age_rf_mad_summary.csv', index=False)

# Combined site MAD plot
fig, ax = plt.subplots(figsize=(8.8, 5.4))
sns.boxplot(data=pred, x='Site', y='MAD', order=site_order, showfliers=False, width=0.46, color='#9ecae1', ax=ax)
sns.stripplot(data=pred, x='Site', y='MAD', order=site_order, size=1.4, alpha=0.18, jitter=0.22, color='#22577a', ax=ax)
ax.axhline(0, linestyle='--', color='#b2182b', linewidth=1.1)
ax.set_xlabel(None)
ax.set_ylabel('MAD (Predicted age - Actual age)')
ax.set_title('Microbiome age deviation by oral site', fontsize=13, fontweight='bold')
ax.tick_params(axis='x', rotation=25)
# annotate median and N
ylim = ax.get_ylim()
y_text = ylim[1] - (ylim[1] - ylim[0]) * 0.06
for i, site in enumerate(site_order):
    row = summary[summary['Site'] == site]
    if row.empty:
        continue
    row = row.iloc[0]
    ax.text(i, y_text, f"N={int(row.N)}\nmed={row.MAD_median:.1f}", ha='center', va='top', fontsize=8)
fig.tight_layout()
fig.savefig(OUT / 'all_sites_MAD_boxplot.pdf')
fig.savefig(OUT / 'all_sites_MAD_boxplot.png', dpi=300)
plt.close(fig)

# Violin + box plot
fig, ax = plt.subplots(figsize=(8.8, 5.4))
sns.violinplot(data=pred, x='Site', y='MAD', order=site_order, inner=None, cut=0, color='#c6dbef', linewidth=0.8, ax=ax)
sns.boxplot(data=pred, x='Site', y='MAD', order=site_order, showfliers=False, width=0.22, color='white', ax=ax)
ax.axhline(0, linestyle='--', color='#b2182b', linewidth=1.1)
ax.set_xlabel(None)
ax.set_ylabel('MAD (Predicted age - Actual age)')
ax.set_title('MAD distribution by oral site', fontsize=13, fontweight='bold')
ax.tick_params(axis='x', rotation=25)
fig.tight_layout()
fig.savefig(OUT / 'all_sites_MAD_violin.pdf')
fig.savefig(OUT / 'all_sites_MAD_violin.png', dpi=300)
plt.close(fig)

# Per-site MAD histogram/density and fold boxplots
for site in site_order:
    df = pred[pred['Site'] == site].copy()
    if df.empty:
        continue
    safe = str(site).replace(' ', '_')
    row = summary[summary['Site'] == site].iloc[0]

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    sns.histplot(df['MAD'].to_numpy(), bins=35, kde=False, color='#4C78A8', alpha=0.55, ax=ax)
    ax.axvline(0, linestyle='--', color='#b2182b', linewidth=1.1)
    ax.axvline(row['MAD_median'], linestyle='-', color='black', linewidth=1.0)
    ax.set_xlabel('MAD (Predicted age - Actual age)')
    ax.set_ylabel('Samples')
    ax.set_title(f'{site}: MAD distribution', fontsize=13, fontweight='bold')
    ax.text(0.98, 0.96, f"N={int(row.N)}\nmedian={row.MAD_median:.2f}\nmean={row.MAD_mean:.2f}", transform=ax.transAxes,
            ha='right', va='top', fontsize=9.5, bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='0.75'))
    fig.tight_layout()
    fig.savefig(OUT / f'{safe}_MAD_distribution.pdf')
    fig.savefig(OUT / f'{safe}_MAD_distribution.png', dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    sns.boxplot(data=df, x='Fold', y='MAD', showfliers=False, color='#b3cde3', ax=ax)
    sns.stripplot(data=df, x='Fold', y='MAD', size=1.3, alpha=0.20, jitter=0.20, color='#22577a', ax=ax)
    ax.axhline(0, linestyle='--', color='#b2182b', linewidth=1.1)
    ax.set_xlabel('Fold')
    ax.set_ylabel('MAD')
    ax.set_title(f'{site}: MAD by BioProject-CV fold', fontsize=13, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUT / f'{safe}_MAD_by_fold.pdf')
    fig.savefig(OUT / f'{safe}_MAD_by_fold.png', dpi=300)
    plt.close(fig)

print(summary.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
print('Output:', OUT)
