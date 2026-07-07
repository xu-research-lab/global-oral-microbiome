#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr, spearmanr

BASE = Path('/beegfs/huangxiaochang/global_oral/new_analysis/age_prediction/grouped_bioproject_cv')
PRED = BASE / 'age_rf_grouped_cv_predictions.csv'
OUT = BASE / 'plots'
OUT.mkdir(parents=True, exist_ok=True)

sns.set_theme(style='whitegrid', context='notebook')

pred = pd.read_csv(PRED)
site_order = ['Tooth', 'Saliva', 'Tongue', 'Oral Mucosa', 'Oral_swab']
site_order = [s for s in site_order if s in set(pred['Site'])] + sorted(set(pred['Site']) - set(site_order))
rows = []

for site in site_order:
    df = pred[pred['Site'] == site].copy()
    if df.empty:
        continue
    y = df['ActualAge'].to_numpy(float)
    p = df['PredictedAge'].to_numpy(float)
    mae = mean_absolute_error(y, p)
    rmse = mean_squared_error(y, p, squared=False)
    r2 = r2_score(y, p) if len(np.unique(y)) > 1 else np.nan
    pr = pearsonr(y, p)[0] if len(np.unique(y)) > 1 and len(np.unique(p)) > 1 else np.nan
    sr = spearmanr(y, p)[0] if len(np.unique(y)) > 1 and len(np.unique(p)) > 1 else np.nan
    rows.append({'Site': site, 'N': len(df), 'BioProjects': df['BioProject'].nunique(), 'MAE': mae, 'RMSE': rmse, 'R2': r2, 'PearsonR': pr, 'SpearmanR': sr})

    min_age = np.floor(min(y.min(), p.min()) / 5) * 5
    max_age = np.ceil(max(y.max(), p.max()) / 5) * 5
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    sns.scatterplot(data=df, x='ActualAge', y='PredictedAge', hue='Fold', palette='tab10', s=18, alpha=0.55, linewidth=0, ax=ax)
    ax.plot([min_age, max_age], [min_age, max_age], linestyle='--', color='black', linewidth=1.1)
    ax.set_xlim(min_age, max_age)
    ax.set_ylim(min_age, max_age)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(site, fontsize=14, fontweight='bold')
    ax.set_xlabel('Actual age')
    ax.set_ylabel('Predicted age')
    label = f'N={len(df)}; BioProjects={df.BioProject.nunique()}\nMAE={mae:.2f}; RMSE={rmse:.2f}\nR²={r2:.2f}; r={pr:.2f}; ρ={sr:.2f}'
    ax.text(0.03, 0.97, label, transform=ax.transAxes, va='top', ha='left', fontsize=9.5,
            bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='0.7', alpha=0.9))
    ax.legend(title='Fold', loc='lower right', fontsize=8, title_fontsize=8, frameon=True)
    fig.tight_layout()
    safe = site.replace(' ', '_')
    fig.savefig(OUT / f'{safe}_actual_vs_predicted_age.pdf')
    fig.savefig(OUT / f'{safe}_actual_vs_predicted_age.png', dpi=300)
    plt.close(fig)

# Combined facet plot
plot_df = pred[pred['Site'].isin(site_order)].copy()
g = sns.FacetGrid(plot_df, col='Site', col_order=site_order, col_wrap=3, height=3.6, sharex=False, sharey=False)
g.map_dataframe(sns.scatterplot, x='ActualAge', y='PredictedAge', s=10, alpha=0.35, linewidth=0, color='#2b6c9c')
for ax, site in zip(g.axes.flat, site_order):
    df = plot_df[plot_df['Site'] == site]
    y = df['ActualAge'].to_numpy(float)
    p = df['PredictedAge'].to_numpy(float)
    min_age = np.floor(min(y.min(), p.min()) / 5) * 5
    max_age = np.ceil(max(y.max(), p.max()) / 5) * 5
    ax.plot([min_age, max_age], [min_age, max_age], linestyle='--', color='black', linewidth=0.9)
    ax.set_xlim(min_age, max_age)
    ax.set_ylim(min_age, max_age)
    ax.set_aspect('equal', adjustable='box')
    stat = next(r for r in rows if r['Site'] == site)
    ax.text(0.04, 0.96, f"MAE={stat['MAE']:.1f}\nR²={stat['R2']:.2f}\nr={stat['PearsonR']:.2f}",
            transform=ax.transAxes, va='top', ha='left', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='0.75', alpha=0.9))
g.set_axis_labels('Actual age', 'Predicted age')
g.set_titles('{col_name}')
g.fig.suptitle('Grouped BioProject CV age prediction', y=1.02, fontsize=14, fontweight='bold')
g.fig.tight_layout()
g.fig.savefig(OUT / 'all_sites_actual_vs_predicted_age.pdf', bbox_inches='tight')
g.fig.savefig(OUT / 'all_sites_actual_vs_predicted_age.png', dpi=300, bbox_inches='tight')
plt.close(g.fig)

summary = pd.DataFrame(rows).sort_values('MAE')
summary.to_csv(OUT / 'age_rf_scatter_plot_metrics.csv', index=False)
print(summary.to_string(index=False, float_format=lambda x: f'{x:.3f}'))
print('Output:', OUT)
