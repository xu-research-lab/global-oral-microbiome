#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

BASE = Path('/beegfs/huangxiaochang/global_oral/new_analysis')
PRED = BASE / 'age_prediction/grouped_bioproject_cv/age_rf_grouped_cv_predictions.csv'
UPSET = BASE / 'maaslin3/disease_shared_plots/upset_set_sizes.csv'
OUT = BASE / 'age_prediction/grouped_bioproject_cv/disease_mad_plots'
OUT.mkdir(parents=True, exist_ok=True)

Q_MIN_N = 5
site_order = ['Saliva', 'Tooth', 'Tongue', 'Oral Mucosa', 'Oral_swab']
label_to_status = {
    'H&N Cancer': 'Head and Neck Cancer',
    'Periodontal Dz': 'Periodontal Disease',
    'Breast Cancer': 'Breast Cancer',
    'Lung Cancer': 'Lung Cancer',
    'OSCC': 'OSCC',
    'HIV': 'HIV',
    'Dental Caries': 'Dental Caries',
    'IBD': 'IBD',
    'RA': 'RA',
    'CRC': 'CRC',
}
status_to_label = {v:k for k,v in label_to_status.items()}

sns.set_theme(style='whitegrid', context='notebook')
pred = pd.read_csv(PRED)
pred['MAD'] = pred['PredictedAge'] - pred['ActualAge']
upset = pd.read_csv(UPSET)
disease_labels = upset['Disease'].tolist()
disease_status = [label_to_status.get(x, x) for x in disease_labels]

keep_status = ['Healthy Control'] + disease_status
plot = pred[pred['Status_Cleaned'].isin(keep_status)].copy()
plot['DiseaseLabel'] = plot['Status_Cleaned'].map(status_to_label).fillna(plot['Status_Cleaned'])
plot['Site'] = plot['Site'].astype(str)
plot = plot[plot['Site'].isin(site_order)].copy()

summary_rows = []
for site in site_order:
    site_df = plot[plot['Site'] == site]
    healthy = site_df[site_df['Status_Cleaned'] == 'Healthy Control']['MAD'].dropna()
    for status,label in zip(disease_status, disease_labels):
        dis = site_df[site_df['Status_Cleaned'] == status]['MAD'].dropna()
        if len(dis) == 0 and len(healthy) == 0:
            continue
        p = np.nan
        if len(dis) >= Q_MIN_N and len(healthy) >= Q_MIN_N:
            try:
                p = mannwhitneyu(dis, healthy, alternative='two-sided').pvalue
            except Exception:
                p = np.nan
        summary_rows.append({
            'Site': site,
            'Disease': label,
            'Status_Cleaned': status,
            'Disease_N': int(len(dis)),
            'Healthy_N': int(len(healthy)),
            'Disease_MAD_median': float(np.median(dis)) if len(dis) else np.nan,
            'Healthy_MAD_median': float(np.median(healthy)) if len(healthy) else np.nan,
            'Delta_median_MAD': float(np.median(dis) - np.median(healthy)) if len(dis) and len(healthy) else np.nan,
            'Disease_MAD_mean': float(np.mean(dis)) if len(dis) else np.nan,
            'Healthy_MAD_mean': float(np.mean(healthy)) if len(healthy) else np.nan,
            'Delta_mean_MAD': float(np.mean(dis) - np.mean(healthy)) if len(dis) and len(healthy) else np.nan,
            'MannWhitney_p': p,
        })
summary = pd.DataFrame(summary_rows)
summary['FDR_BH'] = np.nan
valid = summary['MannWhitney_p'].notna()
if valid.any():
    pvals = summary.loc[valid, 'MannWhitney_p'].to_numpy()
    order = np.argsort(pvals)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(pvals)+1)
    q = pvals * len(pvals) / ranks
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    q_final = np.empty_like(q_sorted)
    q_final[order] = np.minimum(q_sorted, 1)
    summary.loc[valid, 'FDR_BH'] = q_final
summary.to_csv(OUT / 'disease_vs_healthy_same_site_MAD_summary.csv', index=False)

# Long boxplot data: keep only disease-site with disease N >= 5 plus same-site healthy.
box_parts = []
for site in site_order:
    healthy = plot[(plot['Site'] == site) & (plot['Status_Cleaned'] == 'Healthy Control')].copy()
    if len(healthy) == 0:
        continue
    healthy['PlotGroup'] = 'Healthy'
    box_parts.append(healthy)
    for status,label in zip(disease_status, disease_labels):
        dis = plot[(plot['Site'] == site) & (plot['Status_Cleaned'] == status)].copy()
        if len(dis) >= Q_MIN_N:
            dis['PlotGroup'] = label
            box_parts.append(dis)
box = pd.concat(box_parts, ignore_index=True) if box_parts else pd.DataFrame()
box.to_csv(OUT / 'disease_vs_healthy_same_site_MAD_plot_data.csv', index=False)

if not box.empty:
    # Faceted boxplot by site. Disease order within all panels follows upset order.
    group_order = ['Healthy'] + disease_labels
    box['PlotGroup'] = pd.Categorical(box['PlotGroup'], categories=group_order, ordered=True)
    box['Site'] = pd.Categorical(box['Site'], categories=site_order, ordered=True)
    g = sns.catplot(
        data=box, x='PlotGroup', y='MAD', col='Site', col_order=site_order,
        kind='box', showfliers=False, col_wrap=2, height=4.2, aspect=1.45,
        order=group_order, color='#9ecae1', sharex=False, sharey=True
    )
    # Overlay jitter manually per facet, using only groups present in that site.
    for ax, site in zip(g.axes.flat, site_order):
        sdf = box[box['Site'].astype(str) == site]
        present = [x for x in group_order if x in set(sdf['PlotGroup'].astype(str))]
        ax.clear()
        sns.boxplot(data=sdf, x='PlotGroup', y='MAD', order=present, showfliers=False, color='#9ecae1', width=0.45, ax=ax)
        sns.stripplot(data=sdf, x='PlotGroup', y='MAD', order=present, size=1.4, alpha=0.18, jitter=0.18, color='#22577a', ax=ax)
        ax.axhline(0, linestyle='--', color='#b2182b', linewidth=1.0)
        ax.set_title(site, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('MAD' if ax in g.axes.flat[::2] else '')
        ax.tick_params(axis='x', rotation=45, labelsize=8)
    g.fig.suptitle('Disease vs same-site healthy MAD', y=1.02, fontsize=14, fontweight='bold')
    g.fig.tight_layout()
    g.fig.savefig(OUT / 'disease_vs_healthy_same_site_MAD_boxplot.pdf', bbox_inches='tight')
    g.fig.savefig(OUT / 'disease_vs_healthy_same_site_MAD_boxplot.png', dpi=300, bbox_inches='tight')
    plt.close(g.fig)

# Heatmap: disease x site, color = disease median MAD - same-site healthy median MAD.
hm = summary[(summary['Disease_N'] >= Q_MIN_N) & summary['Healthy_N'].ge(Q_MIN_N)].copy()
if not hm.empty:
    mat = hm.pivot(index='Disease', columns='Site', values='Delta_median_MAD').reindex(index=disease_labels, columns=site_order)
    annot = hm.copy()
    annot['label'] = annot.apply(lambda r: f"{r['Delta_median_MAD']:.1f}\nn={int(r['Disease_N'])}" if pd.notna(r['Delta_median_MAD']) else '', axis=1)
    annot_mat = annot.pivot(index='Disease', columns='Site', values='label').reindex(index=disease_labels, columns=site_order).fillna('')
    width = max(6.5, 1.15 * len(site_order) + 2)
    height = max(4.8, 0.42 * len(disease_labels) + 1.8)
    plt.figure(figsize=(width, height))
    sns.heatmap(mat, annot=annot_mat, fmt='', cmap='RdBu_r', center=0, linewidths=0.5,
                cbar_kws={'label':'Median MAD difference\nDisease - same-site healthy'})
    plt.title('Disease-associated microbiome age deviation by site', fontweight='bold')
    plt.xlabel('Sampling site')
    plt.ylabel('Disease')
    plt.tight_layout()
    plt.savefig(OUT / 'disease_vs_healthy_same_site_MAD_delta_heatmap.pdf')
    plt.savefig(OUT / 'disease_vs_healthy_same_site_MAD_delta_heatmap.png', dpi=300)
    plt.close()

print('Diseases from upset:', ', '.join(disease_labels))
print(summary.sort_values(['Disease','Site']).to_string(index=False, float_format=lambda x: f'{x:.3g}'))
print('Output:', OUT)
