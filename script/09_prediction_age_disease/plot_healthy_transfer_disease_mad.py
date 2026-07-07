#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

BASE = Path('/beegfs/huangxiaochang/global_oral/new_analysis')
INDIR = BASE / 'age_prediction' / 'healthy_trained_transfer'
HEALTHY = INDIR / 'healthy_cv_predictions.csv'
DISEASE = INDIR / 'disease_transfer_predictions_sample_mean.csv'
UPSET = BASE / 'maaslin3' / 'disease_shared_plots' / 'upset_set_sizes.csv'
OUT = INDIR / 'disease_mad_plots'
OUT.mkdir(parents=True, exist_ok=True)

MIN_N = 5
SITE_ORDER = ['Saliva', 'Tooth', 'Tongue', 'Oral Mucosa', 'Oral_swab']
LABEL_TO_STATUS = {
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


def sig_stars(q):
    if not np.isfinite(q):
        return ''
    if q < 0.001:
        return '***'
    if q < 0.01:
        return '**'
    if q < 0.05:
        return '*'
    return ''


def bh_fdr(pvals):
    out = np.full(len(pvals), np.nan)
    valid = np.isfinite(pvals)
    if not valid.any():
        return out
    pv = pvals[valid]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.minimum(ranked, 1)
    out[np.where(valid)[0]] = q
    return out


def main():
    sns.set_theme(style='whitegrid', context='notebook')
    disease_labels = pd.read_csv(UPSET)['Disease'].tolist()
    disease_statuses = [LABEL_TO_STATUS.get(x, x) for x in disease_labels]

    healthy = pd.read_csv(HEALTHY)
    disease = pd.read_csv(DISEASE)
    healthy['DiseaseLabel'] = 'Healthy'
    healthy['Status_Cleaned'] = 'Healthy Control'
    disease = disease[disease['Status_Cleaned'].isin(disease_statuses)].copy()

    healthy = healthy[healthy['Site'].isin(SITE_ORDER)].copy()
    disease = disease[disease['Site'].isin(SITE_ORDER)].copy()
    healthy['MAD'] = healthy['PredictedAge'] - healthy['ActualAge']
    disease['MAD'] = disease['PredictedAge'] - disease['ActualAge']

    rows = []
    for site in SITE_ORDER:
        h = healthy.loc[healthy['Site'] == site, 'MAD'].dropna()
        for status, label in zip(disease_statuses, disease_labels):
            d = disease.loc[(disease['Site'] == site) & (disease['Status_Cleaned'] == status), 'MAD'].dropna()
            p = np.nan
            if len(h) >= MIN_N and len(d) >= MIN_N:
                p = mannwhitneyu(d, h, alternative='two-sided').pvalue
            rows.append({
                'Site': site,
                'Disease': label,
                'Status_Cleaned': status,
                'Disease_N': int(len(d)),
                'Healthy_N': int(len(h)),
                'Disease_MAD_median': float(np.median(d)) if len(d) else np.nan,
                'Healthy_MAD_median': float(np.median(h)) if len(h) else np.nan,
                'Delta_median_MAD': float(np.median(d) - np.median(h)) if len(d) and len(h) else np.nan,
                'Disease_MAD_mean': float(np.mean(d)) if len(d) else np.nan,
                'Healthy_MAD_mean': float(np.mean(h)) if len(h) else np.nan,
                'Delta_mean_MAD': float(np.mean(d) - np.mean(h)) if len(d) and len(h) else np.nan,
                'MannWhitney_p': p,
            })
    summary = pd.DataFrame(rows)
    summary['FDR_BH'] = bh_fdr(summary['MannWhitney_p'].to_numpy(float))
    summary['Significance'] = summary['FDR_BH'].apply(sig_stars)
    summary.to_csv(OUT / 'healthy_trained_disease_vs_same_site_healthy_MAD_summary.csv', index=False)

    box_parts = []
    for site in SITE_ORDER:
        h = healthy[healthy['Site'] == site].copy()
        if len(h) >= MIN_N:
            h['PlotGroup'] = 'Healthy'
            box_parts.append(h)
        for status, label in zip(disease_statuses, disease_labels):
            d = disease[(disease['Site'] == site) & (disease['Status_Cleaned'] == status)].copy()
            if len(d) >= MIN_N:
                d['PlotGroup'] = label
                box_parts.append(d)
    box = pd.concat(box_parts, ignore_index=True) if box_parts else pd.DataFrame()
    box.to_csv(OUT / 'healthy_trained_disease_MAD_boxplot_data.csv', index=False)

    if not box.empty:
        group_order = ['Healthy'] + disease_labels
        plot_sites = [s for s in SITE_ORDER if s in set(box['Site'])]
        ncols = 2
        nrows = int(np.ceil(len(plot_sites) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(10.5, 3.45 * nrows), sharey=True)
        axes = np.atleast_1d(axes).ravel()
        for ax, site in zip(axes, plot_sites):
            sdf = box[box['Site'] == site].copy()
            disease_order = (sdf[sdf['PlotGroup'] != 'Healthy']
                             .groupby('PlotGroup')['MAD'].median()
                             .sort_values(ascending=False).index.tolist())
            present = ['Healthy'] + disease_order
            sns.boxplot(data=sdf, x='PlotGroup', y='MAD', order=present, showfliers=False,
                        color='#9ecae1', width=0.42, linewidth=0.8, ax=ax)
            sns.stripplot(data=sdf, x='PlotGroup', y='MAD', order=present, size=1.5,
                          alpha=0.18, jitter=0.18, color='#22577a', ax=ax)
            ax.axhline(0, linestyle='--', color='#b2182b', linewidth=1.0)
            site_summary = summary[summary['Site'] == site].set_index('Disease')
            y_min, y_max = sdf['MAD'].min(), sdf['MAD'].max()
            y_span = max(y_max - y_min, 1)
            ax.set_ylim(y_min - 0.06 * y_span, y_max + 0.18 * y_span)
            y_text = y_max + 0.05 * y_span
            for xpos, group in enumerate(present):
                if group == 'Healthy' or group not in site_summary.index:
                    continue
                star = site_summary.loc[group, 'Significance']
                if star:
                    ax.text(xpos, y_text, star, ha='center', va='bottom',
                            fontsize=9.5, fontweight='bold', color='black')
            ax.set_title(site, fontweight='bold')
            ax.set_xlabel('')
            ax.set_ylabel('Predicted age - actual age' if ax is axes[0] else '')
            ax.tick_params(axis='x', rotation=45, labelsize=7)
            ax.tick_params(axis='y', labelsize=8)
        for ax in axes[len(plot_sites):]:
            ax.axis('off')
        fig.suptitle('Disease age deviation using healthy-trained models', y=0.995, fontsize=12.5, fontweight='bold')
        fig.text(0.5, 0.955, 'Stars: FDR vs same-site healthy controls; no ComBat correction',
                 ha='center', fontsize=8.5)
        fig.tight_layout(rect=[0, 0, 1, 0.90])
        fig.savefig(OUT / 'healthy_trained_disease_vs_same_site_healthy_MAD_boxplot.pdf', bbox_inches='tight')
        fig.savefig(OUT / 'healthy_trained_disease_vs_same_site_healthy_MAD_boxplot.png',
                    dpi=300, bbox_inches='tight')
        plt.close(fig)

    hm = summary[(summary['Disease_N'] >= MIN_N) & (summary['Healthy_N'] >= MIN_N)].copy()
    if not hm.empty:
        mat = hm.pivot(index='Disease', columns='Site', values='Delta_median_MAD').reindex(
            index=disease_labels, columns=SITE_ORDER)
        annot = hm.copy()
        annot['label'] = annot.apply(
            lambda r: f"{r['Delta_median_MAD']:.1f}{r['Significance']}\nn={int(r['Disease_N'])}", axis=1)
        annot_mat = annot.pivot(index='Disease', columns='Site', values='label').reindex(
            index=disease_labels, columns=SITE_ORDER).fillna('')
        plt.figure(figsize=(6.7, max(4.0, 0.34 * len(disease_labels) + 1.45)))
        sns.heatmap(mat, annot=annot_mat, fmt='', cmap='RdBu_r', center=0, linewidths=0.5,
                    annot_kws={'fontsize': 7.5},
                    cbar_kws={'label': 'Median MAD difference\nDisease - same-site healthy',
                              'shrink': 0.72})
        plt.title('Disease age deviation', fontweight='bold')
        plt.xlabel('Sampling site')
        plt.ylabel('Disease')
        plt.tight_layout()
        plt.savefig(OUT / 'healthy_trained_disease_same_site_delta_MAD_heatmap.pdf')
        plt.savefig(OUT / 'healthy_trained_disease_same_site_delta_MAD_heatmap.png', dpi=300)
        plt.close()

    disease_counts = disease.groupby(['DiseaseLabel', 'Site']).size().reset_index(name='N')
    disease_counts.to_csv(OUT / 'healthy_trained_disease_prediction_counts.csv', index=False)
    print(summary.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
    print('Output:', OUT)


if __name__ == '__main__':
    main()
