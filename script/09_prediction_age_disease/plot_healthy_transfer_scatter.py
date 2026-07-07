#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr, spearmanr

BASE = Path('/beegfs/huangxiaochang/global_oral/new_analysis/age_prediction/healthy_trained_transfer')
PRED = BASE / 'healthy_cv_predictions.csv'
OUT = BASE / 'healthy_cv_scatter_plots'
OUT.mkdir(parents=True, exist_ok=True)

SITE_ORDER = ['Tooth', 'Saliva', 'Tongue', 'Oral Mucosa', 'Oral_swab']


def safe_corr(y, p, func):
    if len(y) < 3 or np.std(y) == 0 or np.std(p) == 0:
        return np.nan
    return float(func(y, p)[0])


def calc_metrics(df):
    y = df['ActualAge'].to_numpy(float)
    p = df['PredictedAge'].to_numpy(float)
    return {
        'N': len(df),
        'BioProjects': df['BioProject'].nunique(),
        'MAE': mean_absolute_error(y, p),
        'RMSE': mean_squared_error(y, p, squared=False),
        'R2': r2_score(y, p) if len(np.unique(y)) > 1 else np.nan,
        'PearsonR': safe_corr(y, p, pearsonr),
        'SpearmanR': safe_corr(y, p, spearmanr),
    }


def age_limits(df):
    y = df['ActualAge'].to_numpy(float)
    p = df['PredictedAge'].to_numpy(float)
    lo = np.floor(min(y.min(), p.min()) / 5) * 5
    hi = np.ceil(max(y.max(), p.max()) / 5) * 5
    return lo, hi


def main():
    sns.set_theme(style='whitegrid', context='notebook')
    pred = pd.read_csv(PRED)
    pred = pred[pred['PredictionSet'].eq('Healthy_CV')].copy()
    site_order = [s for s in SITE_ORDER if s in set(pred['Site'])]
    site_order += sorted(set(pred['Site']) - set(site_order))

    rows = []
    for site in site_order:
        df = pred[pred['Site'] == site].copy()
        if df.empty:
            continue
        stat = calc_metrics(df)
        stat['Site'] = site
        rows.append(stat)

        lo, hi = age_limits(df)
        fig, ax = plt.subplots(figsize=(6.0, 5.5))
        sns.scatterplot(data=df, x='ActualAge', y='PredictedAge', hue='Fold',
                        palette='tab10', s=17, alpha=0.50, linewidth=0, ax=ax)
        ax.plot([lo, hi], [lo, hi], linestyle='--', color='black', linewidth=1.0)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect('equal', adjustable='box')
        ax.set_title(f'{site}: healthy controls', fontsize=13, fontweight='bold')
        ax.set_xlabel('Actual age')
        ax.set_ylabel('Predicted age')
        label = (
            f"N={stat['N']}; BioProjects={stat['BioProjects']}\n"
            f"MAE={stat['MAE']:.2f}; RMSE={stat['RMSE']:.2f}\n"
            f"R²={stat['R2']:.2f}; r={stat['PearsonR']:.2f}; ρ={stat['SpearmanR']:.2f}"
        )
        ax.text(0.03, 0.97, label, transform=ax.transAxes, va='top', ha='left', fontsize=9.2,
                bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='0.7', alpha=0.9))
        ax.legend(title='Fold', loc='lower right', fontsize=7.5, title_fontsize=8, frameon=True)
        fig.tight_layout()
        safe = site.replace(' ', '_')
        fig.savefig(OUT / f'{safe}_healthy_actual_vs_predicted_age.pdf')
        fig.savefig(OUT / f'{safe}_healthy_actual_vs_predicted_age.png', dpi=300)
        plt.close(fig)

    summary = pd.DataFrame(rows)
    summary = summary[['Site', 'N', 'BioProjects', 'MAE', 'RMSE', 'R2', 'PearsonR', 'SpearmanR']]
    summary.to_csv(OUT / 'healthy_cv_scatter_plot_metrics.csv', index=False)

    plot_df = pred[pred['Site'].isin(site_order)].copy()
    if not plot_df.empty:
        g = sns.FacetGrid(plot_df, col='Site', col_order=site_order, col_wrap=2,
                          height=3.8, sharex=False, sharey=False)
        g.map_dataframe(sns.scatterplot, x='ActualAge', y='PredictedAge',
                        s=10, alpha=0.32, linewidth=0, color='#2b6c9c')
        stat_by_site = {r['Site']: r for r in rows}
        for ax, site in zip(g.axes.flat, site_order):
            df = plot_df[plot_df['Site'] == site]
            lo, hi = age_limits(df)
            ax.plot([lo, hi], [lo, hi], linestyle='--', color='black', linewidth=0.9)
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
            ax.set_aspect('equal', adjustable='box')
            stat = stat_by_site[site]
            ax.text(0.04, 0.96, f"MAE={stat['MAE']:.1f}\nR²={stat['R2']:.2f}\nr={stat['PearsonR']:.2f}",
                    transform=ax.transAxes, va='top', ha='left', fontsize=8,
                    bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='0.75', alpha=0.9))
        g.set_axis_labels('Actual age', 'Predicted age')
        g.set_titles('{col_name}')
        g.fig.suptitle('Healthy-trained age model performance', y=1.02, fontsize=14, fontweight='bold')
        g.fig.text(0.5, 0.985, 'Healthy controls only; BioProject grouped CV; test age within training range',
                   ha='center', fontsize=9.5)
        g.fig.tight_layout()
        g.fig.savefig(OUT / 'all_sites_healthy_actual_vs_predicted_age.pdf', bbox_inches='tight')
        g.fig.savefig(OUT / 'all_sites_healthy_actual_vs_predicted_age.png', dpi=300, bbox_inches='tight')
        plt.close(g.fig)

    print(summary.sort_values('MAE').to_string(index=False, float_format=lambda x: f'{x:.3f}'))
    print('Output:', OUT)


if __name__ == '__main__':
    main()
