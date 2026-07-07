#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

BASE = Path('/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3')
sns.set_theme(style='whitegrid')


def extract_genus(feature: str) -> str:
    if 'g__' in feature:
        part = feature.split('g__', 1)[1]
        return part.split(';', 1)[0] or feature
    return feature


def collect_result_files(root: Path):
    return sorted(root.glob('*/all_results.tsv'))


def summarize_stratified():
    files = collect_result_files(BASE / 'results')
    if not files:
        return
    frames = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        df['Group'] = f.parent.name
        parts = f.parent.name.split('_', 1)
        df['Site_Cleaned'] = parts[0]
        df['Geography'] = parts[1] if len(parts) > 1 else ''
        frames.append(df)
    all_df = pd.concat(frames, ignore_index=True)
    all_df['Genus'] = all_df['feature'].map(extract_genus)
    all_df.to_csv(BASE / 'stratified_results_all.tsv', sep='\t', index=False)

    sig = all_df[all_df['qval_individual'].lt(0.1) & all_df['qval_individual'].notna()].copy()
    sig.to_csv(BASE / 'stratified_significant_results_q0.1.tsv', sep='\t', index=False)

    counts = sig.groupby(['Group', 'metadata']).size().reset_index(name='n_sig').sort_values('n_sig', ascending=False)
    counts.to_csv(BASE / 'stratified_significant_counts_by_group_metadata.csv', index=False)

    if not counts.empty:
        heat = counts.pivot(index='Group', columns='metadata', values='n_sig').fillna(0)
        plt.figure(figsize=(8, max(4, 0.35 * len(heat))))
        sns.heatmap(heat, annot=True, fmt='.0f', cmap='Blues')
        plt.title('Significant MaAsLin3 Associations per Group')
        plt.tight_layout()
        plt.savefig(BASE / 'stratified_significant_counts_heatmap.png', dpi=300)
        plt.close()

    top_genera = sig.groupby(['Genus', 'metadata']).size().reset_index(name='n').sort_values('n', ascending=False)
    top_genera.to_csv(BASE / 'stratified_top_genera_by_metadata.csv', index=False)
    top20 = top_genera.groupby('Genus', as_index=False)['n'].sum().sort_values('n', ascending=False).head(20)['Genus']
    top_sig = sig[sig['Genus'].isin(top20)].copy()
    if not top_sig.empty:
        hm = top_sig.groupby(['Genus', 'Group']).size().reset_index(name='n')
        mat = hm.pivot(index='Genus', columns='Group', values='n').fillna(0)
        plt.figure(figsize=(max(8, 0.45 * len(mat.columns)), max(6, 0.3 * len(mat))))
        sns.heatmap(mat, cmap='Reds')
        plt.title('Top Significant Genera Across Stratified Groups')
        plt.tight_layout()
        plt.savefig(BASE / 'stratified_top_genera_heatmap.png', dpi=300)
        plt.close()


def summarize_disease():
    files = collect_result_files(BASE / 'disease_results')
    if not files:
        return
    frames = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        df['Disease'] = f.parent.name
        frames.append(df)
    all_df = pd.concat(frames, ignore_index=True)
    all_df['Genus'] = all_df['feature'].map(extract_genus)
    all_df.to_csv(BASE / 'disease_results_all.tsv', sep='\t', index=False)

    sig = all_df[(all_df['metadata'] == 'Disease_Status') & all_df['qval_individual'].lt(0.1) & all_df['qval_individual'].notna()].copy()
    sig.to_csv(BASE / 'disease_significant_results_q0.1.tsv', sep='\t', index=False)

    counts = sig.groupby('Disease').size().reset_index(name='n_sig').sort_values('n_sig', ascending=False)
    counts.to_csv(BASE / 'disease_significant_counts.csv', index=False)
    if not counts.empty:
        plt.figure(figsize=(8, 4))
        sns.barplot(data=counts, x='Disease', y='n_sig', color='#4C78A8')
        plt.xticks(rotation=30, ha='right')
        plt.title('Significant Disease-associated Genera (q < 0.1)')
        plt.tight_layout()
        plt.savefig(BASE / 'disease_significant_counts_barplot.png', dpi=300)
        plt.close()

    if not sig.empty:
        # Multiple MaAsLin3 rows can map to the same genus; collapse before pivoting.
        mat_df = sig.groupby(['Genus', 'Disease'], as_index=False)['coef'].mean()
        mat_df['direction'] = mat_df['coef'].gt(0).map({True: 1, False: -1})
        mat = mat_df.pivot(index='Genus', columns='Disease', values='direction').fillna(0)
        plt.figure(figsize=(max(6, 1.2 * len(mat.columns)), max(6, 0.25 * len(mat))))
        sns.heatmap(mat, cmap='RdBu_r', center=0, vmin=-1, vmax=1)
        plt.title('Direction of Disease-associated Genera')
        plt.tight_layout()
        plt.savefig(BASE / 'disease_significant_direction_heatmap.png', dpi=300)
        plt.close()


def main():
    summarize_stratified()
    summarize_disease()
    print('Done.')


if __name__ == '__main__':
    main()

