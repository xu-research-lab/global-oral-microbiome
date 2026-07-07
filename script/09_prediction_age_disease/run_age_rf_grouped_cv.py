#!/usr/bin/env python3
from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from scipy.stats import pearsonr, spearmanr

BASE = Path('/beegfs/huangxiaochang/global_oral/new_analysis')
ROOT = BASE.parent
OUTDIR = BASE / 'age_prediction' / 'grouped_bioproject_cv'
OUTDIR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(BASE))

from analysis_inputs import load_metadata

GENUS_TABLE = BASE / 'maaslin3' / 'genus_table_relabund.tsv'
METADATA_FILE = BASE / 'metadata_oral_qc_dedup.tsv'
BIOM_FILE = BASE / 'biom_oral_qc_dedup.biom'

SITES = ['Saliva', 'Tooth', 'Tongue', 'Oral Mucosa', 'Oral_swab', 'Gingiva']
MIN_SITE_SAMPLES = 100
MIN_BIOPROJECTS = 5
N_SPLITS = 5
N_TREES = 500
RANDOM_SEED = 20260610
N_JOBS = 8
MIN_PREVALENCE = 0.05
MIN_TOTAL_READS = 1000
AGE_BIN_SIZE = 5
MAX_PER_AGE_BIN = 500
EPS = 1e-6


def safe_corr(func, y_true, y_pred):
    if len(y_true) < 3 or np.std(y_true) == 0 or np.std(y_pred) == 0:
        return np.nan, np.nan
    try:
        r, p = func(y_true, y_pred)
        return float(r), float(p)
    except Exception:
        return np.nan, np.nan


def balance_train_indices(train_meta, rng):
    ages = train_meta['age_numeric'].to_numpy()
    bins = np.floor(ages / AGE_BIN_SIZE).astype(int)
    selected = []
    for b in np.unique(bins):
        idx = np.where(bins == b)[0]
        if len(idx) > MAX_PER_AGE_BIN:
            idx = rng.choice(idx, size=MAX_PER_AGE_BIN, replace=False)
        selected.extend(idx.tolist())
    selected = np.array(selected, dtype=int)
    rng.shuffle(selected)
    return selected


def metrics_for(y_true, y_pred):
    mae = mean_absolute_error(y_true, y_pred)
    rmse = mean_squared_error(y_true, y_pred, squared=False)
    r2 = r2_score(y_true, y_pred) if len(np.unique(y_true)) > 1 else np.nan
    pearson_r, pearson_p = safe_corr(pearsonr, y_true, y_pred)
    spearman_r, spearman_p = safe_corr(spearmanr, y_true, y_pred)
    return dict(MAE=mae, RMSE=rmse, R2=r2, PearsonR=pearson_r, PearsonP=pearson_p,
                SpearmanR=spearman_r, SpearmanP=spearman_p)


def main():
    print('Loading metadata...')
    meta = load_metadata(dedup=True)
    age_col = 'age_cleaned' if 'age_cleaned' in meta.columns else 'age'
    meta['age_numeric'] = pd.to_numeric(meta[age_col], errors='coerce')
    meta = meta[
        meta['SampleID'].notna() & meta['BioProject'].notna() & meta['Site_Cleaned'].notna() &
        meta['age_numeric'].notna() & (meta['age_numeric'] >= 0) & (meta['age_numeric'] <= 110)
    ].copy()
    meta['SampleID'] = meta['SampleID'].astype(str)
    meta['BioProject'] = meta['BioProject'].astype(str)

    if 'total_reads' in meta.columns:
        total_reads = pd.to_numeric(meta['total_reads'], errors='coerce')
        meta = meta[total_reads > MIN_TOTAL_READS].copy()

    print('Loading genus relative abundance table...')
    genus = pd.read_csv(GENUS_TABLE, sep='\t', index_col=0)
    common = [s for s in meta['SampleID'].tolist() if s in genus.columns]
    meta = meta.drop_duplicates('SampleID').set_index('SampleID').loc[common].reset_index()
    genus = genus[common]
    print(f'Common age-valid samples: {len(common)}; genera: {genus.shape[0]}')

    run_config = dict(
        metadata=str(METADATA_FILE), biom=str(BIOM_FILE), genus_table=str(GENUS_TABLE),
        sites=SITES, min_site_samples=MIN_SITE_SAMPLES, min_bioprojects=MIN_BIOPROJECTS,
        n_splits=N_SPLITS, n_trees=N_TREES, random_seed=RANDOM_SEED, n_jobs=N_JOBS,
        min_prevalence=MIN_PREVALENCE, min_total_reads=MIN_TOTAL_READS,
        age_bin_size=AGE_BIN_SIZE, max_per_age_bin=MAX_PER_AGE_BIN,
        age_range_rule='evaluate only test samples with train_min_age <= age <= train_max_age',
        batch_correction='none'
    )
    with open(OUTDIR / 'run_config.json', 'w') as f:
        json.dump(run_config, f, indent=2)

    all_predictions = []
    all_stats = []
    all_importance = []
    fold_reports = []

    for site in SITES:
        site_meta = meta[meta['Site_Cleaned'] == site].copy()
        if len(site_meta) < MIN_SITE_SAMPLES:
            print(f'{site}: skip, samples={len(site_meta)} < {MIN_SITE_SAMPLES}')
            continue
        n_bp = site_meta['BioProject'].nunique()
        if n_bp < MIN_BIOPROJECTS:
            print(f'{site}: skip, BioProjects={n_bp} < {MIN_BIOPROJECTS}')
            continue

        site_samples = site_meta['SampleID'].tolist()
        site_genus = genus[site_samples]
        prevalence = (site_genus > 0).sum(axis=1) / site_genus.shape[1]
        keep_features = prevalence[prevalence > MIN_PREVALENCE].index
        site_genus = site_genus.loc[keep_features]
        if site_genus.shape[0] < 5:
            print(f'{site}: skip, too few genera after prevalence filter')
            continue

        X = np.log10(site_genus.T.to_numpy(dtype=float) + EPS)
        y = site_meta['age_numeric'].to_numpy(dtype=float)
        groups = site_meta['BioProject'].to_numpy()
        feature_names = site_genus.index.to_numpy()
        n_splits = min(N_SPLITS, n_bp)
        gkf = GroupKFold(n_splits=n_splits)
        print(f'\n=== {site}: samples={len(site_meta)}, BioProjects={n_bp}, genera={len(feature_names)}, folds={n_splits} ===')

        rng_site = np.random.default_rng(RANDOM_SEED + abs(hash(site)) % 100000)
        for fold_id, (train_idx, test_idx) in enumerate(gkf.split(X, y, groups), start=1):
            train_meta = site_meta.iloc[train_idx].copy()
            test_meta = site_meta.iloc[test_idx].copy()
            train_min, train_max = float(train_meta['age_numeric'].min()), float(train_meta['age_numeric'].max())
            in_range_mask = (test_meta['age_numeric'] >= train_min) & (test_meta['age_numeric'] <= train_max)
            eval_test_idx = test_idx[in_range_mask.to_numpy()]
            excluded_test_idx = test_idx[~in_range_mask.to_numpy()]

            if len(eval_test_idx) < 5:
                print(f'  fold {fold_id}: skip, in-range test n={len(eval_test_idx)}')
                fold_reports.append(dict(Site=site, Fold=fold_id, TrainN=len(train_idx), TestN=len(test_idx),
                                         EvalN=len(eval_test_idx), ExcludedOutOfRange=len(excluded_test_idx),
                                         TrainMinAge=train_min, TrainMaxAge=train_max, Status='skipped_too_few_eval'))
                continue

            bal_local = balance_train_indices(train_meta, rng_site)
            bal_train_idx = train_idx[bal_local]
            rf = RandomForestRegressor(
                n_estimators=N_TREES,
                max_features='sqrt',
                min_samples_leaf=2,
                n_jobs=N_JOBS,
                random_state=RANDOM_SEED + fold_id,
                oob_score=False
            )
            rf.fit(X[bal_train_idx], y[bal_train_idx])
            pred = rf.predict(X[eval_test_idx])
            y_eval = y[eval_test_idx]
            mad = pred - y_eval
            stats = metrics_for(y_eval, pred)
            stats.update(dict(Site=site, Fold=fold_id, TrainN=len(train_idx), TrainBalancedN=len(bal_train_idx),
                              TestN=len(test_idx), EvalN=len(eval_test_idx), ExcludedOutOfRange=len(excluded_test_idx),
                              TrainBioProjects=train_meta['BioProject'].nunique(),
                              TestBioProjects=test_meta['BioProject'].nunique(),
                              TrainMinAge=train_min, TrainMaxAge=train_max,
                              TestMinAge=float(test_meta['age_numeric'].min()),
                              TestMaxAge=float(test_meta['age_numeric'].max())))
            all_stats.append(stats)
            fold_reports.append({**stats, 'Status':'completed'})

            pred_df = site_meta.iloc[eval_test_idx][['SampleID','BioProject','Site_Cleaned','Geography','Status_Cleaned','Sex','Age_Group','Region_custom']].copy()
            pred_df['Site'] = site
            pred_df['Fold'] = fold_id
            pred_df['ActualAge'] = y_eval
            pred_df['PredictedAge'] = pred
            pred_df['MAD'] = mad
            pred_df['TrainMinAge'] = train_min
            pred_df['TrainMaxAge'] = train_max
            pred_df['AgeInTrainRange'] = True
            all_predictions.append(pred_df)

            if len(excluded_test_idx) > 0:
                excl = site_meta.iloc[excluded_test_idx][['SampleID','BioProject','Site_Cleaned','Geography','Status_Cleaned','Sex','Age_Group','Region_custom']].copy()
                excl['Site'] = site
                excl['Fold'] = fold_id
                excl['ActualAge'] = y[excluded_test_idx]
                excl['TrainMinAge'] = train_min
                excl['TrainMaxAge'] = train_max
                excl['Reason'] = 'age_outside_training_range'
                excl.to_csv(OUTDIR / f'excluded_out_of_range_{site.replace(" ", "_")}_fold{fold_id}.csv', index=False)

            imp = pd.DataFrame({'Feature': feature_names, 'Importance': rf.feature_importances_, 'Site': site, 'Fold': fold_id})
            all_importance.append(imp)
            print(f"  fold {fold_id}: EvalN={len(eval_test_idx)}, excluded={len(excluded_test_idx)}, MAE={stats['MAE']:.2f}, R2={stats['R2']:.3f}")

    if all_predictions:
        pred_all = pd.concat(all_predictions, ignore_index=True)
        pred_all.to_csv(OUTDIR / 'age_rf_grouped_cv_predictions.csv', index=False)
        stats_all = pd.DataFrame(all_stats)
        stats_all.to_csv(OUTDIR / 'age_rf_grouped_cv_fold_stats.csv', index=False)
        summary = stats_all.groupby('Site').agg(
            Folds=('Fold','nunique'),
            EvalN=('EvalN','sum'),
            ExcludedOutOfRange=('ExcludedOutOfRange','sum'),
            MeanMAE=('MAE','mean'), SdMAE=('MAE','std'),
            MeanRMSE=('RMSE','mean'), MeanR2=('R2','mean'),
            MeanPearsonR=('PearsonR','mean'), MeanSpearmanR=('SpearmanR','mean')
        ).reset_index()
        summary.to_csv(OUTDIR / 'age_rf_grouped_cv_summary.csv', index=False)
        imp_all = pd.concat(all_importance, ignore_index=True)
        imp_summary = (imp_all.groupby(['Site','Feature'], as_index=False)['Importance'].mean()
                       .sort_values(['Site','Importance'], ascending=[True,False]))
        imp_summary.to_csv(OUTDIR / 'age_rf_grouped_cv_feature_importance.csv', index=False)
        pd.DataFrame(fold_reports).to_csv(OUTDIR / 'age_rf_grouped_cv_fold_report.csv', index=False)
        print('\nSummary:')
        print(summary.to_string(index=False))
        print(f'\nDone. Output: {OUTDIR}')
    else:
        pd.DataFrame(fold_reports).to_csv(OUTDIR / 'age_rf_grouped_cv_fold_report.csv', index=False)
        print('No predictions generated.')


if __name__ == '__main__':
    main()
