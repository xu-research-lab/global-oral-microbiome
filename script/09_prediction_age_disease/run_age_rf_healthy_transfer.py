#!/usr/bin/env python3
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from scipy.stats import pearsonr, spearmanr

BASE = Path('/beegfs/huangxiaochang/global_oral/new_analysis')
OUTDIR = BASE / 'age_prediction' / 'healthy_trained_transfer'
OUTDIR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(BASE))
from analysis_inputs import load_metadata

GENUS_TABLE = BASE / 'maaslin3' / 'genus_table_relabund.tsv'
UPSET = BASE / 'maaslin3' / 'disease_shared_plots' / 'upset_set_sizes.csv'
SITES = ['Saliva', 'Tooth', 'Tongue', 'Oral Mucosa', 'Oral_swab']
N_SPLITS = 5
N_TREES = 500
N_JOBS = 8
MIN_HEALTHY_SAMPLES = 100
MIN_HEALTHY_BIOPROJECTS = 5
MIN_PREVALENCE = 0.05
MIN_TOTAL_READS = 1000
AGE_BIN_SIZE = 5
MAX_PER_AGE_BIN = 500
EPS = 1e-6
SEED = 20260610
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


def balance_indices(meta, rng):
    ages = meta['age_numeric'].to_numpy()
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


def corr(y, p, func):
    if len(y) < 3 or np.std(y) == 0 or np.std(p) == 0:
        return np.nan
    try:
        return float(func(y, p)[0])
    except Exception:
        return np.nan


def metrics(y, p):
    return {
        'MAE': mean_absolute_error(y, p),
        'RMSE': mean_squared_error(y, p, squared=False),
        'R2': r2_score(y, p) if len(np.unique(y)) > 1 else np.nan,
        'PearsonR': corr(y, p, pearsonr),
        'SpearmanR': corr(y, p, spearmanr),
    }


def main():
    disease_labels = pd.read_csv(UPSET)['Disease'].tolist()
    disease_statuses = [label_to_status.get(x, x) for x in disease_labels]
    keep_status = ['Healthy Control'] + disease_statuses

    print('Loading metadata...')
    meta = load_metadata(dedup=True)
    meta['age_numeric'] = pd.to_numeric(meta.get('age_cleaned', meta.get('age')), errors='coerce')
    meta = meta[
        meta['SampleID'].notna() & meta['BioProject'].notna() & meta['Site_Cleaned'].notna() &
        meta['Status_Cleaned'].isin(keep_status) & meta['age_numeric'].notna() &
        (meta['age_numeric'] >= 0) & (meta['age_numeric'] <= 110)
    ].copy()
    if 'total_reads' in meta.columns:
        meta = meta[pd.to_numeric(meta['total_reads'], errors='coerce') > MIN_TOTAL_READS].copy()
    meta['SampleID'] = meta['SampleID'].astype(str)
    meta['BioProject'] = meta['BioProject'].astype(str)

    print('Loading genus table...')
    genus = pd.read_csv(GENUS_TABLE, sep='\t', index_col=0)
    common = [s for s in meta['SampleID'].tolist() if s in genus.columns]
    meta = meta.drop_duplicates('SampleID').set_index('SampleID').loc[common].reset_index()
    genus = genus[common]
    print(f'Common samples with age/status: {len(meta)}; genera={genus.shape[0]}')

    with open(OUTDIR / 'run_config.json', 'w') as f:
        json.dump({
            'method': 'Healthy Control trained RF, BioProject grouped CV for healthy out-of-fold and disease transfer prediction',
            'batch_correction': 'none', 'sites': SITES, 'diseases': disease_labels,
            'n_splits': N_SPLITS, 'n_trees': N_TREES, 'min_prevalence': MIN_PREVALENCE,
            'age_range_rule': 'disease and healthy test samples evaluated only if age lies within healthy training fold age range'
        }, f, indent=2)

    healthy_preds=[]; disease_preds=[]; fold_stats=[]; importances=[]; excluded=[]
    for site in SITES:
        site_meta = meta[meta['Site_Cleaned'] == site].copy()
        healthy_meta = site_meta[site_meta['Status_Cleaned'] == 'Healthy Control'].copy()
        disease_meta = site_meta[site_meta['Status_Cleaned'].isin(disease_statuses)].copy()
        if len(healthy_meta) < MIN_HEALTHY_SAMPLES or healthy_meta['BioProject'].nunique() < MIN_HEALTHY_BIOPROJECTS:
            print(f'{site}: skip, healthy n={len(healthy_meta)}, healthy BioProjects={healthy_meta.BioProject.nunique()}')
            continue
        samples = site_meta['SampleID'].tolist()
        site_genus = genus[samples]
        prevalence = (site_genus[healthy_meta['SampleID']] > 0).sum(axis=1) / len(healthy_meta)
        keep_features = prevalence[prevalence > MIN_PREVALENCE].index
        site_genus = site_genus.loc[keep_features]
        if site_genus.shape[0] < 5:
            print(f'{site}: skip, too few healthy-prevalent genera')
            continue
        X_all = pd.DataFrame(np.log10(site_genus.T.to_numpy(float) + EPS), index=site_genus.columns, columns=site_genus.index)
        h = healthy_meta.set_index('SampleID')
        Xh = X_all.loc[h.index].to_numpy()
        yh = h['age_numeric'].to_numpy(float)
        groups = h['BioProject'].to_numpy()
        n_splits = min(N_SPLITS, h['BioProject'].nunique())
        gkf = GroupKFold(n_splits=n_splits)
        rng = np.random.default_rng(SEED + abs(hash(site)) % 100000)
        print(f'\n=== {site}: healthy={len(h)}, healthy BP={h.BioProject.nunique()}, disease={len(disease_meta)}, genera={site_genus.shape[0]} ===')
        site_disease_seen = set()
        for fold, (train_idx, test_idx) in enumerate(gkf.split(Xh, yh, groups), start=1):
            train_meta = h.iloc[train_idx].copy()
            test_meta = h.iloc[test_idx].copy()
            train_min, train_max = float(train_meta.age_numeric.min()), float(train_meta.age_numeric.max())
            bal = balance_indices(train_meta, rng)
            bal_train_idx = train_idx[bal]
            rf = RandomForestRegressor(n_estimators=N_TREES, max_features='sqrt', min_samples_leaf=2,
                                       n_jobs=N_JOBS, random_state=SEED + fold)
            rf.fit(Xh[bal_train_idx], yh[bal_train_idx])
            # healthy out-of-fold, age-range controlled
            in_range = (test_meta.age_numeric >= train_min) & (test_meta.age_numeric <= train_max)
            eval_test = test_meta[in_range].copy()
            excl_h = test_meta[~in_range].copy()
            if len(eval_test) > 0:
                ph = rf.predict(X_all.loc[eval_test.index].to_numpy())
                ytrue = eval_test.age_numeric.to_numpy(float)
                tmp = eval_test.reset_index().rename(columns={'index':'SampleID'})
                tmp['Site'] = site; tmp['Fold'] = fold; tmp['ActualAge'] = ytrue; tmp['PredictedAge'] = ph
                tmp['MAD'] = ph - ytrue; tmp['TrainMinAge'] = train_min; tmp['TrainMaxAge'] = train_max
                tmp['PredictionSet'] = 'Healthy_CV'
                healthy_preds.append(tmp)
                st = metrics(ytrue, ph)
            else:
                st = {'MAE':np.nan,'RMSE':np.nan,'R2':np.nan,'PearsonR':np.nan,'SpearmanR':np.nan}
            st.update({'Site':site,'Fold':fold,'HealthyTrainN':len(train_idx),'HealthyTrainBalancedN':len(bal_train_idx),
                       'HealthyTestN':len(test_idx),'HealthyEvalN':len(eval_test),'HealthyExcludedOutOfRange':len(excl_h),
                       'TrainMinAge':train_min,'TrainMaxAge':train_max})
            fold_stats.append(st)
            if len(excl_h):
                e = excl_h.reset_index().rename(columns={'index':'SampleID'})
                e['Site'] = site; e['Fold'] = fold; e['Reason'] = 'healthy_age_outside_training_range'
                excluded.append(e)
            # disease transfer prediction with this healthy fold model
            d = disease_meta.copy()
            d = d[(d.age_numeric >= train_min) & (d.age_numeric <= train_max)]
            if len(d) > 0:
                Xd = X_all.loc[d.SampleID].to_numpy()
                pdisease = rf.predict(Xd)
                td = d.copy()
                td['Site'] = site; td['Fold'] = fold; td['ActualAge'] = td.age_numeric.to_numpy(float)
                td['PredictedAge'] = pdisease; td['MAD'] = td['PredictedAge'] - td['ActualAge']
                td['TrainMinAge'] = train_min; td['TrainMaxAge'] = train_max; td['PredictionSet'] = 'Disease_Transfer'
                td['DiseaseLabel'] = td['Status_Cleaned'].map(status_to_label).fillna(td['Status_Cleaned'])
                disease_preds.append(td)
                site_disease_seen.update(td['Status_Cleaned'].unique().tolist())
            imp = pd.DataFrame({'Site':site,'Fold':fold,'Feature':site_genus.index,'Importance':rf.feature_importances_})
            importances.append(imp)
            print(f"  fold {fold}: healthy_eval={st['HealthyEvalN']}, disease_transfer={len(d)}, MAE={st['MAE']:.2f}")
        print('  disease statuses predicted:', ', '.join(sorted(site_disease_seen)) if site_disease_seen else 'none')

    if healthy_preds:
        healthy = pd.concat(healthy_preds, ignore_index=True)
        healthy.to_csv(OUTDIR / 'healthy_cv_predictions.csv', index=False)
        pd.DataFrame(fold_stats).to_csv(OUTDIR / 'healthy_cv_fold_stats.csv', index=False)
        hs = pd.DataFrame(fold_stats).groupby('Site').agg(
            Folds=('Fold','nunique'), HealthyEvalN=('HealthyEvalN','sum'), HealthyExcludedOutOfRange=('HealthyExcludedOutOfRange','sum'),
            MeanMAE=('MAE','mean'), MeanRMSE=('RMSE','mean'), MeanR2=('R2','mean'), MeanPearsonR=('PearsonR','mean'), MeanSpearmanR=('SpearmanR','mean')
        ).reset_index()
        hs.to_csv(OUTDIR / 'healthy_cv_summary.csv', index=False)
    if disease_preds:
        disease = pd.concat(disease_preds, ignore_index=True)
        disease.to_csv(OUTDIR / 'disease_transfer_predictions_by_fold.csv', index=False)
        # average disease predictions across folds per sample to avoid repeated rows in final MAD summaries
        key_cols = ['SampleID','BioProject','Site_Cleaned','Site','Geography','Status_Cleaned','DiseaseLabel','Sex','Age_Group','Region_custom','age_numeric']
        disease_avg = disease.groupby(key_cols, dropna=False).agg(
            ActualAge=('ActualAge','first'), PredictedAge=('PredictedAge','mean'), MAD=('MAD','mean'),
            TransferFolds=('Fold','nunique'), TrainMinAge=('TrainMinAge','min'), TrainMaxAge=('TrainMaxAge','max')
        ).reset_index()
        disease_avg.to_csv(OUTDIR / 'disease_transfer_predictions_sample_mean.csv', index=False)
    if importances:
        imp = pd.concat(importances, ignore_index=True)
        imp.groupby(['Site','Feature'], as_index=False)['Importance'].mean().sort_values(['Site','Importance'], ascending=[True,False]).to_csv(
            OUTDIR / 'healthy_model_feature_importance.csv', index=False)
    if excluded:
        pd.concat(excluded, ignore_index=True).to_csv(OUTDIR / 'excluded_out_of_range.csv', index=False)
    print('\nDone. Output:', OUTDIR)

if __name__ == '__main__':
    main()
