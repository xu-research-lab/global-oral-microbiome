#!/usr/bin/env python
# =====================================================================
# Re-run the leave-studies-out age RF CV that backs the disease
# microbiome-age-acceleration figure (notebook age_predict.ipynb cell 33).
# For each oral site: pool samples (age 1-80) from bioprojects with >10
# samples, split at the STUDY level with KFold(k, shuffle, seed=42),
# train a Random Forest on HC samples of the train studies only, predict
# all test-study samples, and save per-fold predictions as
#   age_predict_test_{fold}_{site}_genus_range_1_80.csv
# The downstream plot script (plot_age_acceleration_1row.py) reads these
# CSVs, recomputes the per-status delta-age-gap stats and draws the
# figure with all sites in one row.
# =====================================================================
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold

# Genus-level count table (TSV: rows=genera, cols=sample ids). This is the
# local equivalent of the notebook's biom_output/feature-table.biom.
BASE = Path(os.environ.get("GLOBAL_ORAL_ROOT") or Path(__file__).resolve().parents[2])
RESULT = Path(__file__).resolve().parents[1] / "result"

GENUS_TABLE_PATH = (BASE / "data" / "metadata-and-biom-v2-qc" /
                     "genus_table_counts.tsv")
METADATA_PATH = (BASE / "data" / "metadata-and-biom-v2-qc" /
                 "metadata_oral_qc_dedup.tsv")
OUTPUT_DIR = RESULT / "leave_studies_out_rf"

AGE_MIN, AGE_MAX = 1, 80
N_TREES = 500
N_SPLITS = 5
EPS = 1e-6
BAND_SITES = ["saliva", "tooth", "tongue", "oral mucosa"]


def load_genus_table(needed_sids=None):
    """Read the genus count TSV once, return (feature_ids, sample_ids, X)
    where X is a numpy array [n_samples, n_features] of counts.

    If needed_sids (a set of str) is given, only those sample columns are
    read — essential because the TSV has ~104k columns."""
    with open(GENUS_TABLE_PATH) as fh:
        header = fh.readline().rstrip("\n")
    all_sids = header.split("\t")[1:]
    if needed_sids is None:
        keep_cols = all_sids
    else:
        keep_cols = [s for s in all_sids if s in needed_sids]
    # usecols selects only the needed columns (C engine scans header to map)
    gf = pd.read_csv(GENUS_TABLE_PATH, sep="\t", index_col=0, usecols=keep_cols)
    feature_ids = gf.index.values
    sample_ids = np.asarray(gf.columns)
    X = gf.values.T.astype(np.float64)  # [n_samples, n_features]
    return feature_ids, sample_ids, X


def build_xy_pool(table_X, table_sids, md_train, md_test, sid_to_age):
    """Train on md_train HC samples only; feature vocabulary from train;
    test reindexed to train columns. Returns log10(relabund+EPS) features."""
    train_sid = md_train.loc[md_train.status_cleaned == "HC", "sample_id"].astype(str).values
    test_sid = md_test["sample_id"].astype(str).values
    sid_idx = {s: i for i, s in enumerate(table_sids)}
    tr_idx = [sid_idx[s] for s in train_sid if s in sid_idx]
    te_idx = [sid_idx[s] for s in test_sid if s in sid_idx]
    Xtr = table_X[tr_idx].copy()
    Xte = table_X[te_idx].copy()
    # relative abundance (normalize per sample to sum 1)
    Xtr = Xtr / Xtr.sum(axis=1, keepdims=True)
    Xte = Xte / Xte.sum(axis=1, keepdims=True)
    # remove features absent in train (vocabulary from TRAIN only)
    keep = Xtr.sum(axis=0) > 0
    Xtr = Xtr[:, keep]
    Xte = Xte[:, keep]
    ytr = sid_to_age.reindex([table_sids[i] for i in tr_idx]).values.astype(float)
    yte = sid_to_age.reindex([table_sids[i] for i in te_idx]).values.astype(float)
    return np.log10(Xtr + EPS), ytr, np.log10(Xte + EPS), yte


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(METADATA_PATH, sep="\t", low_memory=False)
    metadata = metadata.loc[metadata.age_cleaned > 0]
    metadata = metadata.loc[metadata.site_cleaned.notna()]
    metadata = metadata.loc[metadata.site_cleaned != "other"]
    metadata["sample_id"] = metadata["sample_id"].astype(str)

    # Pre-compute the pool of sample ids actually needed (age 1-80, any of
    # the band sites) so we only read those columns from the huge TSV.
    pool_md = metadata.loc[
        metadata.site_cleaned.isin(BAND_SITES)
        & metadata.age_cleaned.between(AGE_MIN, AGE_MAX)
    ]
    needed = set(pool_md.sample_id.values)
    print(f"Loading genus count table (cols={len(needed)})...", flush=True)
    feature_ids, table_sids, table_X = load_genus_table(needed)
    print(f"  features={len(feature_ids)} samples={len(table_sids)}", flush=True)
    sid_in_table = set(table_sids.tolist())

    # notebook cell 20: keep only samples present in the table
    metadata = metadata.loc[metadata.sample_id.isin(sid_in_table)]
    sid_to_age = metadata.set_index("sample_id")["age_cleaned"]

    for site in BAND_SITES:
        md_site = metadata.loc[
            (metadata.site_cleaned == site)
            & (metadata.age_cleaned >= AGE_MIN)
            & (metadata.age_cleaned <= AGE_MAX)
        ]
        counts = md_site.groupby("bioproject_id")["sample_id"].count()
        selected = counts[counts > 10].index.values
        pool = md_site.loc[md_site.bioproject_id.isin(selected)].copy()
        pool["sample_id"] = pool["sample_id"].astype(str)
        studies = pool.bioproject_id.unique()
        k = min(N_SPLITS, len(studies))
        print(f"[{site}] pool_n={len(pool)} studies={len(studies)} folds={k}", flush=True)

        def subset(units):
            return (pool.loc[pool.bioproject_id.isin(units[0])],
                    pool.loc[pool.bioproject_id.isin(units[1])])

        fold = 1
        for tri, tei in KFold(k, shuffle=True, random_state=42).split(studies):
            md_tr, md_te = subset((studies[tri], studies[tei]))
            if md_tr.shape[0] < 5 or md_te.shape[0] < 2:
                continue

            Xtr, ytr, Xte, yte = build_xy_pool(table_X, table_sids,
                                               md_tr, md_te, sid_to_age)

            rf = RandomForestRegressor(n_estimators=N_TREES, random_state=42, n_jobs=-1,
                                       max_features="sqrt", min_samples_leaf=2)
            rf.fit(Xtr, ytr)
            pred = rf.predict(Xte)
            print(f"[{site}] fold {fold}: n_test={len(yte)} "
                  f"spearman={spearmanr(yte, pred)[0]:+.3f}", flush=True)

            pd.DataFrame({
                "sample_id": md_te.sample_id.values,
                "true_age": yte,
                "pred_age": pred,
            }).to_csv(
                OUTPUT_DIR / f"age_predict_test_{fold}_{site}_genus_range_1_80.csv",
                index=False,
            )
            fold += 1

    print("Done.", flush=True)


if __name__ == "__main__":
    main()