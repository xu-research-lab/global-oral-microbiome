#!/usr/bin/env python
# =====================================================================
# Random (sample-level) split version of the age RF CV.
# Same settings as run_age_acceleration_cv.py (same sample pool: site +
# age 1-80 + bioprojects with >10 samples; RF params identical; train on
# HC only) — ONLY the cross-validation split changes: instead of holding
# out whole bioprojects, samples are shuffled and split with KFold(5),
# so samples from the same study appear in both train and test.
# This isolates the batch-effect penalty in the leave-studies-out figure:
#   random-split  = optimistic, studies mixed
#   leave-studies = honest cross-study generalization
# Per-fold predictions are saved as
#   age_predict_test_{fold}_{site}_random_range_1_80.csv
# to data/age_prediction/random_split_rf/.
# The downstream plot script (plot_random_split_scatter.py) reads these
# CSVs and draws the 1x4 true-vs-predicted scatter figure.
# =====================================================================
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold

BASE = Path(os.environ.get("GLOBAL_ORAL_ROOT") or Path(__file__).resolve().parents[2])
RESULT = Path(__file__).resolve().parents[1] / "result"

GENUS_TABLE_PATH = (BASE / "data" / "metadata-and-biom-v2-qc" /
                     "genus_table_counts.tsv")
METADATA_PATH = (BASE / "data" / "metadata-and-biom-v2-qc" /
                 "metadata_oral_qc_dedup.tsv")
OUTPUT_DIR = RESULT / "random_split_rf"

AGE_MIN, AGE_MAX = 1, 80
N_TREES = 500
N_SPLITS = 5          # same fold count for every site (no study-count cap)
EPS = 1e-6
BAND_SITES = ["saliva", "tooth", "tongue", "oral mucosa"]
SEED = 42


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

    metadata = metadata.loc[metadata.sample_id.isin(sid_in_table)]
    sid_to_age = metadata.set_index("sample_id")["age_cleaned"]

    for site in BAND_SITES:
        md_site = metadata.loc[
            (metadata.site_cleaned == site)
            & (metadata.age_cleaned >= AGE_MIN)
            & (metadata.age_cleaned <= AGE_MAX)
        ]
        # same pool rule as the leave-studies-out run so the two versions
        # are evaluated on the same samples (only the split differs)
        counts = md_site.groupby("bioproject_id")["sample_id"].count()
        selected = counts[counts > 10].index.values
        pool = md_site.loc[md_site.bioproject_id.isin(selected)].copy()
        pool["sample_id"] = pool["sample_id"].astype(str).values
        n_studies = pool.bioproject_id.nunique()
        print(f"[{site}] pool_n={len(pool)} studies={n_studies} "
              f"folds={N_SPLITS} (sample-level random split)", flush=True)

        ids = pool.sample_id.values
        fold = 1
        for tri, tei in KFold(N_SPLITS, shuffle=True, random_state=SEED).split(ids):
            md_tr = pool.loc[pool.sample_id.isin(ids[tri])]
            md_te = pool.loc[pool.sample_id.isin(ids[tei])]
            if md_tr.shape[0] < 5 or md_te.shape[0] < 2:
                continue

            Xtr, ytr, Xte, yte = build_xy_pool(table_X, table_sids,
                                               md_tr, md_te, sid_to_age)

            rf = RandomForestRegressor(n_estimators=N_TREES, random_state=SEED,
                                       n_jobs=-1, max_features="sqrt",
                                       min_samples_leaf=2)
            rf.fit(Xtr, ytr)
            pred = rf.predict(Xte)
            print(f"[{site}] fold {fold}: n_test={len(yte)} "
                  f"spearman={spearmanr(yte, pred)[0]:+.3f}", flush=True)

            pd.DataFrame({
                "sample_id": md_te.sample_id.values,
                "true_age": yte,
                "pred_age": pred,
            }).to_csv(
                OUTPUT_DIR / f"age_predict_test_{fold}_{site}_random_range_1_80.csv",
                index=False,
            )
            fold += 1

    print("Done.", flush=True)


if __name__ == "__main__":
    main()