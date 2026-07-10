#!/usr/bin/env python3
"""
Random Forest v2: Site One-vs-All prediction.
Target: site_cleaned (5 categories: Saliva, Tooth, Tongue, Oral_swab, Oral_Mucosa)
CV: BioProject-stratified fold split
"""
import biom
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

BASE = Path("/beegfs/huangxiaochang/global_oral/new_analysis")
OUTDIR = BASE / "randomforest_v2"
METADATA = BASE / "metadata_oral_qc_dedup.tsv"
BIOM_PATH = BASE / "biom_oral_qc_dedup.biom"

SITE_MAP = {
    "saliva": "Saliva", "tongue": "Tongue", "tooth": "Tooth",
    "oral swab": "Oral_swab", "oral mucosa": "Oral_Mucosa",
    "gingiva": "Gingiva", "other": "Others",
}

TARGET_SITES = ["Saliva", "Tooth", "Tongue", "Oral_swab", "Oral_Mucosa"]
N_REPEATS = 5
MIN_STUDY_SAMPLES = 20
MIN_POS_BP = 3


def stratified_fold_split(bp_pos, bp_neg, n_folds, rng):
    pos, neg = list(bp_pos), list(bp_neg)
    rng.shuffle(pos); rng.shuffle(neg)
    fold_bps = [[] for _ in range(n_folds)]
    for i, bp in enumerate(pos): fold_bps[i % n_folds].append(bp)
    for i, bp in enumerate(neg): fold_bps[i % n_folds].append(bp)
    for i in range(n_folds):
        test_bps = fold_bps[i]
        train_bps = [bp for j, fold in enumerate(fold_bps) if j != i for bp in fold]
        yield train_bps, test_bps


# ============================================================================
# Load metadata
# ============================================================================
print("Loading metadata...")
meta = pd.read_csv(METADATA, sep="\t", low_memory=False)
meta["site_cleaned"] = meta["site_cleaned"].replace(SITE_MAP)

meta = meta[
    meta["site_cleaned"].isin(TARGET_SITES)
    & meta["bioproject_id"].notna()
    & meta["sample_id"].notna()
].copy()
meta["sample_id"] = meta["sample_id"].astype(str)
print(f"Samples in target sites: {len(meta)}")

# Filter to BioProjects with enough samples
study_sizes = (
    meta.groupby(["site_cleaned", "bioproject_id"])
    .size().reset_index(name="n")
)
valid_bps = study_sizes[study_sizes["n"] >= MIN_STUDY_SAMPLES][
    ["bioproject_id"]
].drop_duplicates()
meta = meta[meta["bioproject_id"].isin(valid_bps["bioproject_id"])]
print(f"After filter: {len(meta)} samples from {len(valid_bps)} BioProjects")

# Site distribution
print("\nSite distribution:")
for s in TARGET_SITES:
    n = (meta["site_cleaned"] == s).sum()
    n_bp = meta[meta["site_cleaned"] == s]["bioproject_id"].nunique()
    print(f"  {s}: {n} samples, {n_bp} BioProjects")

# ============================================================================
# Load BIOM
# ============================================================================
print("\nLoading BIOM table...")
table = biom.load_table(str(BIOM_PATH))
table.filter(meta["sample_id"].values, axis="sample", inplace=True)
table.remove_empty(inplace=True)
table = table.norm(axis="sample", inplace=False)

sample_ids = np.array(table.ids(axis="sample"))
sample_to_idx = {sid: i for i, sid in enumerate(sample_ids)}
X_all = table.matrix_data.transpose().tocsr()
print(f"Features: {X_all.shape[1]}")

# ============================================================================
# Run OVA
# ============================================================================
records = []
summary = []

for target_site in TARGET_SITES:
    pos_bp = meta.loc[meta["site_cleaned"] == target_site, "bioproject_id"].dropna().unique()
    neg_bp = meta.loc[meta["site_cleaned"] != target_site, "bioproject_id"].dropna().unique()
    if len(pos_bp) < MIN_POS_BP or len(neg_bp) < MIN_POS_BP:
        print(f"\n{target_site}: skip (pos_BP={len(pos_bp)}, neg_BP={len(neg_bp)})")
        continue

    n_folds = min(5, len(pos_bp), len(neg_bp))
    print(f"\n=== {target_site}: pos_BP={len(pos_bp)}, neg_BP={len(neg_bp)}, folds={n_folds} ===")

    eval_rows = []
    for repeat in range(N_REPEATS):
        rng = np.random.RandomState(42 + repeat)
        for fold_i, (train_bps, test_bps) in enumerate(
            stratified_fold_split(pos_bp, neg_bp, n_folds, rng)):
            train_meta = meta[meta["bioproject_id"].isin(train_bps)]
            test_meta = meta[meta["bioproject_id"].isin(test_bps)]
            if train_meta.empty or test_meta.empty:
                continue

            y_train = (train_meta["site_cleaned"] == target_site).astype(int).to_numpy()
            y_test = (test_meta["site_cleaned"] == target_site).astype(int).to_numpy()
            if np.unique(y_train).size < 2 or np.unique(y_test).size < 2:
                continue

            train_idx = [sample_to_idx[sid] for sid in train_meta["sample_id"] if sid in sample_to_idx]
            test_idx = [sample_to_idx[sid] for sid in test_meta["sample_id"] if sid in sample_to_idx]
            if not train_idx or not test_idx:
                continue

            model = RandomForestClassifier(
                n_estimators=200, random_state=42 + repeat,
                class_weight="balanced", n_jobs=8)
            model.fit(X_all[train_idx], y_train)
            y_prob = model.predict_proba(X_all[test_idx])[:, 1]
            y_pred = model.predict(X_all[test_idx])
            auc = roc_auc_score(y_test, y_prob)
            fold_name = f"repeat{repeat}_fold{fold_i}"
            eval_rows.append((fold_name, auc))

            fold_df = test_meta[["sample_id", "bioproject_id", "site_cleaned"]].copy()
            fold_df["target_site"] = target_site
            fold_df["fold"] = fold_name
            fold_df["true_label"] = y_test
            fold_df["pred_label"] = y_pred
            fold_df["pred_prob"] = y_prob
            records.append(fold_df)

    if eval_rows:
        aucs = np.array([x[1] for x in eval_rows], dtype=float)
        summary.append({
            "target_site": target_site,
            "n_folds": len(eval_rows),
            "mean_auc": float(np.mean(aucs)),
            "median_auc": float(np.median(aucs)),
            "min_auc": float(np.min(aucs)),
            "max_auc": float(np.max(aucs)),
            "positive_bioprojects": len(pos_bp),
            "negative_bioprojects": len(neg_bp),
        })
        print(f"  mean AUC={np.mean(aucs):.4f} over {len(eval_rows)} folds")

if records:
    pd.concat(records, ignore_index=True).to_csv(OUTDIR / "site_ova_predictions.csv", index=False)
if summary:
    pd.DataFrame(summary).sort_values("mean_auc", ascending=False).to_csv(
        OUTDIR / "site_ova_summary.csv", index=False)
print("\nDone.")
