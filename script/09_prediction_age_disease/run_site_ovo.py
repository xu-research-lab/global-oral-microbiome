#!/usr/bin/env python3
"""
Random Forest v2: Site One-vs-One prediction.
Target: site_cleaned pairs among Saliva, Tooth, Tongue, Oral_swab, Oral_Mucosa
CV: BioProject-stratified fold split
"""
import itertools
from pathlib import Path
import biom
import numpy as np
import pandas as pd
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
PAIRS = list(itertools.combinations(TARGET_SITES, 2))
N_REPEATS = 5
MIN_STUDY_SAMPLES = 20
MIN_BP = 3


def stratified_fold_split(bp_1, bp_2, n_folds, rng):
    b1, b2 = list(bp_1), list(bp_2)
    rng.shuffle(b1); rng.shuffle(b2)
    fold_bps = [[] for _ in range(n_folds)]
    for i, bp in enumerate(b1): fold_bps[i % n_folds].append(bp)
    for i, bp in enumerate(b2): fold_bps[i % n_folds].append(bp)
    for i in range(n_folds):
        test_bps = fold_bps[i]
        train_bps = [bp for j, fold in enumerate(fold_bps) if j != i for bp in fold]
        yield train_bps, test_bps


# ============================================================================
# Load
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

study_sizes = (
    meta.groupby(["site_cleaned", "bioproject_id"])
    .size().reset_index(name="n")
)
valid_bps = study_sizes[study_sizes["n"] >= MIN_STUDY_SAMPLES][
    ["bioproject_id"]
].drop_duplicates()
meta = meta[meta["bioproject_id"].isin(valid_bps["bioproject_id"])]
print(f"Samples: {len(meta)}, BioProjects: {len(valid_bps)}")

print("Loading BIOM table...")
table = biom.load_table(str(BIOM_PATH))
table.filter(meta["sample_id"].values, axis="sample", inplace=True)
table.remove_empty(inplace=True)
table = table.norm(axis="sample", inplace=False)

sample_ids = np.array(table.ids(axis="sample"))
sample_to_idx = {sid: i for i, sid in enumerate(sample_ids)}
X_all = table.matrix_data.transpose().tocsr()

# ============================================================================
# Run OVO
# ============================================================================
records = []
summary = []

for s1, s2 in PAIRS:
    pair_meta = meta[meta["site_cleaned"].isin([s1, s2])].copy()
    bp_1 = pair_meta.loc[pair_meta["site_cleaned"] == s1, "bioproject_id"].dropna().unique()
    bp_2 = pair_meta.loc[pair_meta["site_cleaned"] == s2, "bioproject_id"].dropna().unique()
    if len(bp_1) < MIN_BP or len(bp_2) < MIN_BP:
        print(f"  {s1} vs {s2}: skip ({len(bp_1)} vs {len(bp_2)} BioProjects)")
        continue

    n_folds = min(5, len(bp_1), len(bp_2))
    print(f"\n=== {s1} vs {s2}: BP {len(bp_1)}/{len(bp_2)}, folds={n_folds} ===")

    eval_rows = []
    for repeat in range(N_REPEATS):
        rng = np.random.RandomState(42 + repeat)
        for fold_i, (train_bps, test_bps) in enumerate(
            stratified_fold_split(bp_1, bp_2, n_folds, rng)):
            train_meta = pair_meta[pair_meta["bioproject_id"].isin(train_bps)]
            test_meta = pair_meta[pair_meta["bioproject_id"].isin(test_bps)]
            if train_meta.empty or test_meta.empty:
                continue

            y_train = (train_meta["site_cleaned"] == s1).astype(int).to_numpy()
            y_test = (test_meta["site_cleaned"] == s1).astype(int).to_numpy()
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
            fold_df["site_1"] = s1
            fold_df["site_2"] = s2
            fold_df["fold"] = fold_name
            fold_df["true_label"] = y_test
            fold_df["pred_label"] = y_pred
            fold_df["pred_prob"] = y_prob
            records.append(fold_df)

    if eval_rows:
        aucs = np.array([x[1] for x in eval_rows], dtype=float)
        summary.append({
            "site_1": s1, "site_2": s2,
            "n_folds": len(eval_rows),
            "mean_auc": float(np.mean(aucs)), "median_auc": float(np.median(aucs)),
            "min_auc": float(np.min(aucs)), "max_auc": float(np.max(aucs)),
            "n_bioproject_site_1": len(bp_1), "n_bioproject_site_2": len(bp_2),
        })
        print(f"  mean AUC={np.mean(aucs):.4f} over {len(eval_rows)} folds")

if records:
    pd.concat(records, ignore_index=True).to_csv(OUTDIR / "site_ovo_predictions.csv", index=False)
if summary:
    pd.DataFrame(summary).sort_values("mean_auc", ascending=False).to_csv(
        OUTDIR / "site_ovo_summary.csv", index=False)
print("\nDone.")
