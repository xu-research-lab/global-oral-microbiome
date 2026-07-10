#!/usr/bin/env python3
"""
Random Forest v2: Sex Prediction with BioProject Leave-One-Out
Updated to use native metadata column names.
"""
import biom
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

BASE = Path("/beegfs/huangxiaochang/global_oral/new_analysis")
OUTDIR = BASE / "randomforest_v2"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Native metadata column names
METADATA = BASE / "metadata_oral_qc_dedup.tsv"
BIOM_PATH = BASE / "biom_oral_qc_dedup.biom"

# Site map
SITE_MAP = {
    "saliva": "Saliva", "tongue": "Tongue", "tooth": "Tooth",
    "oral swab": "Oral_swab", "oral mucosa": "Oral_Mucosa",
    "gingiva": "Gingiva", "other": "Others",
}

REGIONS = ["North America", "Asia", "Europe", "Oceania", "South America", "Africa"]
SITES = ["Saliva", "Tongue", "Tooth"]
MIN_BP_SAMPLES = 10

# ============================================================================
# Load metadata
# ============================================================================
print("Loading metadata...")
meta = pd.read_csv(METADATA, sep="\t", low_memory=False)

# Apply site map
meta["site_cleaned"] = meta["site_cleaned"].replace(SITE_MAP)

# Filter
meta = meta[
    meta["continent_group"].isin(REGIONS)
    & meta["site_cleaned"].isin(SITES)
    & meta["sex_cleaned"].isin(["F", "M"])
    & meta["bioproject_id"].notna()
    & meta["sample_id"].notna()
].copy()
meta["sample_id"] = meta["sample_id"].astype(str)
print(f"Filtered metadata: {len(meta)} samples")

# ============================================================================
# Load BIOM
# ============================================================================
print("Loading BIOM table...")
table = biom.load_table(str(BIOM_PATH))
table.filter(meta["sample_id"].values, axis="sample", inplace=True)
table.remove_empty(inplace=True)
table = table.norm(axis="sample", inplace=False)

sample_ids = np.array(table.ids(axis="sample"))
sample_to_idx = {sid: i for i, sid in enumerate(sample_ids)}
X_all = table.matrix_data.transpose().tocsr()

# ============================================================================
# Run per site × region
# ============================================================================
records = []
summary = []

for site in SITES:
    site_meta = meta[meta["site_cleaned"] == site].copy()
    site_meta = site_meta[site_meta["sample_id"].isin(sample_to_idx)]
    print(f"\n=== {site} ===")

    for region in REGIONS:
        region_meta = site_meta[site_meta["continent_group"] == region].copy()
        if region_meta.empty:
            continue
        bp_counts = region_meta["bioproject_id"].value_counts()
        valid_bps = []
        for bp, n in bp_counts.items():
            if n < MIN_BP_SAMPLES:
                continue
            bp_meta = region_meta[region_meta["bioproject_id"] == bp]
            if bp_meta["sex_cleaned"].nunique() == 2:
                valid_bps.append(bp)
        if len(valid_bps) < 2:
            print(f"  {region}: skip, <2 valid BioProjects")
            continue

        region_meta = region_meta[region_meta["bioproject_id"].isin(valid_bps)]
        aucs = []
        for test_bp in valid_bps:
            test_meta = region_meta[region_meta["bioproject_id"] == test_bp]
            train_meta = region_meta[region_meta["bioproject_id"] != test_bp]
            y_train = (train_meta["sex_cleaned"] == "F").astype(int).to_numpy()
            y_test = (test_meta["sex_cleaned"] == "F").astype(int).to_numpy()
            if np.unique(y_train).size < 2 or np.unique(y_test).size < 2:
                continue

            train_idx = [sample_to_idx[sid] for sid in train_meta["sample_id"]]
            test_idx = [sample_to_idx[sid] for sid in test_meta["sample_id"]]
            X_train = X_all[train_idx]
            X_test = X_all[test_idx]

            model = RandomForestClassifier(
                n_estimators=200, random_state=42,
                class_weight="balanced", n_jobs=8,
            )
            model.fit(X_train, y_train)
            y_prob = model.predict_proba(X_test)[:, 1]
            y_pred = model.predict(X_test)
            auc = roc_auc_score(y_test, y_prob)
            aucs.append(auc)

            fold_df = test_meta[["sample_id", "bioproject_id", "continent_group", "sex_cleaned"]].copy()
            fold_df["site"] = site
            fold_df["fold"] = test_bp
            fold_df["true_label"] = y_test
            fold_df["pred_label"] = y_pred
            fold_df["pred_prob"] = y_prob
            records.append(fold_df)

        if aucs:
            summary.append({
                "site": site, "region": region,
                "n_bioproject_folds": len(aucs),
                "mean_auc": float(np.mean(aucs)),
                "median_auc": float(np.median(aucs)),
                "min_auc": float(np.min(aucs)),
                "max_auc": float(np.max(aucs)),
            })
            print(f"  {region}: mean AUC={np.mean(aucs):.4f} over {len(aucs)} BioProjects")

if records:
    pd.concat(records, ignore_index=True).to_csv(OUTDIR / "sex_bp_loo_predictions.csv", index=False)
if summary:
    pd.DataFrame(summary).sort_values(["site", "mean_auc"], ascending=[True, False]).to_csv(
        OUTDIR / "sex_bp_loo_summary.csv", index=False)

print("\nDone.")
