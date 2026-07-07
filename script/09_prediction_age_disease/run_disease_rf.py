#!/usr/bin/env python3

import biom
import numpy as np
import pandas as pd
from pathlib import Path
import sys
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

BASE = Path("/beegfs/huangxiaochang/global_oral/new_analysis")
OUTDIR = BASE / "disease_prediction"
OUTDIR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(BASE))

from analysis_inputs import load_biom_table, load_metadata

DISEASES = ["Periodontal Disease", "Dental Caries", "HIV", "Head and Neck Cancer"]

print("Loading metadata and biom table...")
meta = load_metadata(dedup=True)
meta = meta[meta["Status_Cleaned"].isin(["Healthy Control"] + DISEASES)].copy()
meta = meta[meta["BioProject"].notna() & meta["SampleID"].notna()]
meta["SampleID"] = meta["SampleID"].astype(str)

table = load_biom_table(dedup=True)
table.filter(meta["SampleID"].values, axis="sample", inplace=True)
table.remove_empty(inplace=True)
table = table.norm(axis="sample", inplace=False)

sample_ids = np.array(table.ids(axis="sample"))
sample_to_idx = {sid: i for i, sid in enumerate(sample_ids)}
X_all = table.matrix_data.transpose().tocsr()

records = []
summary = []

for disease in DISEASES:
    disease_meta = meta[meta["Status_Cleaned"].isin(["Healthy Control", disease])].copy()
    ctab = pd.crosstab(disease_meta["BioProject"], disease_meta["Status_Cleaned"])
    if disease not in ctab.columns or "Healthy Control" not in ctab.columns:
        continue
    valid_bps = ctab.index[(ctab[disease] > 0) & (ctab["Healthy Control"] > 0)].tolist()
    disease_meta = disease_meta[disease_meta["BioProject"].isin(valid_bps)].copy()
    disease_meta = disease_meta[disease_meta["SampleID"].isin(sample_to_idx)]

    if len(valid_bps) < 3:
        print(f"{disease}: skip, only {len(valid_bps)} matched BioProjects")
        continue

    aucs = []
    print(f"\n=== {disease} ({len(valid_bps)} matched BioProjects) ===")
    for test_bp in valid_bps:
        test_meta = disease_meta[disease_meta["BioProject"] == test_bp]
        train_meta = disease_meta[disease_meta["BioProject"] != test_bp]
        y_train = (train_meta["Status_Cleaned"] == disease).astype(int).to_numpy()
        y_test = (test_meta["Status_Cleaned"] == disease).astype(int).to_numpy()
        if np.unique(y_train).size < 2 or np.unique(y_test).size < 2:
            continue

        train_idx = [sample_to_idx[sid] for sid in train_meta["SampleID"]]
        test_idx = [sample_to_idx[sid] for sid in test_meta["SampleID"]]
        X_train = X_all[train_idx]
        X_test = X_all[test_idx]

        model = RandomForestClassifier(
            n_estimators=300,
            random_state=42,
            class_weight="balanced",
            n_jobs=8,
        )
        model.fit(X_train, y_train)
        y_prob = model.predict_proba(X_test)[:, 1]
        y_pred = model.predict(X_test)
        auc = roc_auc_score(y_test, y_prob)
        aucs.append(auc)

        fold_df = test_meta[["SampleID", "BioProject", "Status_Cleaned", "Site_Cleaned", "Geography"]].copy()
        fold_df["disease"] = disease
        fold_df["fold"] = test_bp
        fold_df["true_label"] = y_test
        fold_df["pred_label"] = y_pred
        fold_df["pred_prob"] = y_prob
        records.append(fold_df)

    if aucs:
        summary.append({
            "disease": disease,
            "matched_bioprojects": len(valid_bps),
            "n_folds": len(aucs),
            "mean_auc": float(np.mean(aucs)),
            "median_auc": float(np.median(aucs)),
            "min_auc": float(np.min(aucs)),
            "max_auc": float(np.max(aucs)),
        })
        print(f"{disease}: mean AUC={np.mean(aucs):.4f} over {len(aucs)} BioProjects")

if records:
    pd.concat(records, ignore_index=True).to_csv(OUTDIR / "disease_rf_predictions.csv", index=False)
if summary:
    pd.DataFrame(summary).sort_values("mean_auc", ascending=False).to_csv(
        OUTDIR / "disease_rf_summary.csv", index=False
    )

print("\nDone.")

