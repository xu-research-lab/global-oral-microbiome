#!/usr/bin/env python3
"""Oral + gut region RF OVA/OVO with BioProject-level 5-fold CV (repeated).

Mirrors the CV design in:
  /beegfs/huangxiaochang/global_oral/new_analysis/randomforest_v2

Key rules:
  - Split by project, never by sample: a project is wholly train or test.
  - Positive/negative projects are stratified separately into folds.
  - n_folds = min(5, n_pos_projects, n_neg_projects).
  - 5 repeats with different project shuffles.
  - Keep only site-region-project groups with >= MIN_STUDY_SAMPLES samples.

Sample-level CV lives in run_oral_gut_sample_cv_repeat5.py.
"""

import os
import json
from pathlib import Path

import biom
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)


ROOT = Path(os.environ.get("GLOBAL_ORAL_ROOT") or Path(__file__).resolve().parents[2])
OUTDIR = (Path(__file__).resolve().parents[1] / "result" / "oral")
ORAL_METADATA = ROOT / "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv"
ORAL_BIOM = ROOT / "data/metadata-and-biom-v2-qc/biom_oral_qc_dedup.biom"
GUT_METADATA = ROOT / "data/metadata_gut_v2.txt"
GUT_BIOM = ROOT / "data/table_filter_reads5000.biom"

REGIONS = [
    "North America",
    "Asia",
    "Europe",
    "Oceania",
    "South America",
    "Africa",
]
ORAL_SITES = ["Saliva", "Tongue", "Tooth"]
N_FOLDS_MAX = 5
N_REPEATS = 5
SEED = 42
N_TREES = 200
N_JOBS = 8
MIN_STUDY_SAMPLES = 20


def stratified_project_fold_split(bp_pos, bp_neg, n_folds, rng):
    """Assign positive/negative projects to folds (project wholly held out)."""
    pos, neg = list(bp_pos), list(bp_neg)
    rng.shuffle(pos)
    rng.shuffle(neg)
    fold_bps = [[] for _ in range(n_folds)]
    for i, bp in enumerate(pos):
        fold_bps[i % n_folds].append(bp)
    for i, bp in enumerate(neg):
        fold_bps[i % n_folds].append(bp)
    for i in range(n_folds):
        test_bps = fold_bps[i]
        train_bps = [
            bp for j, fold in enumerate(fold_bps) if j != i for bp in fold
        ]
        yield train_bps, test_bps


def load_oral():
    meta = pd.read_csv(ORAL_METADATA, sep="\t", low_memory=False)
    required = {"sample_id", "continent_group", "site_cleaned", "bioproject_id"}
    missing = sorted(required.difference(meta.columns))
    if missing:
        raise ValueError("Oral metadata missing columns: " + ", ".join(missing))

    site_map = {"saliva": "Saliva", "tongue": "Tongue", "tooth": "Tooth"}
    meta["sample_id"] = meta["sample_id"].astype("string").str.strip()
    meta["region"] = meta["continent_group"].astype("string").str.strip()
    meta["site"] = (
        meta["site_cleaned"].astype("string").str.strip().replace(site_map)
    )
    meta["project"] = meta["bioproject_id"].astype("string").str.strip()
    meta = meta[
        meta["region"].isin(REGIONS)
        & meta["site"].isin(ORAL_SITES)
        & meta["project"].notna()
        & meta["sample_id"].notna()
    ].copy()
    meta["sample_id"] = meta["sample_id"].astype(str)
    if meta["sample_id"].duplicated().any():
        raise ValueError("Oral metadata contains duplicate sample_id values")

    # Keep site-region-project groups with enough samples (v2 rule).
    study_sizes = (
        meta.groupby(["site", "region", "project"]).size().reset_index(name="n")
    )
    valid = study_sizes[study_sizes["n"] >= MIN_STUDY_SAMPLES][
        ["site", "region", "project"]
    ]
    meta = meta.merge(valid, on=["site", "region", "project"], how="inner")
    print(
        f"Oral keep {len(meta):,} samples from {len(valid):,} "
        f"site-region BioProjects (>= {MIN_STUDY_SAMPLES})"
    )

    table = biom.load_table(str(ORAL_BIOM))
    wanted = set(meta["sample_id"].astype(str))
    shared = [sid for sid in table.ids(axis="sample") if str(sid) in wanted]
    if not shared:
        raise ValueError("No oral sample IDs overlap between metadata and BIOM")
    table = table.filter(shared, axis="sample", inplace=False)
    table.remove_empty(inplace=True)
    table = table.norm(axis="sample", inplace=False)
    sample_ids = np.asarray(table.ids(axis="sample"), dtype=str)
    meta = meta.set_index("sample_id", drop=False).loc[sample_ids].reset_index(drop=True)
    meta["_matrix_index"] = np.arange(len(meta), dtype=int)
    return meta, table.matrix_data.transpose().tocsr()


def load_gut():
    meta = pd.read_csv(GUT_METADATA, sep="\t", low_memory=False)
    required = {"sample_id", "region", "project"}
    missing = sorted(required.difference(meta.columns))
    if missing:
        raise ValueError("Gut metadata missing columns: " + ", ".join(missing))

    meta["sample_id"] = meta["sample_id"].astype("string").str.strip()
    meta["region"] = meta["region"].astype("string").str.strip()
    meta["project"] = meta["project"].astype("string").str.strip()
    meta["site"] = "Gut"
    meta = meta[
        meta["region"].isin(REGIONS)
        & meta["project"].notna()
        & meta["sample_id"].notna()
    ].copy()
    meta["sample_id"] = meta["sample_id"].astype(str)
    if meta["sample_id"].duplicated().any():
        raise ValueError("Gut metadata contains duplicate sample_id values")

    study_sizes = (
        meta.groupby(["site", "region", "project"]).size().reset_index(name="n")
    )
    valid = study_sizes[study_sizes["n"] >= MIN_STUDY_SAMPLES][
        ["site", "region", "project"]
    ]
    meta = meta.merge(valid, on=["site", "region", "project"], how="inner")
    print(
        f"Gut keep {len(meta):,} samples from {len(valid):,} "
        f"region projects (>= {MIN_STUDY_SAMPLES})"
    )

    table = biom.load_table(str(GUT_BIOM))
    wanted = set(meta["sample_id"].astype(str))
    shared = [sid for sid in table.ids(axis="sample") if str(sid) in wanted]
    if not shared:
        raise ValueError("No gut sample IDs overlap between metadata and BIOM")
    table = table.filter(shared, axis="sample", inplace=False)
    table.remove_empty(inplace=True)
    table = table.norm(axis="sample", inplace=False)
    sample_ids = np.asarray(table.ids(axis="sample"), dtype=str)
    meta = meta.set_index("sample_id", drop=False).loc[sample_ids].reset_index(drop=True)
    meta["_matrix_index"] = np.arange(len(meta), dtype=int)
    return meta, table.matrix_data.transpose().tocsr()


def fold_metrics(y_true, y_pred, y_prob):
    return {
        "auc": roc_auc_score(y_true, y_prob),
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "f1": f1_score(y_true, y_pred, zero_division=0),
    }


def summarize(rows, base):
    frame = pd.DataFrame(rows)
    result = dict(base)
    result.update(
        {
            "n_evaluations": len(frame),
            "mean_auc": frame["auc"].mean(),
            "sd_auc": frame["auc"].std(ddof=1),
            "median_auc": frame["auc"].median(),
            "min_auc": frame["auc"].min(),
            "max_auc": frame["auc"].max(),
            "mean_accuracy": frame["accuracy"].mean(),
            "mean_balanced_accuracy": frame["balanced_accuracy"].mean(),
            "mean_precision": frame["precision"].mean(),
            "mean_recall": frame["recall"].mean(),
            "mean_f1": frame["f1"].mean(),
        }
    )
    return result


def append_predictions(frame, path):
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        mode="a",
        header=not path.exists(),
    )


def run_ova(meta, matrix, dataset, prediction_path):
    folds, summaries = [], []
    for site in sorted(
        meta["site"].unique(),
        key=lambda x: ORAL_SITES.index(x) if x in ORAL_SITES else 99,
    ):
        site_meta = meta[meta["site"] == site].copy()
        observed = [
            region for region in REGIONS if region in set(site_meta["region"])
        ]
        print(f"\n=== OVA {dataset}/{site} ===")
        for target in observed:
            pos_bp = (
                site_meta.loc[site_meta["region"] == target, "project"]
                .dropna()
                .unique()
            )
            neg_bp = (
                site_meta.loc[site_meta["region"] != target, "project"]
                .dropna()
                .unique()
            )
            if len(pos_bp) < 2 or len(neg_bp) < 2:
                print(
                    f"SKIP {dataset}/{site}/{target}: "
                    f"projects pos={len(pos_bp)} neg={len(neg_bp)}"
                )
                continue
            n_folds = min(N_FOLDS_MAX, len(pos_bp), len(neg_bp))
            if n_folds < 2:
                continue

            task_rows = []
            for repeat in range(1, N_REPEATS + 1):
                rng = np.random.RandomState(SEED + repeat)
                for fold_i, (train_bps, test_bps) in enumerate(
                    stratified_project_fold_split(pos_bp, neg_bp, n_folds, rng),
                    start=1,
                ):
                    train_meta = site_meta[site_meta["project"].isin(train_bps)]
                    test_meta = site_meta[site_meta["project"].isin(test_bps)]
                    if train_meta.empty or test_meta.empty:
                        continue
                    y_train = (train_meta["region"] == target).astype(int).to_numpy()
                    y_test = (test_meta["region"] == target).astype(int).to_numpy()
                    if np.unique(y_train).size < 2 or np.unique(y_test).size < 2:
                        continue

                    train_idx = train_meta["_matrix_index"].to_numpy(dtype=int)
                    test_idx = test_meta["_matrix_index"].to_numpy(dtype=int)
                    model = RandomForestClassifier(
                        n_estimators=N_TREES,
                        random_state=SEED + repeat,
                        class_weight="balanced",
                        n_jobs=N_JOBS,
                    )
                    model.fit(matrix[train_idx], y_train)
                    y_prob = model.predict_proba(matrix[test_idx])[:, 1]
                    y_pred = model.predict(matrix[test_idx])
                    row = fold_metrics(y_test, y_pred, y_prob)
                    row.update(
                        {
                            "dataset": dataset,
                            "site": site,
                            "target_region": target,
                            "repeat": repeat,
                            "fold": fold_i,
                            "fold_id": f"repeat{repeat}_fold{fold_i}",
                            "n_train": int(len(train_meta)),
                            "n_test": int(len(test_meta)),
                            "n_train_projects": int(train_meta["project"].nunique()),
                            "n_test_projects": int(test_meta["project"].nunique()),
                        }
                    )
                    folds.append(row)
                    task_rows.append(row)

                    frame = test_meta[
                        ["sample_id", "region", "project"]
                    ].copy()
                    frame["dataset"] = dataset
                    frame["site"] = site
                    frame["target_region"] = target
                    frame["repeat"] = repeat
                    frame["fold"] = fold_i
                    frame["fold_id"] = f"repeat{repeat}_fold{fold_i}"
                    frame["true_label"] = y_test
                    frame["pred_label"] = y_pred
                    frame["pred_class"] = np.where(y_pred == 1, target, "Rest")
                    frame["pred_prob_target"] = y_prob
                    append_predictions(frame, prediction_path)

            if not task_rows:
                print(f"SKIP {dataset}/{site}/{target}: no valid folds")
                continue
            summaries.append(
                summarize(
                    task_rows,
                    {
                        "dataset": dataset,
                        "site": site,
                        "target_region": target,
                        "n_positive": int((site_meta["region"] == target).sum()),
                        "n_negative": int((site_meta["region"] != target).sum()),
                        "n_positive_projects": int(len(pos_bp)),
                        "n_negative_projects": int(len(neg_bp)),
                        "n_splits": n_folds,
                        "n_repeats": N_REPEATS,
                    },
                )
            )
            print(
                f"DONE OVA {dataset}/{site}/{target}: "
                f"mean AUC={summaries[-1]['mean_auc']:.4f} "
                f"({len(task_rows)} folds, {n_folds}-fold x {N_REPEATS} repeats)"
            )
    return folds, summaries


def run_ovo(meta, matrix, dataset, prediction_path):
    folds, summaries = [], []
    for site in sorted(
        meta["site"].unique(),
        key=lambda x: ORAL_SITES.index(x) if x in ORAL_SITES else 99,
    ):
        site_meta = meta[meta["site"] == site].copy()
        observed = [
            region for region in REGIONS if region in set(site_meta["region"])
        ]
        print(f"\n=== OVO {dataset}/{site} ===")
        for i, region_1 in enumerate(observed):
            for region_2 in observed[i + 1 :]:
                pair_meta = site_meta[
                    site_meta["region"].isin([region_1, region_2])
                ].copy()
                bp_1 = (
                    pair_meta.loc[pair_meta["region"] == region_1, "project"]
                    .dropna()
                    .unique()
                )
                bp_2 = (
                    pair_meta.loc[pair_meta["region"] == region_2, "project"]
                    .dropna()
                    .unique()
                )
                if len(bp_1) < 2 or len(bp_2) < 2:
                    print(
                        f"SKIP {dataset}/{site}/{region_1} vs {region_2}: "
                        f"projects {[len(bp_1), len(bp_2)]}"
                    )
                    continue
                n_folds = min(N_FOLDS_MAX, len(bp_1), len(bp_2))
                if n_folds < 2:
                    continue

                task_rows = []
                for repeat in range(1, N_REPEATS + 1):
                    rng = np.random.RandomState(SEED + repeat)
                    for fold_i, (train_bps, test_bps) in enumerate(
                        stratified_project_fold_split(bp_1, bp_2, n_folds, rng),
                        start=1,
                    ):
                        train_meta = pair_meta[
                            pair_meta["project"].isin(train_bps)
                        ]
                        test_meta = pair_meta[
                            pair_meta["project"].isin(test_bps)
                        ]
                        if train_meta.empty or test_meta.empty:
                            continue
                        y_train = (
                            train_meta["region"] == region_1
                        ).astype(int).to_numpy()
                        y_test = (
                            test_meta["region"] == region_1
                        ).astype(int).to_numpy()
                        if (
                            np.unique(y_train).size < 2
                            or np.unique(y_test).size < 2
                        ):
                            continue

                        train_idx = train_meta["_matrix_index"].to_numpy(dtype=int)
                        test_idx = test_meta["_matrix_index"].to_numpy(dtype=int)
                        model = RandomForestClassifier(
                            n_estimators=N_TREES,
                            random_state=SEED + repeat,
                            class_weight="balanced",
                            n_jobs=N_JOBS,
                        )
                        model.fit(matrix[train_idx], y_train)
                        y_prob = model.predict_proba(matrix[test_idx])[:, 1]
                        y_pred = model.predict(matrix[test_idx])
                        row = fold_metrics(y_test, y_pred, y_prob)
                        row.update(
                            {
                                "dataset": dataset,
                                "site": site,
                                "region_1": region_1,
                                "region_2": region_2,
                                "repeat": repeat,
                                "fold": fold_i,
                                "fold_id": f"repeat{repeat}_fold{fold_i}",
                                "n_train": int(len(train_meta)),
                                "n_test": int(len(test_meta)),
                                "n_train_projects": int(
                                    train_meta["project"].nunique()
                                ),
                                "n_test_projects": int(
                                    test_meta["project"].nunique()
                                ),
                            }
                        )
                        folds.append(row)
                        task_rows.append(row)

                        frame = test_meta[
                            ["sample_id", "region", "project"]
                        ].copy()
                        frame["dataset"] = dataset
                        frame["site"] = site
                        frame["region_1"] = region_1
                        frame["region_2"] = region_2
                        frame["repeat"] = repeat
                        frame["fold"] = fold_i
                        frame["fold_id"] = f"repeat{repeat}_fold{fold_i}"
                        frame["true_label"] = y_test
                        frame["pred_label"] = y_pred
                        frame["pred_region"] = np.where(
                            y_pred == 1, region_1, region_2
                        )
                        frame["pred_prob_region_1"] = y_prob
                        append_predictions(frame, prediction_path)

                if not task_rows:
                    print(
                        f"SKIP {dataset}/{site}/{region_1} vs {region_2}: "
                        "no valid folds"
                    )
                    continue
                summaries.append(
                    summarize(
                        task_rows,
                        {
                            "dataset": dataset,
                            "site": site,
                            "region_1": region_1,
                            "region_2": region_2,
                            "n_region_1": int(
                                (pair_meta["region"] == region_1).sum()
                            ),
                            "n_region_2": int(
                                (pair_meta["region"] == region_2).sum()
                            ),
                            "n_region_1_projects": int(len(bp_1)),
                            "n_region_2_projects": int(len(bp_2)),
                            "n_splits": n_folds,
                            "n_repeats": N_REPEATS,
                        },
                    )
                )
                print(
                    f"DONE OVO {dataset}/{site}/{region_1} vs {region_2}: "
                    f"mean AUC={summaries[-1]['mean_auc']:.4f} "
                    f"({len(task_rows)} folds, {n_folds}-fold x {N_REPEATS} repeats)"
                )
    return folds, summaries


def main():
    oral_meta, oral_matrix = load_oral()
    print(f"Oral features: {oral_matrix.shape[1]:,}")
    gut_meta, gut_matrix = load_gut()
    print(f"Gut features: {gut_matrix.shape[1]:,}")
    jobs = [
        ("Oral", oral_meta, oral_matrix),
        ("Gut", gut_meta, gut_matrix),
    ]

    print("CV: project-level stratified folds (v2-style); datasets=['oral', 'gut']")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    ova_prediction_path = OUTDIR / "region_ova_predictions.tsv"
    ovo_prediction_path = OUTDIR / "region_ovo_predictions.tsv"
    ova_prediction_path.unlink(missing_ok=True)
    ovo_prediction_path.unlink(missing_ok=True)

    all_ova_folds, all_ova_summaries = [], []
    all_ovo_folds, all_ovo_summaries = [], []
    for dataset, meta, matrix in jobs:
        f, s = run_ova(meta, matrix, dataset, ova_prediction_path)
        all_ova_folds.extend(f)
        all_ova_summaries.extend(s)
        f, s = run_ovo(meta, matrix, dataset, ovo_prediction_path)
        all_ovo_folds.extend(f)
        all_ovo_summaries.extend(s)

    pd.DataFrame(all_ova_folds).to_csv(
        OUTDIR / "region_ova_fold_metrics.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_ova_summaries).to_csv(
        OUTDIR / "region_ova_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_ovo_folds).to_csv(
        OUTDIR / "region_ovo_fold_metrics.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_ovo_summaries).to_csv(
        OUTDIR / "region_ovo_summary.tsv", sep="\t", index=False
    )
    config = {
        "oral_metadata": str(ORAL_METADATA),
        "oral_biom": str(ORAL_BIOM),
        "gut_metadata": str(GUT_METADATA),
        "gut_biom": str(GUT_BIOM),
        "output_dir": str(OUTDIR),
        "datasets": ["oral", "gut"],
        "regions": REGIONS,
        "oral_sites": ORAL_SITES,
        "folds_max": N_FOLDS_MAX,
        "repeats": N_REPEATS,
        "seed": SEED,
        "trees": N_TREES,
        "jobs": N_JOBS,
        "min_study_samples": MIN_STUDY_SAMPLES,
        "cv_strategy": "project-level stratified folds (randomforest_v2 style)",
        "group_column_oral": "bioproject_id",
        "group_column_gut": "project",
        "grouping_used": True,
        "reference": (
            "/beegfs/huangxiaochang/global_oral/new_analysis/randomforest_v2"
        ),
    }
    with (OUTDIR / "run_config.json").open("w") as handle:
        json.dump(config, handle, indent=2, sort_keys=True)
    print(f"\nOutputs written to {OUTDIR}")


if __name__ == "__main__":
    main()
