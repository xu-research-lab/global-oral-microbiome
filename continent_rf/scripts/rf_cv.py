"""Shared binary RF evaluation for OVO and OVA tasks."""

from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, f1_score,
                             precision_score, recall_score, roc_auc_score)

from rf_data import ORAL_SITES, REGIONS


@dataclass(frozen=True)
class CVConfig:
    folds_max: int = 5
    repeats: int = 5
    seed: int = 42
    trees: int = 200
    jobs: int = 8


def stratified_project_fold_split(bp_pos, bp_neg, n_folds, rng):
    """Shuffle classes separately; assign any shared project to one fold only."""
    pos, neg = list(bp_pos), list(bp_neg)
    rng.shuffle(pos)
    rng.shuffle(neg)
    fold_map = {}
    for units in (pos, neg):
        for i, unit in enumerate(units):
            fold_map.setdefault(unit, i % n_folds)
    folds = [[] for _ in range(n_folds)]
    for unit, fold in fold_map.items():
        folds[fold].append(unit)
    for i, test in enumerate(folds):
        yield [unit for j, fold in enumerate(folds) if j != i for unit in fold], test


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
    return {
        **base, "n_evaluations": len(frame),
        "mean_auc": frame["auc"].mean(), "sd_auc": frame["auc"].std(ddof=1),
        "median_auc": frame["auc"].median(), "min_auc": frame["auc"].min(),
        "max_auc": frame["auc"].max(),
        **{f"mean_{name}": frame[name].mean() for name in
           ["accuracy", "balanced_accuracy", "precision", "recall", "f1"]},
    }


def output_columns(task):
    identity = ["target_region"] if task == "ova" else ["region_1", "region_2"]
    common = ["dataset", "site", *identity]
    repeat = ["repeat", "fold", "fold_id"]
    metrics = ["auc", "accuracy", "balanced_accuracy", "precision", "recall", "f1"]
    counts = (["n_positive", "n_negative", "n_positive_units", "n_negative_units"]
              if task == "ova" else
              ["n_region_1", "n_region_2", "n_region_1_units", "n_region_2_units"])
    return {
        "fold_metrics": [*metrics, *common, *repeat, "n_train", "n_test",
                         "n_train_units", "n_test_units"],
        "summary": [*common, *counts, "n_splits", "n_repeats", "n_evaluations",
                    "mean_auc", "sd_auc", "median_auc", "min_auc", "max_auc",
                    *[f"mean_{m}" for m in metrics[1:]]],
        "predictions": ["sample_id", "region", "project", *common, *repeat,
                        "true_label", "pred_label",
                        *( ["pred_class", "pred_prob_target"] if task == "ova" else
                           ["pred_region", "pred_prob_region_1"])],
    }


def run_cv(meta, matrix, dataset, task, prediction_path, config=CVConfig()):
    """Evaluate one task type; preserve historical labels, seeds and TSV schema."""
    if task not in {"ova", "ovo"}:
        raise ValueError(f"Unknown task: {task}")
    folds, summaries = [], []
    for site in sorted(meta["site"].unique(),
                       key=lambda s: ORAL_SITES.index(s) if s in ORAL_SITES else 99):
        site_meta = meta[meta["site"] == site]
        observed = [region for region in REGIONS if region in set(site_meta["region"])]
        pairs = ((r, "Rest") for r in observed) if task == "ova" else combinations(observed, 2)
        for positive, negative in pairs:
            frame = (site_meta if task == "ova" else
                     site_meta[site_meta["region"].isin([positive, negative])])
            is_positive = frame["region"] == positive
            pos_units = frame.loc[is_positive, "unit"].dropna().unique()
            neg_units = frame.loc[~is_positive, "unit"].dropna().unique()
            n_folds = min(config.folds_max, len(pos_units), len(neg_units))
            label = f"{task.upper()} {dataset}/{site}/{positive} vs {negative}"
            if n_folds < 2:
                print(f"SKIP {label}: units pos={len(pos_units)} neg={len(neg_units)}")
                continue
            identity = ({"target_region": positive} if task == "ova" else
                        {"region_1": positive, "region_2": negative})
            task_rows = []
            for repeat in range(1, config.repeats + 1):
                rng = np.random.RandomState(config.seed + repeat)
                for fold, (train_units, test_units) in enumerate(
                    stratified_project_fold_split(pos_units, neg_units, n_folds, rng), 1
                ):
                    train = frame[frame["unit"].isin(train_units)]
                    test = frame[frame["unit"].isin(test_units)]
                    y_train = (train["region"] == positive).astype(int).to_numpy()
                    y_test = (test["region"] == positive).astype(int).to_numpy()
                    if np.unique(y_train).size < 2 or np.unique(y_test).size < 2:
                        continue
                    model = RandomForestClassifier(
                        n_estimators=config.trees, random_state=config.seed + repeat,
                        class_weight="balanced", n_jobs=config.jobs,
                    )
                    model.fit(matrix[train["_matrix_index"].to_numpy(dtype=int)], y_train)
                    test_matrix = matrix[test["_matrix_index"].to_numpy(dtype=int)]
                    y_prob = model.predict_proba(test_matrix)[:, 1]
                    y_pred = model.predict(test_matrix)
                    common = {"dataset": dataset, "site": site, **identity,
                              "repeat": repeat, "fold": fold,
                              "fold_id": f"repeat{repeat}_fold{fold}"}
                    row = {**fold_metrics(y_test, y_pred, y_prob), **common,
                           "n_train": len(train), "n_test": len(test),
                           "n_train_units": train["unit"].nunique(),
                           "n_test_units": test["unit"].nunique()}
                    folds.append(row)
                    task_rows.append(row)
                    predictions = test[["sample_id", "region", "project"]].copy()
                    for key, value in common.items():
                        predictions[key] = value
                    predictions["true_label"] = y_test
                    predictions["pred_label"] = y_pred
                    predictions["pred_class" if task == "ova" else "pred_region"] = np.where(
                        y_pred == 1, positive, negative)
                    predictions["pred_prob_target" if task == "ova" else "pred_prob_region_1"] = y_prob
                    predictions.to_csv(prediction_path, sep="\t", index=False,
                                       mode="a", header=not prediction_path.exists())
            if not task_rows:
                print(f"SKIP {label}: no valid folds")
                continue
            count_names = (["n_positive", "n_negative", "n_positive_units", "n_negative_units"]
                           if task == "ova" else
                           ["n_region_1", "n_region_2", "n_region_1_units", "n_region_2_units"])
            base = {"dataset": dataset, "site": site, **identity,
                    **dict(zip(count_names, [int(is_positive.sum()), int((~is_positive).sum()),
                                            len(pos_units), len(neg_units)])),
                    "n_splits": n_folds, "n_repeats": config.repeats}
            summaries.append(summarize(task_rows, base))
            print(f"DONE {label}: mean AUC={summaries[-1]['mean_auc']:.4f} ({len(task_rows)} folds)")
    return folds, summaries
