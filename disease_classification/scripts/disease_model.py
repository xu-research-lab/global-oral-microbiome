"""Shared RF fitting and AUC calculation for the two LSO genus filters."""
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

from disease_common import (RESULT_DIR, PLOT_SITES, RF_OPTIONS, fold_name, load_labels,
                            predictions_path, read_manifest, split_path)


def fit_fold(train, test, labels, jobs=16):
    train = train.norm(axis="sample", inplace=False).to_dataframe().T
    test = test.norm(axis="sample", inplace=False).to_dataframe().T
    test = test.reindex(columns=train.columns, fill_value=0)
    y_train = labels.loc[train.index]
    if set(y_train.unique()) != {0, 1}:
        raise ValueError("Training split must contain controls and cases")
    model = RandomForestClassifier(**RF_OPTIONS, n_jobs=jobs)
    model.fit(train.values, y_train)
    probabilities = model.predict_proba(test.values)[:, 1]
    return test.index, probabilities, pd.DataFrame({"feature": train.columns, "importance": model.feature_importances_})


def fit_models(scheme, jobs=16):
    import biom

    labels = load_labels()
    predictions, importances = [], []
    for row in read_manifest(scheme).itertuples(index=False):
        name = fold_name(row)
        train = biom.load_table(str(split_path(scheme, "train", name)))
        test = biom.load_table(str(split_path(scheme, "test", name)))
        samples, probabilities, importance = fit_fold(train, test, labels, jobs)
        predictions.append(pd.DataFrame({"sampel_id": samples, "group": name, "pred_prob": probabilities}))
        if scheme == "prevalence":
            importances.append(importance.assign(fold=name))
        print(f"{name}: AUC={safe_auc(labels.loc[samples], probabilities):.4f}", flush=True)
    if not predictions:
        raise ValueError("Empty fold manifest")
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    pd.concat(predictions, ignore_index=True).to_csv(predictions_path(scheme), index=False)
    if importances:
        pd.concat(importances, ignore_index=True).to_csv(RESULT_DIR / "RF_feature_importances_lsoP.csv", index=False)


def safe_auc(truth, scores):
    values = pd.DataFrame({"truth": np.asarray(truth), "scores": np.asarray(scores)}).dropna()
    if values.truth.nunique() < 2:
        return np.nan
    return roc_auc_score(values.truth, values.scores)


def auc_tables():
    """Summarize prevalence-filtered LOSO predictions for the final figures."""
    labels = load_labels()
    predictions = pd.read_csv(predictions_path("prevalence"))
    grouped = dict(tuple(predictions.groupby("group", sort=False)))
    records = []
    for row in read_manifest("prevalence").itertuples(index=False):
        if row.site not in PLOT_SITES:
            continue
        name = fold_name(row)
        if name not in grouped:
            raise ValueError(f"Missing predictions: {name}")
        pred = grouped[name]
        sample_ids = pred.sampel_id.astype(str)
        truth = labels.reindex(sample_ids)
        if truth.isna().any():
            raise ValueError(f"Missing labels in {name}")
        records.append(dict(site=row.site, disease=row.disease, test_study=row.study,
                            auc=safe_auc(truth, pred.pred_prob)))
    result = pd.DataFrame(records)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    for site in PLOT_SITES:
        result.loc[result.site == site, ["disease", "test_study", "auc"]].to_csv(
            RESULT_DIR / f"rf_loso_auc_all_diseases_genus_{site}.csv", index=False)
    return result
