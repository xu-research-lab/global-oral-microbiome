"""Count-table loading, HC feature preparation, and shared age RF fitting."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold

from age_common import GENUS_TABLE_PATH

EPS = 1e-6


@dataclass(frozen=True)
class ModelConfig:
    trees: int = 500
    folds: int = 5
    seed: int = 42
    jobs: int = 16


def load_genus_table(needed_sids=None, path=GENUS_TABLE_PATH):
    """Read selected sample columns together with the genus index column."""
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if len(set(header[1:])) != len(header) - 1:
        raise ValueError("Genus table contains duplicate sample columns")
    columns = [0] + [i for i, sid in enumerate(header[1:], 1)
                     if needed_sids is None or sid in needed_sids]
    if len(columns) == 1:
        raise ValueError("No sample IDs overlap between metadata and genus table")
    table = pd.read_csv(path, sep="\t", usecols=columns, index_col=0)
    matrix = table.to_numpy(dtype=np.float64).T
    if not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError("Genus counts must be finite and nonnegative")
    if (matrix.sum(axis=1) == 0).any():
        raise ValueError("Genus table contains samples with zero total counts")
    return table.index.to_numpy(), table.columns.to_numpy(), matrix


def build_xy_pool(table_X, table_sids, md_train, md_test, sid_to_age):
    """Use HC training samples and a vocabulary selected only from training."""
    sid_idx = {sid: i for i, sid in enumerate(table_sids)}
    train_sid = md_train.loc[md_train.status_cleaned == "HC", "sample_id"].astype(str)
    test_sid = md_test["sample_id"].astype(str)
    missing = (set(train_sid) | set(test_sid)) - set(sid_idx)
    if missing:
        raise ValueError(f"Fold sample IDs absent from genus matrix: {sorted(missing)[:5]}")
    Xtr = table_X[[sid_idx[sid] for sid in train_sid]].copy()
    Xte = table_X[[sid_idx[sid] for sid in test_sid]].copy()
    Xtr /= Xtr.sum(axis=1, keepdims=True)
    Xte /= Xte.sum(axis=1, keepdims=True)
    keep = Xtr.sum(axis=0) > 0
    ytr = sid_to_age.reindex(train_sid).to_numpy(dtype=float)
    yte = sid_to_age.reindex(test_sid).to_numpy(dtype=float)
    return np.log10(Xtr[:, keep] + EPS), ytr, np.log10(Xte[:, keep] + EPS), yte


def split_pool(pool, scheme, config=ModelConfig()):
    if scheme not in {"project", "random"}:
        raise ValueError(f"Unknown split scheme: {scheme}")
    units = pool.bioproject_id.unique() if scheme == "project" else pool.sample_id.to_numpy()
    folds = min(config.folds, len(units)) if scheme == "project" else config.folds
    if folds < 2 or len(units) < folds:
        return
    column = "bioproject_id" if scheme == "project" else "sample_id"
    for train, test in KFold(folds, shuffle=True, random_state=config.seed).split(units):
        yield pool[pool[column].isin(units[train])], pool[pool[column].isin(units[test])]


def fit_fold(table_X, table_sids, md_train, md_test, sid_to_age, config=ModelConfig()):
    if (md_train.status_cleaned == "HC").sum() < 5 or len(md_test) < 2:
        return None
    Xtr, ytr, Xte, yte = build_xy_pool(table_X, table_sids, md_train, md_test, sid_to_age)
    model = RandomForestRegressor(
        n_estimators=config.trees, random_state=config.seed, n_jobs=config.jobs,
        max_features="sqrt", min_samples_leaf=2,
    )
    model.fit(Xtr, ytr)
    return pd.DataFrame({"sample_id": md_test.sample_id.to_numpy(),
                         "true_age": yte, "pred_age": model.predict(Xte)})
