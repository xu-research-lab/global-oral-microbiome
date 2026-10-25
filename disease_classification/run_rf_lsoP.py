"""
RF training + feature-importance extraction for the prevalence-filtered LSO
splits (data/lsoPgenus_* biom tables from gen_lso_splits_prevalence.py).

Mirrors run_rf_genus_fi.py (same make_rf() and name_fn scheme) so the
prevalence-standard results are directly comparable to the mean-based ones:

    result/RF_feature_importances_lsoP.csv   (columns: feature, importance, fold)
    result/RF_results_lsoP.csv               (columns: sampel_id, group, pred_prob)

Run: python run_rf_lsoP.py
"""
import biom
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

N_JOBS = 16


def make_rf():
    # fresh model per split so nothing leaks between folds
    return RandomForestClassifier(
        n_estimators=500, max_depth=5, max_features="log2",
        random_state=42, n_jobs=N_JOBS)


def run_lso():
    metadata = pd.read_csv("meta_genus.tsv", sep="\t", index_col=0)
    run_info = pd.read_csv("run_info_lsoP_genus.csv")
    fi_frames, preds, folds, sids = [], [], [], []
    for i in run_info.index.values:
        row = run_info.loc[i]
        s, d, st = row['site'], row['disease'], row['study']
        fold_name = f"{s}_{d}_{st}"
        train_table = f"data/lsoPgenus_train_{s}_{d}_{st}.biom"
        test_table = f"data/lsoPgenus_test_{s}_{d}_{st}.biom"

        X_train = biom.load_table(train_table).norm(axis='sample', inplace=False).to_dataframe().T
        X_test = biom.load_table(test_table).norm(axis='sample', inplace=False).to_dataframe().T
        test_sid = X_test.index.values

        y_train = metadata.loc[X_train.index.values, "group"]
        y_test = metadata.loc[X_test.index.values, "group"]
        X_test = X_test.reindex(columns=X_train.columns, fill_value=0)

        rf = make_rf()
        rf.fit(X_train.values, y_train)
        y_prob = rf.predict_proba(X_test.values)[:, 1]

        fi_frames.append(pd.DataFrame({
            "feature": X_train.columns,
            "importance": rf.feature_importances_,
            "fold": fold_name,
        }))
        preds += list(y_prob)
        sids += list(test_sid)
        folds += [fold_name] * len(y_prob)

        if len(np.unique(y_test)) == 2:
            print(f"{fold_name}\tn_genus={X_train.shape[1]}\t"
                  f"AUC={roc_auc_score(y_test, y_prob):.4f}", flush=True)
        else:
            print(f"{fold_name}\tn_genus={X_train.shape[1]}\tAUC=NA "
                  f"(single-class test set)", flush=True)

    pd.concat(fi_frames, ignore_index=True).to_csv(
        "result/RF_feature_importances_lsoP.csv", index=False)
    pd.DataFrame({"sampel_id": sids, "group": folds, "pred_prob": preds}).to_csv(
        "result/RF_results_lsoP.csv", index=False)
    print("  -> wrote result/RF_feature_importances_lsoP.csv and "
          "result/RF_results_lsoP.csv", flush=True)


if __name__ == "__main__":
    print("### prevalence-standard LSO RF", flush=True)
    run_lso()