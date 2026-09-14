import biom
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

metadata = pd.read_csv("meta.tsv", sep="\t", index_col=0)


def make_rf():
    # fresh model per split so nothing leaks between folds
    return RandomForestClassifier(
        n_estimators=500, max_depth=5, max_features="log2",
        random_state=42, n_jobs=10)


def run_split(run_info, name_fn, out_csv):
    """
    name_fn(row) -> (train_biom_path, test_biom_path, fold_name)
    Trains one RF per row, collects per-sample predicted probabilities,
    prints per-fold AUC, writes out_csv.
    """
    fold, pred_prob, sample_id = [], [], []
    for i in run_info.index.values:
        row = run_info.loc[i]
        train_table, test_table, fold_name = name_fn(row)

        X_train = biom.load_table(train_table).norm(axis='sample', inplace=False).to_dataframe().T
        X_test  = biom.load_table(test_table).norm(axis='sample', inplace=False).to_dataframe().T
        test_sid = X_test.index.values

        y_train = metadata.loc[X_train.index.values, "group"]
        y_test  = metadata.loc[X_test.index.values,  "group"]

        # align test features to train features (biom filter already does this,
        # but guard against any column mismatch after norm/to_dataframe)
        X_test = X_test.reindex(columns=X_train.columns, fill_value=0)

        rf = make_rf()
        rf.fit(X_train.values, y_train)
        y_prob = rf.predict_proba(X_test.values)[:, 1]

        pred_prob += list(y_prob)
        sample_id += list(test_sid)
        fold      += [fold_name] * len(y_prob)

        # AUC only defined if test has both classes
        if len(np.unique(y_test)) == 2:
            print(f"{fold_name}\tAUC={roc_auc_score(y_test, y_prob):.4f}")
        else:
            print(f"{fold_name}\tAUC=NA (single-class test set)")

    pd.DataFrame({"sampel_id": sample_id, "group": fold, "pred_prob": pred_prob}) \
        .to_csv(out_csv, index=False)
    print(f"  -> wrote {out_csv}\n")


# =========================================================================
# part 1: leave-one-study-out
# =========================================================================
# def name_lso(r):
#     s, d, st = r['site'], r['disease'], r['study']
#     return (f"data/lso_train_{s}_{d}_{st}.biom",
#             f"data/lso_test_{s}_{d}_{st}.biom",
#             f"{s}_{d}_{st}")

# print("### part 1: leave-one-study-out")
# metadata = pd.read_csv("meta.tsv", sep="\t", index_col=0)
# run_split(pd.read_csv("run_info.csv"), name_lso, "result/RF_results_lso.csv")


# =========================================================================
# part 2: five-fold (oral/gut)
# csv columns: site,fold,study  ;  filename order is _{site}_{study}_{fold}
# =========================================================================
def name_5fold(r):
    s, f, st = r['site'], r['fold'], r['study']
    return (f"data/five_fold_train_{s}_{st}_{f}.biom",
            f"data/five_fold_test_{s}_{st}_{f}.biom",
            f"{s}_{st}_{f}")

# print("### part 2: five-fold oral/gut")
# metadata = pd.read_csv("meta_oral_gut.tsv", sep="\t", index_col=0)
# run_split(pd.read_csv("run_info_oral_gut.csv"), name_5fold, "result/RF_results_five_fold.csv")


# =========================================================================
# part 3: continent (within_continent + leave_one_continent_out)
# csv columns: site,disease,continent,test_study,split_type
# within-continent files carry the continent + test_study; LOCO files carry
# only the continent (held-out continent).
# =========================================================================
def name_continent(r):
    s, d = r['site'], r['disease']
    st_type = str(r['split_type'])
    cont, test_study = r['continent'], r['test_study']

    if st_type in ("wc", "within_continent"):
        base = f"wc_{s}_{d}_{cont}_{test_study}"
    else:  # leave_one_continent_out / loco
        base = f"loco_{s}_{d}_{cont}"
    return (f"data/train_{base}.biom",
            f"data/test_{base}.biom",
            base)

# print("### part 3: within-continent + leave-one-continent-out")
# metadata = pd.read_csv("meta.tsv", sep="\t", index_col=0)
# run_split(pd.read_csv("run_info_continent.csv"), name_continent, "result/RF_results_continent.csv")


# =========================================================================
# part 4: cross-disease (train on A, test on B)
# csv columns: site,train_disease,test_disease,n_train,n_test,n_train_case,n_test_case,split_type
# =========================================================================
def name_cross_disease(r):
    s, a, b = r['site'], r['train_disease'], r['test_disease']
    base = f"cd_{s}_train_{a}_test_{b}"
    return (f"data/train_{base}.biom",
            f"data/test_{base}.biom",
            base)

print("### part 4: cross-disease")
metadata = pd.read_csv("meta.tsv", sep="\t", index_col=0)
run_split(pd.read_csv("run_info_cross_disease.csv"), name_cross_disease, "result/RF_results_cross_disease.csv")


# =========================================================================
# part 5: learning curve
# csv columns: site,disease,n_studies,repeat,subset,test_study,split_type
# =========================================================================
def name_lc(r):
    s, d = r['site'], r['disease']
    k, rep, test_study = r['n_studies'], r['repeat'], r['test_study']
    base = f"lc_{s}_{d}_k{k}_r{rep}_test_{test_study}"
    return (f"data/train_{base}.biom",
            f"data/test_{base}.biom",
            base)

print("### part 5: learning curve")
metadata = pd.read_csv("meta.tsv", sep="\t", index_col=0)
run_split(pd.read_csv("run_info_lc.csv"), name_lc, "result/RF_results_lc.csv")