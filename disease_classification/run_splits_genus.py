# Auto-extracted from disease_predict.ipynb (genus version) -- DO NOT EDIT BY HAND
# Part A: split generation (notebook cells except plotting)
import matplotlib
matplotlib.use('Agg')
# ===== notebook cell 2 =====
import sys

import matplotlib as mpl
import os
import math
import biom
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.model_selection import StratifiedKFold, train_test_split

# ===== notebook cell 3 =====
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

# ===== notebook cell 4 =====
BASE = os.environ.get("GLOBAL_ORAL_ROOT", "..")  # repo root (data/ lives there, server-side)
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

# ===== notebook cell 5 =====
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score

# ===== notebook cell 6 =====
from itertools import product

# ===== notebook cell 7 =====
import random
from math import comb
from itertools import combinations

# ===== notebook cell 8 =====
# ---------------------------------------------------------------------------
# Genus-level aggregation helper
# ---------------------------------------------------------------------------
# Taxonomy file for all OTU / feature IDs in the global oral collection.
TAXONOMY_PATH = "/beegfs/huangxiaochang/global_oral/taxonomy.tsv"


def _parse_genus(taxon_str):
    """Extract the genus label (g__...) from a semicolon-separated taxonomy string."""
    if pd.isna(taxon_str):
        return None
    parts = [p.strip() for p in str(taxon_str).split(";")]
    for p in parts:
        if p.startswith("g__"):
            return p
    return None


# Load once and build a global Feature ID -> genus map.
_taxonomy = pd.read_csv(TAXONOMY_PATH, sep="\t")
_taxonomy["genus"] = _taxonomy["Taxon"].apply(_parse_genus)
_taxonomy = _taxonomy.dropna(subset=["genus"])
FID_TO_GENUS = dict(zip(_taxonomy["Feature ID"].astype(str), _taxonomy["genus"]))
print(f"[genus] loaded taxonomy for {len(FID_TO_GENUS)} features")


def aggregate_to_genus(table):
    """Collapse a BIOM feature table to genus level.

    Observations without a genus annotation are dropped; remaining observations
    are summed by their genus label.
    """
    obs_ids = table.ids(axis="observation")
    keep = [fid for fid in obs_ids if fid in FID_TO_GENUS]
    if not keep:
        raise ValueError("No observation IDs in the table could be mapped to genus.")
    table = table.filter(keep, axis="observation", inplace=False)

    def _genus_mapper(id_, md):
        return FID_TO_GENUS[id_]

    return table.collapse(_genus_mapper, axis="observation")

# ===== notebook cell 10 =====
table = biom.load_table(os.path.join(BASE, "data/metadata-and-biom-v2-qc/biom_oral_qc_dedup.biom"))
table = aggregate_to_genus(table)
metadata = pd.read_csv(os.path.join(BASE, "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv"), sep="\t", low_memory=False)
sid = metadata.sample_id.values
table.filter(sid, axis="sample")
table.remove_empty()
# fid = table.ids(axis="observation")
# prevalence = table.nonzero_counts(axis="observation")
# table = table.filter(fid[prevalence > 50], axis='observation', inplace=False)
# table.remove_empty()
sid_reads = table.sum(axis="sample")
sid = table.ids(axis="sample")
table = table.filter(sid[sid_reads > 1000], axis='sample', inplace=False)
table.remove_empty()

sid = table.ids(axis="sample")
metadata = metadata.loc[[i in sid for i in metadata.sample_id.values]]

# ===== notebook cell 11 =====
metadata.columns

# ===== notebook cell 12 =====
metadata.status_cleaned.unique()

# ===== notebook cell 13 =====
table.shape

# ===== notebook cell 14 =====
# ---------------- Parameters ----------------
META_PATH = os.path.join(BASE, "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv")
MIN_CASE_TOTAL   = 30         # minimum pooled case samples across matched BPs
MIN_HC_TOTAL     = 30         # minimum pooled HC samples across matched BPs
MIN_MATCHED_BP   = 2          # minimum number of matched BPs per site
ABUND_THRESH     = 0.01       # keep a genus if mean rel. abundance >= 1% in any cohort
COHORT_COL       = "bioproject_id"
RENORM_AFTER_FILTER = True    # renormalize each sample to 1 after genus filtering

oral_sites = ['saliva', 'tooth', 'tongue', 'oral mucosa']
run_site, run_disease, run_test_study = [], [], []


def genus_keep_by_cohort_abundance(table_rel, meta_sub, thresh=ABUND_THRESH,
                                   cohort_col=COHORT_COL):
    """Keep genera whose mean relative abundance is >= thresh in at least one cohort."""
    fid = np.asarray(table_rel.ids(axis="observation"))
    sid = np.asarray(table_rel.ids(axis="sample"))
    mat = pd.DataFrame(table_rel.matrix_data.toarray(), index=fid, columns=sid)  # genus x sample
    cohort = meta_sub.set_index("sample_id").loc[sid, cohort_col].astype(str).values
    mean_by_cohort = mat.T.groupby(cohort).mean().T
    keep = mean_by_cohort.max(axis=1) >= thresh
    return fid[keep.values], mean_by_cohort


def get_matched_bp(metadata_site, disease):
    """Return matched BPs (both HC and disease present, pooled sizes large enough),
    or [] if the site does not qualify for this disease."""
    counts = (metadata_site.groupby(["bioproject_id", "status_cleaned"])["sample_id"]
                            .count().reset_index(name="n"))
    sub = counts[counts.status_cleaned.isin(["HC", disease])]
    matched_bp = (sub.groupby("bioproject_id")["status_cleaned"].nunique()
                     .pipe(lambda x: x[x == 2]).index.tolist())
    if len(matched_bp) < MIN_MATCHED_BP:
        return []

    sub = sub[sub.bioproject_id.isin(matched_bp)]
    n_case = int(sub.loc[sub.status_cleaned == disease, "n"].sum())
    n_hc   = int(sub.loc[sub.status_cleaned == "HC", "n"].sum())
    if n_case < MIN_CASE_TOTAL or n_hc < MIN_HC_TOTAL:
        return []
    return matched_bp


# ---------------- 0) Load & pre-clean metadata once ----------------
sid_all = table.ids(axis="sample")
metadata_all = pd.read_csv(META_PATH, sep="\t", low_memory=False)
metadata_all = metadata_all.loc[metadata_all.status_cleaned.notna()]
metadata_all = metadata_all.loc[metadata_all.status_cleaned != "Other"]
metadata_all = metadata_all.loc[metadata_all.sample_id.isin(set(sid_all))]
metadata_all = metadata_all.loc[metadata_all.site_cleaned.isin(oral_sites)]

all_diseases = sorted(metadata_all.loc[metadata_all.status_cleaned != "HC",
                                       "status_cleaned"].unique())

# ---------------- 1) Filter by disease first: how many sites does it qualify in ----------------
for disease in all_diseases:
    print(f"\n########## DISEASE: {disease} ##########")

    site_matched_bp = {}   # site -> matched_bp list
    for site in oral_sites:
        metadata_site = metadata_all.loc[metadata_all.site_cleaned == site]
        if metadata_site.empty:
            continue
        matched_bp = get_matched_bp(metadata_site, disease)
        if matched_bp:
            site_matched_bp[site] = matched_bp

    print(f"  Qualifies in {len(site_matched_bp)} site(s): {list(site_matched_bp.keys())}")

    # ---------------- 2) Leave-one-study-out per qualifying site ----------------
    for site, matched_bp in site_matched_bp.items():
        print(f"\n=== Processing: {disease} | {site} ===")

        metadata_site = metadata_all.loc[metadata_all.site_cleaned == site]
        meta_pick = metadata_site.loc[
            metadata_site.bioproject_id.isin(matched_bp)
            & metadata_site.status_cleaned.isin(["HC", disease])
        ].copy()
        sid = meta_pick.sample_id.values
        table_pick = table.filter(sid, axis="sample", inplace=False)

        # Convert to relative abundance, then filter genera
        table_pick.norm(axis="sample", inplace=True)
        keep_fid, mean_by_cohort = genus_keep_by_cohort_abundance(table_pick, meta_pick)
        print(f"  genus: {table_pick.shape[0]} -> {len(keep_fid)} "
              f"(mean relative abundance >= {ABUND_THRESH:.0%} in at least one cohort)")
        if len(keep_fid) == 0:
            print("  No genus passed the abundance filter, skipping...")
            continue

        table_pick = table_pick.filter(keep_fid, axis="observation", inplace=False)
        if RENORM_AFTER_FILTER:
            table_pick.norm(axis="sample", inplace=True)
        if table_pick.shape[1] == 0:
            continue

        # ---- Leave-one-study-out ----
        print(f"  {len(matched_bp)} bioprojects with {len(meta_pick)} samples")
        for j in matched_bp:
            run_disease.append(disease.replace(' ', '_'))
            run_site.append(site.replace(' ', '_'))
            run_test_study.append(j)

            train_sid = meta_pick.loc[meta_pick.bioproject_id != j, "sample_id"].values
            test_sid  = meta_pick.loc[meta_pick.bioproject_id == j, "sample_id"].values
            train_biom = table_pick.filter(train_sid, axis="sample", inplace=False)
            test_biom  = table_pick.filter(test_sid,  axis="sample", inplace=False)
            print(f"    [test={j}] train: {train_biom.shape}, test: {test_biom.shape}")

            tag = f"{site.replace(' ', '_')}_{disease.replace(' ', '_')}_{j}"
            with biom.util.biom_open(f"data/lso_genus_train_{tag}.biom", 'w') as f:
                train_biom.to_hdf5(f, 'example')
            with biom.util.biom_open(f"data/lso_genus_test_{tag}.biom", 'w') as f:
                test_biom.to_hdf5(f, 'example')

# ===== notebook cell 15 =====
run_info = pd.DataFrame({"site": run_site, "disease": run_disease, "study": run_test_study})

# ===== notebook cell 16 =====
run_info.head()

# ===== notebook cell 17 =====
run_info.to_csv("run_info_genus.csv", index=None)

# ===== notebook cell 18 =====
metadata = pd.read_csv(os.path.join(BASE, "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv"), sep="\t", low_memory=False)
metadata = metadata.set_index("sample_id")
metadata['group'] = [0 if i == "HC" else 1 for i in metadata.status_cleaned.values]
metadata.to_csv("meta_genus.tsv", sep="\t")

# ===== notebook cell 22 =====
table = biom.load_table(os.path.join(BASE, "data/oral-gut-pair-v2-qc/biom_oral_gut_pair_v2_qc.biom"))
table = aggregate_to_genus(table)
metadata = pd.read_csv(os.path.join(BASE, "data/oral-gut-pair-v2-qc/metadata_oral_gut_pair_v2_merge_qc.csv"), low_memory=False)
sid = metadata.sample_id.values
table.filter(sid, axis="sample")
table.remove_empty()
# fid = table.ids(axis="observation")
# prevalence = table.nonzero_counts(axis="observation")
# table = table.filter(fid[prevalence > 50], axis='observation', inplace=False)
# table.remove_empty()
sid_reads = table.sum(axis="sample")
sid = table.ids(axis="sample")
table = table.filter(sid[sid_reads > 1000], axis='sample', inplace=False)
table.remove_empty()

sid = table.ids(axis="sample")
metadata = metadata.loc[[i in sid for i in metadata.sample_id.values]]

# ===== notebook cell 23 =====
metadata = metadata.loc[metadata.status_cleaned.notna()]
metadata = metadata.loc[metadata.status_cleaned != "Other"]

# ===== notebook cell 24 =====
metadata = metadata.set_index("sample_id")
metadata['group'] = [0 if i == "HC" else 1 for i in metadata.status_cleaned.values]
metadata.to_csv("meta_oral_gut_genus.tsv", sep="\t")

# ===== notebook cell 25 =====
metadata.head()

# ===== notebook cell 27 =====
metadata = metadata.reset_index()

# ===== notebook cell 28 =====
md_summary = pd.DataFrame(metadata.groupby(["pair_bioproject", "site_group_cleaned", "status_cleaned"])["sample_id"].count())
md_summary = md_summary.loc[md_summary.sample_id > 10]

# ===== notebook cell 29 =====
md_summary

# ===== notebook cell 30 =====
project_id = np.unique([i[0] for i in md_summary.index.values])

# ===== notebook cell 31 =====
N_SPLITS = 5
run_site = []
run_study = []
run_fold = []
for i in project_id:
    if md_summary.loc[i].shape[0] != 4:
        continue
    print(i)
    for site in ["Gut", "Oral"]:
        metadata_pick = metadata.loc[metadata.pair_bioproject == i].copy()
        metadata_pick = metadata_pick.loc[metadata_pick.site_group_cleaned == site].copy()
        metadata_pick["label"] = (metadata_pick["status_cleaned"] != "HC").astype(int)
        sid = metadata_pick.sample_id.values
        y   = metadata_pick["label"].values.astype(int)
        cls_count = np.bincount(y, minlength=2)
        if len(np.unique(y)) < 2 or cls_count.min() < N_SPLITS:
            print(f"[SKIP] study {i} / {site}: "
                  f"n={len(y)} pos={int(y.sum())} cannot run {N_SPLITS}-fold stratified CV")
            continue

        ### filter table
        table_pick = table.filter(sid, axis="sample", inplace=False)
        fid = table_pick.ids(axis="observation")
        prevalence = table_pick.nonzero_counts(axis="observation") 
        table_pick = table_pick.filter(fid[prevalence > 2], axis='observation', inplace=False)
        table_pick.norm(axis='sample')
        if table_pick.shape[1] == 0:
            continue
        sid = table_pick.ids(axis="sample")
        
        skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=42)
        for fold, (tr_idx, te_idx) in enumerate(skf.split(sid, y), start=1):
            train_sid,  Y_all  = sid[tr_idx], y[tr_idx].astype(np.float32)
            test_sid, Y_test = sid[te_idx], y[te_idx].astype(np.float32)
            # Test fold not both-class -> AUC undefined, record NaN
            if len(np.unique(Y_test)) < 2:
                print(f"  study {i}/{site} fold{fold}: test fold not both-class, skip")
                continue
            run_site.append(site)
            run_study.append(i)
            run_fold.append(fold)
            train_biom = table_pick.filter(train_sid, axis="sample", inplace=False)
            test_biom = table_pick.filter(test_sid, axis="sample", inplace=False)

            print(f"train table: {train_biom.shape}, test table: {test_biom.shape}")

            with biom.util.biom_open(f"data/five_fold_genus_train_{site.replace(' ', '_')}_{i}_{fold}.biom", 'w') as f:
                train_biom.to_hdf5(f, 'example')
            with biom.util.biom_open(f"data/five_fold_genus_test_{site.replace(' ', '_')}_{i}_{fold}.biom", 'w') as f:
                test_biom.to_hdf5(f, 'example')

# ===== notebook cell 32 =====
run_info = pd.DataFrame({"site": run_site, "fold": run_fold, "study": run_study})
run_info.to_csv("run_info_oral_gut_genus.csv", index=None)

# ===== notebook cell 37 =====
table = biom.load_table(os.path.join(BASE, "data/metadata-and-biom-v2-qc/biom_oral_qc_dedup.biom"))
table = aggregate_to_genus(table)
metadata = pd.read_csv(os.path.join(BASE, "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv"), sep="\t", low_memory=False)
sid = metadata.sample_id.values
table.filter(sid, axis="sample")
table.remove_empty()
sid_reads = table.sum(axis="sample")
sid = table.ids(axis="sample")
table = table.filter(sid[sid_reads > 1000], axis='sample', inplace=False)
table.remove_empty()

sid = table.ids(axis="sample")
metadata = metadata.loc[[i in sid for i in metadata.sample_id.values]]

# ===== notebook cell 38 =====
import numpy as np
import pandas as pd
import biom

# ---------------- Parameters ----------------
META_PATH = os.path.join(BASE, "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv")
CONTINENT_COL  = "continent_group"   # column representing continent/region in metadata
MIN_CASE_TOTAL = 30                  # minimum pooled case samples across matched BPs
MIN_HC_TOTAL   = 30                  # minimum pooled HC samples across matched BPs
MIN_MATCHED_BP = 2                   # minimum number of matched BPs
ABUND_THRESH   = 0.01                # keep a genus if mean relative abundance >= 1% in any cohort
COHORT_COL     = "bioproject_id"     # definition of "cohort" for the abundance filter
RENORM_AFTER_FILTER = False          # whether to renormalize each sample to 1 after genus filtering

oral_sites = ['saliva', 'tooth', 'tongue', 'oral mucosa']

sid_all = sid          # samples present in the biom table (protect from being overwritten)

# run manifest
run_site, run_disease, run_continent, run_test_study, run_split_type = [], [], [], [], []


def clean(x):
    return str(x).replace(' ', '_')


def prep_table(sample_ids, meta_sub, thresh=ABUND_THRESH, cohort_col=COHORT_COL):
    """Subset table -> convert to relative abundance -> keep genera whose mean
    relative abundance is >= thresh in at least one cohort."""
    tp = table.filter(sample_ids, axis="sample", inplace=False)
    tp.norm(axis="sample", inplace=True)

    fid = np.asarray(tp.ids(axis="observation"))
    sids = np.asarray(tp.ids(axis="sample"))
    mat = pd.DataFrame(tp.matrix_data.toarray(), index=fid, columns=sids)  # genus x sample
    cohort = meta_sub.set_index("sample_id").loc[sids, cohort_col].astype(str).values

    # mean relative abundance, genus x cohort
    mean_by_cohort = mat.T.groupby(cohort).mean().T
    keep_fid = fid[(mean_by_cohort.max(axis=1) >= thresh).values]
    print(f"    genus: {len(fid)} -> {len(keep_fid)} "
          f"(mean relative abundance >= {thresh:.0%} in at least one cohort)")
    if len(keep_fid) == 0:
        return None

    tp = tp.filter(keep_fid, axis="observation", inplace=False)
    if RENORM_AFTER_FILTER:
        tp.norm(axis="sample", inplace=True)
    return tp


def dump(train_biom, test_biom, base):
    with biom.util.biom_open(f"data/genus_train_{base}.biom", 'w') as f:
        train_biom.to_hdf5(f, 'example')
    with biom.util.biom_open(f"data/genus_test_{base}.biom", 'w') as f:
        test_biom.to_hdf5(f, 'example')


for site in oral_sites:
    print(f"\n########## SITE: {site} ##########")

    metadata = pd.read_csv(META_PATH, sep="\t", low_memory=False)
    metadata = metadata.loc[metadata.status_cleaned.notna()]
    metadata = metadata.loc[metadata.status_cleaned != "Other"]
    metadata = metadata.loc[metadata.sample_id.isin(set(sid_all))]
    metadata = metadata.loc[metadata.site_cleaned == site]
    metadata = metadata.loc[metadata[CONTINENT_COL].notna()]

    if metadata.empty:
        print("  No samples for this site, skipping...")
        continue

    # one continent per study (majority vote, so a study is never split)
    proj2cont = metadata.groupby("bioproject_id")[CONTINENT_COL].agg(
        lambda x: x.value_counts().index[0])

    # sample counts per (BP, status); no minimum per group within a BP
    counts = (metadata.groupby(["bioproject_id", "status_cleaned"])["sample_id"]
                      .count().reset_index(name="n"))

    for s in sorted(counts.status_cleaned.unique()):
        if s == "HC":
            continue
        print(f"\n=== Disease: {s} ===")

        # ---- 1) Matched BPs: contain both HC and this disease ----
        sub = counts[counts.status_cleaned.isin(["HC", s])]
        matched_bp = (sub.groupby("bioproject_id")["status_cleaned"].nunique()
                         .pipe(lambda x: x[x == 2]).index.tolist())

        if len(matched_bp) < MIN_MATCHED_BP:
            print(f"  Only {len(matched_bp)} matched BP(s) (< {MIN_MATCHED_BP}), skipping...")
            continue

        # ---- 2) Pooled case / HC sample sizes across all matched BPs ----
        sub = sub[sub.bioproject_id.isin(matched_bp)]
        n_case = int(sub.loc[sub.status_cleaned == s, "n"].sum())
        n_hc   = int(sub.loc[sub.status_cleaned == "HC", "n"].sum())
        print(f"  matched BP: {len(matched_bp)} | pooled case={n_case}, HC={n_hc}")

        # if n_case < MIN_CASE_TOTAL or n_hc < MIN_HC_TOTAL:
        #     print(f"  Pooled sample size too small "
        #           f"(need case>={MIN_CASE_TOTAL}, HC>={MIN_HC_TOTAL}), skipping...")
        #     continue

        # ---- 3) Samples of this matched cohort ----
        meta_s = metadata.loc[
            metadata.bioproject_id.isin(matched_bp)
            & metadata.status_cleaned.isin(["HC", s])
        ].copy()

        proj_cont  = {p: proj2cont[p] for p in matched_bp}
        conts      = sorted(set(proj_cont.values()))
        cont2projs = {c: [p for p in matched_bp if proj_cont[p] == c] for c in conts}
        print("  studies per continent: " +
              ", ".join(f"{c}:{len(v)}" for c, v in cont2projs.items()))

        # ===== (A) within-continent: leave-one-study-out inside each continent =====
        for c, projs in cont2projs.items():
            if len(projs) < 2:
                continue
            print(f"  [within] {c}: {len(projs)} studies")

            meta_wc = meta_s[meta_s.bioproject_id.isin(projs)]
            table_pick = prep_table(meta_wc['sample_id'].values, meta_wc)
            if table_pick is None or table_pick.shape[1] == 0:
                continue
            # generate_normalize_tree(table_pick, f"data/phylogeny_normalize_wc_{clean(site)}_{clean(s)}_{clean(c)}.nwk")

            for j in projs:
                train_study = [p for p in projs if p != j]
                train_sid = meta_wc[meta_wc.bioproject_id.isin(train_study)]['sample_id'].values
                test_sid  = meta_wc[meta_wc.bioproject_id == j]['sample_id'].values

                train_biom = table_pick.filter(train_sid, axis="sample", inplace=False)
                test_biom  = table_pick.filter(test_sid,  axis="sample", inplace=False)
                print(f"    test={j}  train {train_biom.shape}  test {test_biom.shape}")

                dump(train_biom, test_biom, f"wc_{clean(site)}_{clean(s)}_{clean(c)}_{j}")

                run_site.append(clean(site)); run_disease.append(clean(s))
                run_continent.append(clean(c)); run_test_study.append(str(j))
                run_split_type.append("wc")

        # ===== (B) leave-one-continent-out, tested study-by-study =====
        # train = all studies outside continent c; test = ONE study inside c (one fold per study)
        if len(conts) < 2:
            print("  Spans only 1 continent, skipping leave-one-continent-out")
            continue

        table_pick = prep_table(meta_s['sample_id'].values, meta_s)
        if table_pick is None or table_pick.shape[1] == 0:
            continue
        # generate_normalize_tree(table_pick, f"data/phylogeny_normalize_loco_{clean(site)}_{clean(s)}.nwk")

        for c in conts:
            train_projs = [p for p in matched_bp if proj_cont[p] != c]
            if len(train_projs) == 0:
                continue
            train_sid  = meta_s[meta_s.bioproject_id.isin(train_projs)]['sample_id'].values
            train_biom = table_pick.filter(train_sid, axis="sample", inplace=False)

            print(f"  [LOCO] held-out continent={c}: "
                  f"train {len(train_projs)} studies {train_biom.shape}, "
                  f"{len(cont2projs[c])} test study/studies")

            for j in cont2projs[c]:
                test_sid  = meta_s[meta_s.bioproject_id == j]['sample_id'].values
                test_biom = table_pick.filter(test_sid, axis="sample", inplace=False)
                print(f"    test-study={j}  test {test_biom.shape}")

                dump(train_biom, test_biom, f"loco_{clean(site)}_{clean(s)}_{clean(c)}_{j}")

                run_site.append(clean(site)); run_disease.append(clean(s))
                run_continent.append(clean(c)); run_test_study.append(str(j))
                run_split_type.append("loco")

runs = pd.DataFrame({
    "site": run_site, "disease": run_disease, "continent": run_continent,
    "test_study": run_test_study, "split_type": run_split_type,
})

# ===== notebook cell 39 =====
runs.to_csv("run_info_continent_genus.csv", index=False)

# ===== notebook cell 43 =====
oral_sites = ['saliva', 'tooth', 'tongue']

# Typical oral diseases for cross-disease transfer (edit as needed)
ORAL_DISEASES = [
    'Dental Caries', 'OSCC', 'Periodontal Disease',
    'Oral Submucous Fibrosis', 'Oral Mucositis',
    'Oral Verrucous Hyperplasia', 'Aphthous', 'Oral Lichen Planus',
    'Head and Neck Cancer',   # includes oral cavity; remove if you want strictly-oral
]

def clean(x):
    return str(x).replace(' ', '_')

# run manifest (now includes site)
run_site, run_train, run_test = [], [], []
run_ntrain, run_ntest, run_ntrain_case, run_ntest_case = [], [], [], []

for site in oral_sites:
    print(f"\n########## SITE: {site} ##########")

    # ---------- Base filtering, same criteria, now per site ----------
    metadata = pd.read_csv(
        os.path.join(BASE, "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv"),
        sep="\t", low_memory=False)
    metadata = metadata.loc[metadata.status_cleaned.notna()]
    metadata = metadata.loc[metadata.status_cleaned != "Other"]
    metadata = metadata.loc[[i in sid for i in metadata.sample_id]]
    metadata = metadata.loc[metadata.site_cleaned == site]     # restrict to this site
    if metadata.empty:
        print("  no samples at this site, skipping")
        continue

    # Each (bioproject, status) needs > 5 samples -- computed AFTER the site filter
    result = metadata.groupby(["bioproject_id", "status_cleaned"])['sample_id'].count().reset_index()
    result = result.query('sample_id > 5')

    # For each disease, keep only studies where BOTH HC and that disease have > 5 samples at this site
    def valid_studies_for(d):
        pick = result[result['status_cleaned'].isin(["HC", d])]
        two = pick.groupby('bioproject_id')['status_cleaned'].nunique()
        return two[two == 2].index.tolist()

    usable = {}
    for d in ORAL_DISEASES:
        projs = valid_studies_for(d)
        if len(projs) >= 1:
            usable[d] = projs

    if len(usable) < 2:
        print(f"  only {len(usable)} usable disease(s) at this site, "
              f"cannot form a cross-disease pair, skipping")
        continue

    print("  usable diseases: " +
          ", ".join(f"{d}({len(p)} studies)" for d, p in usable.items()))

    # Samples (HC + d) restricted to d's valid studies, within this site
    def disease_samples(d):
        m = metadata[metadata.bioproject_id.isin(usable[d]) &
                     metadata.status_cleaned.isin(["HC", d])]
        return m

    # ---------- Cross-disease: all ordered pairs (train on A, test on B) ----------
    diseases = list(usable.keys())
    for A in diseases:
        mA = disease_samples(A)
        train_sid = mA['sample_id'].values
        for B in diseases:
            if A == B:
                continue
            mB = disease_samples(B)
            test_sid = mB['sample_id'].values

            # Guard against leakage: drop any test sample also present in train
            # (only if a study contributes HC to both A and B at this site)
            overlap = set(train_sid) & set(test_sid)
            if overlap:
                print(f"  [warn] {A}->{B}: {len(overlap)} shared samples removed from test")
                test_sid = np.array([x for x in test_sid if x not in overlap])
            if len(test_sid) == 0:
                print(f"  [skip] {A}->{B}: no test samples left")
                continue
                
            ### filter table
            sid_pick = train_sid.tolist() + test_sid.tolist()
            table_pick = table.filter(sid_pick, axis="sample", inplace=False)
            fid = table_pick.ids(axis="observation")
            prevalence = table_pick.nonzero_counts(axis="observation")
            table_pick = table_pick.filter(fid[prevalence > 2], axis='observation', inplace=False)
            table_pick.norm(axis='sample')
            if table_pick.shape[1] == 0:
                continue
        
            # Build biom tables and align features to the training set
            train_biom = table_pick.filter(train_sid, axis="sample", inplace=False)
            test_biom = table_pick.filter(test_sid, axis="sample", inplace=False)

            n_train_case = int((mA.status_cleaned == A).sum())
            n_test_case  = int((metadata[metadata.sample_id.isin(test_sid)]
                                .status_cleaned == B).sum())

            print(f"  [cross-disease] train={A} (HC+{A}, {train_biom.shape})  "
                  f"test={B} (HC+{B}, {test_biom.shape})")

            base = f"cd_{clean(site)}_train_{clean(A)}_test_{clean(B)}"
            with biom.util.biom_open(f"data/train_{base}.biom", 'w') as f:
                train_biom.to_hdf5(f, 'example')
            with biom.util.biom_open(f"data/test_{base}.biom", 'w') as f:
                test_biom.to_hdf5(f, 'example')

            run_site.append(clean(site))
            run_train.append(clean(A)); run_test.append(clean(B))
            run_ntrain.append(train_biom.shape[1]); run_ntest.append(test_biom.shape[1])
            run_ntrain_case.append(n_train_case); run_ntest_case.append(n_test_case)

# Collect all cross-disease splits for the downstream training loop
runs_cd = pd.DataFrame({
    "site": run_site,
    "train_disease": run_train, "test_disease": run_test,
    "n_train": run_ntrain, "n_test": run_ntest,
    "n_train_case": run_ntrain_case, "n_test_case": run_ntest_case,
    "split_type": "cross_disease_within_site",
})

# ===== notebook cell 44 =====
runs_cd.to_csv("run_info_cross_disease_genus.csv", index=False)

# ===== notebook cell 48 =====
import numpy as np
import pandas as pd
import random
import biom
from math import comb
from itertools import combinations

SEED = 42
N_REPEAT = 3

ABUND_THRESH = 0.01          # keep a genus if mean relative abundance >= 1% in any cohort
COHORT_COL   = "bioproject_id"   # definition of "cohort" for the abundance filter
RENORM_AFTER_FILTER = True   # renormalize each sample to sum 1 after genus filtering

# One analysis per (site, disease): large-sample combos only
LC_TARGETS = [
    ("saliva", "OSCC"),
    ("saliva", "Head and Neck Cancer"),
    ("saliva", "Periodontal Disease"),
    ("tooth",  "Periodontal Disease"),
]

def clean(x):
    return str(x).replace(' ', '_')

def draw_distinct_subsets(items, k, n_want, rng):
    """Return up to n_want distinct k-subsets of items (as sorted tuples)."""
    total = comb(len(items), k)
    n = min(n_want, total)
    if total <= 200:                      # small: enumerate then sample
        allc = list(combinations(sorted(items), k))
        rng.shuffle(allc)
        return allc[:n]
    seen = set()                          # large: rejection-sample
    while len(seen) < n:
        seen.add(tuple(sorted(rng.sample(items, k))))
    return list(seen)

def prep_table(sample_ids, meta_sub, thresh=ABUND_THRESH, cohort_col=COHORT_COL):
    """Subset table -> total-sum normalize -> keep genera whose mean relative
    abundance is >= thresh in at least one cohort -> renormalize (total sum)."""
    tp = table.filter(sample_ids, axis="sample", inplace=False)
    tp.norm(axis="sample", inplace=True)

    fid  = np.asarray(tp.ids(axis="observation"))
    sids = np.asarray(tp.ids(axis="sample"))
    mat = pd.DataFrame(tp.matrix_data.toarray(), index=fid, columns=sids)  # genus x sample
    cohort = meta_sub.drop_duplicates("sample_id").set_index("sample_id") \
                     .loc[sids, cohort_col].astype(str).values

    # mean relative abundance, genus x cohort
    mean_by_cohort = mat.T.groupby(cohort).mean().T
    keep_fid = fid[(mean_by_cohort.max(axis=1) >= thresh).values]
    if len(keep_fid) == 0:
        return None

    tp = tp.filter(keep_fid, axis="observation", inplace=False)
    if RENORM_AFTER_FILTER:
        tp.norm(axis="sample", inplace=True)   # total sum normalization
    return tp

# run manifest
run_site, run_disease, run_k, run_rep, run_subset, run_test = [], [], [], [], [], []

for site, disease in LC_TARGETS:
    print(f"\n########## SITE: {site} | DISEASE: {disease} ##########")

    # ---------- Base filtering, same criteria as before ----------
    metadata = pd.read_csv(
        os.path.join(BASE, "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv"),
        sep="\t", low_memory=False)
    metadata = metadata.loc[metadata.status_cleaned.notna()]
    metadata = metadata.loc[metadata.status_cleaned != "Other"]
    metadata = metadata.loc[[i in sid for i in metadata.sample_id]]
    metadata = metadata.loc[metadata.site_cleaned == site]

    result = metadata.groupby(["bioproject_id", "status_cleaned"])['sample_id'].count().reset_index()
    # Valid studies: within one bioproject, both HC and the disease are present
    pick = result[result['status_cleaned'].isin(["HC", disease])]
    two = pick.groupby('bioproject_id')['status_cleaned'].nunique()
    valid_proj = two[two == 2].index.tolist()
    max_k = len(valid_proj)
    print(f"  valid studies: {max_k} -> {valid_proj}")
    if max_k < 2:
        print("  Fewer than 2 studies, cannot build a learning curve, skipping...")
        continue

    # Restrict to HC + disease samples within valid studies
    meta_d = metadata[metadata.bioproject_id.isin(valid_proj) &
                      metadata.status_cleaned.isin(["HC", disease])]

    rng = random.Random(SEED)

    # ---------- Learning curve: k = 2 .. max_k ----------
    for k in range(2, max_k + 1):
        subsets = draw_distinct_subsets(valid_proj, k, N_REPEAT, rng)
        print(f"  k={k}: {len(subsets)} repeat(s)")
        for rep, subset in enumerate(subsets):
            subset = list(subset)

            # genus filter computed on the k sampled studies
            meta_sub = meta_d[meta_d.bioproject_id.isin(subset)]
            table_pick = prep_table(meta_sub['sample_id'].values, meta_sub)
            if table_pick is None or table_pick.shape[1] == 0:
                print(f"    r{rep}: no genus passed the abundance filter, skipping")
                continue
            print(f"    r{rep}: studies={subset}, genus kept={table_pick.shape[0]}")

            # leave-one-study-out inside the sampled k studies
            for j in subset:
                train_study = [p for p in subset if p != j]
                train_sid = meta_sub[meta_sub.bioproject_id.isin(train_study)]['sample_id'].values
                test_sid  = meta_sub[meta_sub.bioproject_id == j]['sample_id'].values

                train_biom = table_pick.filter(train_sid, axis="sample", inplace=False)
                test_biom  = table_pick.filter(test_sid,  axis="sample", inplace=False)

                base = f"lc_{clean(site)}_{clean(disease)}_k{k}_r{rep}_test_{j}"
                with biom.util.biom_open(f"data/genus_train_{base}.biom", 'w') as f:
                    train_biom.to_hdf5(f, 'example')
                with biom.util.biom_open(f"data/genus_test_{base}.biom", 'w') as f:
                    test_biom.to_hdf5(f, 'example')

                run_site.append(clean(site)); run_disease.append(clean(disease))
                run_k.append(k); run_rep.append(rep)
                run_subset.append("|".join(str(p) for p in subset))
                run_test.append(str(j))

# Collect all learning-curve splits for the downstream training loop
runs_lc = pd.DataFrame({
    "site": run_site, "disease": run_disease,
    "n_studies": run_k, "repeat": run_rep,
    "subset": run_subset, "test_study": run_test,
    "split_type": "learning_curve",
})

# ===== notebook cell 49 =====
runs_lc.head()

# ===== notebook cell 50 =====
runs_lc.to_csv("run_info_lc_genus.csv", index=False)

# ===== notebook cell 52 =====

