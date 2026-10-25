"""Regenerate LSO split tables with a prevalence-based genus filter.

Identical to the LSO part of disease_predict.ipynb / run_splits_genus.py except
that the genus filter is prevalence instead of mean relative abundance:

    original (mean):    keep genus if MEAN relative abundance >= 1% in any one
                        matched bioproject cohort
    here  (prevalence): keep genus if the PROPORTION of samples with a
                        non-zero count >= 1% in any one matched bioproject

Entry gates (matched BP >= 2, pooled case >= 30, pooled HC >= 30) and all other
steps are unchanged; sample sets therefore match run_info_genus.csv exactly.

Outputs (new prefix, existing mean-based files are not touched):
    data/lsoPgenus_train_{site}_{disease}_{study}.biom
    data/lsoPgenus_test_{site}_{disease}_{study}.biom
"""

import biom
import numpy as np
import pandas as pd
import os
BASE = os.environ.get("GLOBAL_ORAL_ROOT", "..")  # repo root (data/ lives there, server-side)

TAXONOMY_PATH = "/beegfs/huangxiaochang/global_oral/taxonomy.tsv"
META_PATH = os.path.join(BASE, "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv")
BIOM_PATH = os.path.join(BASE, "data/metadata-and-biom-v2-qc/biom_oral_qc_dedup.biom")

MIN_CASE_TOTAL = 30
MIN_HC_TOTAL = 30
MIN_MATCHED_BP = 2
PREV_THRESH = 0.01     # non-zero sample proportion >= 1% in any matched cohort
REABUND_AFTER_FILTER = True
ORAL_SITES = ['saliva', 'tooth', 'tongue', 'oral mucosa']


def _parse_genus(taxon_str):
    """Genus label (g__...) from a semicolon-separated taxonomy string.

    Empty GTDB placeholders (a bare 'g__' token, i.e. genus-unassigned) must
    NOT count as a genus: features whose only g__ token is the placeholder are
    treated as genus-unannotated and dropped by the caller, instead of being
    collapsed into one artificial 'g__' bucket."""
    if pd.isna(taxon_str):
        return None
    for p in str(taxon_str).split(";"):
        p = p.strip()
        if p.startswith("g__") and len(p) > 3:
            return p
    return None


_taxonomy = pd.read_csv(TAXONOMY_PATH, sep="\t")
_taxonomy["genus"] = _taxonomy["Taxon"].apply(_parse_genus)
_taxonomy = _taxonomy.dropna(subset=["genus"])
FID_TO_GENUS = dict(zip(_taxonomy["Feature ID"].astype(str), _taxonomy["genus"]))


def aggregate_to_genus(table):
    obs_ids = table.ids(axis="observation")
    keep = [fid for fid in obs_ids if fid in FID_TO_GENUS]
    table = table.filter(keep, axis="observation", inplace=False)
    return table.collapse(lambda id_, md: FID_TO_GENUS[id_], axis="observation")


def get_matched_bp(metadata_site, disease):
    counts = (metadata_site.groupby(["bioproject_id", "status_cleaned"])["sample_id"]
                            .count().reset_index(name="n"))
    sub = counts[counts.status_cleaned.isin(["HC", disease])]
    matched_bp = (sub.groupby("bioproject_id")["status_cleaned"].nunique()
                     .pipe(lambda x: x[x == 2]).index.tolist())
    if len(matched_bp) < MIN_MATCHED_BP:
        return []
    sub = sub[sub.bioproject_id.isin(matched_bp)]
    n_case = int(sub.loc[sub.status_cleaned == disease, "n"].sum())
    n_hc = int(sub.loc[sub.status_cleaned == "HC", "n"].sum())
    if n_case < MIN_CASE_TOTAL or n_hc < MIN_HC_TOTAL:
        return []
    return matched_bp


def genus_keep_by_cohort_prevalence(table_rel, meta_sub, thresh=PREV_THRESH,
                                    cohort_col="bioproject_id"):
    """Keep genera whose non-zero sample proportion is >= thresh in at least
    one cohort (same skeleton as the mean-abundance filter)."""
    fid = np.asarray(table_rel.ids(axis="observation"))
    sid = np.asarray(table_rel.ids(axis="sample"))
    mat = pd.DataFrame(table_rel.matrix_data.toarray(), index=fid, columns=sid)
    cohort = meta_sub.set_index("sample_id").loc[sid, cohort_col].astype(str).values

    bin_mat = (mat > 0).astype(np.float32)          # genus x sample, 1 = detected
    bin_mat = bin_mat.set_axis(cohort, axis=1)
    prop_by_cohort = bin_mat.T.groupby(level=0).mean()   # cohort x genus, prevalence
    keep = (prop_by_cohort.max(axis=0) >= thresh).values   # per-genus max over cohorts
    return fid[keep], prop_by_cohort


def main():
    table = biom.load_table(BIOM_PATH)
    print(f"[genus] aggregating to genus ...", flush=True)
    table = aggregate_to_genus(table)
    metadata = pd.read_csv(META_PATH, sep="\t", low_memory=False)
    sid = metadata.sample_id.values
    table.filter(sid, axis="sample")
    table.remove_empty()
    sid_reads = table.sum(axis="sample")
    sid = table.ids(axis="sample")
    table = table.filter(sid[sid_reads > 1000], axis='sample', inplace=False)
    table.remove_empty()
    sid = table.ids(axis="sample")
    metadata = metadata.loc[[i in sid for i in metadata.sample_id.values]]

    metadata_all = metadata.loc[metadata.status_cleaned.notna()]
    metadata_all = metadata_all.loc[metadata_all.status_cleaned != "Other"]
    metadata_all = metadata_all.loc[metadata_all.sample_id.isin(set(sid))]
    metadata_all = metadata_all.loc[metadata_all.site_cleaned.isin(ORAL_SITES)]

    all_diseases = sorted(metadata_all.loc[metadata_all.status_cleaned != "HC",
                                           "status_cleaned"].unique())
    n_run = []
    clean = lambda x: str(x).replace(' ', '_')   # match run_info_genus.csv convention
    for disease in all_diseases:
        site_matched_bp = {}
        for site in ORAL_SITES:
            metadata_site = metadata_all.loc[metadata_all.site_cleaned == site]
            if metadata_site.empty:
                continue
            matched_bp = get_matched_bp(metadata_site, disease)
            if matched_bp:
                site_matched_bp[site] = matched_bp
        if not site_matched_bp:
            continue
        print(f"\n########## DISEASE: {disease} "
              f"({len(site_matched_bp)} site(s)) ##########", flush=True)

        for site, matched_bp in site_matched_bp.items():
            metadata_site = metadata_all.loc[metadata_all.site_cleaned == site]
            meta_pick = metadata_site.loc[
                metadata_site.bioproject_id.isin(matched_bp)
                & metadata_site.status_cleaned.isin(["HC", disease])
            ].copy()
            sid_pick = meta_pick.sample_id.values
            table_pick = table.filter(sid_pick, axis="sample", inplace=False)
            table_pick.norm(axis="sample", inplace=True)   # zeros preserved by TSS

            keep_fid, prop = genus_keep_by_cohort_prevalence(table_pick, meta_pick)
            print(f"  {disease} | {site}: genus {table_pick.shape[0]} -> "
                  f"{len(keep_fid)} (non-zero proportion >= {PREV_THRESH:.0%} in "
                  f"at least one cohort)", flush=True)
            if len(keep_fid) == 0:
                continue

            table_pick = table_pick.filter(keep_fid, axis="observation", inplace=False)
            if REABUND_AFTER_FILTER:
                table_pick.norm(axis="sample", inplace=True)
            if table_pick.shape[1] == 0:
                continue

            for j in matched_bp:
                train_sid = meta_pick.loc[meta_pick.bioproject_id != j, "sample_id"].values
                test_sid = meta_pick.loc[meta_pick.bioproject_id == j, "sample_id"].values
                train_biom = table_pick.filter(train_sid, axis="sample", inplace=False)
                test_biom = table_pick.filter(test_sid, axis="sample", inplace=False)
                print(f"    [test={j}] train {train_biom.shape}, test {test_biom.shape}",
                      flush=True)
                tag = f"{site.replace(' ', '_')}_{disease.replace(' ', '_')}_{j}"
                with biom.util.biom_open(f"data/lsoPgenus_train_{tag}.biom", 'w') as f:
                    train_biom.to_hdf5(f, 'example')
                with biom.util.biom_open(f"data/lsoPgenus_test_{tag}.biom", 'w') as f:
                    test_biom.to_hdf5(f, 'example')
                n_run.append((clean(site), clean(disease), j))

    print(f"\n[done] {len(n_run)} LSO folds written", flush=True)
    pd.DataFrame(n_run, columns=["site", "disease", "study"]).to_csv(
        "run_info_lsoP_genus.csv", index=False)


if __name__ == "__main__":
    main()