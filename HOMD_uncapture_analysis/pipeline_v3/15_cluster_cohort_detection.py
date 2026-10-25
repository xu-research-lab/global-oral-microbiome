#!/usr/bin/env python
"""Detection of each supported cluster in disease cohorts with >100 subjects.

A cohort is a disease type (status_cleaned) with more than 100 distinct
subjects across all studies that report it.  Its controls are the HC samples
(status_cleaned == 'HC') of the same studies, so the comparison is not
confounded by study, region or amplicon region.

For every cohort and cluster the table records:
  n_samples_disease / n_pos_disease   disease samples positive for the cluster
  n_samples_hc      / n_pos_hc        HC samples positive
  n_samples_absent                    samples (either arm) negative for it

Detection is the cluster's count > 0 in 14_cluster_counts_dedup.tsv.  Cohorts
without HC samples (ICU) are flagged, as is the catch-all status 'Other'; the
plot keeps them out of the main panel.

Output: 15_cluster_cohort_detection.tsv (long: cohort x cluster)
"""
import pandas as pd

PIPE = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
META = f"{PIPE}/14_cluster_sample_metadata.tsv"
COUNTS = f"{PIPE}/14_cluster_counts_dedup.tsv"
OUT = f"{PIPE}/15_cluster_cohort_detection.tsv"

meta = pd.read_csv(META, sep="\t", dtype={"sample_key": str}, low_memory=False)
counts = pd.read_csv(COUNTS, sep="\t", index_col=0)
pos = counts > 0

st = meta["status_cleaned"].fillna("NA")
subj = meta.groupby("status_cleaned")["subject_id"].nunique()
diseases = sorted([s for s in subj.index if s not in ("HC", "NA") and subj[s] > 100])

col_of = {k: i for i, k in enumerate(counts.columns)}
rows = []
for d in diseases:
    d_idx = meta.index[st == d]
    studies = set(meta.loc[d_idx, "bioproject_id"])
    h_idx = meta.index[(st == "HC") & meta["bioproject_id"].isin(studies)]
    d_cols = [col_of[meta.at[i, "sample_key"]] for i in d_idx]
    h_cols = [col_of[meta.at[i, "sample_key"]] for i in h_idx]

    pd_pos = pos.iloc[:, d_cols].sum(axis=1).values
    ph_pos = pos.iloc[:, h_cols].sum(axis=1).values
    nd, nh = len(d_cols), len(h_cols)

    rows.append(pd.DataFrame({
        "cohort": d,
        "cluster_id": counts.index.astype(str),
        "n_samples_disease": nd,
        "n_pos_disease": pd_pos,
        "n_samples_hc": nh,
        "n_pos_hc": ph_pos,
        "n_samples_absent": (nd - pd_pos) + (nh - ph_pos),
        "n_subjects_disease": int(meta.loc[d_idx, "subject_id"].nunique()),
        "n_studies": len(studies),
        "has_hc": nh > 0,
    }))

out = pd.concat(rows, ignore_index=True)
out["n_total"] = out["n_samples_disease"] + out["n_samples_hc"]
# within-arm detection rates: comparable between the two arms
out["rate_disease"] = out["n_pos_disease"] / out["n_samples_disease"]
out["rate_hc"] = out["n_pos_hc"] / out["n_samples_hc"].replace(0, pd.NA)
# shares of the whole cohort (the two arms together): these sum to 1 with
# share_absent and are what the stacked bars are drawn from
out["share_disease"] = out["n_pos_disease"] / out["n_total"]
out["share_hc"] = out["n_pos_hc"] / out["n_total"]
out["share_absent"] = out["n_samples_absent"] / out["n_total"]
out.to_csv(OUT, sep="\t", index=False, float_format="%.8g")

print(f"wrote {OUT}: {len(out)} rows, {out['cohort'].nunique()} cohorts")
print(out.groupby(["cohort", "n_samples_disease", "n_samples_hc", "n_subjects_disease"])
      .size().rename("n_clusters").reset_index().to_string(index=False))
