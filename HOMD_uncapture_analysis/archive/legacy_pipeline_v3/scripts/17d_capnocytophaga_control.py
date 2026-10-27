#!/usr/bin/env python
"""Control figure for cluster 21907: prevalence of every ASV that maps to HOMD at
>=97% identity and is annotated to the genus Capnocytophaga, in the OSCC cohort.

Cluster 21907 is an uncaptured Capnocytophaga (its best HOMD hit is below 97%),
so a reader needs the contrast: what does the genus look like for the part of it
eHOMD does cover?  This script takes all non-chimeric ASVs whose best HOMD hit
is a Capnocytophaga reference at >=97% identity, and asks in what fraction of
the OSCC cases and the healthy controls any of them is detected, and what the
per-ASV prevalences look like.

Cohort, detection and denominators are those of 17b_plot_cluster_prevalence.R:
the OSCC samples plus the HC samples of the same bioprojects, full cohort (no
covariate cut), detected = at least one read.  Per-sample counts come from the
same per-project DADA2 tables that 14_cluster_feature_profile.py aggregates
(observation ids are the ASV sequences).

Output
  17d_capnocytophaga_control_samples.tsv   group-level detection
  17d_capnocytophaga_control_asvs.tsv      per-ASV detection, both arms
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import biom

ROOT = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"
PIPE = os.path.join(ROOT, "pipeline_v3")
HOMD = os.path.join(ROOT, "homd_mapping", "exact_derep_homd_annotated.tsv")
CHIM = os.path.join(ROOT, "chimera_removal", "exact_derep_chimera_annotated.tsv")
TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
META = os.path.join(PIPE, "14_cluster_sample_metadata.tsv")
ASV_ROOT = "/beegfs/dongbiao/oral"
STATUS = "OSCC"
GENUS = "Capnocytophaga"
MIN_IDENT = 97.0

OUT_SAMPLES = os.path.join(PIPE, "17d_capnocytophaga_control_samples.tsv")
OUT_ASVS = os.path.join(PIPE, "17d_capnocytophaga_control_asvs.tsv")
OUT_PERSAMPLE = os.path.join(PIPE, "17d_capnocytophaga_control_persample.tsv")

# ---- 1. HOMD taxonomy: reference -> genus ---------------------------------
tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["genus"] = tax["lineage"].str.split(";").apply(
    lambda x: x[5] if len(x) > 5 else "unclassified")
cap_refs = set(tax.loc[tax["genus"] == GENUS, "hmt_id"])
print(f"[17d] HOMD: {len(cap_refs)} {GENUS} reference sequences")

# ---- 2. dereps whose best HOMD hit is Capnocytophaga at >=97% --------------
keep = []
scanned = 0
for chunk in pd.read_csv(HOMD, sep="\t",
                         usecols=["derep_id", "sequence", "homd_target", "identity"],
                         dtype={"derep_id": str, "sequence": str, "homd_target": str},
                         chunksize=200_000, low_memory=False):
    scanned += len(chunk)
    hmt = chunk["homd_target"].astype(str).str.split("_").str[0]
    ident = pd.to_numeric(chunk["identity"], errors="coerce")
    m = (ident >= MIN_IDENT) & hmt.isin(cap_refs)
    if m.any():
        sub = chunk.loc[m, ["derep_id", "sequence"]].copy()
        sub["hmt_id"] = hmt[m].values
        sub["identity"] = ident[m].values
        keep.append(sub)
hits = pd.concat(keep, ignore_index=True)
print(f"[17d] {len(hits):,} dereps >= {MIN_IDENT:.0f}% to {GENUS} "
      f"(of {scanned:,} scanned)")

# chimeras are not part of the ASV set used everywhere else in the pipeline
chim = []
for chunk in pd.read_csv(CHIM, sep="\t", usecols=["derep_id", "is_chimera"],
                         dtype={"derep_id": str}, chunksize=500_000, low_memory=False):
    chim.append(chunk)
chim = pd.concat(chim, ignore_index=True)
hits = hits.merge(chim, on="derep_id", how="left")
n_chim = int((hits["is_chimera"] == True).sum())
hits = hits[hits["is_chimera"] != True].reset_index(drop=True)
print(f"[17d] {n_chim:,} chimeric dereps dropped, {len(hits):,} remain")

seq_to_derep = dict(zip(hits["sequence"], hits["derep_id"]))
print(f"[17d] {len(seq_to_derep):,} distinct sequences")

# ---- 3. cohort ------------------------------------------------------------
meta = pd.read_csv(META, sep="\t", dtype=str, low_memory=False)
meta["sample_key"] = meta["bioproject_id"] + "|" + meta["sample_id"]
di = meta.index[meta["status_cleaned"] == STATUS]
hi = meta.index[(meta["status_cleaned"] == "HC") &
                meta["bioproject_id"].isin(meta.loc[di, "bioproject_id"])]
coh = meta.loc[sorted(set(di) | set(hi))].copy()
coh["group"] = np.where(coh["status_cleaned"] == "HC", "HC", "Disease")
print(f"[17d] cohort: {int((coh.group=='Disease').sum()):,} disease / "
      f"{int((coh.group=='HC').sum()):,} HC, "
      f"{coh.bioproject_id.nunique()} projects")

# ---- 4. detection per ASV per sample --------------------------------------
samples = coh["sample_key"].tolist()
col_idx = {s: j for j, s in enumerate(samples)}
det = np.zeros((len(hits), len(samples)), dtype=bool)
row_idx = {d: i for i, d in enumerate(hits["derep_id"])}
covered_reads = np.zeros(len(samples), dtype=np.int64)
total_reads = np.zeros(len(samples), dtype=np.int64)

for bp, g in coh.groupby("bioproject_id"):
    path = os.path.join(ASV_ROOT, bp, "results", "DADA2.biom")
    if not os.path.exists(path):
        raise SystemExit(f"missing biom: {path}")
    table = biom.load_table(path)
    obs = [str(x) for x in table.ids(axis="observation")]
    sid = [str(x) for x in table.ids(axis="sample")]

    rows, ridx = [], []
    for r, s in enumerate(obs):
        d = seq_to_derep.get(s)
        if d is not None:
            rows.append(r)
            ridx.append(row_idx[d])
    if not rows:
        print(f"[17d] {bp}: no qualifying ASV in this table")
        continue

    gcols, lcols = [], []
    for k, s in enumerate(sid):
        j = col_idx.get(f"{bp}|{s}")
        if j is not None:
            gcols.append(j)
            lcols.append(k)
    mat = table.matrix_data.tocsr()
    sub = mat[rows][:, lcols].tocoo()
    # any read, not a depth threshold: this is the prevalence definition 17b uses
    det[np.asarray(ridx)[sub.row], np.asarray(gcols)[sub.col]] = True
    # the same denominator 14_cluster_feature_profile.py uses for relative
    # abundance: every read in the sample's ASV table, not the metadata's count
    covered_reads[np.asarray(gcols)] += np.asarray(
        mat[rows][:, lcols].sum(axis=0)).ravel().astype(np.int64)
    total_reads[np.asarray(gcols)] += np.asarray(
        mat[:, lcols].sum(axis=0)).ravel().astype(np.int64)
    print(f"[17d] {bp}: {len(rows):,} ASVs, {len(gcols)} cohort samples")

# ---- 5. summaries ---------------------------------------------------------
is_d = (coh["group"] == "Disease").to_numpy()
n_d, n_h = int(is_d.sum()), int((~is_d).sum())

any_det = det.any(axis=0)
derep_d = det[:, is_d].sum(axis=1)
derep_h = det[:, ~is_d].sum(axis=1)

srow = pd.DataFrame({
    "group": ["Disease", "HC"],
    "n_samples": [n_d, n_h],
    "n_detected": [int(any_det[is_d].sum()), int(any_det[~is_d].sum())],
})
srow["prevalence"] = srow["n_detected"] / srow["n_samples"]
n_asv_per = det.sum(axis=0)
srow["n_asv_per_sample_median"] = [float(np.median(n_asv_per[is_d])),
                                   float(np.median(n_asv_per[~is_d]))]
srow.to_csv(OUT_SAMPLES, sep="\t", index=False, float_format="%.6g")
print(srow.to_string(index=False))

rel = np.divide(covered_reads, total_reads, out=np.zeros(len(samples)),
                where=total_reads > 0)
pd.DataFrame({"sample_key": samples, "group": coh["group"].to_numpy(),
              "bioproject_id": coh["bioproject_id"].to_numpy(),
              "n_asv": n_asv_per,
              "covered_reads": covered_reads,
              "total_reads": total_reads,
              "rel_abund": rel}).to_csv(OUT_PERSAMPLE, sep="\t", index=False,
                                        float_format="%.6g")

asv = hits[["derep_id", "hmt_id", "identity"]].copy()
asv["n_disease"] = derep_d
asv["n_hc"] = derep_h
asv["prev_disease"] = derep_d / n_d
asv["prev_hc"] = derep_h / n_h
asv = asv.sort_values(["n_disease", "n_hc", "derep_id"], ascending=False)
asv.to_csv(OUT_ASVS, sep="\t", index=False, float_format="%.6g")

seen = asv[(asv.n_disease + asv.n_hc) > 0]
print(f"[17d] ASVs detected in >=1 cohort sample: {len(seen):,} / {len(asv):,}")
print(f"[17d]   higher in disease: {int((seen.prev_disease > seen.prev_hc).sum()):,}; "
      f"higher in HC: {int((seen.prev_disease < seen.prev_hc).sum()):,}; equal: "
      f"{int((seen.prev_disease == seen.prev_hc).sum()):,}")
print(f"[17d] wrote {OUT_SAMPLES}\n[17d] wrote {OUT_ASVS}")
