#!/usr/bin/env python
"""Empirical identity density curves for the four alternative panels.

Same clean ASV set as 13r (non-chimeric, non-contamination-flagged,
non-organelle-flagged, mapped to a HOMD phylum).

For each group (phylum or 'All phyla combined') and each metric
(ASV count = 1 per derep ASV, Reads = total_reads,
 Subject prevalence = subject_prevalence) the density in every 0.5%
identity bin over 70-100% is:

    density = (metric mass in bin) / (group metric total) / bin_width

so each group x metric curve integrates to 1.
"""
import pandas as pd
import numpy as np

TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
HOMD = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv"
CHIM = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv"
FLAG = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv"
OUT = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13v_homd_phylum_identity_density.tsv"

BW = 0.5
BINS = np.arange(70.0, 100.0, BW)                  # left edges [70, 99.5]
CENT = BINS + BW / 2

# HOMD taxonomy: target -> phylum
tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["phylum"] = tax["lineage"].str.split(";").apply(
    lambda x: x[1] if len(x) > 1 else "unclassified")
hmt2phylum = dict(zip(tax["hmt_id"], tax["phylum"]))

# clean ASV table (same filter chain as 13r)
df = pd.read_csv(HOMD, sep="\t",
                 usecols=["derep_id", "homd_target", "identity",
                          "total_reads", "subject_prevalence"],
                 low_memory=False)
ch = pd.read_csv(CHIM, sep="\t", usecols=["derep_id", "is_chimera"],
                 low_memory=False)
df = df.merge(ch, on="derep_id", how="left")
df = df[df["is_chimera"] == False].copy()

fl = pd.read_csv(FLAG, sep="\t",
                 usecols=["derep_id", "contamination_flag", "organelle_flag"],
                 low_memory=False)
fl["derep_id"] = fl["derep_id"].astype(str)
df["derep_id"] = df["derep_id"].astype(str)
df = df.merge(fl, on="derep_id", how="left")
df["contamination_flag"] = df["contamination_flag"].fillna(False)
df["organelle_flag"] = df["organelle_flag"].fillna(False)
df = df[~(df["contamination_flag"] | df["organelle_flag"])].copy()

df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
df = df.dropna(subset=["identity"])
df["identity"] = df["identity"].clip(upper=100.0)
df["total_reads"] = pd.to_numeric(df["total_reads"], errors="coerce").fillna(0.0)
df["subject_prevalence"] = pd.to_numeric(
    df["subject_prevalence"], errors="coerce").fillna(0.0)

df["hmt_id"] = df["homd_target"].astype(str).str.split("_").str[0]
df["phylum"] = df["hmt_id"].map(hmt2phylum)
df = df.dropna(subset=["phylum"])
df = df[(df["identity"] >= 70) & (df["identity"] <= 100)]

metrics = {"ASV count": np.ones(len(df)),
           "Reads": df["total_reads"].values,
           "Subject prevalence": df["subject_prevalence"].values}

rows = []
groups = [("phylum", ph) for ph in sorted(df["phylum"].unique())]
groups += [("total", "All phyla combined")]

for kind, name in groups:
    g = df if kind == "total" else df[df["phylum"] == name]
    if len(g) == 0:
        continue
    idx = np.digitize(g["identity"].values, BINS, right=False) - 1
    idx = np.clip(idx, 0, len(CENT) - 1)
    for metric, w in metrics.items():
        ww = w if kind == "total" else g["total_reads"].values
        # rebuild metric weights for this subset
        if metric == "ASV count":
            ww = np.ones(len(g))
        elif metric == "Reads":
            ww = g["total_reads"].values
        else:
            ww = g["subject_prevalence"].values
        tot = ww.sum()
        if tot <= 0:
            continue
        counts = np.bincount(idx, weights=ww, minlength=len(CENT))
        dens = counts / tot / BW
        for i, x in enumerate(CENT):
            rows.append((name, metric, round(float(x), 1),
                         round(float(dens[i]), 8)))

out = pd.DataFrame(rows, columns=["group", "metric", "identity", "density"])
out.to_csv(OUT, sep="\t", index=False, float_format="%.8f")
print(f"ASVs with 70-100% identity: {len(df):,}")
print(f"wrote {OUT}: {len(out)} rows")
