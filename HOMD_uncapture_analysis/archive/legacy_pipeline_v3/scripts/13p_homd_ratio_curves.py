#!/usr/bin/env python
"""Per-phylum and total identity ratio curves from full clean derep ASVs.

For each identity threshold x:
  ratio(metric) = sum(metric of ASVs with identity >= x)
                / sum(metric of ASVs with identity <  x)
Metrics: n_asvs (count), total_reads, subject_prevalence.

Outputs:
  13p_homd_phylum_ratio_curves.tsv  (phylum x 0.1% grid)
  13p_homd_total_ratio_curves.tsv   (0.1% grid, all phyla combined)
"""
import pandas as pd
import numpy as np

TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
HOMD = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv"
CHIM = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv"
FLAG = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv"
OUT_PHY = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13p_homd_phylum_ratio_curves.tsv"
OUT_TOT = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13p_homd_total_ratio_curves.tsv"

# ---- taxonomy: HMT id -> phylum ----
tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["phylum"] = tax["lineage"].str.split(";").apply(lambda x: x[1] if len(x) > 1 else "unclassified")
hmt2phylum = dict(zip(tax["hmt_id"], tax["phylum"]))

# ---- full derep ASV HOMD mapping ----
df = pd.read_csv(HOMD, sep="\t",
                 usecols=["derep_id", "homd_target", "identity",
                          "total_reads", "subject_prevalence"],
                 low_memory=False)
print(f"raw derep ASVs: {len(df):,}")

# ---- chimera filter ----
ch = pd.read_csv(CHIM, sep="\t", usecols=["derep_id", "is_chimera"], low_memory=False)
df = df.merge(ch, on="derep_id", how="left")
df = df[df["is_chimera"] == False].copy()
print(f"after chimera filter: {len(df):,}")

# ---- contamination / organelle filter ----
fl = pd.read_csv(FLAG, sep="\t",
                 usecols=["derep_id", "contamination_flag", "organelle_flag"],
                 low_memory=False)
fl["derep_id"] = fl["derep_id"].astype(str)
df["derep_id"] = df["derep_id"].astype(str)
df = df.merge(fl, on="derep_id", how="left")
df["contamination_flag"] = df["contamination_flag"].fillna(False)
df["organelle_flag"] = df["organelle_flag"].fillna(False)
df = df[~(df["contamination_flag"] | df["organelle_flag"])].copy()
print(f"after contamination/organelle filter: {len(df):,}")

# ---- numeric cleaning ----
df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
df = df.dropna(subset=["identity"])
df["identity"] = df["identity"].clip(upper=100.0)
df["total_reads"] = pd.to_numeric(df["total_reads"], errors="coerce").fillna(0.0)
df["subject_prevalence"] = pd.to_numeric(df["subject_prevalence"], errors="coerce").fillna(0.0)

# ---- phylum mapping ----
df["hmt_id"] = df["homd_target"].astype(str).str.split("_").str[0]
df["phylum"] = df["hmt_id"].map(hmt2phylum)
df = df.dropna(subset=["phylum"])
print(f"clean ASVs with phylum: {len(df):,}")

grid = np.arange(0.0, 100.001, 0.1)

def ratio_curves(g):
    g = g.sort_values("identity")
    ident = g["identity"].values
    n = len(g)
    cum_n = np.concatenate([[0.0], np.cumsum(np.ones(n))])
    cum_r = np.concatenate([[0.0], np.cumsum(g["total_reads"].values)])
    cum_s = np.concatenate([[0.0], np.cumsum(g["subject_prevalence"].values)])
    idx = np.searchsorted(ident, grid, side="left")  # #ASVs with identity < x
    n_lt = cum_n[idx]; n_ge = cum_n[-1] - cum_n[idx]
    r_lt = cum_r[idx]; r_ge = cum_r[-1] - cum_r[idx]
    s_lt = cum_s[idx]; s_ge = cum_s[-1] - cum_s[idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        rn = np.where(n_lt > 0, n_ge / n_lt, np.nan)
        rr = np.where(r_lt > 0, r_ge / r_lt, np.nan)
        rs = np.where(s_lt > 0, s_ge / s_lt, np.nan)
    return grid, rn, rr, rs

# ---- per-phylum ----
rows = []
for ph, g in df.groupby("phylum"):
    x, rn, rr, rs = ratio_curves(g)
    for i in range(len(x)):
        rows.append((ph, round(x[i], 1), round(float(rn[i]), 6),
                     round(float(rr[i]), 6), round(float(rs[i]), 6)))
out_phy = pd.DataFrame(rows, columns=["phylum", "identity",
                                      "ratio_n_asvs", "ratio_reads", "ratio_prevalence"])
out_phy.to_csv(OUT_PHY, sep="\t", index=False, float_format="%.6f")
print(f"wrote {OUT_PHY}: {len(out_phy)} rows, {out_phy['phylum'].nunique()} phyla")

# ---- total (all phyla) ----
x, rn, rr, rs = ratio_curves(df)
out_tot = pd.DataFrame({"identity": np.round(x, 1),
                        "ratio_n_asvs": np.round(rn, 6),
                        "ratio_reads": np.round(rr, 6),
                        "ratio_prevalence": np.round(rs, 6)})
out_tot.to_csv(OUT_TOT, sep="\t", index=False, float_format="%.6f")
print(f"wrote {OUT_TOT}: {len(out_tot)} rows")
