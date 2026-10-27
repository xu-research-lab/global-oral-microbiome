#!/usr/bin/env python
"""Per-phylum and total 'fraction with identity < x' curves.

For each identity threshold x:
  frac_lt(metric) = sum(metric of ASVs with identity < x) / sum(metric of all ASVs)
Metrics: n_asvs (count), total_reads, subject_prevalence.

Outputs:
  08_homd_phylum_lt_fractions.tsv  (phylum x 0.1% grid)
  08_homd_total_lt_fractions.tsv   (0.1% grid, all phyla combined)
"""
import pandas as pd
import numpy as np

from pathlib import Path

PIPE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline")
STEP = PIPE / "08_coverage_analysis"
RESULTS = STEP / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

HOMD = str(PIPE / "02_homd_mapping" / "results" / "02_derep_homd_annotated.tsv")
CHIM = str(PIPE / "03_chimera_removal" / "results" / "03_derep_chimera_annotated.tsv")
FLAG = str(PIPE / "07_contamination_firewall" / "results" / "07_contamination_flagged.tsv")

TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
OUT_PHY = str(RESULTS / "08_homd_phylum_lt_fractions.tsv")
OUT_TOT = str(RESULTS / "08_homd_total_lt_fractions.tsv")

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
# Phylum-unassigned ASVs have no eHOMD hit at all (identity 0): they are part of
# the clean catalogue and, being below 97%, belong in the denominator of the
# total curve -- the data-set flow counts them in its <97% bucket.  Only the
# per-phylum curves have to drop them, since they have no phylum to be counted
# under.  Dropping them here as well used to raise the total >=97% share from
# 51.6% to 52.7%, which disagreed with the flow (1,871,022 clean ASVs =
# 965,488 at >=97% + 905,534 below).
print(f"clean ASVs: {len(df):,}; with phylum: {df['phylum'].notna().sum():,}")

grid = np.arange(0.0, 100.001, 0.1)

def lt_fractions(g):
    g = g.sort_values("identity")
    ident = g["identity"].values
    n = len(g)
    cum_n = np.concatenate([[0.0], np.cumsum(np.ones(n))])
    cum_r = np.concatenate([[0.0], np.cumsum(g["total_reads"].values)])
    cum_s = np.concatenate([[0.0], np.cumsum(g["subject_prevalence"].values)])
    idx = np.searchsorted(ident, grid, side="left")  # #ASVs with identity < x
    fn = cum_n[idx] / cum_n[-1]
    fr = cum_r[idx] / cum_r[-1]
    fs = cum_s[idx] / cum_s[-1]
    return grid, fn, fr, fs

# ---- per-phylum ----
rows = []
for ph, g in df.groupby("phylum"):
    x, fn, fr, fs = lt_fractions(g)
    for i in range(len(x)):
        rows.append((ph, round(x[i], 1), round(float(fn[i]), 6),
                     round(float(fr[i]), 6), round(float(fs[i]), 6)))
out_phy = pd.DataFrame(rows, columns=["phylum", "identity",
                                      "frac_n_lt", "frac_r_lt", "frac_s_lt"])
out_phy.to_csv(OUT_PHY, sep="\t", index=False, float_format="%.6f")
print(f"wrote {OUT_PHY}: {len(out_phy)} rows, {out_phy['phylum'].nunique()} phyla")

# ---- total (all phyla) ----
x, fn, fr, fs = lt_fractions(df)
out_tot = pd.DataFrame({"identity": np.round(x, 1),
                        "frac_n_lt": np.round(fn, 6),
                        "frac_r_lt": np.round(fr, 6),
                        "frac_s_lt": np.round(fs, 6)})
out_tot.to_csv(OUT_TOT, sep="\t", index=False, float_format="%.6f")
print(f"wrote {OUT_TOT}: {len(out_tot)} rows")
