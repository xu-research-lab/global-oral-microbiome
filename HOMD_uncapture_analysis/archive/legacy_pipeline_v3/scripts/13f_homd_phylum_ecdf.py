#!/usr/bin/env python
"""Per-phylum cumulative identity curve (fraction of ASVs with identity >= x)
from full clean (non-chimeric, non-contaminant, non-organelle) derep ASVs."""
import pandas as pd
import numpy as np

TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
HOMD = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv"
CHIM = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv"
FLAG = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv"
OUT = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13f_homd_phylum_ecdf.tsv"

tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["phylum"] = tax["lineage"].str.split(";").apply(lambda x: x[1] if len(x) > 1 else "unclassified")
hmt2phylum = dict(zip(tax["hmt_id"], tax["phylum"]))

df = pd.read_csv(HOMD, sep="\t", usecols=["derep_id", "homd_target", "identity", "has_homd_hit"], low_memory=False)
ch = pd.read_csv(CHIM, sep="\t", usecols=["derep_id", "is_chimera"], low_memory=False)
df = df.merge(ch, on="derep_id", how="left")
df = df[df["is_chimera"] == False].copy()

fl = pd.read_csv(FLAG, sep="\t", usecols=["derep_id", "contamination_flag", "organelle_flag"], low_memory=False)
fl["derep_id"] = fl["derep_id"].astype(str)
df["derep_id"] = df["derep_id"].astype(str)
df = df.merge(fl, on="derep_id", how="left")
df["contamination_flag"] = df["contamination_flag"].fillna(False)
df["organelle_flag"] = df["organelle_flag"].fillna(False)
df = df[~(df["contamination_flag"] | df["organelle_flag"])].copy()

df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
df = df.dropna(subset=["identity"])
df["identity"] = df["identity"].clip(upper=100.0)
df["hmt_id"] = df["homd_target"].astype(str).str.split("_").str[0]
df["phylum"] = df["hmt_id"].map(hmt2phylum)
df = df.dropna(subset=["phylum"])
print(f"clean ASVs with phylum: {len(df):,}")

grid = np.arange(0.0, 100.001, 0.1)
edges = np.concatenate([grid, [100.1]])
rows = []
for ph, g in df.groupby("phylum"):
    n = len(g)
    counts, _ = np.histogram(g["identity"], bins=edges)
    cum_ge = np.cumsum(counts[::-1])[::-1]  # fraction with identity >= grid[i]
    frac_ge = cum_ge / n
    for x, f in zip(grid, frac_ge):
        rows.append((ph, round(x, 1), round(float(f), 6)))
out = pd.DataFrame(rows, columns=["phylum", "identity", "frac_ge"])
out.to_csv(OUT, sep="\t", index=False, float_format="%.6f")
print(f"wrote {OUT}: {len(out)} rows, {out['phylum'].nunique()} phyla")
