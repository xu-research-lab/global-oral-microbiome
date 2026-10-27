#!/usr/bin/env python
"""HOMD genus/species coverage from the FULL derep ASV set (non-chimeric)."""
import pandas as pd
import numpy as np

TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
HOMD = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv"
CHIM = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv"
OUTG = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13a_homd_genus_coverage.tsv"
OUTS = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13a_homd_species_coverage.tsv"

# ---- 1. HOMD taxonomy ----
tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["lineage"] = tax["lineage"].str.split(";")
tax["genus"] = tax["lineage"].apply(lambda x: x[5] if len(x) > 5 else "unclassified")
tax["species"] = tax["lineage"].apply(lambda x: x[6] if len(x) > 6 else "unclassified")
hmt2genus = dict(zip(tax["hmt_id"], tax["genus"]))
hmt2species = dict(zip(tax["hmt_id"], tax["species"]))
print(f"HOMD taxonomy: {len(tax)} seqs, {tax['genus'].nunique()} genera, {tax['species'].nunique()} species")

# ---- 2. Load full HOMD mapping + chimera flag ----
df = pd.read_csv(HOMD, sep="\t", usecols=["derep_id","homd_target","identity","has_homd_hit"], low_memory=False)
ch = pd.read_csv(CHIM, sep="\t", usecols=["derep_id","is_chimera"], low_memory=False)
df = df.merge(ch, on="derep_id", how="left")
df = df[df["is_chimera"] == False].copy()
print(f"Non-chimeric ASVs: {len(df):,}")
df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
df["hmt_id"] = df["homd_target"].astype(str).str.split("_").str[0]
df["genus"] = df["hmt_id"].map(hmt2genus)
df["species"] = df["hmt_id"].map(hmt2species)

# ---- 3. Per-genus aggregation ----
def agg(x):
    n = len(x)
    return pd.Series({
        "n_asvs_total": n,
        "n_asvs_ge97": (x["identity"] >= 97).sum(),
        "n_asvs_ge95": (x["identity"] >= 95).sum(),
        "n_asvs_ge90": (x["identity"] >= 90).sum(),
        "n_asvs_lt90": (x["identity"] < 90).sum(),
        "identity_median": x["identity"].median(),
        "identity_q25": x["identity"].quantile(0.25),
        "identity_q75": x["identity"].quantile(0.75),
    })

g = df.groupby("genus").apply(agg)
g["n_homd_refseq"] = tax.groupby("genus")["hmt_id"].count()
g["frac_ge97"] = g["n_asvs_ge97"] / g["n_asvs_total"]
g["frac_ge95"] = g["n_asvs_ge95"] / g["n_asvs_total"]
g["frac_ge90"] = g["n_asvs_ge90"] / g["n_asvs_total"]
g = g.sort_values("n_asvs_total", ascending=False)
g.to_csv(OUTG, sep="\t", float_format="%.4f")
print(f"Wrote {OUTG}: {len(g)} genera")

# ---- 4. Per-species aggregation ----
s = df.groupby("species").apply(agg)
s["n_homd_refseq"] = tax.groupby("species")["hmt_id"].count()
s["frac_ge97"] = s["n_asvs_ge97"] / s["n_asvs_total"]
s["frac_ge95"] = s["n_asvs_ge95"] / s["n_asvs_total"]
s["frac_ge90"] = s["n_asvs_ge90"] / s["n_asvs_total"]
s = s.sort_values("n_asvs_total", ascending=False)
s.to_csv(OUTS, sep="\t", float_format="%.4f")
print(f"Wrote {OUTS}: {len(s)} species")

# ---- 5. Summary ----
print(f"\nTotal non-chimeric ASVs: {g['n_asvs_total'].sum():,}")
print(f"Total >=97: {g['n_asvs_ge97'].sum():,} ({g['n_asvs_ge97'].sum()/g['n_asvs_total'].sum()*100:.1f}%)")
print(f"Total >=95: {g['n_asvs_ge95'].sum():,}")
print(f"Total >=90: {g['n_asvs_ge90'].sum():,}")
print(f"Genera with frac_ge97==0: {(g['frac_ge97']==0).sum()} / {len(g)}")
print(f"Genera with frac_ge97>=0.5: {(g['frac_ge97']>=0.5).sum()} / {len(g)}")
print(f"Genera with frac_ge97>=0.9: {(g['frac_ge97']>=0.9).sum()} / {len(g)}")
print("\n=== Top 20 genera by ASV count ===")
print(g.head(20).to_string())
