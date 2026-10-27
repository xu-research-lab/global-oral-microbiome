#!/usr/bin/env python
"""HOMD phylum/species coverage from full non-chimeric ASVs,
excluding organelle-flagged and contamination-flagged ASVs."""
import pandas as pd
import numpy as np

TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
HOMD = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv"
CHIM = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv"
FLAG = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv"
OUTP = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13d_homd_phylum_coverage.tsv"
OUTS = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13d_homd_species_coverage.tsv"

# ---- 1. HOMD taxonomy ----
tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["lineage"] = tax["lineage"].str.split(";")
tax["phylum"] = tax["lineage"].apply(lambda x: x[1] if len(x) > 1 else "unclassified")
tax["species"] = tax["lineage"].apply(lambda x: x[6] if len(x) > 6 else "unclassified")
hmt2phylum = dict(zip(tax["hmt_id"], tax["phylum"]))
hmt2species = dict(zip(tax["hmt_id"], tax["species"]))
print(f"HOMD taxonomy: {len(tax)} seqs, {tax['phylum'].nunique()} phyla, {tax['species'].nunique()} species")

# ---- 2. Full non-chimeric ASVs + HOMD mapping ----
df = pd.read_csv(HOMD, sep="\t", usecols=["derep_id","homd_target","identity","has_homd_hit"], low_memory=False)
ch = pd.read_csv(CHIM, sep="\t", usecols=["derep_id","is_chimera"], low_memory=False)
df = df.merge(ch, on="derep_id", how="left")
df = df[df["is_chimera"] == False].copy()
print(f"Non-chimeric ASVs: {len(df):,}")

# ---- 3. Merge contamination + organelle flags (from <97% subset) ----
fl = pd.read_csv(FLAG, sep="\t", usecols=["derep_id","contamination_flag","organelle_flag"], low_memory=False)
fl["derep_id"] = fl["derep_id"].astype(str)
df["derep_id"] = df["derep_id"].astype(str)
df = df.merge(fl, on="derep_id", how="left")
df["contamination_flag"] = df["contamination_flag"].fillna(False)
df["organelle_flag"] = df["organelle_flag"].fillna(False)
n_contam = df["contamination_flag"].sum()
n_org = df["organelle_flag"].sum()
df = df[~(df["contamination_flag"] | df["organelle_flag"])].copy()
print(f"Removed contamination-flagged: {n_contam:,} | organelle-flagged: {n_org:,} | kept: {len(df):,}")

# ---- 4. Map taxonomy ----
df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
df["hmt_id"] = df["homd_target"].astype(str).str.split("_").str[0]
df["phylum"] = df["hmt_id"].map(hmt2phylum)
df["species"] = df["hmt_id"].map(hmt2species)

# ---- 5. Per-phylum aggregation ----
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

p = df.groupby("phylum").apply(agg, include_groups=False)
p["n_homd_refseq"] = tax.groupby("phylum")["hmt_id"].count()
p["frac_ge97"] = p["n_asvs_ge97"] / p["n_asvs_total"]
p["frac_ge95"] = p["n_asvs_ge95"] / p["n_asvs_total"]
p["frac_ge90"] = p["n_asvs_ge90"] / p["n_asvs_total"]
p = p.sort_values("n_asvs_total", ascending=False)
p.to_csv(OUTP, sep="\t", float_format="%.4f")
print(f"Wrote {OUTP}: {len(p)} phyla")

# ---- 6. Per-species aggregation ----
s = df.groupby("species").apply(agg, include_groups=False)
s["n_homd_refseq"] = tax.groupby("species")["hmt_id"].count()
s["frac_ge97"] = s["n_asvs_ge97"] / s["n_asvs_total"]
s["frac_ge95"] = s["n_asvs_ge95"] / s["n_asvs_total"]
s["frac_ge90"] = s["n_asvs_ge90"] / s["n_asvs_total"]
s = s.sort_values("n_asvs_total", ascending=False)
s.to_csv(OUTS, sep="\t", float_format="%.4f")
print(f"Wrote {OUTS}: {len(s)} species")

# ---- 7. Summary ----
print(f"\nTotal clean ASVs: {p['n_asvs_total'].sum():,}")
print(f"Total >=97: {p['n_asvs_ge97'].sum():,} ({p['n_asvs_ge97'].sum()/p['n_asvs_total'].sum()*100:.1f}%)")
print("\n=== Phylum table ===")
print(p.to_string())
