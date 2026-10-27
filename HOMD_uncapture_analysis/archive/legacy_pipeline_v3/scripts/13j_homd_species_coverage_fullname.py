#!/usr/bin/env python
"""Per-species (full Genus species name) HOMD coverage from full clean
(non-chimeric, non-contaminant, non-organelle) derep ASVs."""
import pandas as pd
import numpy as np

TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
HOMD = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv"
CHIM = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv"
FLAG = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv"
OUTS = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13j_homd_species_coverage_fullname.tsv"

tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["lineage"] = tax["lineage"].str.split(";")
tax["genus"] = tax["lineage"].apply(lambda x: x[5] if len(x) > 5 else "unclassified")
tax["species"] = tax["lineage"].apply(lambda x: x[6] if len(x) > 6 else "unclassified")
tax["full"] = tax["genus"] + " " + tax["species"]
hmt2full = dict(zip(tax["hmt_id"], tax["full"]))
print(f"HOMD taxonomy: {len(tax)} seqs, {tax['full'].nunique()} full species names")

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
df["hmt_id"] = df["homd_target"].astype(str).str.split("_").str[0]
df["species"] = df["hmt_id"].map(hmt2full)
df = df.dropna(subset=["species"])
print(f"clean ASVs with species: {len(df):,}")

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

s = df.groupby("species").apply(agg, include_groups=False)
s["n_homd_refseq"] = tax.groupby("full")["hmt_id"].count()
s["frac_ge97"] = s["n_asvs_ge97"] / s["n_asvs_total"]
s["frac_ge95"] = s["n_asvs_ge95"] / s["n_asvs_total"]
s["frac_ge90"] = s["n_asvs_ge90"] / s["n_asvs_total"]
s = s.sort_values("n_asvs_total", ascending=False)
s.to_csv(OUTS, sep="\t", float_format="%.4f")
print(f"Wrote {OUTS}: {len(s)} species")
print(f"Total: {s['n_asvs_total'].sum():,} | >=97: {s['n_asvs_ge97'].sum():,}")
