#!/usr/bin/env python
"""Per-species fraction of ASVs with best HOMD identity < 79%,
from the same clean (non-chimeric, non-contaminant, non-organelle)
derep ASV set used for the species coverage table (13j).
Output adds n_asvs_lt79 / frac_lt79 alongside n_homd_refseq."""
import pandas as pd

TAX  = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
HOMD = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv"
CHIM = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv"
FLAG = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv"
J13  = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13j_homd_species_coverage_fullname.tsv"
OUT  = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13o_homd_species_lt79.tsv"

tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["lineage"] = tax["lineage"].str.split(";")
tax["genus"]  = tax["lineage"].apply(lambda x: x[5] if len(x) > 5 else "unclassified")
tax["species"]= tax["lineage"].apply(lambda x: x[6] if len(x) > 6 else "unclassified")
tax["full"]   = tax["genus"] + " " + tax["species"]
hmt2full = dict(zip(tax["hmt_id"], tax["full"]))
n_refseq = tax.groupby("full")["hmt_id"].count()

df = pd.read_csv(HOMD, sep="\t", usecols=["derep_id", "homd_target", "identity"], low_memory=False)
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

g = df.groupby("species")["identity"]
s = pd.DataFrame({"n_asvs_total": g.size(),
                  "n_asvs_lt79": g.apply(lambda x: (x < 79).sum())})
s["frac_lt79"] = s["n_asvs_lt79"] / s["n_asvs_total"]
s["n_homd_refseq"] = n_refseq
s = s.sort_values("n_asvs_total", ascending=False)
s.to_csv(OUT, sep="\t", float_format="%.6f")
print(f"Wrote {OUT}: {len(s)} species")
print(f"Total ASVs: {s['n_asvs_total'].sum():,} | lt79: {s['n_asvs_lt79'].sum():,} ({100*s['n_asvs_lt79'].sum()/s['n_asvs_total'].sum():.2f}%)")
# cross-check with 13j
j = pd.read_csv(J13, sep="\t", usecols=["species", "n_asvs_total"])
m = s.reset_index()[["species", "n_asvs_total"]].merge(j, on="species", suffixes=("_new", "_13j"))
print("species match:", (m["n_asvs_total_new"] == m["n_asvs_total_13j"]).all())
