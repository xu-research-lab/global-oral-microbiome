#!/usr/bin/env python
"""Count full clean ASVs whose best HOMD identity is 70 <= identity < 97,
by the phylum of the best-hit HOMD reference species.

Uses exactly the same filtered HOMD best-hit table that feeds 13f / 13i
(non-chimeric, non-contamination-flagged, non-organelle-flagged ASVs).
Output is consumed by 13i_plot_homd_phylum_identity_facet.R for the
70-97% stacked-bar panel.
"""
import pandas as pd

TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
HOMD = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv"
CHIM = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv"
FLAG = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv"
OUT = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13u_homd_lt97_ge70_phylum.tsv"
OUT_GT97 = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13u_homd_gt97_phylum.tsv"

# ---- HOMD taxonomy: best-hit reference -> phylum ----
tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["phylum"] = tax["lineage"].str.split(";").apply(lambda x: x[1] if len(x) > 1 else "unclassified")
hmt2phylum = dict(zip(tax["hmt_id"], tax["phylum"]))

# ---- full clean ASVs (same filter chain as 13d/13f) ----
df = pd.read_csv(HOMD, sep="\t",
                 usecols=["derep_id", "homd_target", "identity", "has_homd_hit"],
                 low_memory=False)
ch = pd.read_csv(CHIM, sep="\t", usecols=["derep_id", "is_chimera"], low_memory=False)
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

# ---- 70 <= identity < 97, phylum of best-hit species ----
sub = df[(df["identity"] >= 70) & (df["identity"] < 97)].copy()
sub["hmt_id"] = sub["homd_target"].astype(str).str.split("_").str[0]
sub["phylum"] = sub["hmt_id"].map(hmt2phylum)
sub = sub.dropna(subset=["phylum"])

cnt = sub.groupby("phylum", observed=True)["derep_id"].count().reset_index()
cnt.columns = ["phylum", "n_asvs_70_97"]
cnt = cnt.sort_values("n_asvs_70_97", ascending=False)

print(f"ASVs with 70 <= identity < 97 and mapped phylum: {cnt['n_asvs_70_97'].sum():,}")
print(cnt.to_string(index=False))
cnt.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT}: {len(cnt)} phyla")

# ---- identity > 97% (complementary stacked bar) ----
gt97 = df[(df["identity"] > 97) & (df["identity"] <= 100)].copy()
gt97["hmt_id"] = gt97["homd_target"].astype(str).str.split("_").str[0]
gt97["phylum"] = gt97["hmt_id"].map(hmt2phylum)
gt97 = gt97.dropna(subset=["phylum"])

cnt_gt = gt97.groupby("phylum", observed=True)["derep_id"].count().reset_index()
cnt_gt.columns = ["phylum", "n_asvs_gt97"]
cnt_gt = cnt_gt.sort_values("n_asvs_gt97", ascending=False)
print(f"ASVs with identity > 97% and mapped phylum: {cnt_gt['n_asvs_gt97'].sum():,}")
cnt_gt.to_csv(OUT_GT97, sep="\t", index=False)
print(f"wrote {OUT_GT97}: {len(cnt_gt)} phyla")
