#!/usr/bin/env python
"""Per-species coverage of each eHOMD species, under the three weightings.

For every eHOMD species that at least one clean ASV hits, the share of that
species' own signal whose best-hit identity reaches 97%:

  frac_asv_ge97    ASVs with identity >= 97 / all ASVs of the species
  frac_reads_ge97  their reads / all reads of the species
  frac_prev_ge97   their subject-level detections / all detections

The same filters and species assignment as 13j_homd_species_coverage_fullname.py,
so the two tables describe the same 801 species; this one adds the read and
prevalence weightings and the species' phylum (for the heatmap track beside the
ranking in 13m_plot_homd_species_metric_dotline.R).
"""
import pandas as pd

TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
HOMD = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv"
CHIM = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv"
FLAG = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv"
OUT = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13m_homd_species_metric_coverage.tsv"

# ---- HOMD taxonomy: reference -> full species name, phylum, reference count --
tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["lineage"] = tax["lineage"].str.split(";")
tax["phylum"] = tax["lineage"].apply(lambda x: x[1] if len(x) > 1 else "unclassified")
tax["genus"] = tax["lineage"].apply(lambda x: x[5] if len(x) > 5 else "unclassified")
tax["species"] = tax["lineage"].apply(lambda x: x[6] if len(x) > 6 else "unclassified")
tax["full"] = tax["genus"] + " " + tax["species"]
hmt2full = dict(zip(tax["hmt_id"], tax["full"]))
sp_phylum = tax.drop_duplicates("full").set_index("full")["phylum"]
sp_refs = tax.groupby("full")["hmt_id"].count()

# ---- clean ASV set, as in 13j ----
df = pd.read_csv(HOMD, sep="\t",
                 usecols=["derep_id", "homd_target", "identity",
                          "total_reads", "subject_prevalence"],
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
df["total_reads"] = pd.to_numeric(df["total_reads"], errors="coerce").fillna(0.0)
df["subject_prevalence"] = pd.to_numeric(df["subject_prevalence"], errors="coerce").fillna(0.0)
df["hmt_id"] = df["homd_target"].astype(str).str.split("_").str[0]
df["species"] = df["hmt_id"].map(hmt2full)
df = df.dropna(subset=["species"])
print(f"clean ASVs with a species: {len(df):,}")

# ---- aggregate per species ----
ge = df["identity"] >= 97
df["ge97"] = ge
g = df.groupby("species")
s = pd.DataFrame({
    "n_asvs_total": g.size(),
    "n_asvs_ge97": g["ge97"].sum(),
    "reads_total": g["total_reads"].sum(),
    "reads_ge97": g.apply(lambda x: x.loc[x["ge97"], "total_reads"].sum()),
    "prev_total": g["subject_prevalence"].sum(),
    "prev_ge97": g.apply(lambda x: x.loc[x["ge97"], "subject_prevalence"].sum()),
})
s["frac_asv_ge97"] = s["n_asvs_ge97"] / s["n_asvs_total"]
s["frac_reads_ge97"] = s["reads_ge97"] / s["reads_total"].replace(0, pd.NA)
s["frac_prev_ge97"] = s["prev_ge97"] / s["prev_total"].replace(0, pd.NA)
s["phylum"] = s.index.map(sp_phylum)
s["n_homd_refseq"] = s.index.map(sp_refs)
s = s.sort_values("frac_asv_ge97", ascending=False)
s.index.name = "species"
s.to_csv(OUT, sep="\t", float_format="%.6f")
print(f"wrote {OUT}: {len(s)} species")
print(f"  ASVs   >=97: {s['n_asvs_ge97'].sum():,} / {s['n_asvs_total'].sum():,}")
print(f"  reads  >=97: {100 * s['reads_ge97'].sum() / s['reads_total'].sum():.2f}%")
print(f"  prev   >=97: {100 * s['prev_ge97'].sum() / s['prev_total'].sum():.2f}%")
print(f"  species with zero reads: {(s['reads_total'] == 0).sum()}; "
      f"zero prevalence: {(s['prev_total'] == 0).sum()}")
