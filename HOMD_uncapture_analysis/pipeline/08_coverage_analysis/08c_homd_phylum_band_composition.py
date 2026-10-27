#!/usr/bin/env python
"""Phylum composition of each identity band, under all three weightings.

The two identity bands are the ones the coverage panel already uses:

    <97%  band : 70 <= identity < 97
    >=97% band : identity > 97 and identity <= 100

Within a band the composition is the share of that band's own total, so the
three weightings give three different pictures of who is under-covered:

    ASVs       - one per non-redundant ASV sequence
    Reads      - sum of total_reads over those ASVs
    Prevalence - sum of subject_prevalence over those ASVs

The same filtered HOMD best-hit table feeds the coverage panel and curves: non-chimeric,
non-contamination-flagged, non-organelle-flagged ASVs only, best-hit reference
phylum taken from eHOMD V16.03 taxonomy field 2.

Only the five phyla with the most ASVs overall are named; the rest are pooled
into "Other phyla".  The named set is decided once, on ASV counts, so that the
same five phyla appear in all three facets and the facets stay comparable.
"""
import pandas as pd

from pathlib import Path

PIPE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline")
STEP = PIPE / "08_coverage_analysis"
RESULTS = STEP / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

HOMD = str(PIPE / "02_homd_mapping" / "results" / "02_derep_homd_annotated.tsv")
CHIM = str(PIPE / "03_chimera_removal" / "results" / "03_derep_chimera_annotated.tsv")
FLAG = str(PIPE / "07_contamination_firewall" / "results" / "07_contamination_flagged.tsv")

TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
OUT = str(RESULTS / "08_homd_phylum_band_composition.tsv")

N_NAMED = 5

# ---- HOMD taxonomy: best-hit reference -> phylum ----
tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["phylum"] = tax["lineage"].str.split(";").apply(lambda x: x[1] if len(x) > 1 else "unclassified")
hmt2phylum = dict(zip(tax["hmt_id"], tax["phylum"]))

# ---- the same filter chain as 13d / 13r / 13u ----
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
df["identity"] = df["identity"].clip(upper=100.0)
df["total_reads"] = pd.to_numeric(df["total_reads"], errors="coerce").fillna(0.0)
df["subject_prevalence"] = pd.to_numeric(df["subject_prevalence"], errors="coerce").fillna(0.0)

df["hmt_id"] = df["homd_target"].astype(str).str.split("_").str[0]
df["phylum"] = df["hmt_id"].map(hmt2phylum)
df = df.dropna(subset=["phylum"]).copy()
print(f"clean ASVs with an assignable phylum: {len(df):,}")

# ---- band membership -------------------------------------------------------
df["band"] = pd.NA
df.loc[(df["identity"] >= 70) & (df["identity"] < 97), "band"] = "<97% identity"
df.loc[(df["identity"] > 97) & (df["identity"] <= 100), "band"] = ">=97% identity"
banded = df.dropna(subset=["band"]).copy()

# the panel's bands deliberately leave two slivers out; report them so the
# numbers in the figure are never read as the whole catalogue
lo = int((df["identity"] < 70).sum())
eq = int((df["identity"] == 97).sum())
print(f"outside both bands: identity < 70 -> {lo:,} | identity == 97 -> {eq:,}")
for b in ("<97% identity", ">=97% identity"):
    n = int((banded["band"] == b).sum())
    print(f"  {b:<15} {n:,} ASVs")

# ---- the five largest phyla, decided on ASV counts, once ----
by_rank = banded.groupby("phylum", observed=True).size().sort_values(ascending=False)
NAMED = list(by_rank.index[:N_NAMED])
print(f"named phyla: {', '.join(NAMED)}")
print(f"pooled into Other phyla: {by_rank.index.size - N_NAMED} phyla, "
      f"{int(by_rank.iloc[N_NAMED:].sum()):,} ASVs")

banded["phylum_grp"] = banded["phylum"].where(banded["phylum"].isin(NAMED), "Other phyla")

# ---- composition per metric per band ---------------------------------------
# the ASV weighting is a row count: carry it as an explicit 1.0 column rather
# than summing derep_id, which is a string
banded["one"] = 1.0
LONG = banded.melt(id_vars=["band", "phylum_grp"],
                   value_vars=["one", "total_reads", "subject_prevalence"],
                   var_name="metric", value_name="value")
LONG["value"] = pd.to_numeric(LONG["value"], errors="coerce").fillna(0.0)
LONG["metric"] = LONG["metric"].map({"one": "ASVs", "total_reads": "Reads",
                                     "subject_prevalence": "Prevalence"})

comp = (LONG.groupby(["metric", "band", "phylum_grp"], observed=True)["value"]
            .sum().reset_index())
comp["share"] = comp["value"] / comp.groupby(["metric", "band"])["value"].transform("sum")

# a fixed stack order, held across every bar: ranked by overall ASV counts, so
# the bars can be read against each other without re-finding each phylum
order = NAMED + ["Other phyla"]
comp["phylum_grp"] = pd.Categorical(comp["phylum_grp"], categories=order, ordered=True)
comp["band"] = pd.Categorical(comp["band"], categories=["<97% identity", ">=97% identity"],
                              ordered=True)
comp = comp.sort_values(["metric", "band", "phylum_grp"])

print("\nshare of each band, %:")
print(comp.assign(pct=(100 * comp["share"]).round(2))
          .pivot_table(index=["metric", "band"], columns="phylum_grp",
                       values="pct", observed=True).to_string())

comp.to_csv(OUT, sep="\t", index=False, float_format="%.6f")
print(f"\nwrote {OUT}: {len(comp)} rows")
