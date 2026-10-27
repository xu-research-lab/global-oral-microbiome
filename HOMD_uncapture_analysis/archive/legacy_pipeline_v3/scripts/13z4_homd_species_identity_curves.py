#!/usr/bin/env python
"""Per-species identity curves to eHOMD, the species-level counterpart of 13z3.

Only ASVs whose best hit reaches 97 % identity -- the species-level boundary --
are counted.  For each eHOMD species (genus + epithet of the best-hit reference)
and each identity x in [97, 100]:

    frac(x) = share of the species' ASVs with a best hit >= 97 % whose best hit
              has identity >= x

Every curve therefore starts at 100 % at x = 97.  The denominator is the
species' own >=97 % matches, so the curves compare how the ASVs that identify
to that species at all are spread across identity -- the sub-species (strain)
structure of the matches -- rather than how many matches the species has.
The panel starts at 97 % because that is the species boundary itself; a best
hit below it is not accepted as the species, so a curve over 79-97 % would be
counting matches that do not identify the species.

The 99 % line marks the strain-level threshold: an ASV at 99-100 % identity to
a reference is conventionally taken to be the same strain, whereas 97-99 %
separates species.

Output: 13z4_homd_species_identity_curves.tsv
        species  identity  frac  n_asvs
"""
import numpy as np
import pandas as pd

TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"
HOMD = ("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/"
        "HOMD_uncapture_analysis/homd_mapping/exact_derep_homd_annotated.tsv")
CHIM = ("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/"
        "HOMD_uncapture_analysis/chimera_removal/exact_derep_chimera_annotated.tsv")
FLAG = ("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/"
        "HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv")
OUT = ("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/"
       "HOMD_uncapture_analysis/pipeline_v3/13z4_homd_species_identity_curves.tsv")

MIN_ID = 97.0

# ---- HOMD taxonomy: HMT id -> species (genus + epithet) --------------------
tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
sp = tax["lineage"].str.split(";")
tax["species"] = [f"{x[5]} {x[6]}" if len(x) > 6 and x[6].strip() else None
                  for x in sp]
hmt2species = dict(zip(tax["hmt_id"], tax["species"]))

# ---- ASVs: best hit and chimera flag ---------------------------------------
df = pd.read_csv(HOMD, sep="\t",
                 usecols=["derep_id", "homd_target", "identity"], low_memory=False)
ch = pd.read_csv(CHIM, sep="\t", usecols=["derep_id", "is_chimera"],
                 low_memory=False)
df = df.merge(ch, on="derep_id", how="left")
df = df[df["is_chimera"] == False].copy()

# same clean catalogue as the rest of the figure
fl = pd.read_csv(FLAG, sep="\t",
                 usecols=["derep_id", "contamination_flag", "organelle_flag"],
                 low_memory=False)
fl["derep_id"] = fl["derep_id"].astype(str)
df["derep_id"] = df["derep_id"].astype(str)
df = df.merge(fl, on="derep_id", how="left")
df["contamination_flag"] = df["contamination_flag"].fillna(False)
df["organelle_flag"] = df["organelle_flag"].fillna(False)
df = df[~(df["contamination_flag"] | df["organelle_flag"])].copy()
print(f"[13z4] clean ASVs: {len(df):,}")

df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
df["species"] = df["homd_target"].astype(str).str.split("_").str[0].map(hmt2species)
df = df.dropna(subset=["species", "identity"])
n_all = len(df)
df = df[df["identity"] >= MIN_ID].copy()
print(f"[13z4] {len(df):,} of {n_all:,} ASVs reach {MIN_ID:.0f}% identity, "
      f"{df['species'].nunique()} species")

# ---- cumulative share at each identity -------------------------------------
xs = np.round(np.arange(MIN_ID, 100.0001, 0.1), 1)
rows = []
for species, sub in df.groupby("species"):
    ident = sub["identity"].to_numpy()
    n = len(ident)
    for x in xs:
        rows.append((species, x, float((ident >= x).mean()), n))
out = pd.DataFrame(rows, columns=["species", "identity", "frac", "n_asvs"])
out.to_csv(OUT, sep="\t", index=False, float_format="%.6f")
print(f"[13z4] wrote {OUT}: {out['species'].nunique()} species x {len(xs)} identities")
