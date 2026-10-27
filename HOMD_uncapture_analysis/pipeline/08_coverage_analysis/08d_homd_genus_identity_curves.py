#!/usr/bin/env python
"""Per-genus identity curves to eHOMD.

The top-left panel of panel B shows how much of the catalogue as a whole is
matched at or above each best-hit identity.  The same curve can be drawn for
every genus, which says where the coverage gap sits taxonomically rather than
in aggregate.

Only ASVs whose best hit reaches 95 % identity are counted.  For each genus and
each identity x in [95, 100]:

    frac(x) = share of the genus's ASVs with a best hit >= 95 % whose best hit
              has identity >= x

Every curve therefore starts at 100 % at x = 95.  The curves compare how the
ASVs that are matched at all are spread across identity, not how many of them
there are; the genera with no ASV at or above 95 % drop out of the file.

Output: 08_homd_genus_identity_curves.tsv
        genus  identity  frac  n_asvs
"""
import numpy as np
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
OUT = str(RESULTS / "08_homd_genus_identity_curves.tsv")

# ---- HOMD taxonomy ----------------------------------------------------------
tax = pd.read_csv(TAX, sep="\t", header=None, names=["hmt_id", "lineage"])
tax["hmt_id"] = tax["hmt_id"].str.split("_").str[0]
tax["genus"] = tax["lineage"].str.split(";").apply(
    lambda x: x[5] if len(x) > 5 else "unclassified")
hmt2genus = dict(zip(tax["hmt_id"], tax["genus"]))

# ---- ASVs: best hit and chimera flag ---------------------------------------
df = pd.read_csv(HOMD, sep="\t",
                 usecols=["derep_id", "homd_target", "identity"], low_memory=False)
ch = pd.read_csv(CHIM, sep="\t", usecols=["derep_id", "is_chimera"],
                 low_memory=False)
df = df.merge(ch, on="derep_id", how="left")
df = df[df["is_chimera"] == False].copy()

# The clean catalogue of the rest of the figure -- chimeric, contamination- and
# organelle-flagged ASVs removed -- is what panel B's counts refer to; without
# this filter the curves are drawn on a different denominator and their values
# cannot be read against it.
fl = pd.read_csv(FLAG, sep="\t",
                 usecols=["derep_id", "contamination_flag", "organelle_flag"],
                 low_memory=False)
fl["derep_id"] = fl["derep_id"].astype(str)
df["derep_id"] = df["derep_id"].astype(str)
df = df.merge(fl, on="derep_id", how="left")
df["contamination_flag"] = df["contamination_flag"].fillna(False)
df["organelle_flag"] = df["organelle_flag"].fillna(False)
df = df[~(df["contamination_flag"] | df["organelle_flag"])].copy()
print(f"[08d] after contamination/organelle filter: {len(df):,}")

df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
df["genus"] = df["homd_target"].astype(str).str.split("_").str[0].map(hmt2genus)
df = df.dropna(subset=["genus", "identity"])
n_all = len(df)
# the denominator of every curve: only ASVs that reach 95 % identity
df = df[df["identity"] >= 95].copy()
print(f"[08d] {len(df):,} of {n_all:,} non-chimeric ASVs reach 95% identity, "
      f"{df['genus'].nunique()} genera")

# ---- cumulative share at each identity -------------------------------------
xs = np.round(np.arange(95, 100.0001, 0.1), 1)
rows = []
for genus, sub in df.groupby("genus"):
    ident = sub["identity"].to_numpy()
    n = len(ident)
    for x in xs:
        rows.append((genus, x, float((ident >= x).mean()), n))
out = pd.DataFrame(rows, columns=["genus", "identity", "frac", "n_asvs"])
out.to_csv(OUT, sep="\t", index=False, float_format="%.6f")
print(f"[08d] wrote {OUT}: {out['genus'].nunique()} genera x {len(xs)} identities")
