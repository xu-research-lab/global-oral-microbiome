#!/usr/bin/env python
"""Annotate organelle (mitochondrion/chloroplast) flags in 02_contamination_flagged.tsv.

Reads vsearch hits (13b_organelle.hits.tsv) against the GG2 organelle reference,
classifies each hit via GG2 taxonomy, keeps best hit per ASV (identity, then coverage),
and adds organelle_flag / organelle_target / organelle_identity columns.
"""
import pandas as pd
import os

PIPE_V2 = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2"
PIPE_V3 = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
GG2_TAX = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"
HITS = os.path.join(PIPE_V3, "13b_organelle.hits.tsv")
SRC = os.path.join(PIPE_V2, "02_contamination_flagged.tsv")

log = lambda s: print(f"[13c_organelle] {s}", flush=True)

# ---- 1. GG2 taxonomy: ref -> organelle type ----
tax = {}
with open(GG2_TAX) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) == 2:
            tax[p[0]] = p[1]

def org_type(ref):
    t = tax.get(ref, "")
    if "f__Mitochondria" in t:
        return "mitochondrion"
    if "c__Chloroplast" in t:
        return "chloroplast"
    return "unknown"

# ---- 2. Parse vsearch hits (query, target, id, alnlen, ...) ----
hits = pd.read_csv(HITS, sep="\t", header=None,
                   names=["derep_id", "ref", "identity", "alnlen", "mism", "qcov", "tcov"])
hits["identity"] = pd.to_numeric(hits["identity"], errors="coerce")
hits["qcov"] = pd.to_numeric(hits["qcov"], errors="coerce")
hits["organelle_target"] = hits["ref"].map(org_type)
log(f"raw hits: {len(hits):,} | "
    f"chloroplast={(hits['organelle_target']=='chloroplast').sum():,} | "
    f"mitochondrion={(hits['organelle_target']=='mitochondrion').sum():,} | "
    f"unknown={(hits['organelle_target']=='unknown').sum():,}")

# best hit per ASV: highest identity, then coverage
best = (hits.sort_values(["identity", "qcov"], ascending=False)
             .drop_duplicates("derep_id"))
log(f"best-hit ASVs: {len(best):,} | "
    f"chloroplast={(best['organelle_target']=='chloroplast').sum():,} | "
    f"mitochondrion={(best['organelle_target']=='mitochondrion').sum():,}")

# ---- 3. Merge into 02_contamination_flagged.tsv ----
df = pd.read_csv(SRC, sep="\t", low_memory=False)
df["derep_id"] = df["derep_id"].astype(str)
best = best[["derep_id", "organelle_target", "identity"]].rename(
    columns={"organelle_target": "_org_target", "identity": "_org_identity"})
df = df.merge(best, on="derep_id", how="left")
df["organelle_flag"] = df["_org_target"].notna()
df["organelle_target"] = df["_org_target"].fillna("")
df["organelle_identity"] = df["_org_identity"].fillna("")
df = df.drop(columns=["_org_target", "_org_identity"])

# append ";organelle" to flag_type for organelle-flagged ASVs
org_mask = df["organelle_flag"].fillna(False)
ft = df["flag_type"].fillna("clean").astype(str)
df.loc[org_mask, "flag_type"] = ft[org_mask].map(
    lambda s: s if "organelle" in s else f"{s};organelle")

df.to_csv(SRC, sep="\t", index=False)
log(f"updated {SRC}: {len(df):,} rows | "
    f"organelle-flagged: {int(org_mask.sum()):,} "
    f"(chloroplast={int((df['organelle_target']=='chloroplast').sum()):,} | "
    f"mitochondrion={int((df['organelle_target']=='mitochondrion').sum()):,})")
log("Done!")
