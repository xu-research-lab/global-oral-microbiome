#!/usr/bin/env python3
"""
Step 05 — Organelle removal, part 2: classify the organelle-reference hits.

Reads the vsearch hits against the GG2 organelle reference, classifies each
reference sequence through the GG2 taxonomy (mitochondrion / chloroplast /
unclassified), and keeps the best hit per candidate ASV (identity, then query
coverage).

A candidate with any hit at the search threshold counts as organelle-flagged,
including hits to references the taxonomy does not place.  This is the
coverage screen: step 08 removes these ASVs from the HOMD coverage statistics.
It is deliberately separate from the cluster-gate organelle mapping used in
step 12, which has its own frozen input.

Outputs (in 05_organelle_removal/results/):
  05_organelle_hits.tsv      raw vsearch hits (from 05b)
  05_organelle_flagged.tsv   best hit per ASV: derep_id, organelle_target, identity
  05_organelle_excluded.tsv  the same rows with sequence and hit identity
  05_organelle_summary.json  counts
"""
import json
import os
from pathlib import Path

import pandas as pd

PIPE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline")
STEP = PIPE / "05_organelle_removal"
RESULTS = STEP / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

GG2_TAX = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"
HITS = RESULTS / "05_organelle_hits.tsv"
CANDIDATES = PIPE / "04_candidate_selection" / "results" / "04_candidates.tsv"

log = lambda s: print(f"[05_organelle] {s}", flush=True)

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

# ---- 3. Write the flag table used downstream ----
flagged = best[["derep_id", "organelle_target", "identity"]].copy()
flagged.to_csv(RESULTS / "05_organelle_flagged.tsv", sep="\t", index=False)

# ---- 4. Write the excluded ASVs with their sequence ----
cand = pd.read_csv(CANDIDATES, sep="\t", usecols=["derep_id", "sequence"], low_memory=False)
excluded = flagged.merge(cand, on="derep_id", how="left")[
    ["derep_id", "organelle_target", "identity", "sequence"]]
excluded = excluded.rename(columns={"identity": "organelle_identity"})
excluded.to_csv(RESULTS / "05_organelle_excluded.tsv", sep="\t", index=False)

# ---- 5. Summary ----
summary = {
    "raw_hits": int(len(hits)),
    "flagged_asvs": int(len(flagged)),
    "chloroplast": int((flagged["organelle_target"] == "chloroplast").sum()),
    "mitochondrion": int((flagged["organelle_target"] == "mitochondrion").sum()),
    "unknown_target": int((flagged["organelle_target"] == "unknown").sum()),
    "missing_sequence": int(excluded["sequence"].isna().sum()),
}
with open(RESULTS / "05_organelle_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

log(f"wrote 05_organelle_flagged.tsv: {len(flagged):,} ASVs")
log(f"wrote 05_organelle_excluded.tsv: {len(excluded):,} ASVs")
log(f"summary: {RESULTS / '05_organelle_summary.json'}")
log("Done!")
