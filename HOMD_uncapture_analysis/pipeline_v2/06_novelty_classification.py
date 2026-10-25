#!/usr/bin/env python3
"""
pipeline_v2 Step 06: Novelty Classification — placement-clade + HOMD catalog.

Rewritten 2026-08-13. Novelty-to-HOMD is a BINARY question: does the nearest
reference species appear in HOMD? Since 2026-08-14 the species is the 05d
PER-ASV TREE-DISTANCE ANCHOR (rep ASV's nearest-tip species at D<=0.03), not
the 05c clade-majority species of the placement edge (which misclassified novel
ASVs embedded in HOMD-dominated clades).

  homd_captured    — phylo_species ∈ HOMD species set  OR  ≥97% identity to an
                     HOMD reference sequence (safety net for low-LWR placements)
                     → the phylotype's nearest species IS in the oral catalog.
                     ALSO: a resolved GTDB species whose epithet (species part)
                     maps to exactly ONE HOMD species — GTDB 2022.10 reclassified
                     many genera, so the same organism is often named differently
                     (e.g. 'Eubacterium yurii' = HOMD 'Peptoanaerobacter yurii'
                     HMT-106). These are captured too and flagged by the new
                     `renamed_in_homd` column. Epithets that hit multiple HOMD
                     species are ambiguous and stay uncaptured.
  homd_uncaptured  — otherwise → the phylotype's nearest species is NOT in HOMD
                     (a named non-HOMD species, or unresolved) → candidate new
                     lineage

The auxiliary `oral_clade` column (phylo_genus_base ∈ HOMD genus set) records
whether the phylotype still falls in an HOMD genus (a divergent variant of a
known oral genus) or outside it.

No comparison against Greengenes2 / SILVA sequence databases is performed:
the GG2 reference is already in play via the SEPP placement, and the decision
that matters is HOMD membership of the nearest species.

Species names are matched to HOMD headers ("Genus species", e.g. HMT-432
Treponema vincentii) after stripping GTDB genus tokens (e.g.
"Treponema_B_986142 vincentii" → "Treponema vincentii").

Subject-level prevalence (distinct (bioproject, subject_id) individuals) drives
the candidate gate, replacing the sample-based count.

Inputs:
  - 05_phylotypes_all.tsv               : phylotypes (+ max_subject_prevalence)
  - 05c_phylotype_phylo_annotation.tsv  : placement-clade genus/species per phylotype
  - 06_novelty/hits_homd_iddef2.tsv     : vsearch vs HOMD (safety net; cached)
  - HOMD_16S_rRNA_RefSeq_V16.03_full.fasta

Outputs:
  - 06_phylotypes_novelty_classified.tsv
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
NOVELTY_DIR = PIPE_V2 / "06_novelty"
os.makedirs(NOVELTY_DIR, exist_ok=True)
log = lambda s: print(f"[06_novelty] {s}", flush=True)

VSEARCH = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"
THREADS = 26
ID_THRESHOLD = 97.0
IDDEF = "2"

HOMD_FASTA = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta"

# ============================================================
# 1. HOMD species / genus sets (from reference headers)
# ============================================================
log("Loading HOMD catalog (genus + species sets)...")
homd_genus = set()
homd_species = set()
homd_epithets = {}   # species epithet -> set of HOMD species carrying it
n_homd_refs = 0
with open(HOMD_FASTA) as f:
    for line in f:
        if line.startswith(">"):
            parts = line.split()
            if len(parts) >= 3:
                n_homd_refs += 1
                homd_genus.add(parts[2].strip())
                if len(parts) >= 4:
                    sp = f"{parts[2].strip()} {parts[3].strip()}"
                    homd_species.add(sp)
                    homd_epithets.setdefault(parts[3].strip(), set()).add(sp)
log(f"  HOMD refs: {n_homd_refs:,} | genera: {len(homd_genus):,} | species: {len(homd_species):,} | epithets: {len(homd_epithets):,}")


def clean_species(s):
    """GG2 s__ ('Treponema_B_986142 vincentii') → 'Genus species', or None.

    Placeholder species (spNNN, e.g. 'Dialister sp001553355') are NOT names →
    return None (treated as unresolved → novel)."""
    if not isinstance(s, str) or not s:
        return None
    parts = s.split(" ")
    if len(parts) < 2:
        return None
    genus = parts[0].split("_")[0]
    sp = " ".join(parts[1:]).strip()
    if re.fullmatch(r"sp\d+", sp):
        return None
    return f"{genus} {sp}".strip()


def clean_genus(s):
    """GG2 g__ ('Leptotrichia_A_993758') → 'Leptotrichia' (strip GTDB strain token)."""
    if not isinstance(s, str) or not s:
        return None
    return s.split("_")[0]


# ============================================================
# 2. Load phylotypes + placement-clade annotation
# ============================================================
log("Loading phylotypes...")
phylo = pd.read_csv(PIPE_V2 / "05_phylotypes_all.tsv", sep="\t", low_memory=False)
log(f"  {len(phylo):,} phylotypes")

ann = pd.read_csv(PIPE_V2 / "05c_phylotype_phylo_annotation.tsv", sep="\t", low_memory=False)
log(f"  placement annotation: {len(ann):,} phylotypes")
phylo = phylo.merge(ann[["phylotype_id", "phylo_genus", "phylo_species", "median_lwr"]],
                    on="phylotype_id", how="left")
# NOVELTY uses the 05d per-ASV species anchor (rep ASV's nearest-tip annotation),
# NOT the 05c clade-majority species (2026-08-14, 05e re-clustering).
# Species-level phylotypes carry the clean species name; novel phylotypes
# (D>0.03, or unnamed nearest tip) carry a BLANK anchor, so they can never be
# spuriously captured just because the placement-edge clade is HOMD-dominated.
# (05c phylo_genus/phylo_species are still kept as label columns for downstream.)
phylo["phylo_species_clean"] = phylo["anchor_species"].replace("", np.nan)
phylo["phylo_genus_base"] = phylo["anchor_genus"]
log(f"  with resolvable species name (05d anchor): {int(phylo['phylo_species_clean'].notna().sum()):,}")

# ============================================================
# 3. HOMD sequence search (safety net for low-LWR placements) — cached
# ============================================================
log("\nHOMD vsearch (iddef 2)...")
hits_path = NOVELTY_DIR / "hits_homd_iddef2.tsv"
if not hits_path.exists():
    log("  running vsearch vs HOMD...")
    passed = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False)
    seq_map = dict(zip(passed["derep_id"].astype(str), passed["sequence"].astype(str)))
    rep_fa = NOVELTY_DIR / "phylotype_reps.fasta"
    with open(rep_fa, "w") as f:
        for _, row in phylo.iterrows():
            seq = seq_map.get(str(row["rep_derep_id"]), "")
            if seq:
                f.write(f">{row['rep_derep_id']}\n{seq}\n")
    subprocess.run([
        VSEARCH, "--usearch_global", str(rep_fa),
        "--db", HOMD_FASTA,
        "--id", "0.0", "--iddef", IDDEF,
        "--maxaccepts", "1", "--maxrejects", "32",
        "--strand", "both",
        "--userout", str(hits_path),
        "--userfields", "query+target+id+alnlen+qcov",
        "--threads", str(THREADS), "--quiet"
    ], check=True)
    log("  HOMD search done")
else:
    log(f"  using cached: {hits_path}")

homd_hits = {}
with open(hits_path) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
        qid, target, identity = parts[0], parts[1], float(parts[2])
        if qid not in homd_hits:
            homd_hits[qid] = {"id": identity, "target": target}
log(f"  {len(homd_hits):,} phylotypes with an HOMD hit")

# ============================================================
# 4. Assign novelty classification
# ============================================================
log("\nAssigning novelty (placement-clade species vs HOMD)...")

def species_renamed_in_homd(sp):
    """GTDB species 'G epithet' vs HOMD: GTDB 2022.10 reclassified many genera, so
    the same organism can be listed in HOMD under a different genus name (e.g.
    'Eubacterium yurii' = HOMD 'Peptoanaerobacter yurii' HMT-106). If the species
    epithet maps to EXACTLY ONE HOMD species, treat the GTDB name as that HOMD
    species → homd_captured. Epithets hitting multiple HOMD species (e.g.
    'hongkongensis' → Leptotrichia/Schaalia) are ambiguous → NOT auto-captured."""
    if not isinstance(sp, str):
        return False
    parts = sp.split(" ")
    if len(parts) < 2:
        return False
    epithet = parts[-1].strip()
    return len(homd_epithets.get(epithet, ())) == 1


def classify_novelty(row):
    rid = row["rep_derep_id"]
    sp = row["phylo_species_clean"]
    hid = homd_hits.get(rid, {}).get("id", np.nan)
    hid = float(hid) if pd.notna(hid) else 0.0

    if sp is not None and sp in homd_species:
        cls = "homd_captured"          # nearest species on placement edge ∈ HOMD
        renamed = False
    elif hid >= ID_THRESHOLD:
        cls = "homd_captured"          # ≥97% to HOMD ref sequence (low-LWR safety net)
        renamed = False
    elif sp is not None and species_renamed_in_homd(sp):
        cls = "homd_captured"          # same organism, GTDB genus name ≠ HOMD genus name
        renamed = True
    else:
        cls = "homd_uncaptured"        # nearest species NOT in HOMD (or unresolved)
        renamed = False

    return pd.Series({
        "homd_best_id": round(hid, 2) if pd.notna(hid) else np.nan,
        "homd_target":  homd_hits.get(rid, {}).get("target", ""),
        "novelty_class": cls,
        "oral_clade":   bool(row["phylo_genus_base"] in homd_genus),
        "renamed_in_homd": renamed,
    })

nov = phylo.apply(classify_novelty, axis=1)
phylo_novel = pd.concat([phylo, nov], axis=1)

# ============================================================
# 5. Candidate novel gates (subject-based prevalence)
# ============================================================
log("\nApplying candidate novel gates...")
prev_col = "max_subject_prevalence"
if prev_col in phylo_novel.columns:
    phylo_novel[prev_col] = pd.to_numeric(phylo_novel[prev_col], errors="coerce").fillna(0)
    subjects_gate = phylo_novel[prev_col] >= 10
    log(f"  subject gate: {prev_col} >= 10 (distinct (bioproject,subject_id))")
else:
    subjects_gate = pd.Series(True, index=phylo_novel.index)
    log("  WARNING: max_subject_prevalence column missing; subject gate skipped")

novelty_gate = phylo_novel["novelty_class"] == "homd_uncaptured"
phylo_novel["candidate_novel"] = subjects_gate & novelty_gate
log(f"  candidate novel phylotypes: {int(phylo_novel['candidate_novel'].sum()):,} / {len(phylo_novel):,}")

# ============================================================
# 6. Summary
# ============================================================
log(f"\n=== Novelty Classification Summary ===")
log(f"Total phylotypes: {len(phylo_novel):,}")
for cls in ["homd_captured", "homd_uncaptured"]:
    n = (phylo_novel["novelty_class"] == cls).sum()
    log(f"  {cls}: {n:,}")
n_renamed = int(phylo_novel["renamed_in_homd"].sum())
log(f"  (of homd_captured, via GTDB→HOMD genus-synonym epithet: {n_renamed:,})")
log(f"\nnovelty_class x oral_clade:")
print(pd.crosstab(phylo_novel["novelty_class"], phylo_novel["oral_clade"]).to_string())

cn = phylo_novel[phylo_novel["candidate_novel"]]
log(f"\nCandidate novel ({len(cn):,}):")
log(f"  in HOMD genus clade (oral_clade): {int(cn['oral_clade'].sum()):,} | outside HOMD genus: {int((~cn['oral_clade']).sum()):,}")
log(f"  with resolved species name (non-HOMD): {int(cn['phylo_species_clean'].notna().sum()):,}")

# ============================================================
# 7. Save
# ============================================================
out = PIPE_V2 / "06_phylotypes_novelty_classified.tsv"
phylo_novel.to_csv(out, sep="\t", index=False)
log(f"Saved: {out} ({len(phylo_novel):,} rows)")

summary = {
    "total_phylotypes": int(len(phylo_novel)),
    "by_class": {c: int((phylo_novel["novelty_class"] == c).sum())
                 for c in ["homd_captured", "homd_uncaptured"]},
    "candidate_novel": int(phylo_novel["candidate_novel"].sum()),
    "subject_gate": f"{prev_col} >= 10",
    "method": "05d tree-distance species anchor vs HOMD species set (+ HOMD sequence safety net)",
}
with open(PIPE_V2 / "06_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
log("Saved: 06_summary.json")
log("Done!")
