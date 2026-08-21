#!/usr/bin/env python3
"""
pipeline_v2 Step 06b: Authenticity Tiering.

Replaces the strict residence hard-gate with a tiered confidence scheme that works
across ALL 334 BioProjects, because only 20 BioProjects contain blank controls and
~86% of V2 phylotypes are never testable by the oral-vs-blank residence test.

Tiers (priority order):
  Tier 1  residence-validated   — >=1 member ASV passed the oral-vs-blank residence
                                  test (reuses pipeline/04_oral_vs_blank results).
                                  Strongest evidence: empirically enriched in oral
                                  vs blank within blank-bearing studies.
  Tier 3  concern               — tested-but-failed residence (member ASVs appeared
                                  in blank-bearing studies but none passed the oral-vs-
                                  blank test) OR flagged by the contamination firewall.
  Tier 2  reproducibility-validated — untested for residence, not flagged, and passes
                                  the authenticity gate below (blank-free criteria).
  below   below threshold       — untested and fails the authenticity gate (weak tail).

Authenticity gate (Tier 2, blank-free; rationale: kit contaminants are lab-specific
and do not independently recur across sequencing centers, Salter 2014):
  not any_contamination
  AND n_bioprojects_total >= 2            (>=2 independent BioProjects)
  AND >=2 of {n_bioprojects_total >= 3,   (independent labs/centers)
              max_prevalence    >= 50,    (independent subjects)
              n_continents      >= 2}     (independent continents/geography)

A Tier-1 phylotype whose residence-passing member is ITSELF flagged is marked
tier1_flag_conflict=1 for transparency (evidence conflict), but stays Tier 1
(empirical oral enrichment is treated as the stronger signal).

Outputs:
  - 06b_phylotypes_tiered.tsv   : all phylotypes with tier columns
  - 06b_summary.json            : tier distribution + candidate breakdown + sensitivity
"""
import json
import os
import sys
from pathlib import Path

import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
STRICT_RES = BASE / "pipeline" / "04_oral_vs_blank"
log = lambda s: print(f"[06b_tier] {s}", flush=True)

# ============================================================
# 1. Load V2 phylotypes
# ============================================================
log("Loading V2 novelty-classified phylotypes...")
phy = pd.read_csv(PIPE_V2 / "06_phylotypes_novelty_classified.tsv", sep="\t", low_memory=False)
log(f"  {len(phy):,} phylotypes")

# ============================================================
# 2. Load residence test data (strict pipeline)
# ============================================================
log("Loading residence test data (strict pipeline)...")
ab = pd.read_csv(STRICT_RES / "abundance_long.tsv", sep="\t", low_memory=False)
eval_set = set(ab["feature_id"].astype(str))
res = pd.read_csv(STRICT_RES / "oral_vs_blank_enrichment.tsv", sep="\t", low_memory=False)
res_pass = set(res.loc[res["residence_pass"] == True, "feature_id"].astype(str))  # noqa: E712
log(f"  evaluated ASVs: {len(eval_set):,} | residence-pass ASVs: {len(res_pass):,}")

# ============================================================
# 3. Load contamination-flagged ASVs (V2 step 02)
# ============================================================
flagged_set = set()
fl_path = PIPE_V2 / "02_contamination_flagged.tsv"
if fl_path.exists():
    fl = pd.read_csv(fl_path, sep="\t", low_memory=False, dtype=str)
    # file contains ALL screened candidates; keep only rows actually flagged
    flag_col = "contamination_flag" if "contamination_flag" in fl.columns else None
    if flag_col:
        flagged_set = set(fl.loc[fl[flag_col].astype(str).str.lower().isin(["true", "1", "t"]), "derep_id"].astype(str))
    else:
        flagged_set = set(fl["derep_id"].astype(str))
    log(f"  V2 contamination-flagged ASVs: {len(flagged_set):,}")
else:
    log(f"  WARNING: {fl_path} not found; using phylotype-level any_contamination only")

# ============================================================
# 4. Assign tiers
# ============================================================
log("\nAssigning tiers...")

# numeric coercion
for c in ["n_bioprojects_total", "max_prevalence", "n_continents"]:
    if c in phy.columns:
        phy[c] = pd.to_numeric(phy[c], errors="coerce").fillna(0)

rows = []
for _, row in phy.iterrows():
    members = [m for m in str(row.get("member_derep_ids", "")).split(";") if m]
    n_tested = sum(1 for m in members if m in eval_set)
    n_resid = sum(1 for m in members if m in res_pass)
    # Tier-1 conflict: the residence-passing member is itself contamination-flagged
    flag_conflict = 0
    if n_resid >= 1:
        resid_members = [m for m in members if m in res_pass]
        if any(m in flagged_set for m in resid_members):
            flag_conflict = 1

    # contamination signal at phylotype level (V2 step 02)
    any_contam = bool(row.get("any_contamination", False))
    if not any_contam and any(m in flagged_set for m in members):
        any_contam = True

    bp = int(row.get("n_bioprojects_total", 0) or 0)
    subj = int(row.get("max_prevalence", 0) or 0)
    cont = int(row.get("n_continents", 0) or 0)

    if n_resid >= 1:
        tier, reason = 1, "residence_pass_member"
    elif n_tested >= 1:
        tier, reason = 3, "tested_failed_residence"
    elif any_contam:
        tier, reason = 3, "contamination_flagged"
    elif bp >= 2 and (int(bp >= 3) + int(subj >= 50) + int(cont >= 2)) >= 2:
        tier, reason = 2, "reproducibility_validated"
    else:
        tier, reason = 0, "below_threshold"

    # evidence annotation: the path that validated/excluded each phylotype
    evidence = {1: "residence", 2: "reproducibility",
                3: "failed_residence", 0: "insufficient_support"}[tier]

    rows.append({
        "phylotype_id": row["phylotype_id"],
        "tier": tier,
        "evidence": evidence,
        "tier_label": {1: "Tier 1 residence-validated",
                       2: "Tier 2 reproducibility-validated",
                       3: "Tier 3 concern",
                       0: "below threshold"}[tier],
        "tier_reason": reason,
        "tier1_flag_conflict": flag_conflict,
        "n_members": len(members),
        "n_members_tested": n_tested,
        "n_members_residence_pass": n_resid,
        "candidate_novel": bool(row.get("candidate_novel", False)),
        "novelty_class": row.get("novelty_class", ""),
        "n_bioprojects_total": bp,
        "max_prevalence": subj,
        "n_continents": cont,
        "any_contamination": any_contam,
    })

t = pd.DataFrame(rows)
t["final_novel"] = t["candidate_novel"] & t["tier"].isin([1, 2])

# ============================================================
# 5. Summary
# ============================================================
log("\n=== Tier distribution (all phylotypes) ===")
for k in [1, 2, 3, 0]:
    n = (t["tier"] == k).sum()
    log(f"  {k} ({t.loc[t['tier']==k,'tier_label'].iloc[0] if n else ''}): {n:,}")

log("\n=== Tier distribution (candidate novel, 7343) ===")
cn = t[t["candidate_novel"]]
for k in [1, 2, 3, 0]:
    n = (cn["tier"] == k).sum()
    log(f"  Tier {k}: {n:,} ({n/len(cn)*100:.1f}%)")

log(f"\n  final_novel (candidate & tier in {1,2}): {int(t['final_novel'].sum()):,}")
t1_conf = int(t[(t["tier"] == 1) & (t["tier1_flag_conflict"] == 1)]["candidate_novel"].sum())
log(f"  Tier-1 candidates with flag-conflict: {t1_conf}")

# breakdown by novelty class within final novel
fn = t[t["final_novel"]]
log("\n=== final_novel by novelty_class ===")
print(fn["novelty_class"].value_counts().to_string())

# ============================================================
# 6. Save
# ============================================================
out = PIPE_V2 / "06b_phylotypes_tiered.tsv"
t.to_csv(out, sep="\t", index=False)
log(f"\nSaved: {out} ({os.path.getsize(out):,} bytes)")

summary = {
    "total_phylotypes": int(len(t)),
    "tier_distribution_all": {int(k): int(n) for k, n in t["tier"].value_counts().items()},
    "candidate_novel": int(len(cn)),
    "tier_distribution_candidate": {int(k): int(n) for k, n in cn["tier"].value_counts().items()},
    "final_novel": int(t["final_novel"].sum()),
    "final_novel_by_novelty": {k: int(v) for k, v in fn["novelty_class"].value_counts().items()},
    "tier1_flag_conflict_candidates": int(t1_conf),
    "residence_coverage": {
        "blank_bearing_bioprojects": 20,
        "phylotypes_tested": int((t["n_members_tested"] > 0).sum()),
        "phylotypes_residence_pass": int((t["n_members_residence_pass"] > 0).sum()),
    },
    "gates": {
        "tier2_bioprojects_min": 2,
        "tier2_ge3_bp": 3,
        "tier2_subjects_min": 50,
        "tier2_continents_min": 2,
        "tier2_requires_2of3": True,
    },
}
with open(PIPE_V2 / "06b_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
log(f"Saved: {PIPE_V2 / '06b_summary.json'}")

log("\nDone!")
