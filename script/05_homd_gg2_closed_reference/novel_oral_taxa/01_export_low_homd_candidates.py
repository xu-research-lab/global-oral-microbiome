#!/usr/bin/env python3
"""Phase 1: export low-HOMD candidate ASVs, keeping the rare tail.

Reads the dereplicated ASV x HOMD table, applies the low-HOMD definition explicitly
(no hit OR HOMD identity < HOMD_ID_CAPTURE OR query coverage < QCOV_MIN, identity
query-global), applies a RELAXED support floor (>= MIN_BIOPROJECTS, >= MIN_TOTAL_READS)
so rare / population-specific taxa are retained, and LABELS each candidate with a support
tier instead of deleting the low-support tail. Emits a per-region candidate table and a
FASTA of representative sequences for the chimera screen (step 02) and placement (step 05).

Expected input columns (superset; names overridable):
  derep_id, sequence, has_homd_hit, identity (HOMD % global), query_cov (%),
  n_bioprojects, total_reads, subject_prevalence, sequencing_regions (';'-joined)
Optional passthrough if present: continents, sites, gg2_genus, gg2_family, gg2_phylum.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402


def split_regions(x):
    if pd.isna(x):
        return []
    return [v.strip() for v in str(x).split(";") if v.strip() and v.strip() not in {"nan", "NA"}]


def support_tier(subjects: float) -> str:
    if subjects >= C.TIER_HIGH_SUBJECTS:
        return "high"
    if subjects >= C.TIER_MED_SUBJECTS:
        return "medium"
    return "low_rare_tail"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", default=C.LOW_HOMD_ANNOTATED, help="dereplicated ASV x HOMD table (TSV)")
    ap.add_argument("--outdir", default=os.path.join(C.OUT_ROOT, "01_candidates"))
    ap.add_argument("--homd-id-max", type=float, default=C.HOMD_ID_CAPTURE)
    ap.add_argument("--qcov-min", type=float, default=C.QCOV_MIN)
    ap.add_argument("--min-bioprojects", type=int, default=C.MIN_BIOPROJECTS)
    ap.add_argument("--min-reads", type=int, default=C.MIN_TOTAL_READS)
    ap.add_argument("--id-col", default="identity")
    ap.add_argument("--qcov-col", default="query_cov")
    ap.add_argument("--hit-col", default="has_homd_hit")
    ap.add_argument("--regions-col", default="sequencing_regions")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.table, sep="\t", low_memory=False)
    n0 = len(df)

    # --- low-HOMD definition (explicit; identity is query-global) --------------
    ident = pd.to_numeric(df.get(args.id_col), errors="coerce")
    qcov = pd.to_numeric(df.get(args.qcov_col), errors="coerce")
    has_hit = df.get(args.hit_col)
    if has_hit is not None:
        has_hit = has_hit.astype(str).str.lower().isin({"true", "1", "yes", "t"})
    else:
        has_hit = ident.notna()
    is_low_homd = (~has_hit) | (ident < args.homd_id_max) | (qcov < args.qcov_min)
    df = df[is_low_homd].copy()
    n_low = len(df)

    # --- relaxed support floor (keep the rare tail) ---------------------------
    nbp = pd.to_numeric(df.get("n_bioprojects"), errors="coerce").fillna(0)
    reads = pd.to_numeric(df.get("total_reads"), errors="coerce").fillna(0)
    df = df[(nbp >= args.min_bioprojects) & (reads >= args.min_reads)].copy()
    n_supported = len(df)

    subjects = pd.to_numeric(df.get("subject_prevalence"), errors="coerce").fillna(0)
    df["support_tier"] = [support_tier(s) for s in subjects]

    # --- write candidate table + per-region tables + FASTA --------------------
    df.to_csv(os.path.join(args.outdir, "candidates_all.tsv"), sep="\t", index=False)

    per_region_counts = {}
    df["_regions"] = df[args.regions_col].map(split_regions) if args.regions_col in df.columns else [[]] * len(df)
    all_regions = sorted({r for rs in df["_regions"] for r in rs})
    for reg in all_regions:
        sub = df[df["_regions"].map(lambda rs, reg=reg: reg in rs)]
        safe = reg.replace("/", "_").replace(" ", "")
        sub.drop(columns=["_regions"]).to_csv(
            os.path.join(args.outdir, f"candidates_region_{safe}.tsv"), sep="\t", index=False)
        per_region_counts[reg] = len(sub)

    if "sequence" in df.columns:
        with open(os.path.join(args.outdir, "candidates.fasta"), "w") as fh:
            for row in df.itertuples(index=False):
                seq = getattr(row, "sequence", "")
                if not isinstance(seq, str) or not seq:
                    continue
                fh.write(f">{row.derep_id}\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")

    tier_counts = df["support_tier"].value_counts().to_dict()
    summ = pd.DataFrame(
        [{"metric": "input_rows", "value": n0},
         {"metric": "low_homd", "value": n_low},
         {"metric": "after_relaxed_support_floor", "value": n_supported},
         *[{"metric": f"tier_{k}", "value": v} for k, v in tier_counts.items()],
         *[{"metric": f"region_{r}", "value": n} for r, n in per_region_counts.items()]])
    summ.to_csv(os.path.join(args.outdir, "candidate_summary.tsv"), sep="\t", index=False)

    print(f"input rows: {n0}")
    print(f"low-HOMD: {n_low}")
    print(f"after relaxed support floor (>= {args.min_bioprojects} BP, >= {args.min_reads} reads): {n_supported}")
    print(f"support tiers: {tier_counts}")
    print(f"regions: {per_region_counts}")
    print(f"wrote candidates to {args.outdir}")


if __name__ == "__main__":
    main()
