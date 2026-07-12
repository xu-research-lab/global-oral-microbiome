#!/usr/bin/env python3
"""Phase 4: define candidate phylotypes from SEPP placements and apply gates.

The reporting unit is the phylotype (candidate taxon), not the raw ASV: many error /
hypervariable / cross-region ASV variants collapse to one taxon, and placement lets ASVs
from different amplicon regions contribute to the same phylotype. A candidate novel oral
phylotype must sit in a HOMD-free neighborhood on the GG2 tree, clear the support floor,
include a chimera-passing member, and include a member that passed the oral-vs-blank
residence test.

Inputs
------
--candidates   candidates_all.tsv from step 01 (derep_id, support_tier, subject_prevalence,
               total_reads, n_bioprojects, sequencing_regions)
--placements   TSV: query_id, nearest_ref, nearest_ref_is_homd (bool), pendant_length,
               and optionally 'cluster' (e.g. from `gappa examine cluster`). Phylotypes are
               grouped by 'cluster' if present, else by nearest_ref.
--nonchimeric  optional: text file of derep_ids that passed 02_chimera_screen (one per line)
--residence    optional: 04 output (feature_id, residence_pass) to gate/annotate
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402


def as_bool(s):
    return s.astype(str).str.lower().isin({"true", "1", "yes", "t", "y"})


def split_regions(x):
    if pd.isna(x):
        return []
    return [v.strip() for v in str(x).split(";") if v.strip() and v.strip() not in {"nan", "NA"}]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", default=os.path.join(C.OUT_ROOT, "01_candidates", "candidates_all.tsv"))
    ap.add_argument("--placements", required=True)
    ap.add_argument("--nonchimeric", default=None)
    ap.add_argument("--residence", default=None)
    ap.add_argument("--outdir", default=os.path.join(C.OUT_ROOT, "06_phylotypes"))
    ap.add_argument("--min-subjects", type=int, default=C.TIER_MED_SUBJECTS)
    ap.add_argument("--min-studies", type=int, default=C.MIN_BIOPROJECTS)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    cand = pd.read_csv(args.candidates, sep="\t", low_memory=False)
    plc = pd.read_csv(args.placements, sep="\t", low_memory=False)
    plc["nearest_ref_is_homd"] = as_bool(plc["nearest_ref_is_homd"])

    nonchim = None
    if args.nonchimeric and os.path.exists(args.nonchimeric):
        nonchim = {l.strip() for l in open(args.nonchimeric) if l.strip()}
    residence = None
    if args.residence and os.path.exists(args.residence):
        rdf = pd.read_csv(args.residence, sep="\t")
        residence = set(rdf.loc[as_bool(rdf["residence_pass"]), "feature_id"])

    df = cand.merge(plc, left_on="derep_id", right_on="query_id", how="inner")
    df["_regions"] = df.get("sequencing_regions", "").map(split_regions)
    df["_nonchim"] = True if nonchim is None else df["derep_id"].isin(nonchim)
    df["_residence"] = pd.NA if residence is None else df["derep_id"].isin(residence)
    group_col = "cluster" if "cluster" in df.columns else "nearest_ref"

    rows = []
    for gid, sub in df.groupby(group_col):
        members = sub["derep_id"].tolist()
        subjects = pd.to_numeric(sub["subject_prevalence"], errors="coerce").fillna(0)
        rep = sub.loc[subjects.idxmax(), "derep_id"] if len(subjects) else members[0]
        in_homd_gap = not sub["nearest_ref_is_homd"].any()   # no HOMD reference nearby
        n_resid = int(sub["_residence"].fillna(False).sum()) if residence is not None else None
        rows.append(dict(
            phylotype_id=str(gid), rep_derep_id=rep, n_asvs=len(members),
            member_derep_ids=";".join(members),
            max_subject_prevalence=int(subjects.max()),
            sum_total_reads=int(pd.to_numeric(sub["total_reads"], errors="coerce").fillna(0).sum()),
            max_n_bioprojects=int(pd.to_numeric(sub["n_bioprojects"], errors="coerce").fillna(0).max()),
            regions=";".join(sorted({r for rs in sub["_regions"] for r in rs})),
            in_homd_gap=in_homd_gap,
            n_nonchimeric=int(sub["_nonchim"].sum()),
            n_residence_pass=n_resid,
        ))
    ph = pd.DataFrame(rows)

    gate = (ph["in_homd_gap"]
            & (ph["max_subject_prevalence"] >= args.min_subjects)
            & (ph["max_n_bioprojects"] >= args.min_studies)
            & (ph["n_nonchimeric"] >= 1))
    if residence is not None:
        gate &= ph["n_residence_pass"].fillna(0) >= 1
    ph["candidate_novel"] = gate

    ph.to_csv(os.path.join(args.outdir, "phylotypes_all.tsv"), sep="\t", index=False)
    ph[ph["candidate_novel"]].to_csv(
        os.path.join(args.outdir, "phylotypes_candidate_novel.tsv"), sep="\t", index=False)
    print(f"phylotypes: {len(ph)}  candidate-novel (HOMD-gap + support + non-chimeric"
          f"{' + residence' if residence is not None else ''}): {int(gate.sum())}")
    print(f"wrote {args.outdir}")


if __name__ == "__main__":
    main()
