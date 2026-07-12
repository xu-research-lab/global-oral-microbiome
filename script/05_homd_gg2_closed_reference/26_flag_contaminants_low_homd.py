#!/usr/bin/env python3
"""Screen low-HOMD candidate ASVs against a curated contaminant / oral-core reference.

Fast-path decontamination for the "novel oral taxa not in HOMD" list: applies the
`contaminant_reference_genera.tsv` blacklist/whitelist to the GG2-annotated low-HOMD
table and labels each ASV keep / flag / remove with a reason. Oral-core genera are
protected first (never removed), so genuine oral taxa cannot be lost to a genus that is
also a common contaminant elsewhere.

IMPORTANT — this is the *weak* screen the user opted into (detection-based + blacklist).
It removes reagent/marine/gut/plant genera by name, but it does NOT establish oral
residence. It will NOT catch: gut/marine/plant lineages whose genus is absent from the
reference list, novel contaminants, or GG2-Unassigned junk. For a defensible result add
blank-based decontam (see 25_negative_control_audit.py for which BioProjects have blanks)
and an oral-vs-gut enrichment test. Treat the "kept" set as candidates, not conclusions.

Inputs
------
  --table      GG2-annotated low-HOMD table (cluster intermediate from script 15),
               e.g. low_identity_lt97_subject_prevalence_filtered_gg2_annotated.tsv
               Expected taxonomy columns: gg2_genus, gg2_family, gg2_phylum, gg2_taxonomy
  --reference  contaminant_reference_genera.tsv (defaults to the copy beside this script)

Outputs (under --outdir)
------------------------
  low_homd_candidates_screened.tsv   all ASVs + columns: screen_action, screen_reason,
                                     matched_taxon, matched_category
  low_homd_candidates_kept.tsv       screen_action == keep
  low_homd_candidates_removed.tsv    screen_action == remove
  low_homd_screen_summary.tsv        counts by action and category

Not executed on real data here (table is a cluster intermediate); logic smoke-tested on
synthetic input.
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import Counter
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT_REF = HERE / "contaminant_reference_genera.tsv"


def norm(name) -> str:
    """Normalize a GG2 taxon label to a bare lowercase genus/phylum token.

    'g__Vreelandella_640984' -> 'vreelandella'; 'g__Coprococcus_A_187866' -> 'coprococcus';
    'p__Cyanobacteriota' -> 'cyanobacteriota'.
    """
    if not isinstance(name, str) or not name:
        return ""
    for pref in ("g__", "p__", "f__", "c__", "o__", "s__", "d__"):
        if name.startswith(pref):
            name = name[len(pref):]
            break
    return name.split("_")[0].strip().lower()


def load_reference(path: Path):
    ref = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("taxon\t") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            taxon, rank, classification, category = parts[0], parts[1], parts[2], parts[3]
            ref[taxon.strip().lower()] = (classification, category)
    return ref


def screen_row(row, ref: dict, gcol, fcol, pcol, taxcol):
    genus = norm(row.get(gcol, ""))
    family = norm(row.get(fcol, ""))
    phylum = norm(row.get(pcol, ""))
    taxstr = str(row.get(taxcol, "")).lower()

    # 1. oral-core protection wins over everything
    for lvl in (genus, family, phylum):
        if lvl in ref and ref[lvl][0] == "oral_core":
            return "keep", "oral_core", lvl, ref[lvl][1]
    # 2. explicit plastid/host strings anywhere in the taxonomy
    for kw in ("chloroplast", "mitochondria"):
        if kw in taxstr:
            return "remove", "dietary_plastid", kw, "dietary_plastid"
    # 3. auto_remove / flag_review by genus, then family, then phylum
    for lvl in (genus, family, phylum):
        if lvl in ref:
            cls, cat = ref[lvl]
            if cls == "auto_remove":
                return "remove", cat, lvl, cat
            if cls == "flag_review":
                return "flag", cat, lvl, cat
    # 4. unmatched -> retained as candidate novel (but unproven oral)
    return "keep", "unmatched_candidate", "", ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="GG2-annotated low-HOMD table (TSV)")
    ap.add_argument("--reference", type=Path, default=DEFAULT_REF)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--genus-col", default="gg2_genus")
    ap.add_argument("--family-col", default="gg2_family")
    ap.add_argument("--phylum-col", default="gg2_phylum")
    ap.add_argument("--taxonomy-col", default="gg2_taxonomy")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ref = load_reference(args.reference)
    df = pd.read_csv(args.table, sep="\t", low_memory=False)

    actions, reasons, matched, cats = [], [], [], []
    for _, row in df.iterrows():
        a, r, m, c = screen_row(row, ref, args.genus_col, args.family_col,
                                args.phylum_col, args.taxonomy_col)
        actions.append(a); reasons.append(r); matched.append(m); cats.append(c)
    df["screen_action"] = actions
    df["screen_reason"] = reasons
    df["matched_taxon"] = matched
    df["matched_category"] = cats

    df.to_csv(os.path.join(args.outdir, "low_homd_candidates_screened.tsv"), sep="\t", index=False)
    df[df["screen_action"] == "keep"].to_csv(
        os.path.join(args.outdir, "low_homd_candidates_kept.tsv"), sep="\t", index=False)
    df[df["screen_action"] == "remove"].to_csv(
        os.path.join(args.outdir, "low_homd_candidates_removed.tsv"), sep="\t", index=False)

    summ = (df.groupby(["screen_action", "screen_reason"]).size()
              .reset_index(name="n_asvs").sort_values(["screen_action", "n_asvs"], ascending=[True, False]))
    summ.to_csv(os.path.join(args.outdir, "low_homd_screen_summary.tsv"), sep="\t", index=False)

    n = len(df)
    counts = Counter(actions)
    print(f"screened {n} low-HOMD ASVs")
    for a in ("keep", "flag", "remove"):
        print(f"  {a:7s}: {counts.get(a,0)} ({counts.get(a,0)/n*100:.1f}%)")
    print("\nby reason:")
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
