#!/usr/bin/env python3
"""Phase 4: novelty verdict for candidate phylotypes against multiple references.

The database used for placement (GG2) should not be the sole arbiter of novelty. This
checks each phylotype's representative sequence against HOMD, GG2, and SILVA/NCBI
independently and classifies:
  in_homd            best HOMD identity >= --id-known  (should be rare post-filter; audit)
  known_not_in_homd  no HOMD match but a GG2/SILVA match >= --id-known
                     -> candidate ADDITION to the oral catalog (known bug, new to oral)
  novel_all_refs     no reference >= --id-known anywhere
                     -> candidate new species/lineage

Two modes:
  (a) --hits <tsv>  : classify from a precomputed best-identity table
        columns: phylotype_id, homd_best_id, gg2_best_id, silva_best_id  (percent)
  (b) --rep-fasta   : run vsearch against each DB to build the hits table, then classify
                      (cluster mode; needs vsearch + DB paths)
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402


def vsearch_best_ids(vsearch, fasta, db, threads):
    """Return {query_id: best percent identity} via vsearch --usearch_global."""
    if not db or not os.path.exists(db):
        print(f"[novelty] reference not found, skipping: {db}")
        return {}
    with tempfile.NamedTemporaryFile("w+", suffix=".tsv", delete=False) as tf:
        out = tf.name
    cmd = [vsearch, "--usearch_global", fasta, "--db", db, "--id", "0.5",
           "--maxaccepts", "1", "--maxrejects", "0", "--strand", "both",
           "--userfields", "query+id", "--userout", out, "--threads", str(threads), "--quiet"]
    subprocess.run(cmd, check=True)
    best = {}
    for line in open(out):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            q, pid = parts[0], float(parts[1])
            best[q] = max(best.get(q, 0.0), pid)
    os.unlink(out)
    return best


def classify(homd, gg2, silva, id_known):
    best_other = max([x for x in (gg2, silva) if pd.notna(x)] or [0.0])
    if pd.notna(homd) and homd >= id_known:
        return "in_homd"
    if best_other >= id_known:
        return "known_not_in_homd"
    return "novel_all_refs"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phylotypes", default=os.path.join(C.OUT_ROOT, "06_phylotypes", "phylotypes_candidate_novel.tsv"))
    ap.add_argument("--hits", default=None, help="precomputed phylotype_id,homd_best_id,gg2_best_id,silva_best_id TSV")
    ap.add_argument("--rep-fasta", default=None, help="representative sequences (cluster mode)")
    ap.add_argument("--homd", default=C.HOMD_FASTA)
    ap.add_argument("--gg2", default=C.GG2_FASTA)
    ap.add_argument("--silva", default=C.SILVA_FASTA)
    ap.add_argument("--vsearch", default=C.VSEARCH)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--id-known", type=float, default=97.0)
    ap.add_argument("--outdir", default=os.path.join(C.OUT_ROOT, "07_novelty"))
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ph = pd.read_csv(args.phylotypes, sep="\t", low_memory=False)

    if args.hits:
        hits = pd.read_csv(args.hits, sep="\t")
    elif args.rep_fasta:
        homd = vsearch_best_ids(args.vsearch, args.rep_fasta, args.homd, args.threads)
        gg2 = vsearch_best_ids(args.vsearch, args.rep_fasta, args.gg2, args.threads)
        silva = vsearch_best_ids(args.vsearch, args.rep_fasta, args.silva, args.threads)
        # rep-fasta headers are rep_derep_id; map back to phylotype_id
        rep2ph = dict(zip(ph["rep_derep_id"], ph["phylotype_id"]))
        recs = []
        for rid, pid in rep2ph.items():
            recs.append(dict(phylotype_id=pid, homd_best_id=homd.get(rid),
                             gg2_best_id=gg2.get(rid), silva_best_id=silva.get(rid)))
        hits = pd.DataFrame(recs)
    else:
        raise SystemExit("provide --hits or --rep-fasta")

    m = ph.merge(hits, on="phylotype_id", how="left")
    m["novelty_class"] = [
        classify(r.get("homd_best_id"), r.get("gg2_best_id"), r.get("silva_best_id"), args.id_known)
        for _, r in m.iterrows()
    ]
    out = os.path.join(args.outdir, "phylotypes_novelty_classified.tsv")
    m.to_csv(out, sep="\t", index=False)
    print(m["novelty_class"].value_counts().to_string())
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
