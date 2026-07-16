#!/usr/bin/env python3
"""Phase 1: reference-based chimera screen for HOMD-uncaptured candidates.

DADA2 removes only within-run chimeras; cross-study / cross-taxon chimeras survive and
can masquerade as novel taxa (a chimera partial-matches HOMD, so under GLOBAL identity it
also looks HOMD-uncaptured). This runs `vsearch --uchime_ref` against a combined HOMD + GG2
(+ SILVA) reference and splits candidates into non-chimeric / chimeric / borderline. Must
run before any candidate is promoted to a novel phylotype.

Cluster step (invokes vsearch). Use --dry-run to print the commands without executing.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402


def log(m):
    print(f"[chimera] {m}", flush=True)


def build_combined_db(refs, out_path, dry):
    log(f"combining {len(refs)} reference FASTA(s) -> {out_path}")
    if dry:
        return
    with open(out_path, "w") as out:
        for r in refs:
            if not r or not os.path.exists(r):
                log(f"WARNING: reference not found, skipping: {r}")
                continue
            with open(r) as fh:
                for line in fh:
                    out.write(line)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", default=os.path.join(C.OUT_ROOT, "01_candidates", "candidates.fasta"))
    ap.add_argument("--outdir", default=os.path.join(C.OUT_ROOT, "02_chimera"))
    ap.add_argument("--vsearch", default=C.VSEARCH)
    ap.add_argument("--homd", default=C.HOMD_FASTA)
    ap.add_argument("--gg2", default=C.GG2_FASTA)
    ap.add_argument("--silva", default=C.SILVA_FASTA)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    combined = os.path.join(args.outdir, "combined_reference.fasta")
    build_combined_db([args.homd, args.gg2, args.silva], combined, args.dry_run)

    nonchim = os.path.join(args.outdir, "candidates_nonchimeric.fasta")
    chim = os.path.join(args.outdir, "candidates_chimeric.fasta")
    border = os.path.join(args.outdir, "candidates_borderline.fasta")
    tab = os.path.join(args.outdir, "uchime_ref.tsv")
    cmd = [
        args.vsearch, "--uchime_ref", args.candidates, "--db", combined,
        "--nonchimeras", nonchim, "--chimeras", chim, "--borderline", border,
        "--uchimeout", tab, "--threads", str(args.threads),
    ]
    log("command: " + " ".join(cmd))
    if args.dry_run:
        log("dry-run: not executing")
        return
    subprocess.run(cmd, check=True)

    def count_headers(p):
        return sum(1 for line in open(p) if line.startswith(">")) if os.path.exists(p) else 0
    log(f"non-chimeric: {count_headers(nonchim)}  chimeric: {count_headers(chim)}  "
        f"borderline: {count_headers(border)}")
    log(f"non-chimeric candidates written to {nonchim}")


if __name__ == "__main__":
    main()
