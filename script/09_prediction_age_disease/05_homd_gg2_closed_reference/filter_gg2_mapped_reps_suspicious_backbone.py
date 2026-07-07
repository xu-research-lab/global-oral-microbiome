#!/usr/bin/env python3
"""Filter GG2-mapped containment representatives using GG2 backbone taxonomy.

Rule:
- Re-align already GG2-mapped representative sequences to GG2 backbone.
- For each query, keep all final best hits by highest identity, then highest
  query coverage, then longest alignment.
- If any final best-hit taxonomy indicates mitochondria, chloroplast, or the
  eukaryotic genus g__GWA2-37-10, mark the query suspicious and remove it.
- Queries without a backbone hit are retained, because they already passed the
  earlier closed-reference mapping and this step is only for suspicious filtering.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, Set, Tuple

import pandas as pd


DEFAULT_BASE = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep"
DEFAULT_QUERY = f"{DEFAULT_BASE}/gg2_closed_reference_id99/containment_longest_vs_gg2_id99_matched.fa"
DEFAULT_HITS = f"{DEFAULT_BASE}/gg2_backbone_suspicious_filter/containment_longest_matched_vs_gg2_backbone.top_hits.tsv"
DEFAULT_TAXONOMY = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"
DEFAULT_OUTDIR = f"{DEFAULT_BASE}/gg2_backbone_suspicious_filter"


def log(message: str) -> None:
    print(f"[suspicious-filter] {message}", flush=True)


def read_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    header = None
    chunks = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def write_fasta_record(handle, header: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{header}\n")
    for i in range(0, len(sequence), width):
        handle.write(sequence[i : i + width] + "\n")


def suspicious_taxonomy(taxon: str) -> bool:
    t = str(taxon).lower()
    return (
        "mitochondria" in t
        or "mitochondr" in t
        or "chloroplast" in t
        or "chloro" in t
        or "g__gwa2-37-10" in t
    )


def suspicious_reason(taxon: str) -> str:
    t = str(taxon).lower()
    reasons = []
    if "mitochondria" in t or "mitochondr" in t:
        reasons.append("mitochondria")
    if "chloroplast" in t or "chloro" in t:
        reasons.append("chloroplast")
    if "g__gwa2-37-10" in t:
        reasons.append("g__GWA2-37-10")
    return ";".join(reasons) if reasons else "NA"


def load_taxonomy(path: Path) -> Dict[str, str]:
    tax = pd.read_csv(path, sep="\t", dtype=str)
    tax = tax.rename(columns={"Feature ID": "target", "Taxon": "taxonomy"})
    return dict(zip(tax["target"], tax["taxonomy"].fillna("")))


def final_best_hits(hits_path: Path, taxonomy: Dict[str, str], out_best_hits: Path, out_suspicious_hits: Path):
    cols = [
        "query",
        "target",
        "identity",
        "alnlen",
        "mismatches",
        "gaps",
        "qlo",
        "qhi",
        "tlo",
        "thi",
        "query_length",
        "target_length",
        "query_cov",
        "target_cov",
    ]
    hits = pd.read_csv(hits_path, sep="\t", names=cols)
    if hits.empty:
        return set(), set(), 0, 0
    for col in ["identity", "alnlen", "query_cov"]:
        hits[col] = pd.to_numeric(hits[col], errors="coerce").fillna(0)

    best_rows = []
    suspicious_queries: Set[str] = set()
    suspicious_rows = []

    for query, sub in hits.groupby("query", sort=False):
        max_id = sub["identity"].max()
        sub = sub[sub["identity"] == max_id]
        max_qcov = sub["query_cov"].max()
        sub = sub[sub["query_cov"] == max_qcov]
        max_aln = sub["alnlen"].max()
        sub = sub[sub["alnlen"] == max_aln].copy()
        sub["taxonomy"] = sub["target"].map(lambda x: taxonomy.get(x, "NO_TAXONOMY"))
        sub["is_suspicious_hit"] = sub["taxonomy"].map(suspicious_taxonomy)
        sub["suspicious_reason"] = sub["taxonomy"].map(suspicious_reason)
        if sub["is_suspicious_hit"].any():
            suspicious_queries.add(query)
            suspicious_rows.append(sub[sub["is_suspicious_hit"]])
        best_rows.append(sub)

    best = pd.concat(best_rows, ignore_index=True) if best_rows else pd.DataFrame()
    best.to_csv(out_best_hits, sep="\t", index=False)
    if suspicious_rows:
        pd.concat(suspicious_rows, ignore_index=True).to_csv(out_suspicious_hits, sep="\t", index=False)
    else:
        pd.DataFrame(columns=list(best.columns)).to_csv(out_suspicious_hits, sep="\t", index=False)
    return set(best["query"]), suspicious_queries, len(hits), len(best)


def split_fasta(query_fasta: Path, hit_queries: Set[str], suspicious_queries: Set[str], outdir: Path):
    retained_fa = outdir / "containment_longest_vs_gg2_id99_matched.backbone_filtered_retained.fa"
    suspicious_fa = outdir / "containment_longest_vs_gg2_id99_matched.backbone_suspicious_removed.fa"
    no_backbone_fa = outdir / "containment_longest_vs_gg2_id99_matched.no_backbone_hit_retained.fa"

    retained = 0
    suspicious = 0
    no_hit = 0
    total = 0
    with retained_fa.open("w") as kept, suspicious_fa.open("w") as bad, no_backbone_fa.open("w") as nohit:
        for header, seq in read_fasta(query_fasta):
            query = header.split()[0]
            total += 1
            if query in suspicious_queries:
                write_fasta_record(bad, header, seq)
                suspicious += 1
            else:
                write_fasta_record(kept, header, seq)
                retained += 1
                if query not in hit_queries:
                    write_fasta_record(nohit, header, seq)
                    no_hit += 1
    return {
        "query_fasta_records": total,
        "retained_records": retained,
        "suspicious_removed_records": suspicious,
        "no_backbone_hit_retained_records": no_hit,
        "retained_fasta": str(retained_fa),
        "suspicious_removed_fasta": str(suspicious_fa),
        "no_backbone_hit_retained_fasta": str(no_backbone_fa),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-fasta", default=DEFAULT_QUERY)
    parser.add_argument("--hits", default=DEFAULT_HITS)
    parser.add_argument("--taxonomy", default=DEFAULT_TAXONOMY)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    best_hits_path = outdir / "containment_longest_matched_vs_gg2_backbone.final_best_hits_with_taxonomy.tsv"
    suspicious_hits_path = outdir / "containment_longest_matched_vs_gg2_backbone.suspicious_best_hits.tsv"
    suspicious_ids_path = outdir / "containment_longest_vs_gg2_id99_matched.suspicious_query_ids.txt"
    retained_ids_path = outdir / "containment_longest_vs_gg2_id99_matched.backbone_filtered_retained_ids.txt"
    stats_path = outdir / "containment_longest_vs_gg2_id99_matched.backbone_suspicious_filter.stats.json"

    log("loading taxonomy")
    taxonomy = load_taxonomy(Path(args.taxonomy))
    log(f"taxonomy records: {len(taxonomy):,}")

    log("selecting final best hits and suspicious queries")
    hit_queries, suspicious_queries, candidate_hit_rows, final_best_hit_rows = final_best_hits(
        Path(args.hits), taxonomy, best_hits_path, suspicious_hits_path
    )
    log(f"queries with backbone hits: {len(hit_queries):,}")
    log(f"suspicious queries: {len(suspicious_queries):,}")

    log("splitting FASTA")
    fasta_stats = split_fasta(Path(args.query_fasta), hit_queries, suspicious_queries, outdir)

    with suspicious_ids_path.open("w") as handle:
        for q in sorted(suspicious_queries):
            handle.write(q + "\n")

    all_retained = []
    for header, _ in read_fasta(Path(args.query_fasta)):
        q = header.split()[0]
        if q not in suspicious_queries:
            all_retained.append(q)
    with retained_ids_path.open("w") as handle:
        for q in all_retained:
            handle.write(q + "\n")

    stats = {
        **fasta_stats,
        "queries_with_backbone_hit": len(hit_queries),
        "queries_without_backbone_hit_retained": fasta_stats["no_backbone_hit_retained_records"],
        "candidate_hit_rows_from_vsearch": candidate_hit_rows,
        "final_best_hit_rows": final_best_hit_rows,
        "suspicious_query_ids": str(suspicious_ids_path),
        "retained_query_ids": str(retained_ids_path),
        "final_best_hits_with_taxonomy": str(best_hits_path),
        "suspicious_best_hits": str(suspicious_hits_path),
        "filter_rule": "remove query if any final best backbone hit taxonomy contains mitochondria, chloroplast, or g__GWA2-37-10",
    }
    with stats_path.open("w") as handle:
        json.dump(stats, handle, indent=2)
    log(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
