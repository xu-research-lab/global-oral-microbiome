#!/usr/bin/env python3
"""
Build a unified GG2 closed-reference abundance profile across all retained datasets.

Mapping chain:
  original project-local ASV
    -> exact derep sequence ID
    -> longest containment representative
    -> GG2 reference ID from closed-reference mapping

The primary output is a sparse long-format table:
  bioproject, sample_id, sample_key, gg2_target, count

This avoids creating an extremely wide dense matrix by default. A per-feature
summary and per-sample summary are also written for QC and downstream filtering.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Tuple

import pandas as pd


DEFAULT_ROOT = "/beegfs/huangxiaochang/global_oral"
DEFAULT_ASV_ROOT = "/beegfs/dongbiao/oral"
DEFAULT_METADATA = f"{DEFAULT_ROOT}/metadata_v2_merge/metadata-v2-merge-oral.csv"
DEFAULT_ASV_META = f"{DEFAULT_ROOT}/ASV_vs_HOMD_analysis/region_aware_pipeline/asv_abundance_metadata.tsv"
DEFAULT_EXACT_SUMMARY = f"{DEFAULT_ROOT}/ASV_vs_HOMD_analysis/region_aware_pipeline/exact_derep_summary.tsv"
DEFAULT_CONTAINMENT_MAP = f"{DEFAULT_ROOT}/ASV_vs_HOMD_analysis/containment_derep/containment_derep_mapping.tsv"
DEFAULT_GG2_HITS = (
    f"{DEFAULT_ROOT}/ASV_vs_HOMD_analysis/containment_derep/"
    "gg2_closed_reference_id99/containment_longest_vs_gg2_id99_best_hits.tsv"
)
DEFAULT_OUTDIR = f"{DEFAULT_ROOT}/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile"


def log(message: str) -> None:
    print(f"[profile] {message}", flush=True)


def load_retained_bioprojects(metadata_path: Path) -> list[str]:
    md = pd.read_csv(metadata_path, low_memory=False, usecols=["bioproject"])
    bps = sorted(md["bioproject"].dropna().astype(str).str.strip().loc[lambda x: x.ne("")].unique())
    return bps


def load_best_gg2_hits(path: Path) -> Dict[str, dict]:
    """Load one deterministic best GG2 target per containment representative."""
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
    best: Dict[str, dict] = {}
    ties = Counter()
    with path.open() as handle:
        reader = csv.DictReader(handle, fieldnames=cols, delimiter="\t")
        for row in reader:
            q = row["query"]
            row["identity"] = float(row["identity"])
            row["query_cov"] = float(row["query_cov"])
            row["alnlen"] = int(row["alnlen"])
            key = (row["identity"], row["query_cov"], row["alnlen"], row["target"])
            if q not in best:
                row["_key"] = key
                best[q] = row
            else:
                ties[q] += 1
                old = best[q]
                # Highest identity/qcov/alnlen; lexicographically smallest target for stable tie-break.
                old_rank = (old["identity"], old["query_cov"], old["alnlen"], old["target"])
                new_rank = (row["identity"], row["query_cov"], row["alnlen"], row["target"])
                if new_rank[:3] > old_rank[:3] or (new_rank[:3] == old_rank[:3] and row["target"] < old["target"]):
                    row["_key"] = key
                    best[q] = row
    for row in best.values():
        row.pop("_key", None)
    log(f"GG2 mapped containment representatives: {len(best):,}; representatives with extra hit rows: {len(ties):,}")
    return best


def load_containment_to_gg2(containment_map: Path, gg2_hits: Dict[str, dict]) -> Dict[str, dict]:
    """Return exact_derep_id -> GG2 assignment for mapped representatives only."""
    out: Dict[str, dict] = {}
    with containment_map.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rep = row["representative_derep_id"]
            hit = gg2_hits.get(rep)
            if hit is None:
                continue
            member = row["member_derep_id"]
            out[member] = {
                "containment_representative": rep,
                "gg2_target": hit["target"],
                "gg2_identity": hit["identity"],
                "gg2_query_cov": hit["query_cov"],
                "gg2_alnlen": hit["alnlen"],
                "assignment_source": "direct_representative" if member == rep else "inherited_containment",
            }
    log(f"exact derep sequences with GG2 assignment after containment propagation: {len(out):,}")
    return out


def load_sequence_to_derep(exact_summary: Path) -> Dict[str, str]:
    seq_to_derep: Dict[str, str] = {}
    with exact_summary.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            seq_to_derep[row["sequence"]] = row["derep_id"]
    log(f"sequence -> exact derep map loaded: {len(seq_to_derep):,}")
    return seq_to_derep


def build_original_asv_assignment(
    asv_meta: Path,
    seq_to_derep: Dict[str, str],
    derep_to_gg2: Dict[str, dict],
    out_path: Path,
    chunksize: int,
) -> tuple[Dict[str, Dict[str, dict]], dict]:
    """Build BioProject -> original ASV ID -> assignment and write audit table."""
    bp_to_asv: Dict[str, Dict[str, dict]] = defaultdict(dict)
    stats = Counter()
    fields = [
        "bioproject",
        "original_asv_id",
        "renamed_asv_id",
        "exact_derep_id",
        "containment_representative",
        "gg2_target",
        "gg2_identity",
        "gg2_query_cov",
        "gg2_alnlen",
        "assignment_source",
    ]

    usecols = ["renamed_asv_id", "bioproject", "original_asv_id", "sequence"]
    with out_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for i, chunk in enumerate(pd.read_csv(asv_meta, sep="\t", usecols=usecols, chunksize=chunksize), start=1):
            for row in chunk.itertuples(index=False):
                stats["original_asvs_total"] += 1
                derep = seq_to_derep.get(row.sequence)
                if derep is None:
                    stats["missing_exact_derep"] += 1
                    continue
                hit = derep_to_gg2.get(derep)
                if hit is None:
                    stats["unmapped_to_gg2"] += 1
                    continue
                record = {
                    "bioproject": row.bioproject,
                    "original_asv_id": row.original_asv_id,
                    "renamed_asv_id": row.renamed_asv_id,
                    "exact_derep_id": derep,
                    **hit,
                }
                writer.writerow(record)
                bp_to_asv[str(row.bioproject)][str(row.original_asv_id)] = record
                stats["original_asvs_with_gg2_assignment"] += 1
            if i % 10 == 0:
                log(f"processed ASV metadata chunks={i}; assigned ASVs={stats['original_asvs_with_gg2_assignment']:,}")
    log(f"original ASV assignment complete: {dict(stats)}")
    return bp_to_asv, dict(stats)


def aggregate_project_counts(
    bp: str,
    counts_path: Path,
    asv_to_assignment: Dict[str, dict],
    profile_handle,
    sample_summary_rows: list[dict],
    feature_counter: Counter,
    feature_sample_counter: Counter,
) -> dict:
    if not counts_path.exists():
        return {"bioproject": bp, "status": "missing_counts", "input_asvs": 0, "mapped_asvs": 0}

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts.index = counts.index.astype(str)
    mapped_asvs = [asv for asv in counts.index if asv in asv_to_assignment]
    if not mapped_asvs:
        return {
            "bioproject": bp,
            "status": "no_mapped_asvs",
            "input_asvs": int(counts.shape[0]),
            "mapped_asvs": 0,
            "samples": int(counts.shape[1]),
        }

    sub = counts.loc[mapped_asvs].apply(pd.to_numeric, errors="coerce").fillna(0)
    targets = pd.Series({asv: asv_to_assignment[asv]["gg2_target"] for asv in mapped_asvs}, name="gg2_target")
    grouped = sub.groupby(targets, sort=False).sum()

    total_reads = float(counts.apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy().sum())
    mapped_reads = float(grouped.to_numpy().sum())

    for sample_id in grouped.columns:
        col = grouped[sample_id]
        nonzero = col[col > 0]
        sample_key = f"{bp}|{sample_id}"
        sample_summary_rows.append(
            {
                "bioproject": bp,
                "sample_id": sample_id,
                "sample_key": sample_key,
                "closed_reference_reads": f"{float(nonzero.sum()):.12g}",
                "closed_reference_features": int(nonzero.shape[0]),
            }
        )
        for gg2_target, value in nonzero.items():
            count = float(value)
            profile_handle.write(f"{bp}\t{sample_id}\t{sample_key}\t{gg2_target}\t{count:.12g}\n")
            feature_counter[gg2_target] += count
            feature_sample_counter[gg2_target] += 1

    return {
        "bioproject": bp,
        "status": "ok",
        "input_asvs": int(counts.shape[0]),
        "mapped_asvs": int(len(mapped_asvs)),
        "samples": int(counts.shape[1]),
        "gg2_features": int(grouped.shape[0]),
        "total_reads": f"{total_reads:.12g}",
        "closed_reference_reads": f"{mapped_reads:.12g}",
        "closed_reference_read_fraction": f"{(mapped_reads / total_reads) if total_reads else 0:.8g}",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default=DEFAULT_METADATA)
    parser.add_argument("--asv-root", default=DEFAULT_ASV_ROOT)
    parser.add_argument("--asv-meta", default=DEFAULT_ASV_META)
    parser.add_argument("--exact-summary", default=DEFAULT_EXACT_SUMMARY)
    parser.add_argument("--containment-map", default=DEFAULT_CONTAINMENT_MAP)
    parser.add_argument("--gg2-hits", default=DEFAULT_GG2_HITS)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--chunksize", type=int, default=200000)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    assignment_path = outdir / "original_asv_to_gg2_closed_reference_assignment.tsv"
    sparse_profile_path = outdir / "all_samples_gg2_closed_reference_profile.long.tsv"
    sample_summary_path = outdir / "all_samples_gg2_closed_reference_sample_summary.tsv"
    feature_summary_path = outdir / "all_samples_gg2_closed_reference_feature_summary.tsv"
    project_summary_path = outdir / "all_samples_gg2_closed_reference_project_summary.tsv"
    stats_path = outdir / "all_samples_gg2_closed_reference_profile.stats.json"

    bps = load_retained_bioprojects(Path(args.metadata))
    log(f"retained BioProjects: {len(bps):,}")
    gg2_hits = load_best_gg2_hits(Path(args.gg2_hits))
    derep_to_gg2 = load_containment_to_gg2(Path(args.containment_map), gg2_hits)
    seq_to_derep = load_sequence_to_derep(Path(args.exact_summary))
    bp_to_asv, assignment_stats = build_original_asv_assignment(
        Path(args.asv_meta),
        seq_to_derep,
        derep_to_gg2,
        assignment_path,
        args.chunksize,
    )

    sample_summary_rows: list[dict] = []
    project_rows = []
    feature_counter: Counter = Counter()
    feature_sample_counter: Counter = Counter()

    with sparse_profile_path.open("w") as profile:
        profile.write("bioproject\tsample_id\tsample_key\tgg2_target\tcount\n")
        for i, bp in enumerate(bps, start=1):
            counts_path = Path(args.asv_root) / bp / "results" / "ASVs_counts.tsv"
            row = aggregate_project_counts(
                bp,
                counts_path,
                bp_to_asv.get(bp, {}),
                profile,
                sample_summary_rows,
                feature_counter,
                feature_sample_counter,
            )
            project_rows.append(row)
            if i % 10 == 0:
                log(f"aggregated {i}/{len(bps)} BioProjects")

    pd.DataFrame(sample_summary_rows).to_csv(sample_summary_path, sep="\t", index=False)
    pd.DataFrame(project_rows).to_csv(project_summary_path, sep="\t", index=False)

    feature_rows = [
        {
            "gg2_target": target,
            "closed_reference_reads": f"{reads:.12g}",
            "sample_prevalence": int(feature_sample_counter[target]),
        }
        for target, reads in feature_counter.items()
    ]
    feature_rows.sort(key=lambda r: (-float(r["closed_reference_reads"]), r["gg2_target"]))
    pd.DataFrame(feature_rows).to_csv(feature_summary_path, sep="\t", index=False)

    stats = {
        "assignment_stats": assignment_stats,
        "retained_bioprojects": len(bps),
        "projects_with_profile": sum(1 for r in project_rows if r.get("status") == "ok"),
        "samples_with_closed_reference_reads": len(sample_summary_rows),
        "gg2_features_observed": len(feature_rows),
        "outputs": {
            "assignment": str(assignment_path),
            "sparse_profile": str(sparse_profile_path),
            "sample_summary": str(sample_summary_path),
            "feature_summary": str(feature_summary_path),
            "project_summary": str(project_summary_path),
        },
    }
    with stats_path.open("w") as handle:
        json.dump(stats, handle, indent=2)
    log(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
