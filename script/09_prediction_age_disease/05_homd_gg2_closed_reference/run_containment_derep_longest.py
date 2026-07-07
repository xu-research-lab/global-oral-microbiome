#!/usr/bin/env python3
"""
Collapse exact-dereplicated ASV sequences by exact containment.

The intent is computational: if a shorter exact-derep sequence is fully contained
in a longer exact-derep sequence at 100% identity and 100% query coverage, the
longest sequence is used as the representative for downstream closed-reference
mapping. The original exact-derep IDs and source ASV abundance metadata are kept
in mapping tables so assignments can be propagated or audited later.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Tuple


DEFAULT_IN_FASTA = (
    "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/"
    "region_aware_pipeline/exact_derep_representatives.fa"
)
DEFAULT_IN_SUMMARY = (
    "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/"
    "region_aware_pipeline/exact_derep_summary.tsv"
)
DEFAULT_OUTDIR = (
    "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/"
    "containment_derep"
)
DEFAULT_VSEARCH = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"


def log(message: str) -> None:
    print(f"[containment] {message}", flush=True)


def read_fasta_lengths(path: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    header = None
    total = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    lengths[header] = total
                header = line[1:].split()[0]
                total = 0
            else:
                total += len(line)
        if header is not None:
            lengths[header] = total
    return lengths


def iter_fasta(path: Path) -> Iterable[Tuple[str, str]]:
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
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def run_vsearch(args: argparse.Namespace, uc_path: Path, centroid_path: Path) -> None:
    cmd = [
        args.vsearch,
        "--cluster_fast",
        args.input_fasta,
        "--id",
        "1.0",
        "--query_cov",
        "1.0",
        "--iddef",
        "0",
        "--strand",
        args.strand,
        "--centroids",
        str(centroid_path),
        "--uc",
        str(uc_path),
        "--threads",
        str(args.threads),
        "--minseqlength",
        "1",
    ]
    log("running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_uc(uc_path: Path, lengths: Dict[str, int], mapping_path: Path) -> Tuple[Dict[str, str], Dict[str, dict]]:
    member_to_rep: Dict[str, str] = {}
    rep_stats: Dict[str, dict] = {}

    with uc_path.open() as inp, mapping_path.open("w", newline="") as out:
        fields = [
            "member_derep_id",
            "representative_derep_id",
            "member_is_representative",
            "member_length",
            "representative_length",
            "vsearch_identity",
            "strand",
            "uc_record_type",
        ]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for line in inp:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            record_type = parts[0]
            if record_type == "S":
                rep = parts[8]
                member = rep
                identity = "100.0"
                strand = "+"
                is_rep = 1
            elif record_type == "H":
                member = parts[8]
                rep = parts[9]
                identity = parts[3]
                strand = parts[4]
                is_rep = 0
            else:
                continue

            member_to_rep[member] = rep
            if rep not in rep_stats:
                rep_stats[rep] = {
                    "contained_exact_derep_sequences": 0,
                    "contained_shorter_sequences": 0,
                    "representative_length": lengths.get(rep, 0),
                    "min_member_length": lengths.get(member, 0),
                    "max_member_length": lengths.get(member, 0),
                }
            rep_stats[rep]["contained_exact_derep_sequences"] += 1
            if member != rep:
                rep_stats[rep]["contained_shorter_sequences"] += 1
            member_len = lengths.get(member, 0)
            rep_stats[rep]["min_member_length"] = min(rep_stats[rep]["min_member_length"], member_len)
            rep_stats[rep]["max_member_length"] = max(rep_stats[rep]["max_member_length"], member_len)

            writer.writerow(
                {
                    "member_derep_id": member,
                    "representative_derep_id": rep,
                    "member_is_representative": is_rep,
                    "member_length": member_len,
                    "representative_length": lengths.get(rep, 0),
                    "vsearch_identity": identity,
                    "strand": strand,
                    "uc_record_type": record_type,
                }
            )

    return member_to_rep, rep_stats


def load_exact_summary(summary_path: Path) -> Dict[str, dict]:
    summary: Dict[str, dict] = {}
    with summary_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            derep_id = row["derep_id"]
            summary[derep_id] = row
    return summary


def write_representative_fasta(input_fasta: Path, reps: set[str], output_fasta: Path) -> None:
    with output_fasta.open("w") as out:
        for derep_id, sequence in iter_fasta(input_fasta):
            if derep_id in reps:
                write_fasta_record(out, derep_id, sequence)


def safe_float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return 0.0


def safe_int_from_float(value: str) -> int:
    return int(round(safe_float(value)))


def write_summary(
    rep_stats: Dict[str, dict],
    member_to_rep: Dict[str, str],
    exact_summary: Dict[str, dict],
    summary_path: Path,
) -> None:
    rep_to_members: Dict[str, list[str]] = defaultdict(list)
    for member, rep in member_to_rep.items():
        rep_to_members[rep].append(member)

    fields = [
        "representative_derep_id",
        "representative_length",
        "contained_exact_derep_sequences",
        "contained_shorter_sequences",
        "min_member_length",
        "max_member_length",
        "source_asvs",
        "n_bioprojects",
        "total_reads",
        "prevalence_samples",
        "representative_source_asv_id",
        "bioprojects",
        "sequencing_regions",
        "platforms",
        "countries",
        "continents",
        "sites",
        "disease_studies",
        "statuses",
    ]

    meta_set_fields = [
        "bioprojects",
        "sequencing_regions",
        "platforms",
        "countries",
        "continents",
        "sites",
        "disease_studies",
        "statuses",
    ]

    rows = []
    for rep, members in rep_to_members.items():
        source_asvs = 0
        total_reads = 0.0
        prevalence_samples = 0
        merged_sets = {field: set() for field in meta_set_fields}

        for member in members:
            row = exact_summary[member]
            source_asvs += safe_int_from_float(row.get("n_asvs", "0"))
            total_reads += safe_float(row.get("total_reads", "0"))
            prevalence_samples += safe_int_from_float(row.get("prevalence_samples", "0"))
            for field in meta_set_fields:
                for value in row.get(field, "").split(";"):
                    value = value.strip()
                    if value and value not in {"NA", "nan"}:
                        merged_sets[field].add(value)

        rep_exact = exact_summary[rep]
        stat = rep_stats[rep]
        out = {
            "representative_derep_id": rep,
            "representative_length": stat["representative_length"],
            "contained_exact_derep_sequences": stat["contained_exact_derep_sequences"],
            "contained_shorter_sequences": stat["contained_shorter_sequences"],
            "min_member_length": stat["min_member_length"],
            "max_member_length": stat["max_member_length"],
            "source_asvs": source_asvs,
            "n_bioprojects": len(merged_sets["bioprojects"]),
            "total_reads": f"{total_reads:.6g}",
            "prevalence_samples": prevalence_samples,
            "representative_source_asv_id": rep_exact.get("representative_asv_id", ""),
        }
        for field in meta_set_fields:
            out[field] = ";".join(sorted(merged_sets[field])) if merged_sets[field] else "NA"
        rows.append(out)

    rows.sort(
        key=lambda r: (
            -safe_float(str(r["total_reads"])),
            -int(r["contained_exact_derep_sequences"]),
            r["representative_derep_id"],
        )
    )
    with summary_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", default=DEFAULT_IN_FASTA)
    parser.add_argument("--exact-summary", default=DEFAULT_IN_SUMMARY)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--vsearch", default=DEFAULT_VSEARCH)
    parser.add_argument("--threads", type=int, default=28)
    parser.add_argument("--strand", choices=["plus", "both"], default="both")
    parser.add_argument("--skip-vsearch", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_fasta = Path(args.input_fasta)
    exact_summary_path = Path(args.exact_summary)
    uc_path = outdir / "containment_clusters_id100_qcov100.uc"
    centroid_path = outdir / "vsearch_centroids_id100_qcov100.fa"
    mapping_path = outdir / "containment_derep_mapping.tsv"
    representative_fasta = outdir / "containment_longest_representatives.fa"
    summary_path = outdir / "containment_representative_summary.tsv"
    stats_path = outdir / "containment_derep_stats.json"

    log("reading input FASTA lengths")
    lengths = read_fasta_lengths(input_fasta)
    log(f"input exact-derep sequences: {len(lengths):,}")

    if not args.skip_vsearch:
        run_vsearch(args, uc_path, centroid_path)
    elif not uc_path.exists():
        raise FileNotFoundError(f"--skip-vsearch was used but UC file is missing: {uc_path}")

    log("parsing UC mapping")
    member_to_rep, rep_stats = parse_uc(uc_path, lengths, mapping_path)
    log(f"containment representatives: {len(rep_stats):,}")

    log("writing representative FASTA")
    write_representative_fasta(input_fasta, set(rep_stats), representative_fasta)

    log("loading exact-derep summary")
    exact_summary = load_exact_summary(exact_summary_path)
    log("writing representative summary")
    write_summary(rep_stats, member_to_rep, exact_summary, summary_path)

    singleton_reps = sum(1 for stat in rep_stats.values() if stat["contained_exact_derep_sequences"] == 1)
    stats = {
        "input_exact_derep_sequences": len(lengths),
        "mapped_exact_derep_sequences": len(member_to_rep),
        "containment_representatives": len(rep_stats),
        "singleton_representatives": singleton_reps,
        "collapsed_nonrepresentative_sequences": len(member_to_rep) - len(rep_stats),
        "vsearch_identity": 1.0,
        "vsearch_query_coverage": 1.0,
        "vsearch_strand": args.strand,
        "outputs": {
            "mapping": str(mapping_path),
            "representative_fasta": str(representative_fasta),
            "representative_summary": str(summary_path),
            "uc": str(uc_path),
            "vsearch_centroids": str(centroid_path),
        },
    }
    with stats_path.open("w") as handle:
        json.dump(stats, handle, indent=2)
    log(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
