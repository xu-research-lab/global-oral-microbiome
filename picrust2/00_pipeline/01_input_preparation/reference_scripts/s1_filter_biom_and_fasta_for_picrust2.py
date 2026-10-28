#!/usr/bin/env python3
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from biom import load_table
from biom.util import biom_open


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fasta(path):
    order = []
    seq_by_id = {}
    current_id = None
    chunks = []

    with open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    seq_by_id[current_id] = "".join(chunks).upper()

                current_id = line[1:].strip().split()[0]
                order.append(current_id)
                chunks = []
            else:
                chunks.append(line.strip())

        if current_id is not None:
            seq_by_id[current_id] = "".join(chunks).upper()

    return order, seq_by_id


def wrap_seq(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser(
        description="Filter BIOM samples/features and matching FASTA for PICRUSt2."
    )
    ap.add_argument("--biom", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", default="picrust2_filtered")

    # min-feature-count 1 means remove zero-count features.
    ap.add_argument(
        "--min-feature-count",
        type=float,
        default=1,
        help="Keep features with total count >= this value. Default=1 removes zero-count features."
    )
    ap.add_argument(
        "--min-sample-count",
        type=float,
        default=1,
        help="Keep samples with total count >= this value. Default=1 removes zero-count samples only."
    )

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    table = load_table(args.biom)

    original_sample_ids = list(table.ids(axis="sample"))
    original_feature_ids = list(table.ids(axis="observation"))

    sample_sums = np.asarray(table.sum(axis="sample"), dtype=float)
    feature_sums = np.asarray(table.sum(axis="observation"), dtype=float)

    keep_samples = [
        sid for sid, s in zip(original_sample_ids, sample_sums)
        if s >= args.min_sample_count
    ]

    table_sample_filtered = table.filter(
        keep_samples,
        axis="sample",
        inplace=False
    )

    # Recalculate feature sums after sample filtering.
    feature_ids_after_sample = list(table_sample_filtered.ids(axis="observation"))
    feature_sums_after_sample = np.asarray(
        table_sample_filtered.sum(axis="observation"),
        dtype=float
    )

    keep_features = [
        fid for fid, s in zip(feature_ids_after_sample, feature_sums_after_sample)
        if s >= args.min_feature_count
    ]

    filtered_table = table_sample_filtered.filter(
        keep_features,
        axis="observation",
        inplace=False
    )

    final_sample_ids = list(filtered_table.ids(axis="sample"))
    final_feature_ids = list(filtered_table.ids(axis="observation"))
    final_feature_set = set(final_feature_ids)

    # Write BIOM
    out_biom = outdir / f"{args.prefix}.feature_table.biom"
    with biom_open(str(out_biom), "w") as f:
        filtered_table.to_hdf5(f, "filtered_for_picrust2")

    # Write FASTA filtered to final BIOM feature IDs.
    _, seq_by_id = read_fasta(args.fasta)

    missing_fasta = sorted(final_feature_set - set(seq_by_id))
    if missing_fasta:
        missing_file = outdir / f"{args.prefix}.missing_fasta_ids.txt"
        with open(missing_file, "w") as w:
            for x in missing_fasta:
                w.write(x + "\n")
        raise SystemExit(
            f"ERROR: {len(missing_fasta)} final BIOM feature IDs missing from FASTA. "
            f"See {missing_file}"
        )

    out_fasta = outdir / f"{args.prefix}.rep_seqs.fna"
    with open(out_fasta, "w") as w:
        for fid in final_feature_ids:
            w.write(f">{fid}\n")
            w.write(wrap_seq(seq_by_id[fid]) + "\n")

    # Write dropped lists
    dropped_samples = sorted(set(original_sample_ids) - set(final_sample_ids))
    dropped_features = sorted(set(original_feature_ids) - set(final_feature_ids))

    dropped_samples_file = outdir / f"{args.prefix}.dropped_samples.txt"
    dropped_features_file = outdir / f"{args.prefix}.dropped_features.txt"

    with open(dropped_samples_file, "w") as w:
        for x in dropped_samples:
            w.write(x + "\n")

    with open(dropped_features_file, "w") as w:
        for x in dropped_features:
            w.write(x + "\n")

    summary = {
        "input_biom": args.biom,
        "input_fasta": args.fasta,
        "output_biom": str(out_biom),
        "output_fasta": str(out_fasta),
        "min_sample_count": args.min_sample_count,
        "min_feature_count": args.min_feature_count,

        "n_samples_before": len(original_sample_ids),
        "n_samples_after": len(final_sample_ids),
        "n_samples_dropped": len(dropped_samples),

        "n_features_before": len(original_feature_ids),
        "n_features_after": len(final_feature_ids),
        "n_features_dropped": len(dropped_features),

        "total_reads_before": float(np.sum(sample_sums)),
        "total_reads_after": float(np.sum(filtered_table.sum(axis="sample"))),
    }

    summary_file = outdir / f"{args.prefix}.summary.json"
    with open(summary_file, "w") as w:
        json.dump(summary, w, indent=2)
        w.write("\n")

    print(json.dumps(summary, indent=2))
    print("")
    print("DONE. Use these files for PICRUSt2:")
    print(f"  -i {out_biom}")
    print(f"  -s {out_fasta}")


if __name__ == "__main__":
    main()
