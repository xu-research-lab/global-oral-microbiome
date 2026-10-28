#!/usr/bin/env python3
import argparse
import json
import numpy as np
from biom import load_table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--biom", required=True)
    ap.add_argument(
        "--thresholds",
        default="1,100,500,1000,5000,10000,20000",
        help="Comma-separated sample count thresholds."
    )
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    table = load_table(args.biom)
    sample_sums = np.asarray(table.sum(axis="sample"), dtype=float)
    total_reads = float(sample_sums.sum())
    thresholds = [float(x) for x in args.thresholds.split(",") if x.strip()]

    rows = []
    for t in thresholds:
        keep = sample_sums >= t
        rows.append({
            "min_sample_count": t,
            "n_samples_keep": int(keep.sum()),
            "n_samples_drop": int((~keep).sum()),
            "frac_samples_keep": float(keep.mean()),
            "reads_keep": float(sample_sums[keep].sum()),
            "frac_reads_keep": float(sample_sums[keep].sum() / total_reads) if total_reads > 0 else 0.0,
        })

    with open(args.out, "w") as w:
        w.write(
            "min_sample_count\tn_samples_keep\tn_samples_drop\tfrac_samples_keep\treads_keep\tfrac_reads_keep\n"
        )
        for r in rows:
            w.write(
                f"{r['min_sample_count']}\t"
                f"{r['n_samples_keep']}\t"
                f"{r['n_samples_drop']}\t"
                f"{r['frac_samples_keep']}\t"
                f"{r['reads_keep']}\t"
                f"{r['frac_reads_keep']}\n"
            )

    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
