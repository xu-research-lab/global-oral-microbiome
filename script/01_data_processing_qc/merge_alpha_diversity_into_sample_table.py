#!/usr/bin/env python3
import csv
from pathlib import Path


ROOT = Path("/beegfs/huangxiaochang/global_oral/new_analysis")
ALPHA_PATH = ROOT / "alpha_diversity_data.csv"
SAMPLE_PATH = ROOT / "sample_reads_asv_chim_perc.tsv"


def main() -> None:
    with ALPHA_PATH.open(newline="") as handle:
        reader = csv.DictReader(handle)
        alpha_fieldnames = reader.fieldnames or []
        merge_columns = [c for c in alpha_fieldnames if c not in {"bioproject", "sample_id"}]
        alpha_by_key = {
            (row["bioproject"], row["sample_id"]): {col: row[col] for col in merge_columns}
            for row in reader
        }

    with SAMPLE_PATH.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        sample_fieldnames = reader.fieldnames or []
        merged_fieldnames = [c for c in sample_fieldnames if c not in merge_columns] + merge_columns
        rows = []
        matched = 0
        unmatched = 0
        for row in reader:
            key = (row["bioproject_id"], row["sample_id"])
            alpha_row = alpha_by_key.get(key)
            if alpha_row is None:
                unmatched += 1
                for col in merge_columns:
                    row[col] = "NA"
            else:
                matched += 1
                row.update(alpha_row)
            rows.append(row)

    with SAMPLE_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=merged_fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"matched={matched}")
    print(f"unmatched={unmatched}")
    print(f"added_columns={len(merge_columns)}")


if __name__ == "__main__":
    main()
