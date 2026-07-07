#!/usr/bin/env python3
"""Post-process wc -l output into sample-level raw read counts."""
import os
import sys
from collections import defaultdict

LINES_FILE = "/beegfs/huangxiaochang/global_oral/new_analysis/file_lines_raw.txt"
OUTFILE = "/beegfs/huangxiaochang/global_oral/new_analysis/raw_reads_stats.tsv"

sample_data = defaultdict(lambda: {"has_paired": False, "lines_paired": 0, "lines_single": 0})

with open(LINES_FILE) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.endswith(" total"):
            continue
        try:
            count_str, filepath = line.rsplit(" ", 1)
            count = int(count_str.strip())
        except ValueError:
            continue

        parts = filepath.split("/")
        project = ""
        for i, p in enumerate(parts):
            if p == "rawdata" and i > 0:
                project = parts[i - 1]
                break
        if not project:
            continue

        fname = os.path.basename(filepath)
        sample = fname
        is_paired = False
        for suffix in ("_1.fastq", "_2.fastq", ".fastq"):
            if sample.endswith(suffix):
                sample = sample[:-len(suffix)]
                if suffix == "_1.fastq":
                    is_paired = True
                break

        key = (project, sample)
        if is_paired:
            sample_data[key]["has_paired"] = True
            sample_data[key]["lines_paired"] += count
        else:
            # single-end file: only use if no paired data exists for this sample
            sample_data[key]["lines_single"] += count

with open(OUTFILE, "w") as out:
    out.write("Project\tSample\tTotal_Reads\n")
    for (project, sample), info in sorted(sample_data.items()):
        if info["has_paired"]:
            reads = info["lines_paired"] // 4 * 2
        else:
            reads = info["lines_single"] // 4
        out.write(f"{project}\t{sample}\t{reads}\n")

print(f"Done: {OUTFILE} with {len(sample_data)} samples")
