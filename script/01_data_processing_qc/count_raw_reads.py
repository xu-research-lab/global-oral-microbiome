#!/usr/bin/env python3
"""
Efficiently count raw reads in all FASTQ files under /beegfs/dongbiao/oral/*/rawdata/.
Uses parallel chunked wc -l to minimize process overhead.
"""

import os
import sys
import glob
import subprocess
import multiprocessing
from collections import defaultdict

BASE_DIR = "/beegfs/dongbiao/oral"
OUTDIR = "/beegfs/huangxiaochang/global_oral/new_analysis"
OUTFILE = os.path.join(OUTDIR, "raw_reads_stats.tsv")
LOGFILE = os.path.join(OUTDIR, "count_raw_reads.log")


def get_sample_info(filepath):
    """Extract project and sample base name from filepath."""
    parts = filepath.split('/')
    project = ""
    for i, p in enumerate(parts):
        if p == "rawdata" and i > 0:
            project = parts[i - 1]
            break
    filename = os.path.basename(filepath)
    sample = filename
    for suffix in ('_1.fastq', '_2.fastq', '.fastq'):
        if sample.endswith(suffix):
            sample = sample[:-len(suffix)]
            break
    return project, sample


def count_lines_chunk(filepaths):
    """Run wc -l on a chunk of files and return list of (filepath, lines)."""
    try:
        result = subprocess.run(
            ['wc', '-l'] + filepaths,
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        # Fallback: count individually if batch fails (e.g., too many args)
        return [(fp, -1) for fp in filepaths]

    results = []
    for line in result.stdout.strip().split('\n'):
        line = line.strip()
        if not line or line.endswith(' total'):
            continue
        # wc -l output: "  12345 /path/to/file"
        # Split from the right once to handle paths with spaces
        try:
            count_str, filepath = line.rsplit(' ', 1)
            count = int(count_str.strip())
            results.append((filepath, count))
        except ValueError:
            continue
    return results


def log(msg):
    print(msg, flush=True)
    with open(LOGFILE, 'a') as f:
        f.write(msg + '\n')


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    with open(LOGFILE, 'w'):
        pass

    log("Scanning for FASTQ files...")
    files = []
    for rawdata in glob.glob(os.path.join(BASE_DIR, "*/rawdata/")):
        for f in glob.glob(os.path.join(rawdata, "*.fastq")):
            files.append(f)

    total_files = len(files)
    log(f"Found {total_files} FASTQ files")
    if total_files == 0:
        sys.exit("No FASTQ files found.")

    # Chunk size: balance between minimizing wc invocations and avoiding ARG_MAX
    chunk_size = 1000
    chunks = [files[i:i + chunk_size] for i in range(0, total_files, chunk_size)]
    log(f"Split into {len(chunks)} chunks of ~{chunk_size} files each")

    sample_lines = defaultdict(int)
    processed = 0
    errors = 0

    n_workers = min(multiprocessing.cpu_count(), 28)
    log(f"Using {n_workers} parallel workers")

    with multiprocessing.Pool(processes=n_workers) as pool:
        for chunk_idx, chunk_results in enumerate(
            pool.imap_unordered(count_lines_chunk, chunks), start=1
        ):
            for filepath, count in chunk_results:
                if count < 0:
                    errors += 1
                    continue
                project, sample = get_sample_info(filepath)
                if not project:
                    errors += 1
                    continue
                sample_lines[(project, sample)] += count

            processed += len(chunk_results)
            if chunk_idx % 10 == 0 or chunk_idx == len(chunks):
                log(f"Progress: {processed}/{total_files} files processed "
                    f"({100*processed/total_files:.1f}%), "
                    f"{len(sample_lines)} samples found, "
                    f"{errors} errors")

    log(f"Writing results to {OUTFILE}...")
    with open(OUTFILE, 'w') as f:
        f.write("Project\tSample\tTotal_Reads\n")
        for (project, sample), lines in sorted(sample_lines.items()):
            reads = lines // 4
            f.write(f"{project}\t{sample}\t{reads}\n")

    log(f"Done. Total samples: {len(sample_lines)}, Total files: {total_files}")


if __name__ == '__main__':
    main()
