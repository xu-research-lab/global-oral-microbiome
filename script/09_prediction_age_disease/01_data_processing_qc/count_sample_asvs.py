#!/usr/bin/env python3
import os
import sys
import h5py
import numpy as np


def bioproject_id(path):
    marker = '/beegfs/dongbiao/oral/'
    if path.startswith(marker) and path.endswith('/results/DADA2.biom'):
        return path[len(marker): -len('/results/DADA2.biom')]
    return os.path.basename(os.path.dirname(os.path.dirname(path)))


def count_file(path):
    study = bioproject_id(path)
    with h5py.File(path, 'r') as f:
        sample_ids = [x.decode() if isinstance(x, bytes) else x
                      for x in f['sample/ids'][:]]
        n_samples = len(sample_ids)

        indices = f['observation/matrix/indices'][:]
        indptr = f['observation/matrix/indptr'][:]

        counts = np.zeros(n_samples, dtype=int)
        for i in range(len(indptr) - 1):
            start, end = indptr[i], indptr[i + 1]
            if start < end:
                cols = indices[start:end]
                for col in cols:
                    if col < n_samples:
                        counts[col] += 1

        samples = sample_ids
    return study, samples, counts


def main():
    if len(sys.argv) != 3:
        raise SystemExit('Usage: count_sample_asvs.py <file_list> <output_tsv>')

    file_list, output_tsv = sys.argv[1:]
    with open(file_list, 'rt', encoding='utf-8') as handle:
        paths = [line.strip() for line in handle if line.strip()]

    with open(output_tsv, 'wt', encoding='utf-8', newline='') as out:
        col_sep = chr(9)
        nl = chr(10)
        out.write(col_sep.join(['bioproject_id', 'sample_id', 'asv_count']) + nl)
        for idx, path in enumerate(paths, start=1):
            study, samples, counts = count_file(path)
            for sample, count in zip(samples, counts):
                out.write(col_sep.join([study, str(sample), str(count)]) + nl)
            out.flush()
            print(f'[{idx}/{len(paths)}] {study}: {len(samples)} samples', file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
