#!/usr/bin/env python
"""
Precompute DeepPhylo inputs on a CPU node (PCA embedding is CPU-bound).

Outputs (written to --output_dir):
  - otu_ids.npy            : fixed OTU order (1-D array of feature ids)
  - phy_emb_global_hs{HS}.npy : phylogenetic embedding, shape (n_otus, hidden_size)

The GPU training job then loads these three files and skips PCA / normalization entirely.
"""

import sys
sys.path.insert(0, '/home/dongbiao/software/DeepPhylo')

import os
import time
import argparse
import numpy as np
import pandas as pd
import biom

from deepphylo.pre_dataset import set_seed, reducer


def parse_args():
    p = argparse.ArgumentParser(description="Precompute DeepPhylo embedding + abundance matrix")
    p.add_argument("--biom",
                   default="../../data/oral-gut-pair-v2-qc/biom_oral_gut_pair_v2_qc.biom")
    p.add_argument("--metadata",
                   default="../../data/oral-gut-pair-v2-qc/metadata_oral_gut_pair_v2_merge_qc.csv")
    p.add_argument("--cophenetic", default="data/cophenetic_matrix.csv",
                   help="OTU x OTU cophenetic distance matrix (index = feature ids)")
    p.add_argument("--output_dir", default="data/")
    p.add_argument("-hs", "--hidden_size", type=int, default=80,
                   help="PCA embedding dimension (must be < number of OTUs)")
    p.add_argument("--prevalence", type=int, default=1,
                   help="keep OTUs present in more than this many samples")
    p.add_argument("--min_reads", type=int, default=100,
                   help="keep samples with more than this many reads")
    return p.parse_args()


def qc_table(biom_path, metadata_path, prevalence, min_reads):
    """Reproduce the QC pipeline; return the filtered table."""
    table = biom.load_table(biom_path)
    metadata = pd.read_csv(metadata_path, sep=",", low_memory=False)

    # Keep only samples present in BOTH metadata and table (avoids KeyError)
    table_samples = set(table.ids(axis="sample"))
    sid = [s for s in metadata.sample_id.values if s in table_samples]
    table.filter(sid, axis="sample")
    table.remove_empty()

    # Prevalence filter on OTUs
    fid = table.ids(axis="observation")
    prev = table.nonzero_counts(axis="observation")
    table = table.filter(fid[prev > prevalence], axis="observation", inplace=False)
    table.remove_empty()

    # Read-depth filter on samples
    sid_reads = table.sum(axis="sample")
    sid = table.ids(axis="sample")
    table = table.filter(sid[sid_reads > min_reads], axis="sample", inplace=False)
    table.remove_empty()

    return table


def main():
    args = parse_args()
    set_seed(1234)
    os.makedirs(args.output_dir, exist_ok=True)

    t0 = time.time()
    print("[1/4] Loading + QC biom table ...", flush=True)
    table = qc_table(args.biom, args.metadata, args.prevalence, args.min_reads)
    print(f"      table: {table.shape[0]} OTUs x {table.shape[1]} samples "
          f"({time.time()-t0:.1f}s)", flush=True)

    t1 = time.time()
    print("[2/4] Loading cophenetic matrix ...", flush=True)
    C_full = pd.read_csv(args.cophenetic, index_col=0)
    C_full.columns = C_full.index  # ensure square / consistent labels
    print(f"      C: {C_full.shape} ({time.time()-t1:.1f}s)", flush=True)

    # Fixed global OTU set = QC OTUs that also exist in C
    otu_ids = [o for o in table.ids(axis="observation") if o in C_full.index]
    n_missing = table.shape[0] - len(otu_ids)
    if n_missing:
        print(f"      [WARN] {n_missing} OTUs absent from C, dropped", flush=True)
    if len(otu_ids) <= args.hidden_size:
        raise ValueError(f"Only {len(otu_ids)} OTUs but hidden_size={args.hidden_size}")
    print(f"      global OTUs: {len(otu_ids)}", flush=True)

    np.save(os.path.join(args.output_dir, "otu_ids.npy"),
            np.array(otu_ids, dtype=object))

    t2 = time.time()
    print(f"[3/4] Computing PCA embedding (hs={args.hidden_size}) ...", flush=True)
    emb_path = os.path.join(args.output_dir, f"phy_emb_global_hs{args.hidden_size}.npy")

    C_global = C_full.loc[otu_ids, otu_ids].values.astype(np.float32)
    phy_emb = reducer(C_global, "pca", args.hidden_size, whiten=True)
    np.save(emb_path, phy_emb)
    print(f"      embedding {phy_emb.shape} ({time.time()-t2:.1f}s)", flush=True)

    t3 = time.time()

    print(f"\nDONE in {time.time()-t0:.1f}s. Artifacts in {args.output_dir}:")
    print(f"  otu_ids.npy                      ({len(otu_ids)} OTUs)")
    print(f"  phy_emb_global_hs{args.hidden_size}.npy   {phy_emb.shape}")

if __name__ == "__main__":
    main()
