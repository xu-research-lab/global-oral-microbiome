#!/usr/bin/env python
"""Cluster x sample feature profile for the 486 supported clusters.

Rows are the supported clusters of 08_supported_clusters.tsv; columns are the
samples of the deduplicated oral QC metadata (metadata_oral_qc_dedup.tsv,
51,186 samples, 288 studies).  A cluster's value in a sample is the number of
reads carried by its member ASVs, obtained from each project's DADA2.biom,
whose observation ids are the ASV sequences themselves.  The cluster members
are exact-derep sequences (08 -> 07 -> sequence), so the match to the biom is
exact.

Relative abundance is the cluster's reads divided by the sample's total reads
over all ASVs in that table -- not over the member ASVs only -- so the value
is the cluster's share of the whole community.

Outputs (all with the same row and column order):
  14_cluster_counts_dedup.tsv        integer read counts
  14_cluster_relabund_dedup.tsv     counts / sample total reads
  14_cluster_sample_metadata.tsv    sample_key first, then the QC metadata
  14_supported_cluster_annotation.tsv  cluster-level annotation for the rows
"""
from __future__ import annotations

import os
import sys
import json

import numpy as np
import pandas as pd
from scipy import sparse
import biom

PIPE = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
CLUSTERS = os.path.join(PIPE, "08_supported_clusters.tsv")
MASTER = os.path.join(PIPE, "07_master_asv_table.tsv")
META = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv"
ASV_ROOT = "/beegfs/dongbiao/oral"

OUT_COUNT = os.path.join(PIPE, "14_cluster_counts_dedup.tsv")
OUT_REL = os.path.join(PIPE, "14_cluster_relabund_dedup.tsv")
OUT_META = os.path.join(PIPE, "14_cluster_sample_metadata.tsv")
OUT_ANN = os.path.join(PIPE, "14_supported_cluster_annotation.tsv")
OUT_SUM = os.path.join(PIPE, "14_summary.json")


def load_members():
    """cluster_id -> list of member derep ids, in the order of the cluster file."""
    cl = pd.read_csv(CLUSTERS, sep="\t", dtype={"cluster_id": str}, low_memory=False)
    members = {}
    for cid, ids in zip(cl["cluster_id"], cl["member_derep_ids"].astype(str)):
        members[cid] = [x for x in ids.split(";") if x]
    return cl, members


def load_member_sequences(members):
    """derep id -> sequence, for the member dereps only."""
    want = {d for ids in members.values() for d in ids}
    seq = {}
    for chunk in pd.read_csv(MASTER, sep="\t", usecols=["derep_id", "sequence"],
                             dtype=str, chunksize=200_000, low_memory=False):
        sub = chunk[chunk["derep_id"].isin(want)]
        seq.update(dict(zip(sub["derep_id"], sub["sequence"])))
    missing = want - set(seq)
    if missing:
        raise SystemExit(f"{len(missing)} member dereps have no sequence, e.g. {sorted(missing)[:5]}")
    return seq


def main():
    # a project limit exists only to smoke-test the script on a few tables
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    suffix = ".test" if limit else ""

    cl, members = load_members()
    derep_seq = load_member_sequences(members)
    print(f"[14] {len(members)} clusters, {sum(len(v) for v in members.values())} member dereps")

    # member sequence -> cluster.  Exact derep means one sequence per derep, so a
    # sequence maps to a single cluster.
    seq_cluster = {}
    for cid, ids in members.items():
        for d in ids:
            s = derep_seq[d]
            if s in seq_cluster and seq_cluster[s] != cid:
                raise SystemExit(f"sequence shared by clusters {seq_cluster[s]} and {cid}")
            seq_cluster[s] = cid

    meta = pd.read_csv(META, sep="\t", dtype=str, low_memory=False)
    meta.insert(0, "sample_key", meta["bioproject_id"] + "|" + meta["sample_id"])
    # the metadata's own total_reads is the read count entering closed-reference
    # mapping, not the ASV table's; the denominator used here is written out
    # explicitly so it cannot be confused with it

    if meta["sample_key"].duplicated().any():
        raise SystemExit("duplicate sample_key in metadata")
    col_index = {k: j for j, k in enumerate(meta["sample_key"])}
    n_samples = len(meta)
    print(f"[14] metadata: {n_samples} samples, {meta['bioproject_id'].nunique()} projects")

    cluster_ids = list(members)
    row_index = {c: i for i, c in enumerate(cluster_ids)}
    counts = np.zeros((len(cluster_ids), n_samples), dtype=np.int64)
    totals = np.zeros(n_samples, dtype=np.int64)
    seen = np.zeros(n_samples, dtype=bool)

    projects = sorted(meta["bioproject_id"].unique())
    if limit:
        projects = projects[:limit]
    for i, bp in enumerate(projects, 1):
        path = os.path.join(ASV_ROOT, bp, "results", "DADA2.biom")
        if not os.path.exists(path):
            raise SystemExit(f"missing biom: {path}")
        table = biom.load_table(path)
        sid = [str(x) for x in table.ids(axis="sample")]
        obs = [str(x) for x in table.ids(axis="observation")]

        # samples of this project that survive into the metadata
        gcols, lcols = [], []
        for k, s in enumerate(sid):
            j = col_index.get(f"{bp}|{s}")
            if j is not None and not seen[j]:
                gcols.append(j)
                lcols.append(k)
        if not gcols:
            continue
        gcols = np.asarray(gcols, dtype=np.int64)
        lcols = np.asarray(lcols, dtype=np.int64)

        mat = table.matrix_data.tocsr()
        sub = mat[:, lcols]                       # all ASVs x selected samples

        totals[gcols] = np.asarray(sub.sum(axis=0)).ravel().astype(np.int64)
        seen[gcols] = True

        # rows of this biom that are members of a supported cluster
        rows, aidx = [], []
        for r, s in enumerate(obs):
            c = seq_cluster.get(s)
            if c is not None:
                rows.append(r)
                aidx.append(row_index[c])
        if rows:
            A = sparse.csr_matrix((np.ones(len(rows)), (aidx, rows)),
                                  shape=(len(cluster_ids), len(obs)))
            C = (A @ sub).toarray()
            counts[:, gcols] += np.rint(C).astype(np.int64)

        if i % 25 == 0 or i == len(projects):
            print(f"[14] {i}/{len(projects)} projects", flush=True)

    if not seen.all():
        print(f"[14] WARNING: {int((~seen).sum())} metadata samples not found in any biom",
              file=sys.stderr)

    rel = np.divide(counts, totals[None, :],
                    out=np.zeros_like(counts, dtype=np.float64),
                    where=totals[None, :] > 0)

    cols = meta["sample_key"].tolist()
    print(f"[14] matrix {counts.shape}; non-zero cells "
          f"{int((counts > 0).sum()):,} / {counts.size:,}; "
          f"clusters with no reads {int((counts.sum(axis=1) == 0).sum())}")

    pd.DataFrame(counts, index=cluster_ids, columns=cols).to_csv(
        OUT_COUNT.replace(".tsv", suffix + ".tsv"), sep="\t", index_label="cluster_id")
    pd.DataFrame(rel, index=cluster_ids, columns=cols).to_csv(
        OUT_REL.replace(".tsv", suffix + ".tsv"), sep="\t", index_label="cluster_id", float_format="%.6g")
    meta["total_reads_asv_table"] = totals
    meta["sum_cluster_reads"] = counts.sum(axis=0)
    meta.to_csv(OUT_META.replace(".tsv", suffix + ".tsv"), sep="\t", index=False)

    ann_cols = ["cluster_id", "n_asvs", "n_refs_in_cluster", "annotation",
                "annotation_level", "annotation_basis", "annotation_lineage",
                "cluster_novelty_class", "candidate_novel", "tree_tier", "headline",
                "max_subject_prevalence", "n_bioprojects_total", "n_continents",
                "total_reads_sum", "rep_homd_identity", "rep_homd_target"]
    ann = cl[[c for c in ann_cols if c in cl.columns]].copy()
    ann["phylum"] = ann["annotation_lineage"].fillna("").str.split(";").apply(
        lambda x: x[1].strip() if len(x) > 1 else "")
    ann["cluster_total_reads_dedup"] = counts.sum(axis=1)
    ann["n_samples_detected"] = (counts > 0).sum(axis=1)
    ann["prevalence_dedup"] = ann["n_samples_detected"] / n_samples
    ann.to_csv(OUT_ANN.replace(".tsv", suffix + ".tsv"), sep="\t", index=False)

    summary = {
        "n_clusters": len(cluster_ids),
        "n_member_dereps": sum(len(v) for v in members.values()),
        "n_samples": int(n_samples),
        "n_projects": int(len(projects)),
        "samples_not_found": int((~seen).sum()),
        "nonzero_cells": int((counts > 0).sum()),
        "clusters_without_reads": int((counts.sum(axis=1) == 0).sum()),
        "total_reads_all_asvs": int(totals.sum()),
        "total_cluster_reads": int(counts.sum()),
    }
    with open(OUT_SUM.replace(".json", suffix + ".json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("[14] " + json.dumps(summary))


if __name__ == "__main__":
    main()
