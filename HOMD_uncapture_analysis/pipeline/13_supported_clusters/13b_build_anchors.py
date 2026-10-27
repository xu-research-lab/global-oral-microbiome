#!/usr/bin/env python3
"""Step 13b - anchor table for the supported catalog.

One row per supported-catalog cluster: representative tip, cluster metrics,
display label, oral_clade and phylum.  It is the per-cluster input of the
profile panel.

oral_clade = annotation's first token is a genus in the HOMD genus set AND the
cluster is HOMD-uncaptured (same concept as v2: novel lineages inside vs
outside known oral genera).
"""
import os
from pathlib import Path

import pandas as pd

PIPE = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline"
PIPE_V3 = PIPE
HOMD_FASTA = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta"
log = lambda s: print(f"[13_anchors] {s}", flush=True)


def load_homd_genus():
    g = set()
    with open(HOMD_FASTA) as f:
        for line in f:
            if line.startswith(">"):
                parts = line.split()
                if len(parts) >= 3:
                    g.add(parts[2].strip())
    return g


def main():
    sup = pd.read_csv(os.path.join(PIPE_V3, "13_supported_clusters/results/13_supported_clusters.tsv"), sep="\t",
                      low_memory=False)
    homd_genus = load_homd_genus()
    log(f"supported clusters: {len(sup):,}")

    phy = {}
    with open(os.path.join(PIPE_V3, "11_geography/results/11_candidate_novel_with_phylum.tsv")) as f:
        for r in pd.read_csv(f, sep="\t", usecols=["cluster_id", "gg2_phylum"],
                              dtype={"cluster_id": str}).itertuples(index=False):
            phy[r.cluster_id] = r.gg2_phylum

    def oral(row):
        first = str(row["annotation"]).split(" ")[0]
        return bool(first in homd_genus and row["cluster_novelty_class"] == "homd_uncaptured")

    sup["oral_clade"] = sup.apply(oral, axis=1)
    sup["phylum"] = sup["cluster_id"].astype(str).map(phy).fillna("unknown")
    sup["tree_label"] = sup["annotation"]

    out = pd.DataFrame({
        "anchor_tip": sup["rep_derep_id"],
        "cluster_id": sup["cluster_id"],
        "n_members": sup["n_asvs"],
        "max_subject_prevalence": sup["max_subject_prevalence"],
        "n_bioprojects_total": sup["n_bioprojects_total"],
        "n_continents": sup["n_continents"],
        "n_regions": sup["n_regions"],
        "oral_clade": sup["oral_clade"],
        "tree_label": sup["tree_label"],
        "phylum": sup["phylum"],
        "rep_nearest_ref_dist": sup["rep_nearest_ref_dist"],
        "headline": sup["headline"],
    })
    out = out.sort_values("max_subject_prevalence", ascending=False)
    out.to_csv(os.path.join(PIPE_V3, "13_supported_clusters/results/13_anchors_supported.tsv"), sep="\t", index=False)
    log(f"saved 13_anchors_supported.tsv ({len(out):,} rows)")


if __name__ == "__main__":
    main()
