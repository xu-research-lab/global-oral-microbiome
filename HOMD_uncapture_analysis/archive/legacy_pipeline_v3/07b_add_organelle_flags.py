#!/usr/bin/env python3
"""
pipeline_v3 Step 07b: annotate organelle (mitochondrion / chloroplast) flags.

ASVs were mapped with vsearch against the Greengenes2 organelle reference set
(gg2_organelle_ref.fasta, extracted via GG2 taxonomy f__Mitochondria /
c__Chloroplast). This step re-attaches the organelle evidence:

  - 02_clusters_aggregated.tsv  : + n_organelles, organelle_targets, any_organelle
  - 07_master_asv_table.tsv     : + organelle_flag, organelle_target, organelle_identity
                                  and flag_type gets ";organelle" appended

Outputs are written in place.
"""
import os
from pathlib import Path

import pandas as pd

PIPE_V3 = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
GG2_TAX = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"
log = lambda s: print(f"[07b_organelle] {s}", flush=True)


def load_mapping():
    """organelle_mapping.tsv: derep_id, ref, identity, coverage (vsearch).
    Keep the best hit per ASV; classify the ref via GG2 taxonomy."""
    om = pd.read_csv(os.path.join(PIPE_V3, "organelle_mapping.tsv"), sep="\t", header=None,
                     names=["derep_id", "ref", "identity", "coverage"])
    om["derep_id"] = om["derep_id"].astype(str)

    tax = {}
    with open(GG2_TAX) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == 2:
                tax[p[0]] = p[1]

    def org_type(ref):
        t = tax.get(ref, "")
        if "f__Mitochondria" in t:
            return "mitochondrion"
        if "c__Chloroplast" in t:
            return "chloroplast"
        return "unknown"

    om["organelle_target"] = om["ref"].map(org_type)
    # best hit per ASV: highest identity, then coverage
    om = om.sort_values(["identity", "coverage"], ascending=False).drop_duplicates("derep_id")
    return om


def main():
    omap = load_mapping()
    log(f"organelle-mapped ASVs: {len(omap):,} "
        f"(chloroplast={(omap['organelle_target']=='chloroplast').sum():,} | "
        f"mitochondrion={(omap['organelle_target']=='mitochondrion').sum():,})")

    # ---- 1. cluster-level table: n_organelles / targets / any_organelle ----
    agg = pd.read_csv(os.path.join(PIPE_V3, "02_clusters_aggregated.tsv"), sep="\t",
                      low_memory=False)
    cl = agg[["cluster_id", "rep_derep_id", "member_derep_ids"]].copy()
    cl_members = cl.assign(members=cl["member_derep_ids"].str.split(";")).explode("members")
    cl_org = cl_members.merge(omap[["derep_id", "organelle_target"]],
                              left_on="members", right_on="derep_id", how="inner")
    cl_org = cl_org.groupby("cluster_id").agg(
        n_organelles=("members", "size"),
        organelle_targets=("organelle_target", lambda x: ";".join(sorted(set(x)))),
    ).reset_index()
    cl_org["any_organelle"] = cl_org["n_organelles"] > 0
    agg = agg.merge(cl_org, on="cluster_id", how="left")
    agg["n_organelles"] = agg["n_organelles"].fillna(0).astype(int)
    agg["organelle_targets"] = agg["organelle_targets"].fillna("")
    agg["any_organelle"] = agg["any_organelle"].fillna(False)
    if "organelle_targets" in agg.columns and agg["organelle_targets"].duplicated().any():
        # keep existing column position (append was avoided if already present)
        pass
    agg.to_csv(os.path.join(PIPE_V3, "02_clusters_aggregated.tsv"), sep="\t", index=False)
    log(f"02_clusters_aggregated.tsv updated "
        f"({len(agg):,} clusters, {(agg['any_organelle']).sum():,} organelle-flagged)")

    # ---- 2. master ASV table: organelle_flag / target / identity + flag_type ----
    d = pd.read_csv(os.path.join(PIPE_V3, "07_master_asv_table.tsv"), sep="\t",
                    low_memory=False)
    if "organelle_flag" not in d.columns:
        d = d.merge(omap[["derep_id", "organelle_target", "identity"]].rename(
            columns={"organelle_target": "_org_target", "identity": "_org_identity"}),
            on="derep_id", how="left")
        d["organelle_flag"] = d["_org_target"].notna()
        d["organelle_target"] = d["_org_target"].fillna("")
        d["organelle_identity"] = d["_org_identity"].fillna("")
        d = d.drop(columns=["_org_target", "_org_identity"])

        # append ";organelle" to flag_type for organelle-flagged ASVs
        org_mask = d["organelle_flag"].fillna(False)
        ft = d["flag_type"].fillna("clean").astype(str)
        d.loc[org_mask, "flag_type"] = ft[org_mask].map(
            lambda s: s if "organelle" in s else f"{s};organelle")
    else:
        log("master table already has organelle columns - skipped")
    d.to_csv(os.path.join(PIPE_V3, "07_master_asv_table.tsv"), sep="\t", index=False)
    log(f"07_master_asv_table.tsv updated ({len(d):,} ASVs, "
        f"{int(d['organelle_flag'].fillna(False).sum()):,} organelle-flagged)")


if __name__ == "__main__":
    main()
