#!/usr/bin/env python3
"""
Step 12 - master per-ASV supplementary table.

One row per placed ASV (54,277). Column blocks:
  A. sequence + support      : asv_id, sequence, n_asvs, n_bioprojects,
                               total_reads, prevalence_samples, subject_prevalence
  B. TreeCluster cluster     : cluster_id, cluster_size, n_refs_in_cluster,
                               cluster_novelty_class, candidate_novel, clean_cluster
  C. tree annotation         : annotation, annotation_level, annotation_basis,
                               nearest_ref, nearest_ref_species, nearest_ref_genus,
                               nearest_ref_dist, nearest_ref_in_homd,
                               nearest_ref_renamed, annotation_resolution
  D. HOMD status             : covered, nearest_homd_ref, nearest_homd_ref_species,
                               nearest_homd_ref_dist, homd_identity, homd_target,
                               has_homd_hit
  E. geography / regions     : bioprojects, amplicon_regions, sequencing_regions,
                               continents, countries
  F. contamination flags     : contamination_flag, flag_type

Output: 12_master_asv_table.tsv
"""
import json
from pathlib import Path

import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE = BASE / "pipeline"
PIPE_V3 = PIPE
PIPE_V2 = PIPE
log = lambda s: print(f"[12_master_table] {s}", flush=True)


def main():
    log("Loading 10_asv_clusters.tsv ...")
    asv = pd.read_csv(PIPE_V3 / "10_clustering/results/10_asv_clusters.tsv", sep="\t", low_memory=False)
    asv["asv_id"] = asv["asv_id"].astype(str)
    log(f"  {len(asv):,} ASVs")

    log("Loading support-floor metadata (step 09)...")
    base = pd.read_csv(PIPE_V2 / "09_support_floor/results/09_support_floor_passed.tsv", sep="\t", low_memory=False)
    base["derep_id"] = base["derep_id"].astype(str)
    base = base[["derep_id", "sequence", "n_bioprojects", "total_reads",
                 "prevalence_samples", "contamination_flag", "flag_type",
                 "bioprojects", "amplicon_regions", "sequencing_regions",
                 "continents", "countries"]]
    log(f"  {len(base):,} ASVs")

    log("Loading subject prevalence (step 09)...")
    subj = pd.read_csv(PIPE_V2 / "09_support_floor/results/09_subject_prevalence.tsv", sep="\t",
                       low_memory=False, usecols=["derep_id", "subject_prevalence"])
    subj["derep_id"] = subj["derep_id"].astype(str)
    base = base.merge(subj, on="derep_id", how="left")

    log("Loading HOMD mapping (step 04)...")
    h01 = pd.read_csv(PIPE_V2 / "04_candidate_selection/results/04_candidates.tsv", sep="\t", low_memory=False,
                      usecols=["derep_id", "n_asvs", "homd_target", "identity",
                               "has_homd_hit"])
    h01["derep_id"] = h01["derep_id"].astype(str)
    h01 = h01.rename(columns={"identity": "homd_identity"})
    base = base.merge(h01, on="derep_id", how="left")

    log("Loading cluster-level evidence (step 10)...")
    cl = pd.read_csv(PIPE_V3 / "10_clustering/results/10_clusters_aggregated.tsv", sep="\t", low_memory=False)
    cl_keep = ["cluster_id", "cluster_novelty_class", "candidate_novel", "clean_cluster",
               "n_asvs", "max_subject_prevalence", "n_bioprojects_total", "n_regions",
               "n_continents", "regions", "total_reads_sum", "support_score"]
    cl = cl[cl_keep].rename(columns={
        "n_asvs": "cluster_n_asvs",
        "max_subject_prevalence": "cluster_max_subject_prevalence",
        "n_bioprojects_total": "cluster_n_bioprojects",
        "n_regions": "cluster_n_regions",
        "n_continents": "cluster_n_continents",
        "regions": "cluster_regions",
        "total_reads_sum": "cluster_total_reads",
        "support_score": "cluster_support_score"})

    df = asv.merge(base, left_on="asv_id", right_on="derep_id", how="left")
    df = df.merge(cl, on="cluster_id", how="left")

    first = ["asv_id", "sequence",
             "n_asvs", "n_bioprojects", "total_reads", "prevalence_samples",
             "subject_prevalence",
             "cluster_id", "cluster_size", "n_refs_in_cluster",
             "cluster_novelty_class", "candidate_novel", "clean_cluster",
             "cluster_n_asvs", "cluster_max_subject_prevalence",
             "cluster_n_bioprojects", "cluster_n_regions", "cluster_n_continents",
             "cluster_regions", "cluster_total_reads", "cluster_support_score",
             "annotation", "annotation_level", "annotation_basis",
             "nearest_ref", "nearest_ref_species", "nearest_ref_genus",
             "nearest_ref_dist", "nearest_ref_in_homd", "nearest_ref_renamed",
             "annotation_resolution",
             "covered", "nearest_homd_ref", "nearest_homd_ref_species",
             "nearest_homd_ref_dist", "homd_identity", "homd_target", "has_homd_hit",
             "bioprojects", "amplicon_regions", "sequencing_regions",
             "continents", "countries",
             "contamination_flag", "flag_type"]
    rest = [c for c in df.columns if c not in first]
    df = df[first + rest]

    out = PIPE_V3 / "12_master_table/results/12_master_asv_table.tsv"
    df.to_csv(out, sep="\t", index=False)
    log(f"  saved {out} ({len(df):,} rows x {len(df.columns)} cols, "
        f"{out.stat().st_size / 1e6:.1f} MB)")

    summary = {
        "n_rows": int(len(df)),
        "n_columns": int(len(df.columns)),
        "n_with_cluster": int(df["cluster_id"].notna().sum()),
        "n_covered": int(df["covered"].astype(str).str.lower().isin(["true", "1"]).sum()),
        "n_novel": int((~df["covered"].astype(str).str.lower().isin(["true", "1"])).sum()),
        "n_contamination_flagged": int(df["contamination_flag"].astype(str)
                                       .str.lower().isin(["true", "1"]).sum()),
        "n_in_candidate_novel_clusters": int(
            df["candidate_novel"].fillna(False).astype(bool).sum()),
    }
    with open(PIPE_V3 / "12_master_table/results/12_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log("  saved 12_summary.json")
    log("Done!")


if __name__ == "__main__":
    main()
