#!/usr/bin/env python3
"""
pipeline_v2 Step 03: 97% self-clustering of clean candidate ASVs.
Groups ASVs that are ≥97% identical → each cluster ≈ a tentative phylotype.
This merges cross-region/cross-BP variants before support floor filtering.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import os, sys, subprocess, json
from collections import defaultdict

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
os.makedirs(PIPE_V2, exist_ok=True)
log = lambda s: print(f"[03_cluster] {s}", flush=True)

# ============================================================
# 1. Load clean ASVs
# ============================================================
log("Loading clean candidates...")
df = pd.read_csv(PIPE_V2 / "02_contamination_flagged.tsv", sep="\t", low_memory=False)
clean = df[~df["contamination_flag"]].copy()
log(f"  Clean ASVs: {len(clean):,}")

for c in ["n_bioprojects", "total_reads", "prevalence_samples"]:
    if c in clean.columns:
        clean[c] = pd.to_numeric(clean[c], errors="coerce").fillna(0)

# ============================================================
# 2. Write FASTA, run vsearch --cluster_fast --id 0.97
# ============================================================
cluster_dir = PIPE_V2 / "tmp_self_97cluster"
os.makedirs(cluster_dir, exist_ok=True)
fa = cluster_dir / "clean_candidates.fasta"
uc_out = cluster_dir / "clusters.uc"

if not uc_out.exists():
    log(f"Writing FASTA: {len(clean):,} sequences...")
    with open(fa, "w") as f:
        for _, row in clean.iterrows():
            f.write(f">{row['derep_id']}\n{row['sequence']}\n")

    log("Running vsearch --cluster_fast --id 0.97...")
    vsearch_bin = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"
    subprocess.run([
        vsearch_bin, "--cluster_fast", str(fa),
        "--id", "0.97", "--uc", str(uc_out),
        "--threads", "16", "--quiet"
    ], check=True)
    log("  vsearch done")
else:
    log(f"  Using existing UC: {uc_out}")

# ============================================================
# 3. Parse UC → cluster-level summary
# ============================================================
log("Parsing clusters...")

cluster_members = defaultdict(list)
with open(uc_out) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        rec_type, cluster_id, label = parts[0], parts[1], parts[8]
        if rec_type in ("H", "S"):
            cluster_members[cluster_id].append(label)

n_clusters = len(cluster_members)
log(f"  Clusters: {n_clusters:,}")
log(f"  Avg ASVs/cluster: {len(clean)/n_clusters:.1f}")

# Map derep_id → cluster_id
derep_to_cluster = {}
for cid, members in cluster_members.items():
    for m in members:
        derep_to_cluster[m] = cid

# Aggregate per cluster
clean["cluster_id"] = clean["derep_id"].map(derep_to_cluster)

cluster_stats = clean.groupby("cluster_id").agg(
    n_asvs=("derep_id", "count"),
    n_bioprojects_total=("n_bioprojects", "sum"),
    total_reads=("total_reads", "sum"),
    prevalence_samples=("prevalence_samples", "sum"),
    rep_derep=("derep_id", "first"),
    rep_sequence=("sequence", "first"),
    gg2_genus=("gg2_genus", lambda x: x.value_counts().index[0] if len(x) > 0 else "Unknown"),
).reset_index()

log(f"\n=== Cluster-level summary ===")
log(f"Total clusters: {len(cluster_stats):,}")

for lo, hi, label in [(1,2,"1 ASV"),(2,5,"2-4"),(5,10,"5-9"),(10,100,"10-99"),(100,10000,"100+")]:
    n = len(cluster_stats[(cluster_stats["n_asvs"]>=lo)&(cluster_stats["n_asvs"]<hi)])
    log(f"  {label}: {n:,} clusters")

log(f"\nBy cluster-level n_bioprojects (sum):")
for lo, hi, label in [(0,2,"1"),(2,5,"2-4"),(5,10,"5-9"),(10,20,"10-19"),(20,50,"20-49"),(50,1000,"50+")]:
    n = len(cluster_stats[(cluster_stats["n_bioprojects_total"]>=lo)&(cluster_stats["n_bioprojects_total"]<hi)])
    log(f"  {label} BPs: {n:,} clusters")

# Cluster level support floors
log(f"\n=== Effect of cluster-level support floor ===")

for bp_min, rd_min, label in [(2,50,"≥2 BP, ≥50 reads"),
                                (2,20,"≥2 BP, ≥20 reads"),
                                (2,1,"≥2 BP only"),
                                (3,50,"≥3 BP, ≥50 reads"),
                                (5,50,"≥5 BP, ≥50 reads")]:
    n_pass = len(cluster_stats[(cluster_stats["n_bioprojects_total"]>=bp_min)&(cluster_stats["total_reads"]>=rd_min)])
    bp_check = len(cluster_stats[cluster_stats["n_bioprojects_total"]>=bp_min])
    n_asvs_pass = cluster_stats[(cluster_stats["n_bioprojects_total"]>=bp_min)&(cluster_stats["total_reads"]>=rd_min)]["n_asvs"].sum()
    log(f"  {label}: {n_pass:,} clusters ({n_asvs_pass:,} ASVs)")

# Save
outpath = PIPE_V2 / "03_cluster_stats.tsv"
cluster_stats.to_csv(outpath, sep="\t", index=False)
log(f"\nSaved: {outpath} ({len(cluster_stats):,} rows)")

# Add cluster_id to clean table and save
clean_out = PIPE_V2 / "03_clean_with_clusters.tsv"
clean.to_csv(clean_out, sep="\t", index=False)
log(f"Saved: {clean_out}")

summary = {
    "clean_asvs": int(len(clean)),
    "clusters": int(n_clusters),
    "avg_asvs_per_cluster": round(len(clean)/n_clusters, 1),
    "clusters_2bp_50reads": int(len(cluster_stats[(cluster_stats["n_bioprojects_total"]>=2)&(cluster_stats["total_reads"]>=50)])),
    "clusters_2bp_20reads": int(len(cluster_stats[(cluster_stats["n_bioprojects_total"]>=2)&(cluster_stats["total_reads"]>=20)])),
    "clusters_2bp_any": int(len(cluster_stats[cluster_stats["n_bioprojects_total"]>=2])),
}
with open(PIPE_V2 / "03_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

log("Done!")
