#!/usr/bin/env python3
"""
pipeline_v2 Step 02: Two-Stage Contamination Firewall

Stage 1: 97% OTU clustering with control HOMD-uncaptured ASVs
  - No taxonomy needed, pure sequence-based

Stage 2: GG2 genus annotation (via taxonomy.tsv) + contaminant genus list
  - Only for ASVs NOT flagged in Stage 1
  - vsearch vs GG2 backbone → lookup genus in taxonomy.tsv
  - Flag if genus in contaminant list (Salter 2014 + Weyrich 2019 + controls)
"""
import pandas as pd
import numpy as np
import re as re_mod
import os, sys, subprocess, json
from pathlib import Path
from collections import defaultdict, Counter

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
os.makedirs(PIPE_V2, exist_ok=True)
log = lambda s: print(f"[02_contam] {s}", flush=True)

# ============================================================
# 1. Load candidates
# ============================================================
log("Loading candidates from step 01...")
candidates = pd.read_csv(PIPE_V2 / "01_chimera_filtered.tsv", sep="\t", low_memory=False)
log(f"  {len(candidates):,} candidate ASVs (HOMD-uncaptured + non-chimeric)")

# ============================================================
# 2. Control ASVs — ALL control-catalog ASVs (no HOMD filter)
#    (decision 2026-08-05: include every control-sourced ASV in Stage 1
#     clustering, not only the HOMD-uncaptured subset — common environmental
#     contaminants are often HOMD-captured, so they must also flag candidates)
# ============================================================
log("\nLoading control ASV data (all catalog ASVs)...")
catalog = pd.read_csv(BASE / "contaminant_catalog/contaminant_asv_catalog.tsv", sep="\t", low_memory=False)
ctrl_uncaptured_seqs = catalog[["asv_id", "sequence"]].drop_duplicates(subset="sequence")
log(f"  Control ASVs (all): {len(ctrl_uncaptured_seqs):,}")

# ============================================================
# Stage 1: 97% clustering
# ============================================================
log("\n=== Stage 1: 97% Clustering with ALL Control ASVs ===")
cluster_dir = PIPE_V2 / "tmp_97cluster"
os.makedirs(cluster_dir, exist_ok=True)
combined_fa = cluster_dir / "combined_all_ctrl_candidates.fasta"
uc_out = cluster_dir / "clusters_all_ctrl.uc"

if not uc_out.exists():
    log("Building combined FASTA...")
    with open(combined_fa, "w") as f:
        for _, row in ctrl_uncaptured_seqs.iterrows():
            f.write(f">CTRL|{row['asv_id']}\n{row['sequence']}\n")
        for _, row in candidates.iterrows():
            f.write(f">CAND|{row['derep_id']}\n{row['sequence']}\n")
    n_ctrl, n_cand = len(ctrl_uncaptured_seqs), len(candidates)
    log(f"  Combined FASTA: {n_ctrl:,} control + {n_cand:,} candidate = {n_ctrl+n_cand:,} total")

    log("Running vsearch --cluster_fast --id 0.97...")
    vsearch_bin = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"
    subprocess.run([
        vsearch_bin, "--cluster_fast", str(combined_fa),
        "--id", "0.97", "--uc", str(uc_out),
        "--threads", "26", "--quiet"
    ], check=True)
    log("  vsearch done")
else:
    log(f"  Using existing: {uc_out}")

# Parse UC
log("Parsing UC clusters...")
cluster_members = defaultdict(list)
with open(uc_out) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 9: continue
        rec_type, cluster_id, label = parts[0], parts[1], parts[8]
        if rec_type in ("H", "S"):
            is_ctrl = label.startswith("CTRL|")
            cluster_members[cluster_id].append((label.split("|", 1)[1], is_ctrl))

contaminated_clusters = {cid for cid, m in cluster_members.items()
                         if any(is_ctrl for _, is_ctrl in m)}

cluster97_contaminated = set()
for cid in contaminated_clusters:
    for derep_id, is_ctrl in cluster_members[cid]:
        if not is_ctrl:
            cluster97_contaminated.add(derep_id)

log(f"  Total clusters: {len(cluster_members):,}")
log(f"  Clusters with control HOMD-uncaptured: {len(contaminated_clusters):,}")
log(f"  Candidate ASVs flagged: {len(cluster97_contaminated):,}")

candidates["stage1_cluster97_flag"] = candidates["derep_id"].isin(cluster97_contaminated)
n_stage1 = candidates["stage1_cluster97_flag"].sum()
log(f"  Stage 1 flagged: {n_stage1:,} / {len(candidates):,} ({n_stage1/len(candidates)*100:.1f}%)")

# ============================================================
# Stage 2: GG2 genus via taxonomy.tsv
# ============================================================
log("\n=== Stage 2: GG2 Genus Annotation + Contaminant List ===")

# Load GG2 taxonomy → genus lookup
log("Loading GG2 taxonomy.tsv...")
gg2_tax_path = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"
tax_df = pd.read_csv(gg2_tax_path, sep="\t", low_memory=False, dtype=str)
gg2_tax_map = {}
for _, row in tax_df.iterrows():
    fid = str(row["Feature ID"]).replace("-", "_")  # normalize: hyphens → underscores
    g_match = re_mod.search(r'g__([^;]+)', str(row["Taxon"]))
    if g_match:
        gg2_tax_map[fid] = g_match.group(1)
log(f"  Genus-mapped GG2 IDs: {len(gg2_tax_map):,} / {len(tax_df):,}")

# Only annotate Stage 1 clean ASVs
to_annotate = candidates[~candidates["stage1_cluster97_flag"]]
log(f"  ASVs to annotate (Stage 1 clean): {len(to_annotate):,}")

# Run vsearch vs GG2 backbone
gg2_fa = cluster_dir / "stage2_to_annotate.fasta"
with open(gg2_fa, "w") as f:
    for _, row in to_annotate.iterrows():
        f.write(f">{row['derep_id']}\n{row['sequence']}\n")

gg2_db = "/beegfs/dongbiao/greengene2/backbone_fasta/dna-sequences.fasta"
gg2_hits = cluster_dir / "gg2_annotations_v2.tsv"

if not gg2_hits.exists():
    log(f"  Running vsearch vs GG2 ({len(to_annotate):,} queries)...")
    vsearch_bin = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"
    subprocess.run([
        vsearch_bin, "--usearch_global", str(gg2_fa),
        "--db", gg2_db, "--id", "0.0", "--iddef", "1",
        "--maxaccepts", "1", "--maxrejects", "32", "--strand", "both",
        "--userout", str(gg2_hits),
        "--userfields", "query+target+id+alnlen",
        "--threads", "26", "--quiet"
    ], check=True)
    log("  vsearch done")
else:
    log(f"  Using existing: {gg2_hits}")

# Parse hits → lookup genus
log("Mapping GG2 hits to genera via taxonomy.tsv...")
derep_to_genus, derep_to_gg2_name, derep_to_gg2_id = {}, {}, {}

with open(gg2_hits) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 4: continue
        derep_id, target, identity = parts[0], parts[1], float(parts[2])
        if derep_id not in derep_to_genus:  # first hit = best
            derep_to_genus[derep_id] = gg2_tax_map.get(target, "Unknown")
            derep_to_gg2_name[derep_id] = target
            derep_to_gg2_id[derep_id] = identity

n_annotated = len(derep_to_genus)
n_unknown = sum(1 for g in derep_to_genus.values() if g == "Unknown")
log(f"  Annotated: {n_annotated:,}")
log(f"  Unknown genus: {n_unknown:,} ({n_unknown/max(n_annotated,1)*100:.1f}%)")

genus_counts = Counter(derep_to_genus.values())
log(f"  Unique genera: {len(genus_counts):,}")
log("  Top 20 genera:")
for genus, count in genus_counts.most_common(20):
    log(f"    {genus}: {count:,}")

# Add annotations
candidates["gg2_genus"] = candidates["derep_id"].map(derep_to_genus).fillna("Unannotated")
candidates["gg2_name"] = candidates["derep_id"].map(derep_to_gg2_name).fillna("")
candidates["gg2_identity"] = candidates["derep_id"].map(derep_to_gg2_id).fillna(0)

# ============================================================
# Stage 2 flag: genus in contaminant list
# ============================================================
log("\nChecking against contaminant genus list...")
contam_genus_df = pd.read_csv(PIPE_V2 / "contaminant_genera.tsv", sep="\t", comment="#")
contam_genera = set(contam_genus_df["genus"].str.strip().str.lower())
log(f"  {len(contam_genera)} contaminant genera")

candidates["stage2_genus_flag"] = candidates["gg2_genus"].apply(
    lambda g: str(g).strip().lower() in contam_genera
    if pd.notna(g) and g not in ("Unknown", "Unannotated") else False
)
n_stage2 = candidates["stage2_genus_flag"].sum()
log(f"  Stage 2 flagged: {n_stage2:,}")

# Combined
candidates["contamination_flag"] = candidates["stage1_cluster97_flag"] | candidates["stage2_genus_flag"]
n_any = candidates["contamination_flag"].sum()

def flag_type(row):
    if row["stage1_cluster97_flag"] and row["stage2_genus_flag"]: return "both"
    elif row["stage1_cluster97_flag"]: return "cluster97"
    elif row["stage2_genus_flag"]: return "genus_list"
    else: return "clean"

candidates["flag_type"] = candidates.apply(flag_type, axis=1)

# ============================================================
# Summary
# ============================================================
log(f"\n=== Contamination Firewall Summary ===")
log(f"Total candidates: {len(candidates):,}")
log(f"  Stage 1 (97% cluster): {n_stage1:,}")
log(f"  Stage 2 (genus list):  {n_stage2:,}")
log(f"  Any flag: {n_any:,} ({n_any/len(candidates)*100:.1f}%)")
log(f"  Clean: {len(candidates)-n_any:,} ({(len(candidates)-n_any)/len(candidates)*100:.1f}%)")

for ft in candidates["flag_type"].value_counts().index:
    n = (candidates["flag_type"] == ft).sum()
    log(f"    {ft}: {n:,}")

if "n_bioprojects" in candidates.columns:
    candidates["n_bioprojects"] = pd.to_numeric(candidates["n_bioprojects"], errors="coerce").fillna(0)
    log(f"\nBy n_bioprojects:")
    for lo, hi, label in [(0,2,"1 BP"),(2,5,"2-4"),(5,10,"5-9"),(10,20,"10-19"),(20,50,"20-49"),(50,1000,"50+")]:
        s = candidates[(candidates["n_bioprojects"]>=lo)&(candidates["n_bioprojects"]<hi)]
        nf = s["contamination_flag"].sum()
        log(f"  {label}: {len(s):,} total, {nf:,} flagged ({nf/max(len(s),1)*100:.1f}%)")

# ============================================================
# Save
# ============================================================
outpath = PIPE_V2 / "02_contamination_flagged.tsv"
candidates.to_csv(outpath, sep="\t", index=False)
log(f"\nSaved: {outpath} ({os.path.getsize(outpath):,} bytes)")

summary = {
    "total_candidates": int(len(candidates)),
    "stage1_cluster97_flagged": int(n_stage1),
    "stage2_genus_flagged": int(n_stage2),
    "any_flag": int(n_any),
    "clean": int(len(candidates)-n_any),
    "flag_type_distribution": {k: int(v) for k,v in candidates["flag_type"].value_counts().items()},
    "annotated_with_gg2": int(sum(candidates["derep_id"].isin(derep_to_genus))),
    "unique_genera_found": len(genus_counts),
}
with open(PIPE_V2 / "02_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

log(f"Saved summary: {PIPE_V2 / '02_summary.json'}")
log("Done!")
