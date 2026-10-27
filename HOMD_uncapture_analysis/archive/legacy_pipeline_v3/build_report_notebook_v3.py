#!/usr/bin/env python3
"""
Build the reproducible analysis-report Jupyter notebook for pipeline_v3
(novel_oral_taxa_report_v3.ipynb): TreeCluster OTU-like clustering + HOMD
coverage annotation + SUPPORTED CATALOG.

Run: python build_report_notebook_v3.py
"""
import nbformat as nbf
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "analysis_report" / "novel_oral_taxa_report_v3.ipynb"
PIPE_V3 = HERE

nb = nbf.v4.new_notebook()
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}


def md(src):
    return nbf.v4.new_markdown_cell(src)


def code(src):
    return nbf.v4.new_code_cell(src)


def fig(name, caption):
    return md(f"![{caption}](figs/{name})")


cells = []

cells.append(md("""# Identification of Oral Taxa Absent from the HOMD Database — v3

**pipeline_v3 · TreeCluster OTU-like clustering scheme** — reproducible notebook ·
SEPP placement (2022.10 backbone) · TreeCluster `max_clade` t=0.03 ·
per-ASV tree-distance HOMD coverage annotation · cluster-level support ·
SUPPORTED CATALOG · 2026-08-25

> Design decisions (2026-08-25):
> 1. Cluster = TreeCluster `max_clade` at t=0.03 on the GG2 2022.10 SEPP
>    insertion tree (backbone reference tips kept as cluster members).
> 2. Per-ASV annotation: nearest reference tip within 0.03 whose species is in
>    HOMD -> covered; otherwise cluster-reference LCA / nearest-reference
>    taxonomy (capped at genus for singletons).
> 3. Cluster novelty: `homd_captured` if any member ASV is covered (or rep HOMD
>    identity >=97%); otherwise `homd_uncaptured`.
> 4. SUPPORTED CATALOG = clean clusters passing subjects > 50 & BioProjects >= 2
>    & continents >= 2. HOMD capture status is shown as a color column in the
>    prevalence profile, which is restricted to species-annotated supported
>    clusters (no rep-tree-distance filter).
> 5. Singleton clusters (size-1 TreeCluster units = unclustered ASVs) are their
>    own clusters and are INCLUDED in the supported catalog.
> 6. Organelle-flagged clusters (GG2 annotation mitochondrion/chloroplast) are
>    excluded from the supported catalog.

## Pipeline overview

| # | Stage | Script |
|---|-------|--------|
| 01 | TreeCluster clustering + HOMD coverage annotation | `01_treecluster_cluster.py` |
| 02 | Cluster aggregation + novelty classification | `02_cluster_aggregate.py` |
| 05 | Phylum distribution of candidate novel clusters | `05_cluster_geography.py` |
| 06 | ALL-candidate collapsed tree | `06_build_cluster_tree.py` |
| 07 | Master per-ASV table | `07_master_table.py` |
| 08 | Supported-cluster catalog | `08_supported_clusters.py` |
| 09 | Summary figures | `09_make_summary_figs.py` |
| 12e | Supported-catalog circular tree | `12e_ggtree_circles_supported_v3.R` |
| 12f | Supported-catalog prevalence profile | `12f_plot_headline_prevalence_profile_v3.R` |
"""))

cells.append(md("## Setup"))
cells.append(code("""import pandas as pd, numpy as np, json
from pathlib import Path
BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V3 = BASE / "pipeline_v3"
FIG = PIPE_V3 / "analysis_report" / "figs"
pd.set_option("display.max_columns", 40)
pd.set_option("display.width", 180)
print("OK")"""))

cells.append(md("## Step 01 — TreeCluster clustering + HOMD coverage annotation"))
cells.append(md("""The 54,277 support-floor ASVs were placed on the GG2 2022.10 tree with SEPP.
TreeCluster `max_clade` at t=0.03 grouped all tips (ASVs + backbone references)
into OTU-like clades. Each ASV was annotated by its nearest reference tip and
classified as `covered` if the nearest HOMD-species reference is within 0.03
(species boundary)."""))
cells.append(code("""# heavy step (exact code that ran, outputs cached):
# python3 pipeline_v3/01_treecluster_cluster.py \\
#   --tree pipeline/05_sepp/tree_export/tree.nwk \\
#   --threshold 0.03 --method max_clade --outdir pipeline_v3

asv = pd.read_csv(PIPE_V3/"asv_clusters.tsv", sep="\\t", low_memory=False)
print(f"ASVs: {len(asv):,}")
print(f"clusters: {asv['cluster_id'].nunique():,}")
print(f"covered: {int(asv['covered'].astype(str).str.lower().isin(['true','1']).sum()):,} "
      f"| novel: {int((~asv['covered'].astype(str).str.lower().isin(['true','1'])).sum()):,}")
print("annotation_level:")
print(asv["annotation_level"].value_counts().to_string())"""))

cells.append(md("## Step 02 — Cluster aggregation + novelty classification"))
cells.append(md("""Per-cluster support statistics were aggregated from the ASV-level metadata.
A cluster is `homd_captured` if any member ASV is covered (or the rep ASV has
>=97% HOMD identity); otherwise `homd_uncaptured`."""))
cells.append(code("""agg = pd.read_csv(PIPE_V3/"02_clusters_aggregated.tsv", sep="\\t", low_memory=False)
print(f"clusters: {len(agg):,}")
print(agg["cluster_novelty_class"].value_counts().rename("novelty").to_string())
print(f"clean_clusters: {int(agg['clean_cluster'].sum()):,}")
print("annotation_level:")
print(agg["annotation_level"].value_counts().to_string())"""))

cells.append(md("## Step 05 — Phylum distribution of candidate novel clusters"))
cells.append(md("""Each candidate novel cluster was mapped to a phylum/class via its member
ASVs' nearest-reference taxonomy (GG2 2022.10). The figure splits the phylum
distribution by whether the annotation genus is in the HOMD genus set (oral
clade)."""))
cells.append(code("""cn = pd.read_csv(PIPE_V3/"05_candidate_novel_with_phylum.tsv", sep="\\t", low_memory=False)
print(f"candidate novel clusters: {len(cn):,}")
print("top phyla:")
print(cn["gg2_phylum"].value_counts().head(15).to_string())
print(f"oral clade (in HOMD genus): {int(cn['is_oral_clade'].sum()):,} | "
      f"outside: {int((~cn['is_oral_clade']).sum()):,}")"""))
cells.append(fig("fig_cluster_phylum.png",
                 "Phylum distribution of candidate novel clusters (split by oral clade)"))

cells.append(md("## Step 06 — ALL-candidate collapsed tree"))
cells.append(md("""Because candidate novel clusters span most of the 385K-tip tree, a local
subtree is impossible. Every clade with no candidate member is collapsed to a
`CLADE_<n>` tip; all candidate rep ASV tips are retained in real phylogenetic
position."""))
cells.append(code("""d2c = pd.read_csv(PIPE_V3/"analysis_report/derep_to_cluster_allcand.tsv", sep="\\t", low_memory=False)
print(f"candidate rep tips: {len(d2c):,}")
print(d2c[["tip","cluster_id","annotation","annotation_level","n_members","subjects",
           "bps","continents"]].head(10).to_string(index=False))"""))

cells.append(md("## Step 07 — Master per-ASV table"))
cells.append(md("""One row per placed ASV (54,277) with sequence/support metadata, TreeCluster
cluster assignment, tree annotation, HOMD status, geography, and contamination
flags."""))
cells.append(code("""master = pd.read_csv(PIPE_V3/"07_master_asv_table.tsv", sep="\\t", low_memory=False)
print(f"rows: {len(master):,} | cols: {len(master.columns):,}")
print(f"with cluster: {int(master['cluster_id'].notna().sum()):,} | covered: {int(master['covered'].astype(str).str.lower().isin(['true','1']).sum()):,}")
print(f"contamination-flagged: {int(master['contamination_flag'].astype(str).str.lower().isin(['true','1']).sum()):,}")
print(f"organelle-flagged: {int(master['organelle_flag'].astype(str).str.lower().isin(['true','1']).sum()):,}")"""))

cells.append(md("## Step 08 — Supported-cluster catalog"))
cells.append(md("""SUPPORTED CATALOG = clean clusters passing subjects > 50 & BioProjects >= 2 &
continents >= 2. HOMD capture status is shown as a color column in the
prevalence profile (no rep-tree-distance / HOMD-uncaptured filter)."""))
cells.append(code("""sup = pd.read_csv(PIPE_V3/"08_supported_clusters.tsv", sep="\\t", low_memory=False)
print(f"supported clusters: {len(sup):,}")
print("tree tier:")
print(sup["tree_tier"].value_counts().to_string())
print("novelty:")
print(sup["cluster_novelty_class"].value_counts().to_string())
print(sup[["cluster_id","annotation","annotation_level",
      "max_subject_prevalence","n_bioprojects_total","n_continents","n_regions",
      "rep_nearest_ref_dist","rep_homd_identity","cluster_novelty_class"]].head(15).to_string(index=False))"""))
cells.append(fig("fig_supported.png",
                 "Supported catalog — tree distance tier x HOMD novelty composition"))
cells.append(fig("fig_tree_circles_supported_v3.png",
                 "Supported catalog — circular cluster-circle tree (size = member ASVs, fill = oral clade; rings: prevalence / BioProjects / Continents / Regions / Phylum)"))
cells.append(fig("fig_v3_headline_prevalence_profile.png",
                 "Supported catalog — prevalence profile (species-annotated supported clusters, 26; HOMD-captured/uncaptured color column; tracks: No. derep ASVs, BioProjects, reads, HOMD id, Continents, Regions; boxplot+jitter of per-ASV subject prevalence)"))

cells.append(md("## Step 09 — Summary figures"))
cells.append(md("""The summary figures show the pipeline funnel, HOMD novelty classes + oral
clade split, and continent coverage."""))
cells.append(code("""import subprocess, sys
r = subprocess.run([sys.executable, str(PIPE_V3/"09_make_summary_figs.py")],
                   capture_output=True, text=True)
print(r.stdout[-1500:])
print("stderr:", r.stderr[-500:] if r.stderr else "none")"""))
cells.append(fig("fig_funnel.png", "Pipeline funnel — ASV/cluster counts at each stage"))
cells.append(fig("fig_novelty.png", "HOMD novelty of TreeCluster clusters + oral clade split"))

cells.append(md("## Data funnel (summary)"))
cells.append(code("""d1  = len(pd.read_csv(BASE/"pipeline_v2/01_chimera_filtered.tsv", sep="\\t", low_memory=False, usecols=["derep_id"]))
d4  = len(pd.read_csv(BASE/"pipeline_v2/04_support_floor_passed.tsv", sep="\\t", low_memory=False, usecols=["derep_id"]))
s01 = json.load(open(PIPE_V3/"summary.json"))
s02 = json.load(open(PIPE_V3/"02_summary.json"))
s08 = json.load(open(PIPE_V3/"08_supported_summary.json"))
funnel = pd.DataFrame({
    "Stage": ["ASVs screened (334 BioProjects)",
              "HOMD-uncaptured & non-chimeric",
              "Support floor (>=2 BP, >=50 reads)",
              "ASVs placed on GG2 tree",
              "TreeCluster clusters (t=0.03)",
              "Supported catalog (clean)"],
    "Count": [3396776, d1, d4, s01["n_asv"], s02["n_clusters"], s08["n_supported_clean"]],
})
print(funnel.to_string(index=False))"""))

cells.append(md("## Limitations"))
cells.append(md("""1. **TreeCluster t=0.03 is a tree-distance threshold**, not a sequence-identity
   threshold; branch lengths on the SEPP placement tree are approximate.
2. **Cluster annotation is reference-anchored**: ASVs far from any GG2 reference
   are capped at genus or higher, so many clusters have no species-level name.
3. **99%/97% identity is computed on short amplicon fragments**; species calls
   from a single region are provisional.
4. **No oral-vs-environment adjudication**: some supported clusters have very low
   HOMD identity and could be environmental/transient; retained by user decision.

---
*Report v3 built 2026-08-29 (singleton ASV clusters included; organelle-flagged clusters excluded).*"""))

nb["cells"] = cells
OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as f:
    nbf.write(nb, f)
print("wrote", OUT)
