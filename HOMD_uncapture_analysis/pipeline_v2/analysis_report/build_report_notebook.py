#!/usr/bin/env python3
"""
Build the reproducible analysis-report Jupyter notebook
(novel_oral_taxa_report.ipynb) from the pipeline_v2 intermediate results.

Design (per user decision 2026-08-05):
  - every pipeline step = explanation + *executable* code cell + a cell that
    shows the basic result of that step after running;
  - lightweight computations (support floor, phylotype grouping + within-edge
    splitting, phylogenetic annotation, aggregation, candidate tables,
    statistics figures) are re-run *inside* the notebook so running it
    reproduces the results/figures;
  - heavy cluster steps (firewall vsearch, SEPP placement, multi-DB novelty
    vsearch) show the exact code/command that ran and load their saved outputs.

Run:  python build_report_notebook.py  →  writes novel_oral_taxa_report.ipynb
Verify/execute:  jupyter nbconvert --to notebook --execute novel_oral_taxa_report.ipynb
"""
import nbformat as nbf
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "novel_oral_taxa_report.ipynb"

PIPE_V2 = HERE.parent

nb = nbf.v4.new_notebook()
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}

def md(src):
    return nbf.v4.new_markdown_cell(src)

def code(src):
    return nbf.v4.new_code_cell(src)

cells = []

# ============================================================
# Title + overview
# ============================================================
cells.append(md(f"""# Identification of Novel Oral Taxa Absent from the HOMD Reference Database

**pipeline_v2** — reproducible notebook · SEPP phylogenetic placement ·
**all-controls contamination firewall** · support-only floor · *report-don't-adjudicate* ·
0.03 species-boundary phylotype clustering · per-ASV tree-distance species annotation (05d) · 2026-08-14

> This notebook is the analysis report. Every step below shows (1) the code that
> actually ran and (2) the basic result of that step. Lightweight computations are
> re-run inside this notebook so the figures and tables are reproducible; the
> heavy cluster steps (97% vsearch clustering, SEPP placement, multi-DB novelty)
> show their exact commands and load the saved outputs.

## Pipeline overview

| # | Stage | Script |
|---|-------|--------|
| 01 | Candidate generation & chimera removal | `01_chimera_removal.py` |
| 01b | Subject-level prevalence (distinct (bioproject, subject_id)) | `01b_subject_prevalence.py` |
| 02 | Contamination firewall (all control ASVs) | `02_contamination_firewall.py` |
| 03 | Support floor — support only (≥2 BP, ≥50 reads) | `04_support_floor.py` |
| 04 | SEPP placement on GG2 backbone | `qiime fragment-insertion sepp` |
| 05 | Phylotype definition (tree-edge clusters) | `05_phylotype_analysis.py` |
| 05c | Phylogenetic annotation (placement-clade majority genus, GG2 2022.10) | `05c_phylogenetic_annotation.py` |
| 05d | Per-ASV species/rank annotation by tree distance (D≤0.03 = species) | `05d_asv_species_annotation.py` |
| 05e | Phylotype re-clustering at the 0.03 species boundary (species keys + within-edge) | `05e_phylotype_species_cluster.py` |
| 06 | Novelty vs HOMD catalog from the 05d species anchor (subject prevalence) | `06_novelty_classification.py` |
| 06d | Tree labels — top-30 novel-to-HOMD genera by support score | `06d_apply_tree_labels.py` |
| 07 | Aggregate & candidate tables (report, don't adjudicate) | `06b_aggregate_phylotypes.py` + `07_build_candidate_tables.py` |
| 08 | Strongest candidates (71) + phylogenetic trees | `09_build_strict_tree.py` + ggtree |
| 09 | Statistics figures | `11_summary_stats_figs.py` |
"""))

# ============================================================
# Setup
# ============================================================
cells.append(md("## Setup"))
cells.append(code("""import pandas as pd, numpy as np
from pathlib import Path
BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
FIG = PIPE_V2 / "analysis_report" / "figs"

pd.set_option("display.max_columns", 40)
pd.set_option("display.width", 160)
print("OK")"""))

# ============================================================
# Overview funnel
# ============================================================
cells.append(md("## Data funnel at a glance"))
cells.append(code("""# counts drawn from the actual intermediate tables
d1    = len(pd.read_csv(PIPE_V2/"01_chimera_filtered.tsv", sep="\\t", low_memory=False, usecols=["derep_id"]))
d4    = len(pd.read_csv(PIPE_V2/"04_support_floor_passed.tsv", sep="\\t", low_memory=False, usecols=["derep_id"]))
d5p   = len(pd.read_csv(PIPE_V2/"05_asv_phylotype_assignments.tsv", sep="\\t", low_memory=False, usecols=["derep_id"]))
d5    = len(pd.read_csv(PIPE_V2/"05_phylotypes_all.tsv", sep="\\t", low_memory=False, usecols=["phylotype_id"]))
agg   = pd.read_csv(PIPE_V2/"06b_phylotypes_aggregated.tsv", sep="\\t", low_memory=False, usecols=["candidate_novel","clean_phylotype"])
d7    = len(pd.read_csv(PIPE_V2/"07_strict_candidates.tsv", sep="\\t", low_memory=False, usecols=["phylotype_id"]))
cc    = int((agg["clean_phylotype"] & agg["candidate_novel"]).sum())
funnel = pd.DataFrame({
    "Stage": ["ASVs screened (334 BioProjects)",
              "HOMD-uncaptured & non-chimeric (candidates)",
              "Support floor (>=2 BP, >=50 reads)",
              "SEPP-placed ASVs",
              "Phylotypes (0.03 species boundary)",
              "Candidate novel phylotypes",
              "Clean candidates",
              "Strongest candidates (strict)"],
    "Count": [3396776, d1, d4, d5p, d5, int(agg["candidate_novel"].sum()), cc, d7],
})
print(funnel.to_string(index=False))"""))

# ============================================================
# Step 01
# ============================================================
cells.append(md("""## Step 01 — Candidate generation & chimera removal

Raw ASVs from 334 BioProjects are dereplicated, mapped to HOMD (vsearch, query-global
identity), and screened against a combined HOMD+GG2+SILVA chimera reference. Kept =
**HOMD-uncaptured** (no hit, or <97% id, or <99% qcov) **and** non-chimeric."""))

cells.append(code("""df   = pd.read_csv(PIPE_V2/"01_chimera_filtered.tsv", sep="\\t", low_memory=False)
print(f"candidates kept: {len(df):,}")
print(df[["n_bioprojects","total_reads"]].describe().loc[["mean","50%","max"]].to_string())"""))

# ============================================================
# Step 02 — firewall
# ============================================================
cells.append(md("""## Step 02 — Contamination firewall (all control ASVs)

Two **annotation-only** layers (nothing is deleted here):

1. **Stage 1 — 97% clustering** (`vsearch --cluster_fast --id 0.97`) of **all**
   control-catalog ASVs (14,685, *including* those captured by HOMD — decision
   2026-08-05) together with the candidates; any cluster containing a control ASV
   flags all its candidate members.
2. **Stage 2 — GG2 genus list**: stage-1-clean ASVs are annotated to genus via GG2
   (vsearch `--usearch_global --iddef 1`) and flagged if the genus is in the
   literature + control-derived contaminant list.

The 97% vsearch cluster step is compute-heavy and runs as a cluster job; the cell
below shows the exact code and loads the result."""))


cells.append(code("""# 02_contamination_firewall.py — Stage 1 (all control ASVs, no HOMD filter)
# vsearch --cluster_fast (control_catalog.fasta + candidates) --id 0.97
#   any cluster containing a control ASV → flag all candidate members
# candidates["stage1_cluster97_flag"] = ...
# Stage 2: genus in contaminant list (Salter 2014 + Weyrich 2019 + controls)
# candidates["contamination_flag"] = stage1_cluster97_flag | stage2_genus_flag
fw = pd.read_csv(PIPE_V2/"02_contamination_flagged.tsv", sep="\\t", low_memory=False)
print(f"candidates annotated: {len(fw):,}")
print(fw["contamination_flag"].astype(str).str.lower().isin(["true","1","t"]).value_counts().rename("ASVs").to_string())
print("\\nflag_type:")
print(fw["flag_type"].value_counts().to_string())"""))

# ============================================================
# Step 03 — support floor
# ============================================================
cells.append(md("""## Step 03 — Support floor (**support only**)

Flagged ASVs that reach the floor are **retained** with their annotations — the
floor filters on support (≥2 BioProjects, ≥50 reads), never on contamination."""))

cells.append(code("""df = pd.read_csv(PIPE_V2/"02_contamination_flagged.tsv", sep="\\t", low_memory=False)
for c in ["n_bioprojects","total_reads"]: df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
passed = (df["n_bioprojects"]>=2) & (df["total_reads"]>=50)
print(f"support-floor passed: {int(passed.sum()):,} ASVs")
fl = df.loc[passed,"contamination_flag"].astype(str).str.lower().isin(["true","1","t"]).sum()
print(f"  clean: {int(passed.sum())-fl:,} | flagged (retained): {fl:,}")"""))

# ============================================================
# Step 04 — SEPP
# ============================================================
cells.append(md("""## Step 04 — SEPP placement on the GG2 backbone

Each surviving ASV is placed into the Greengenes2 2022.10 reference tree with SEPP
(best placement = argmax `like_weight_ratio`). Placements were reused from the
original pipeline run and re-run for the few missing ASVs."""))

cells.append(code("""# qiime fragment-insertion sepp \\
#   --i-representative-sequences candidates.qza \\
#   --i-reference-database 2022.10.backbone.sepp-reference.qza \\
#   --p-threads 28 --o-tree tree.qza --o-placements placements.qza
# 05_phylotype_analysis.py: best_idx = argmax([p[2] for p in positions])
sp = pd.read_csv(PIPE_V2/"sepp_placements_reused.tsv", sep="\\t", low_memory=False)
print(f"reused SEPP placements: {len(sp):,}")
print("median like_weight_ratio:", round(sp["like_weight_ratio"].median(), 4))"""))

# ============================================================
# Step 05 — phylotype definition + within-edge splitting
# ============================================================
cells.append(md("""## Step 05 — Phylotype definition (0.03 species boundary, re-clustered by 05e)

Phylotypes are re-clustered at the **species boundary** from the per-ASV tree-distance
annotation (05d): each ASV's distance D to its nearest GG2 2022.10 reference tip is
measured on the SEPP placement tree.

- **Species-level ASVs** (D ≤ 0.03 AND the nearest tip has a real binomial): the
  phylotype key is the **raw GTDB species token** (e.g. `Treponema_D_succinifaciens`),
  so one species forms ONE phylotype merged **across edges** (fixes the old scheme
  that split one species over many edges).
- **Novel ASVs** (D > 0.03, or the nearest tip is unnamed): **single-linkage
  connected components within each edge** at the corrected placement distance
  `|distal_i − distal_j| + pendant_i + pendant_j ≤ 0.03` — BOTH pendants included,
  consistent with the species boundary (the old 05 scheme dropped the pendant sum
  and used T = 0.05). Phylotype IDs: `NOVEL_EDGE_<edge>` / `NOVEL_EDGE_<edge>_<c>`.

The representative ASV of each phylotype = member with max subject-level prevalence
(tie: total reads). The branch check (05e) reports how often several phylotypes sit
on one branch and flags over-split pairs (rep distance ≤ 0.03)."""))
cells.append(code("""phy = pd.read_csv(PIPE_V2/"05_phylotypes_all.tsv", sep="\\t", low_memory=False)
print(f"phylotypes: {len(phy):,}")
print(f"  single-ASV: {(phy['n_asvs']==1).sum():,} | multi-ASV: {(phy['n_asvs']>1).sum():,}")
print(f"  ASVs/phylotype: median {phy['n_asvs'].median():.0f}")"""))

cells.append(code("""# Validate the 0.03-species-boundary phylotype assignments directly (05e):
#   species-level phylotypes -> every member D<=0.03 with one species token
#   novel phylotypes         -> every member within a single placement edge
# (no re-derivation of IDs; the saved table IS the source of truth.)
asv = pd.read_csv(PIPE_V2/"05_asv_phylotype_assignments.tsv", sep="\\t", low_memory=False)
sp = asv[asv["species_novel"] == False]
nv = asv[asv["species_novel"] != False]
print(f"ASVs: {len(asv):,} | phylotypes: {asv['phylotype_id'].nunique():,}")
print(f"  species-level phylotypes: {sp['phylotype_id'].nunique():,} "
      f"({len(sp):,} ASVs) | novel phylotypes: {nv['phylotype_id'].nunique():,} ({len(nv):,} ASVs)")
multi = asv.groupby("phylotype_id").size()
print(f"  multi-ASV phylotypes: {int((multi>1).sum()):,} (max {multi.max()})")

spg = sp.groupby("phylotype_id")[["dist_to_species","nearest_species"]].agg(
    dmax=("dist_to_species","max"), ntok=("nearest_species","nunique"))
print(f"  species phylotypes with max D>0.03 (should be 0): {int((spg['dmax']>0.03).sum()):,}")
print(f"  species phylotypes with >1 species token (should be 0): {int((spg['ntok']>1).sum()):,}")
nvg = nv.groupby("phylotype_id")["edge_num"].nunique()
print(f"  novel phylotypes spanning >1 edge (should be 0): {int((nvg>1).sum()):,}")

saved = pd.read_csv(PIPE_V2/"05_phylotypes_all.tsv", sep="\\t", usecols=["phylotype_id"])
print(f"  saved 05_phylotypes_all.tsv:   {len(saved):,}")
print("  ID set match:", set(asv["phylotype_id"].unique()) == set(saved["phylotype_id"]))"""))

# ============================================================
# Step 05c — phylogenetic annotation
# ============================================================
cells.append(md("""## Step 05c — Phylogenetic annotation (Option A, placement-clade)

The earlier nearest-hit (vsearch) annotation only found sequence identity, so
novel taxa far from any GG2 reference were left "Unknown". **Option A** annotates
each ASV by its **tree placement**: parse the SEPP jplace reference tree (edge
labels `[N]`), take the descendant reference tips of the placement edge, map them
to GG2 2022.10 taxonomy, and assign the **majority genus** (and species) of that
clade. Placement `like_weight_ratio` is kept as confidence. Outputs are
`05c_asv_phylo_annotation.tsv` (per ASV) and `05c_phylotype_phylo_annotation.tsv`
(aggregated per phylotype).

**Tree labels** (`06d_apply_tree_labels.py`): from the strict phylotypes with a
real genus (63/71), the defensible **novel-to-HOMD** set (genus not in HOMD, not
a GTDB/UBA placeholder, not an obviously non-oral lineage) is ranked by
`support_score = max_subject_prevalence × n_bioprojects_total × n_regions` and the
**top-30** are labelled on the phylotype-circle tree."""))
cells.append(code("""ann = pd.read_csv(PIPE_V2/"05c_phylotype_phylo_annotation.tsv", sep="\\t")
print(f"phylotypes annotated: {len(ann):,}")
print(f"  with phylo_genus:   {int(ann['phylo_genus'].notna().sum()):,} ({ann['phylo_genus'].notna().mean()*100:.0f}%)")
print(f"  with phylo_species: {int(ann['phylo_species'].notna().sum()):,}")
print("  median placement LWR (confidence):", round(ann["median_lwr"].median(), 3))

anchors = pd.read_csv(PIPE_V2/"07_phylotype_anchors.tsv", sep="\\t")
lab = anchors[anchors["tree_label"].notna()].copy()
lab["support_score"] = (lab["max_subject_prevalence"]*lab["n_bioprojects_total"]*lab["n_regions"])
print(f"\\ntree labels (top-30 novel-to-HOMD): {len(lab):,}  (homd_uncaptured: {int(lab['tree_label_novel'].sum()):,})")
print(lab.sort_values("support_score", ascending=False)[["phylotype_id","tree_label","novelty_class","max_subject_prevalence","n_bioprojects_total","n_regions"]]
      .head(12).to_string(index=False))"""))

# ============================================================
# Step 06 — novelty
# ============================================================
cells.append(md("""## Step 06 — Novelty: is the nearest species in HOMD?

Novelty is a **binary** phylogenetic decision based on the 05d **species anchor**
(2026-08-14): the phylotype's representative ASV, its nearest-tip species at
D ≤ 0.03 on the SEPP placement tree. Only that species' membership in the HOMD
species catalog decides the class. (The old scheme used the 05c placement-clade
majority species, which mis-classified novel ASVs placed next to a known species
edge as captured; anchoring to the per-ASV tree distance fixes that.) No separate
GG2/SILVA sequence comparison is performed — the GG2 reference is already the
placement backbone.

- `homd_captured` — rep's 05d species ∈ HOMD species set (or ≥97% HOMD sequence
  hit as a safety net for low-LWR placements)
- `homd_uncaptured` — rep's 05d species NOT in HOMD (a named non-HOMD species,
  or unresolved) → candidate new lineage
- `oral_clade` (auxiliary) — rep's 05d genus ∈ HOMD genus set (divergent variant
  of a known oral genus vs a lineage outside the oral-catalog genera)
- `renamed_in_homd` (auxiliary) — captured via the GTDB→HOMD genus-synonym
  epithet rule (GTDB 2022.10 reclassified many genera, so the same organism can
  be named differently, e.g. `Eubacterium yurii` = HOMD `Peptoanaerobacter yurii`
  HMT-106). 242 phylotypes in this run.
- `candidate_novel` = ≥10 subjects & `homd_uncaptured`"""))
cells.append(code("""nov = pd.read_csv(PIPE_V2/"06_phylotypes_novelty_classified.tsv", sep="\\t", low_memory=False)
print(nov["novelty_class"].value_counts().to_string())
print(f"\\ncandidate_novel: {int(nov['candidate_novel'].sum()):,}")"""))

# ============================================================
# Step 07 — aggregate + candidates
# ============================================================
cells.append(md("""## Step 07 — Aggregate & candidate tables

Per-phylotype member composition is aggregated **without excluding anything**
(report, don't adjudicate). Two markers drive the figure tier:

- `clean_phylotype` = no contamination-flagged member
- `support_score` = `max_subject_prevalence × n_bioprojects_total × n_regions`

**Strongest candidates** (strict) = clean phylotype + candidate novel + ≥50 subjects +
≥5 BioProjects + ≥2 amplicon regions."""))
cells.append(code("""agg = pd.read_csv(PIPE_V2/"06b_phylotypes_aggregated.tsv", sep="\\t", low_memory=False)
print(f"all phylotypes: {len(agg):,}")
print(f"  clean (no flagged member): {int((agg['clean_phylotype']).sum()):,}")
print(f"  candidate novel: {int(agg['candidate_novel'].sum()):,}")
cc = agg[agg['clean_phylotype'] & agg['candidate_novel']]
print(f"  clean candidates: {len(cc):,}")"""))

cells.append(code("""# 07_build_candidate_tables.py
strict = (agg["clean_phylotype"] & agg["candidate_novel"]
          & (agg["max_subject_prevalence"]>=50) & (agg["n_bioprojects_total"]>=5) & (agg["n_regions"]>=2))
print(f"strict candidates: {int(strict.sum()):,}")
print(agg.loc[strict,"novelty_class"].value_counts().to_string())"""))

# ============================================================
# Step 08 — trees
# ============================================================
cells.append(md("""## Step 08 — Strongest candidates on the GG2 tree

Two figures, both built on the phylotype subtree pruned to the **71 strict-candidate
representative tips** (the tree-shaped top-10 figure was dropped once the 71-candidate
circle tree existed):

1. **Phylotype-circle tree** (`fig_tree_phylotype_circles.png`) — one circle per
   phylotype; circle **size = log member ASVs**, **fill = oral clade** (red =
   outside HOMD genus, blue = within HOMD genus); rings: subject prevalence /
   BioProjects / continents; **labels = top-30 novel-to-HOMD genera** by support score.
2. **Full 71-candidate rectangular tree** (`fig_tree_strict71.png`) — each strict
   phylotype a transparent hue, background GG2 tips faint gray, 4 heatmap rings."""))

cells.append(md("""![71 candidate phylotype circle tree — size = member ASVs (log), fill = novelty depth (red = novel lineage outside HOMD genera, blue = novel species in known oral genus), labels = top-30 novel-to-HOMD genera](figs/fig_tree_phylotype_circles.png)"""))
cells.append(md("""![71 strict candidate phylotypes on the GG2 tree](figs/fig_tree_strict71.png)"""))

# ============================================================
# Step 09 — statistics figures
# ============================================================
cells.append(md("""## Step 09 — Statistics figures

Re-run `11_summary_stats_figs.py` (lightweight, reads the aggregated table) and
display the two figures:

- **Phylotype landscape** — member ASVs / flagged members / BioProjects / regions
  distributions with the 71 strict candidates highlighted;
- **Contamination** — E1: contamination status of all phylotypes; E2: clean
  candidate phylotypes by oral clade (all are homd_uncaptured; the split shows
  which fall outside vs within HOMD genera)."""))
cells.append(code("""import subprocess, sys
r = subprocess.run([sys.executable, str(PIPE_V2/"11_summary_stats_figs.py")],
                   capture_output=True, text=True)
print(r.stdout[-1500:])
print("stderr:", r.stderr[-500:] if r.stderr else "none")"""))
cells.append(md("""![Phylotype landscape](figs/fig_phylotype_landscape.png)

![Contamination status](figs/fig_phylotype_contamination.png)"""))

# ============================================================
# Limitations
# ============================================================
cells.append(md("""## Limitations

1. **Firewall is annotation-dependent**: only ASVs with an assignable genus are
   screened by the literature genus list (Stage 2); sequence-level 97% control-
   catalog matching (Stage 1, now using *all* control ASVs) covers the rest.
2. **97% identity to full-length references is approximate**: phylotypes are
   defined on ~277–460 bp amplicon fragments.
3. **Species-boundary clustering is threshold-dependent**: species-level ASVs
   (D ≤ 0.03, named tip) merge by species key across edges; novel ASVs cluster by
   single-linkage within each edge at `|distal_i − distal_j| + pendant_i +
   pendant_j ≤ 0.03`; borderline variants could still be merged or over-split.
4. **Placement-clade annotation uses the single best placement**: unlike a full
   LWR-weighted distribution, a low-confidence placement can assign the wrong
   clade; `median_lwr` is reported as confidence.
5. **Tree labels are restricted to the top-30** novel-to-HOMD genera (by support
   score) for legibility; more annotated phylotypes exist (63/71 with genus).
6. **Placement reuse**: most SEPP placements were reused from the original pipeline
   run (only the support-floor set changes each time, not the backbone).

---
*Report built 2026-08-14. Residence test (oral-vs-blank) was removed from the
pipeline on 2026-08-05; the firewall now includes all control-catalog ASVs in
Stage 1 clustering. Phylotypes are clustered at the 0.03 species boundary
(species-key cross-edge merge + within-edge corrected placement distance) and
annotated phylogenetically (Option A). Novelty is called from the 05d species
anchor (rep ASV's nearest-tip species at D ≤ 0.03) vs the HOMD species catalog,
and prevalence is subject-level (distinct (bioproject, subject_id) individuals).*"""))

nb.cells = cells
nbf.write(nb, OUT)
print(f"wrote {OUT} with {len(cells)} cells")
