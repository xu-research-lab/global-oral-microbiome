# HOMD comprehensiveness of the global oral 16S ASV catalogue

This repository contains the analysis behind one main figure and three
supplementary figures. It asks two questions of a global collection of oral
16S samples:

1. how much of the observed ASV catalogue eHOMD represents at species level;
2. which lineages fall below that boundary, whether they form coherent
   sequence clusters, and how those clusters are distributed across projects,
   subjects, geography and disease.

## Where things are

```
pipeline/     the final analysis, ordered by step; start at pipeline/README.md
archive/      superseded analysis schemes, kept for provenance, not runnable
```

`pipeline/README.md` documents the 19 steps, the figure → script → table
chain, the frozen inputs, and how to run each stage. `pipeline/run_order.sh`
lists the commands in order.

## Deliverables

All four PDFs are in `pipeline/19_figures/final/`:

- `fig_v3_master_alt.pdf` — main figure. A data-set flow, B eHOMD coverage,
  C per-genus identity curves, D multi-ASV profile of a supported cluster.
- `fig_v3_supported_singleton_top40_profile.pdf` — the 40 most prevalent
  supported singleton clusters.
- `fig_head_and_neck_cancer_ctrl_only.pdf` — head-and-neck cancer, control
  cohorts only.
- `fig_periodontal_disease_ctrl_only.pdf` — periodontal disease, control
  cohorts only.

## Terminology

- **ASV** — an amplicon sequence variant called by DADA2 within one project.
- **derep sequence** — an exact-dereplicated sequence pooled over all
  projects; `derep_id` (`DEREP_<n>`) names it.
- **eHOMD-captured / uncaptured** — captured means an eHOMD best hit at
  ≥ 97 % identity over ≥ 99 % of the query.
- **cluster** — a TreeCluster `max_clade` unit at t = 0.03 on the reference
  tree, made of placed ASVs and the reference tips inside the same clade.
  The final id is `c<number>` for a cluster of two or more ASVs and
  `s<number>` for a single ASV.
- **supported cluster** — a cluster that passes the prevalence, project,
  geography and cleanliness criteria in step 13.

## Raw data

Raw reads and the per-project DADA2 outputs are not duplicated here; step 01
reads them from their original location on the cluster, which is recorded in
the script headers. The tables needed to rebuild every downstream result are
in each step's `results/` directory.
