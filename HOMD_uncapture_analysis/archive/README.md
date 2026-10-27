# Superseded analysis schemes

These directories are kept for provenance. They are **not runnable**: the
final analysis lives in `../pipeline/`, and the scripts here still contain the
paths they had when they were written.

| Directory | What it is |
|---|---|
| `legacy_pipeline/` | first pipeline: candidate integration, SEPP placement, phylotype trees, report figures |
| `legacy_pipeline_v2/` | second pipeline: chimera filter, contamination firewall, support floor, SEPP, phylotype analysis and its figures |
| `legacy_pipeline_v3/` | third pipeline: TreeCluster clustering, eHOMD coverage variants, supported clusters, disease models, and the report that produced the final figures |
| `legacy_root_scripts/` | scripts that lived at the repository root before the restructure: the per-project DADA2 entry point, eHOMD mapping, chimera removal, the control-catalog chain (`C01`–`C09`, including the `C08b` sensitivity that was not used), and assorted decontamination and tree-plotting scripts |
| `legacy_figures/` | earlier figure set, before the four deliverables were fixed |
| `legacy_cluster_id_backup/` | tables as they stood before cluster ids were renamed |

Notes for anyone reading the history:

- The three pipeline generations overlap in content. The same logical step may
  appear in more than one of them, at different filenames and different points
  in the order.
- `legacy_pipeline_v3/analysis_report/` holds the figure scripts and the
  delivered PDFs as they were at submission time.
- The control catalogue used downstream came from `C08` (in
  `legacy_root_scripts/`), not from the later `C08b` sensitivity check.

The mapping from every script in `../pipeline/` to its original filename and
directory is in `../pipeline/unrenamed_backup/MAPPING.tsv`.

Reference FASTA files and vsearch scratch directories inside these trees are
on disk but excluded by `../.gitignore` because they are too large for the
repository; the code and the small result tables are tracked.
