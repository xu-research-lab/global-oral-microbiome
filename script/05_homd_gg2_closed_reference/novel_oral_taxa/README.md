# Novel-oral-taxa (not-in-HOMD) pipeline

Implements the plan in `results/novel_oral_taxa_pipeline_plan.md`: identify candidate oral
taxa absent from HOMD, with a contamination firewall and phylotype-level (not per-ASV)
reporting. Fixes the failure mode of the legacy HOMD-uncaptured funnel, whose recurrent "novel
oral" set was dominated by reagent/marine/gut contaminants (recurrence across BioProjects
is also the signature of contamination).

Upstream QC-gate changes are a **separate workstream** and are not part of this pipeline.

## Steps

| # | Script | Phase | Runs where |
|---|---|---|---|
| 01 | `01_export_homd_uncaptured_candidates.py` | 1 candidate generation | local/cluster |
| 02 | `02_chimera_screen.py` | 1 chimera removal | cluster (vsearch) |
| 03 | `03_decontam_within_blanks.R` | 2 contamination firewall | cluster (R/decontam) |
| 04 | `04_oral_vs_blank_enrichment.py` | 2/3 residence test | local/cluster |
| 05 | `05_sepp_place.sh` | 4 phylotype placement | cluster (QIIME2/SEPP) |
| 06 | `06_define_phylotypes.py` | 4 phylotype definition | local/cluster |
| 07 | `07_multidb_novelty.py` | 4 novelty verdict | cluster (vsearch) + local |

Shared paths/thresholds are in `config.py` (cluster defaults; override via flags).

## Key design decisions (see the plan doc for rationale)

- **Keep the rare tail.** Candidate generation uses a relaxed support floor
  (>=2 BioProjects, >=50 reads) and *labels* support tiers rather than deleting rare taxa
  at subject-prevalence >= 50 — the rare/population-specific taxa are the target.
- **Global identity.** HOMD-uncaptured = no hit OR identity < 97% OR qcov < 99%, with identity
  computed query-globally; the qcov clause is then largely redundant. Trim primers before
  mapping. Confirm vsearch `--iddef` is query-global.
- **Chimera removal precedes novelty.** A low-coverage ASV is ambiguous (chimera vs
  divergent); resolve with reference-based `--uchime_ref` before promoting to a candidate.
- **Contamination via blanks, residence via blanks — not gut.** Oral-gut co-abundant taxa
  are real biology, so oral-vs-gut is *not* a residence test; use oral-vs-blank enrichment
  plus within-oral ecology. 23 BioProjects carry blanks (see `25_negative_control_audit.py`).
- **Phylotypes via GG2 SEPP placement.** GG2 is SEPP-native and unifies amplicon regions on
  one tree; HOMD is overlaid as labels. Novelty verdict cross-checks HOMD + GG2 + SILVA/NCBI.

## Status

All seven steps drafted. 01, 04, 06, 07 have real logic and are smoke-tested on synthetic
input. 02 (vsearch) and 05 (QIIME2/SEPP) are tool wrappers, syntax/dry-run verified. 03
(R/decontam) is syntax-checked; needs R + decontam on the cluster. None have been run on
the real cluster data yet — thresholds (support, residence, novelty identity) are defaults
to tune per data.
