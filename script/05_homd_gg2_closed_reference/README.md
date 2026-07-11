# Longest Containment Dereplication

This directory contains a supplementary dereplication workflow for closed-reference
profile generation.

Purpose:

- Start from exact-dereplicated ASVs from the region-aware HOMD pipeline.
- Collapse a shorter sequence into a longer representative only when the shorter
  sequence is fully contained in the longer sequence at 100% identity and 100%
  query coverage.
- Keep mapping tables so abundance can still be traced back to the original
  exact-derep and source ASV records.

Main command:

```bash
sbatch ASV_vs_HOMD_analysis/containment_derep/run_containment_derep_longest.sbatch
```

Main outputs:

- `containment_longest_representatives.fa`: longest representative sequences for
  downstream closed-reference mapping.
- `containment_derep_mapping.tsv`: mapping from each exact-derep sequence to its
  longest containment representative.
- `containment_representative_summary.tsv`: abundance and metadata summarized by
  containment representative.
- `containment_derep_stats.json`: run-level statistics.

Important interpretation:

This is not intended to replace exact-ASV analysis of HOMD representation. It is
intended as a containment-assisted preprocessing step for downstream
closed-reference profile construction.

---

## Geographic reference-bias analyses (scripts 23-24)

Two analyses assess whether the GG2 closed-reference mapping rate reflects a
reference-representation bias against non-Western populations, rather than pure
data quality.

- `23_cr_mapping_by_continent_amplicon_adjusted.py` — **runnable locally.** Plots
  per-sample and per-study GG2 mapping rate by continent (all amplicons and V4-only),
  and fits a study-level weighted regression of `logit(mapping rate) ~ continent +
  amplicon region` so the continent effect is quantified net of amplicon and without
  sample-level pseudoreplication. Consumes the committed QC metadata
  (`data/metadata-and-biom-v2-qc/metadata_all_for_qc_raw.tsv`); writes figures to
  `results/figures/05_homd_gg2_closed_reference_profiling/` and tables to
  `results/tables/05_homd_gg2_closed_reference/`. Key result: after amplicon
  adjustment, Asia and Africa map significantly worse than North America while Europe
  (Western) does not.

- `24_geography_reference_gap_bacterial_asvs.py` — **cluster-only (not yet run).**
  Turns the sample-level signal into a per-taxon, contamination-guarded statement:
  which *bacterial* low-HOMD/GG2 ASVs (already screened to `d__Bacteria` past the
  chloroplast/mito filter by scripts 14-15) drive the reference gap, and whether that
  gap concentrates in non-Western continents. Consumes
  `low_identity_lt97_subject_prevalence_filtered_gg2_bacterial_only.tsv` (cluster
  intermediate) + oral metadata; normalizes burden per 1000 subjects and reports the
  top GG2 genera per continent plus the ASVs seen only in non-Western regions. Logic
  is smoke-tested on synthetic input; adjust `--*-col` flags if the upstream schema
  differs.
