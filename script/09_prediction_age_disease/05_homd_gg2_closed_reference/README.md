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
