# Region-Aware ASV vs HOMD Pipeline

This pipeline evaluates whether observed oral ASVs are captured by HOMD at high
identity while respecting sequencing region, platform, geography, disease/status,
and ASV abundance/prevalence from `ASVs_counts.tsv`.

## Steps

1. Collect retained metadata bioprojects, rename ASVs as `bioproject|ASV_ID`, and
   map `ASVs_counts.tsv` abundance by ASV ID:

```bash
python /beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline/01_collect_asvs_and_abundance.py
```

2. Align exact-dereplicated ASV sequences to HOMD and summarize coverage:

```bash
python /beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline/02_align_homd_and_summarize.py --threads 16
```

## Main Outputs

- `asv_abundance_metadata.tsv`: original ASV-level sequence, abundance, prevalence,
  and project metadata.
- `exact_derep_summary.tsv`: exact sequence representatives with source ASV counts,
  bioproject recurrence, total reads, prevalence, and metadata strata.
- `homd_coverage_by_derep_sequence.tsv`: HOMD best hit and capture flags at 100%,
  99.5%, 99%, and 97% identity with query coverage filtering.
- `homd_coverage_grouped_summary.tsv`: coverage summaries by region, platform,
  continent, country, site, disease-study flag, and status.
- `reproducible_abundant_not_high_identity_homd.tsv`: exact ASV sequences recurring
  across projects and abundant/prevalent but not captured by HOMD at the primary
  high-identity threshold.

Mitochondria/chloroplast/non-bacterial screening should be applied after this step
to the low-HOMD-match candidates, so off-target reads are removed from the final
biological denominator without biasing the collection and dereplication stage.
