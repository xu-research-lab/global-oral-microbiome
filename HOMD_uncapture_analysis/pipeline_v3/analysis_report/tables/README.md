# Supplementary tables: phylotypes

Terminology: what the pipeline calls a *cluster* is a **phylotype** in these
tables. The id keeps its form: `c<number>` for a phylotype of two or more ASVs,
`s<number>` for a single ASV that clustered with nothing else. Two single-ASV
phylotypes that contain reference sequences carry a plain number (`7846`,
`32458`): their treecluster id is positive, so there is no sign to replace.

## Files

| File | Rows | What it is |
|---|---|---|
| `TableS_supported_phylotypes.tsv` / `.xlsx` | 486 | the supported catalog, one row per phylotype |
| `TableS_differential_phylotypes_HNC_PD.tsv` | 49 | one row per cohort × phylotype called differential |
| `TableS_differential_phylotypes_HNC_PD_members.tsv` | 200 | one row per member ASV of those phylotypes, with its full 16S sequence |
| `TableS_differential_phylotypes_HNC_PD.xlsx` | — | both differential tables as worksheets |

Each workbook also carries a `column_definitions` worksheet. Builders:
`scripts/19_build_supported_phylotype_table.py` and
`scripts/18_build_disease_feature_table.py`.

## The supported catalog

486 phylotypes passing the support gate — detected in more than 50 subjects,
in at least two bioprojects and on at least two continents, with no
contamination-flagged member and no organelle member. 443 are single ASVs and
43 are multi-ASV (the largest is `c21907`, 112 ASVs). 421 are eHOMD-uncaptured
and 65 eHOMD-captured. Rows keep the source order, support descending.

## The differential tables

Head and neck cancer: 36 phylotypes (3 multi-ASV, 33 single ASV).
Periodontal disease: 13 phylotypes (2 multi-ASV, 11 single ASV).
`c21907` is differential in both cohorts, so it appears twice; 48 distinct
phylotypes in all.

The multi-ASV phylotypes are `c21907` (112 ASVs, *Capnocytophaga*), `c6615`
(40 ASVs, *Selenomonas sputigena*), `c29573` (2 ASVs, *Helicobacter pylori*)
and `c36553` (2 ASVs, *Phnomibacter ginsenosidimutans*).

## Definitions worth keeping in view

- **Distance to the annotation** (`annotation_distance_*`, member table
  `annotation_distance`) is the phylogenetic distance, in substitutions per
  site, from an ASV tip up to the reference tip or ancestor node that carries
  the annotation. For a species-level annotation it is the distance to that
  species; for a genus- or family-level annotation it is the distance to that
  node. In the differential table the representative's value is given with the
  minimum, median and maximum over the member ASVs; for a single ASV all four
  are the same number. In the catalog table only the minimum and maximum are
  available.
- **Taxonomy** is the reference-tree lineage cut at its deepest named rank.
  Sixty-nine catalog phylotypes and seven differential phylotypes
  (`s6635, s6973, s9612, s21076, s23151, s26286, s27943`) stop above the genus
  rank: the pipeline's `annotation` is empty for them, so `assigned_taxon` and
  `taxon_rank` are the deepest rank the lineage actually names (order or
  family), not a genus name.
- **`annotation_source`** says where the annotation comes from: reference
  sequences within the phylotype, the ancestor node of a phylotype that holds
  no reference, or the nearest reference-bearing ancestor of a single ASV.
- **Representative ASV** is the sequence the phylotype was aggregated on
  (`rep_derep_id`); for a single-ASV phylotype it is the ASV itself.
- **`novelty_class`** is eHOMD-captured when a member is annotated as an eHOMD
  species, eHOMD-uncaptured otherwise. `n_members_annotation_in_homd` and
  `genus_in_homd` give the finer counts.
- **`nearest_reference_identity_tier`** bins the distance from the
  representative ASV to its nearest reference sequence: ≥99%, 97–99%, <97%.
  `homd_identity_pct` is the identity to the nearest **eHOMD** reference and is
  a different comparison.
- **Statistics** (differential table) come from the per-cohort MaAsLin3 model
  (disease arm plus the healthy controls of the same bioprojects; covariates
  age, sex, site, region and continent, bioproject as random intercept).
  `qval_joint` is the joint q value across the prevalence and abundance models;
  phylotypes were retained at q < 0.25. `higher_in` says which group has the
  higher value (`Disease` or `HC`) and `coef` is the disease-term coefficient.
- **`support_score`** = `max_n_subjects` × `n_bioprojects` ×
  `n_amplicon_regions`.

## Source tables

- `08_supported_clusters.tsv` — supported catalog
- `14_supported_cluster_annotation.tsv` — taxonomy, phylum, member count, support
- `17_maaslin3_significant.tsv` — tested phylotypes and statistics
- `07_master_asv_table.tsv` — member ASVs, sequences, per-ASV annotation distance and eHOMD match
