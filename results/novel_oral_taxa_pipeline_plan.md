# Pipeline plan: identifying oral microbes not represented in HOMD

Draft version: 2026-07-12

A planning document for the analysis that identifies candidate **novel oral taxa absent
from HOMD**, and connects them to the geographic reference-bias result. It records the
design decisions reached during review, the rationale for each, and the concrete steps.

Scope note: upstream QC-gate changes (replacing the ASV-count filter with depth /
coverage / chimera / closed-reference gates) are a **separate workstream** and are not
part of this pipeline.

## Motivation

Expanding the oral reference beyond HOMD is an explicit project aim. The existing
low-HOMD funnel (2.90M dereplicated ASVs -> 1.90M low-HOMD -> 26k -> 1,416 -> 736
bacterial) is not fit for that claim as built: its most prevalent "novel oral" ASVs are
dominated by contaminants (marine *Vreelandella*, kit *Methylobacterium*/*Shewanella*,
gut *Ruminococcus*/*Coprococcus*/*Enterocloster*, dietary *Cyanobacteriota*). The root
cause is structural — recurrence across many BioProjects is the study's confidence signal
**and** the signature of reagent contamination, so the filter selects for cosmopolitan
contaminants while the aggressive prevalence cut discards the rare, population-specific
novel taxa that are the actual discovery. This plan fixes that.

## Definitions

**Low-HOMD (candidate not-in-HOMD) ASV** = no HOMD hit, OR HOMD identity < 97%, OR
query coverage (qcov) < 99%, evaluated **per amplicon region**.

- Identity is **global** (computed over the full query), so identity <= qcov and an
  identity >= 97% already implies qcov >= 97%. Under this definition the `qcov < 99%`
  clause is largely redundant with the identity clause and mainly trims high-identity
  matches that are ~97-99% covered; global identity — not qcov — is what guards against
  partial/chimeric matches. (Confirm the vsearch `--iddef` used is query-global, e.g.
  `--iddef 1`; the default `--iddef 2` ignores terminal gaps, in which case qcov still
  does independent work at the sequence ends.)
- **Primers must be trimmed before mapping**, otherwise primer remnants deflate qcov and
  identity and inflate false novelty.

**Phylotype** = a candidate taxon defined by phylogenetic placement (below), the reporting
unit of the analysis. Raw ASVs are *not* the unit: thousands of error / hypervariable /
cross-region variants collapse to one taxon.

## Phase 1 — Candidate generation

1. Collect dereplicated ASVs with per-ASV read and subject support, **per amplicon region**.
2. Map to HOMD (vsearch `--usearch_global`, query-global identity). Flag low-HOMD ASVs by
   the definition above. **Keep the rare tail** — do not pre-filter at subject-prevalence
   >= 50; rare and population-specific taxa are the target, not noise to be removed.
3. **Reference-based chimera removal** against a combined HOMD + GG2 + SILVA database
   (DADA2 removes only within-run chimeras; cross-study/cross-taxon chimeras remain).
   This must precede any promotion of a low-coverage ASV to a novel candidate, because a
   low-qcov ASV is ambiguous — chimera vs genuinely divergent — and only chimera screening
   resolves it.

## Phase 2 — Contamination firewall (before trusting recurrence)

The audit (`25_negative_control_audit.py`) found **23 BioProjects with candidate blanks/
NTC/mock controls (17 with explicit extraction-blank/NTC)**, skewed toward the UCSD/
Knight-lab consortium. Decontamination is therefore a hybrid, not confined to those 23 and
not a naive global blacklist:

4. Run `decontam` (prevalence mode) **within** each of the 23 blank-bearing studies —
   gold standard where blanks exist.
5. **Learn a cross-lab contaminant catalog** from those blanks (taxa recurrent in blanks
   across multiple independent labs = high-confidence cosmopolitan contaminants). This
   extends the curated literature blacklist (`contaminant_reference_genera.tsv`;
   Salter 2014 + Eisenhofer 2019 + oral-core protection).
6. **Broadcast the catalog to all studies as a flag, not an automatic delete**
   (`26_flag_contaminants_low_homd.py`), because contaminant profiles are kit/lab-specific
   and the 23 are consortium-biased. Supplement with blank-free contamination signals:
   within-study abundance ceiling (contaminants stay trace) and the
   high-prevalence/low-abundance-everywhere pattern; use frequency-based decontam where any
   DNA-quantitation metadata exists.
7. **Batch-stratify the discovery itself**: stop trusting an ASV because it recurs across
   many studies (contaminant signature). Require real **within-study abundance**, and where
   a study has a blank, enrichment over that study's own blank.
8. Record per-ASV which evidence flagged it (in-blank, blacklist, low-abundance-cosmopolitan)
   so removals are auditable; state the residual limitation for blank-free studies.

## Phase 3 — Oral-residence test

Do **not** use oral-vs-gut enrichment. Genuine oral taxa can be abundant in the gut
(*Fusobacterium nucleatum*, oral streptococci, *Veillonella*, *Klebsiella* translocate and
colonize the gut, especially in IBD/CRC), so an "oral > gut" filter would wrongly reject
real oral microbes; oral-gut co-abundance is if anything evidence *for* being a genuine
host-associated microbe.

9. Test residence against **blanks, not gut**: an oral resident is far more abundant/
   prevalent in oral samples than in blanks; a contaminant is similar in both. This targets
   contamination without penalizing oral-gut shared taxa.
10. Add within-oral **ecological structure** as positive evidence: reaches high relative
    abundance in a subset of hosts, shows site-specificity (tongue / tooth / plaque), and
    co-occurs with core oral taxa. Contaminants are low-abundance and structureless.
11. (Optional, separate question) Reserve oral-vs-gut comparison for translocation /
    oral-specificity biology — never as a residence gate.

## Phase 4 — Phylotype-level analysis

12. **Phylogenetic placement** of surviving ASVs onto the **Greengenes2** backbone via SEPP
    (`q2-fragment-insertion`). GG2 is chosen as the scaffold because it is SEPP-native
    (built for fragment insertion), unifies multiple amplicon regions on one tree (a V4 and
    a V1-V3 ASV of the same organism cannot be sequence-clustered but *can* be co-placed),
    and keeps the analysis in the same coordinate system as the rest of the pipeline.
    **HOMD taxa are overlaid as labels** on the GG2 tree (HOMD alone is too sparse to be a
    placement backbone).
13. **Define candidate phylotypes** from the placements (e.g. `gappa examine cluster` or
    monophyletic groups): clusters of placed ASVs that (a) group together on the tree,
    (b) have no HOMD reference in their genus/species-level neighborhood, (c) are supported
    by >= N independent subjects across >= M studies and reach an abundance threshold, and
    (d) survive Phases 2-3.
14. **Novelty verdict via multiple independent references** — check phylotype representative
    sequences against **HOMD + GG2 + SILVA/NCBI 16S** (independence matters: the DB used for
    placement should not be the sole arbiter of novelty). Split into two publishable bins:
    - *known bacterium, absent from HOMD* (match in SILVA/NCBI but not HOMD) -> candidate
      additions to the oral catalog;
    - *novel to all references* -> candidate new species/lineages.

## Phase 5 — Characterization and connection to reference bias

15. Profile each phylotype by oral site, geography, prevalence, and co-occurrence.
16. Cross-link to the geographic reference-bias result: phylotypes that are novel **and**
    restricted to under-referenced populations (Asia, Africa — the continents that map
    significantly worse to GG2 after amplicon adjustment) are the headline — concrete,
    named diversity the reference is missing in exactly the populations shown to be
    under-represented.

## Existing assets (branch `worktree-cr-reference-bias`)

| Item | Path | Status |
|---|---|---|
| Reference-bias regression + figures | `script/05_homd_gg2_closed_reference/23_cr_mapping_by_continent_amplicon_adjusted.py` | runs locally |
| Per-taxon geographic reference gap | `script/05_homd_gg2_closed_reference/24_geography_reference_gap_bacterial_asvs.py` | cluster-only, smoke-tested |
| Negative-control audit | `script/05_homd_gg2_closed_reference/25_negative_control_audit.py` | runs locally (23 projects) |
| Contaminant genus reference | `script/05_homd_gg2_closed_reference/contaminant_reference_genera.tsv` | curated |
| Contaminant screen (Phase 2 core) | `script/05_homd_gg2_closed_reference/26_flag_contaminants_low_homd.py` | cluster input, smoke-tested |
| Reference-bias writeup | `results/geographic_reference_bias.md` | draft |

## To build next

- Phase 1 chimera-removal + per-region low-HOMD candidate export.
- Phase 2 `decontam` within the 23 blank studies + cross-lab catalog builder + oral-vs-blank
  enrichment (extends scripts 25/26).
- Phase 4 GG2 SEPP placement + phylotype definition + multi-DB novelty classification.

## Open items to confirm

1. vsearch `--iddef` actually used (query-global vs terminal-gap-excluding default).
2. HOMD version (V16.03 vs newer / eHOMD) for the "not in HOMD" determination.
3. Abundance and support thresholds (N subjects, M studies, min relative abundance) for
   phylotype calling — set from data, not a priori.
