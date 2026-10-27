# Identification of Novel Oral Taxa Absent from the HOMD Reference Database

**Analysis Report** — 2026-07-21

---

## Executive Summary

We screened **3,396,776 ASVs** from **334** oral microbiome BioProjects
against HOMD V16.03, with rigorous contamination controls and phylotype-level clustering.

| Phase | Key Step | Result |
|-------|----------|--------|
| **1. Candidate Generation** | HOMD mapping → chimera removal → relaxed support floor | **54,290** non-chimeric HOMD-uncaptured candidates |
| **2. Contamination Firewall** | Control QC → decontam → oral-vs-blank residence | **47,203** clean, **14** cross-lab contaminants identified |
| **4. Phylotype Analysis** | GG2 nearest-neighbor → phylotype clustering → multi-DB novelty | **68** candidate novel oral phylotypes |

### Final Classification
- **32 phylotypes** — **novel to ALL references** (HOMD + GG2 + SILVA) → candidate new species/lineages
- **36 phylotypes** — **known bacteria, absent from HOMD** → candidate additions to the oral catalog

---

## 1. Phase 1 — Candidate Generation

**Input**: 3,396,776 ASVs from 334 BioProjects → 2,844,174 dereplicated sequences → vsearch global vs HOMD V16.03

| HOMD Mapping | Count | % |
|-------------|-------|---|
| Captured (id>=97% & qcov>=99%) | 997,603 | 35.1 |
| Uncaptured | 1,846,571 | 64.9 |
| — Low identity (<97%) | 1,804,155 | 97.7 |
| — No HOMD hit | 39,146 | 2.1 |
| — Low qcov (<99%) | 3,270 | 0.2 |

**Chimera removal** (vsearch --uchime_ref vs HOMD+GG2+SILVA): 744,715 chimeras removed (26.18% rate).
**Relaxed support floor** (>=2 BioProjects, >=50 reads): 69,858 candidates retained (rare tail preserved).

![Pipeline Flow](fig1_pipeline.png)
![HOMD Mapping](fig2_homd.png)

---

## 2. Phase 2 — Contamination Firewall

### Control Sample QC

1,469 candidate control samples (extraction blanks, NTCs, PCR blanks) from 23 BioProjects were
re-analyzed via DADA2 and filtered with multi-dimensional criteria:

| Criterion | Threshold | Rationale |
|-----------|-----------|-----------|
| Total reads | 1 – 5,000 | Genuine blanks have minimal DNA |
| ASV count | <= 100 | True blanks have low diversity |
| Oral genera fraction | < 20% | High oral signal = likely carryover |
| Per-study | >=3 controls + >=5 oral | Must have sufficient power for decontam |

**Result**: **967** reliable controls retained across **15** BioProjects;
**502** excluded (incl. 5 entire studies where all controls were high-read).

### Decontam & Residence

- **Decontam** (prevalence mode, 10 BioProjects): 655 features flagged within single studies; **14** cross-lab contaminants (>=2 independent BioProjects)
- **Oral-vs-Blank Residence**: 79 / 6,294 features pass (oral enrichment over blanks)

### Contamination Concern Summary

| Level | Count | % |
|-------|-------|---|
| Clean | 47,203 | 86.9 |
| Residence fail | 4,468 | 8.2 |
| In contaminant catalog | 1,978 | 3.6 |
| Decontam (single-study) | 627 | 1.2 |
| Cross-lab contaminant | **14** | **0.0** |

![Contamination Firewall](fig3_contamination.png)

---

## 3. Phase 4 — Phylotype-Level Analysis

### Phylogenetic Placement

Candidate ASVs mapped to the Greengenes2 backbone (337,506 full-length 16S sequences) via
vsearch global alignment (`--iddef 2`). Each candidate assigned to nearest GG2 reference.
GG2 sequences matching HOMD at >=97% identity + >=98% qcov marked as "HOMD-proximal".

- **Candidates placed**: 104,262
- **HOMD-gap** (nearest ref NOT HOMD): 101,923 (97.8%)
- **Unique GG2 clusters**: 16,731
- **Median identity to GG2 ref**: 90.6%

*True SEPP phylogenetic placement via `q2-fragment-insertion sepp` currently running (SLURM 3051311).*

### Phylotype Definition Gates

Phylotypes = groups of candidates sharing the same nearest GG2 reference. A phylotype is
a **candidate novel oral phylotype** if it:
1. Has no HOMD reference in its GG2 neighborhood
2. Has >=10 independent subjects
3. Has >=2 independent BioProjects
4. Has >=1 member passing oral-vs-blank residence

### Multi-Database Novelty

| Classification | Count | % |
|----------------|-------|---|
| **novel_all_refs** (new to HOMD + GG2 + SILVA) | **32** | **47.1** |
| **known_not_in_homd** (in GG2/SILVA but not HOMD) | **36** | **52.9** |

### Top Phylotypes

| Phylotype | Subjects | BPs | Reads | Novelty |
|-----------|----------|-----|-------|---------|
| GG2_G001821375 | 3,497 | 17 | 3,611,242 | known not in homd |
| GG2_RS-GCF-000210395.1-NC-017175.1--13 | 2,893 | 33 | 35,375,427 | known not in homd |
| GG2_GB-GCA-009724635.1-RMVO01000037.1 | 2,547 | 27 | 5,450,170 | novel all refs |
| GG2_LS992100.88603424.88604816 | 2,547 | 27 | 11,307,097 | known not in homd |
| GG2_RS-GCF-001984285.1-NZ-MUJR01000009.1--2 | 2,462 | 32 | 12,562,970 | known not in homd |
| GG2_RS-GCF-013293565.1-NZ-JAAKAJ010000083.1 | 1,962 | 31 | 60,689,397 | novel all refs |
| GG2_RS-GCF-013304685.1-NZ-JAAITH010000167.1 | 1,962 | 13 | 2,131,426 | known not in homd |
| GG2_RS-GCF-003260125.1-NZ-QHKO01000005.1 | 1,877 | 27 | 13,174,634 | known not in homd |
| GG2_G000013145 | 1,867 | 17 | 654,880 | novel all refs |
| GG2_G001793335 | 1,867 | 17 | 654,219 | novel all refs |

Median support: 502 subjects, 18 BioProjects, 215,998 reads.

![Phylotype Results](fig4_phylotypes.png)
![Characterization](fig5_characterization.png)

---

## 4. Output Files

| File | Description |
|------|-------------|
| `pipeline/05_integrated_candidates/integrated_candidate_annotations.tsv` | Full ASV annotation (54K x 50 cols) |
| `pipeline/06_phylotypes/phylotypes_candidate_novel.tsv` | 68 candidate novel phylotypes |
| `pipeline/07_novelty/phylotypes_novelty_classified.tsv` | Novelty-classified phylotypes |
| `control_quality_viz/genuine_controls_keep_v2.csv` | 967 reliable control samples |

---

## 5. Methods

| Step | Method | Tool |
|------|--------|------|
| HOMD mapping | Global alignment (`--iddef 1`) | VSEARCH |
| Chimera removal | Reference-based (`--uchime_ref`) vs HOMD+GG2+SILVA | VSEARCH |
| Control QC | Multi-dimensional filtering | Custom Python |
| Decontam | Prevalence mode (threshold=0.5) | R/decontam |
| Residence | Oral-vs-blank prevalence & enrichment | Custom Python |
| Phylogenetic placement | Nearest-neighbor to GG2 backbone | VSEARCH (`--iddef 2`) |
| Phylotype definition | Nearest-ref clustering + 4-gate filter | Custom Python |
| Multi-DB novelty | Global alignment vs HOMD, GG2, SILVA | VSEARCH |

*SEPP placement via `q2-fragment-insertion sepp` is in progress and will be updated.*
