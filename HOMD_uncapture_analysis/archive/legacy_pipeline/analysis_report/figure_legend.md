# Figure Legend: HOMD Coverage Analysis

## Composite Figure: Assessment of HOMD Reference Database Coverage

**A** | HOMD capture overview. Of 2,844,174 dereplicated 16S rRNA ASV sequences collected from 334 oral microbiome BioProjects, 997,603 (35.1%) were captured by HOMD V16.03 (identity ≥ 97%), while 1,846,571 (64.9%) remained uncaptured.

**B** | HOMD capture rate by 16S rRNA amplicon region. Capture rates vary substantially across sequencing regions: V1-V3 achieved the highest capture rate (51.7%), while V4 (26.6%) and V6-V8 (16.6%) had the lowest. Regions are ordered by total ASV count. Only regions with >50,000 ASVs are shown.

**C** | Reasons for HOMD non-capture. Among 1,846,571 uncaptured ASVs, the overwhelming majority (97.7%) failed due to low sequence identity (<97%), while 2.1% had no HOMD hit and 0.2% had insufficient query coverage (<99%).

**D** | Candidate novel phylotype landscape. Each point represents one of 56 candidate novel oral phylotypes identified after rigorous quality control and contamination filtering. x-axis: number of independent BioProjects; y-axis: number of independent subjects; point size: total read count; color: novelty classification. Labeled phylotypes have the highest subject or bioproject counts.

**E** | Novelty classification of candidate phylotypes. Of 56 candidate novel oral phylotypes, 14 (25%) are novel to all reference databases (HOMD, Greengenes2, SILVA), representing potential new species or lineages. The remaining 42 (75%) are known bacteria absent from HOMD, representing candidates for addition to the oral reference catalog.

---

## Methods Summary

Dereplicated 16S rRNA ASV sequences from 334 BioProjects were mapped against HOMD V16.03 using VSEARCH global alignment (identity threshold: 97%). Uncaptured sequences underwent chimera removal (vsearch --uchime_ref vs HOMD+GG2+SILVA), contamination filtering (decontam within-study, cross-lab contaminant catalog, literature blacklist, oral-vs-blank residence test), and phylogenetic placement into the Greengenes2 backbone tree via SEPP. Phylotypes were defined as groups of candidates sharing the same nearest GG2 reference, filtered for ≥10 subjects, ≥2 BioProjects, and ≥1 member passing oral residence.
