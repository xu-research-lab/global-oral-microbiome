# Residual Limitations for Blank-Free Studies

## Background
Of 334 BioProjects in this analysis, 316 (94.6%) have NO 
negative control samples (blanks, NTCs, extraction blanks). For these studies, 
blank-based decontamination and oral-vs-blank residence testing cannot be performed.

## Contamination filtering available for ALL studies (blank-free included)

1. **Literature blacklist** (`contaminant_reference_genera.tsv`): 43 auto_remove + 19 flag_review genera
   - Applied to ALL candidates regardless of study
   - Sources: Salter 2014, Weyrich 2019 (Eisenhofer), cross-verified

2. **Genus-level pooled decontam**: Fisher exact test on ALL 1,469 controls vs ALL 7,281 oral samples
   - Cross-study contamination signal
   - 659 genera flagged (p<0.5)

3. **Cross-lab contaminant ASV catalog**: 14,685 ASV sequences from 1,469 controls
   - Sequence-level matching

## Not available for blank-free studies

4. **Within-study decontam** (Strategy A): requires blanks within same study
   - 17/20 blank-bearing studies had sufficient blanks

5. **Oral-vs-blank residence test**: requires paired blanks
   - The residence pass gate (10,924→56 phylotypes) relies on oral-vs-blank comparison

## Residual Risk

- Contaminants specific to a single blank-free study's reagent kit cannot be detected
- Novel taxa restricted to blank-free studies have weaker oral-residence evidence
- For the 316 blank-free studies, contamination flags rely on cross-study signals only

## Recommended Mitigation

- Flag all candidates from blank-free studies as "weaker evidence"
- Require >=3 independent BioProjects for candidates from blank-free studies
- Manual review of any candidate found primarily in blank-free studies
