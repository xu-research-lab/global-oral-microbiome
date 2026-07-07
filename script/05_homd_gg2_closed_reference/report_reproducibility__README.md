# HOMD ASV representation report bundle

This directory contains the final English HTML report, final figures, selected summary tables, and scripts used for the HOMD/ASV representation analysis.

Main report:
- HOMD_ASV_representation_report.html

Subdirectories:
- figures/: final PNG/PDF figures used in the report
- tables/: selected small summary tables used in the report
- code/: analysis and plotting scripts required to reproduce the workflow

Large raw/intermediate files are not copied into this bundle. Their absolute paths are recorded in manifest.json.

Primary capture threshold: HOMD identity >=99% and query coverage >=99%.
Supplementary relaxed threshold: HOMD identity >=97% and query coverage >=99%.

Additional documentation:
- ANALYSIS_METHODS_AND_FILTERS.md: detailed ASV processing steps, filtering thresholds, mapping criteria, and reproduction order.

Report update:
- HOMD_ASV_representation_report.html is now organized by analysis module, with methods, filtering criteria, figures, interpretations, and related tables kept together.

Subject-level discovery update:
- Added fig15/fig15b and subject-ASV presence tables. Module 5 now uses subject-level discovery curves as the main discovery analysis.
