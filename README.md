# Global Oral Microbiome

This repository brings together curated 16S rRNA sequencing metadata, microbial
feature tables, and analysis workflows for studying the oral microbiome across
geographies, sampling sites, host characteristics, and disease cohorts. It
includes data-curation records and analyses of diversity, taxonomic
composition, disease associations, prediction, microbial networks, and
reference-database coverage.

The repository is a research workspace, not a single installable software
package. Analyses are maintained as scripts, notebooks, and workflow-specific
reports; inputs, dependencies, and execution requirements vary between them.

## Repository layout

| Path | Contents |
| --- | --- |
| `data/` | Curated project metadata, QC metadata and BIOM tables, oral-gut paired data, and per-study metadata-curation files. |
| `script/` | Shared and staged workflows for data processing, QC, sampling overview, diversity, reference profiling, host and technical factors, geography, disease, and prediction. |
| `HOMD_uncapture_analysis/` | Analysis of oral ASVs not represented in eHOMD, including control-based contamination screening and supported phylotype analyses. See its [README](HOMD_uncapture_analysis/README.md). |
| `overview_*/`, `diversity_*/`, `disease_*`, `differential_abundance/`, `age_*`, `continent_rf/`, and `network/` | Topic-specific analysis scripts, notebooks, results, and reports. |
| `results/` | Selected analysis outputs and reports. |

## Data

The `data/` directory contains the central metadata and feature-table resources.
The QC release documentation describes the available table versions, filtering
criteria, and matched metadata/BIOM pairs:

- [Metadata and BIOM QC release](data/metadata-and-biom-v2-qc/README.md)
- [Oral-gut paired QC release](data/oral-gut-pair-v2-qc/README.md)
- [Metadata cleaning mappings](data/mapping_sumary.md)

For analyses using the QC release, use the documented metadata and feature-table
pair together. The QC README also explains when the deduplicated oral subset is
appropriate and how to handle repeated measures. Individual studies are
organized under `data/projects/`; their curation notebooks and notes describe
study-specific metadata decisions.

## Analyses and reproducibility

Start with the README or report in the relevant analysis directory. In
particular, [`script/01_data_processing_qc/metadata_curation_v2.md`](script/01_data_processing_qc/metadata_curation_v2.md)
documents the metadata-curation conventions, while
[`HOMD_uncapture_analysis/README.md`](HOMD_uncapture_analysis/README.md)
describes the order and dependencies of that multi-stage workflow.

There is no repository-wide setup command. Workflows use a mix of Python and R,
and some require specialized tools such as QIIME 2, DADA2, VSEARCH, or SLURM
for high-performance computing. Scripts and notebooks may rely on large
reference databases, intermediate files, or compute-cluster paths that are not
included or configured for a local machine. Check the workflow-specific
documentation and input paths before running an analysis; do not assume a
notebook or script can be run directly from the repository root.

## Data use

This repository incorporates data and metadata from multiple studies. Before
redistributing or reusing any source data, check the terms and attribution
requirements of the original study and data repository.
