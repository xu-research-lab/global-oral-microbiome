R Notebook: Metadata Processing Workflow
================

# Step 1: Obtain Initial Metadata from Each Project

To obtain the initial metadata: 1. Search for the BioProject ID (e.g.,
“PRJEBXXXX”) on the [NCBI SRA
website](https://www.ncbi.nlm.nih.gov/Traces/study/). 2. In the “Select”
table, click on “Metadata” and download the file. 3. The downloaded file
is typically named “SraRunTable.csv” and contains metadata for the
project.

``` r
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(cowplot)

# Set working directory (update this path to your local directory)
setwd("D:/metadata/PRJEB51725")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries non-oral samples
#initial_meta <- subset(initial_meta,isolate %in% c("saliva","tongue coat"))

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen     Bases BioProject     BioSample     Bytes
    ## 1 ERR9321399      OTHER        501 276512103 PRJEB51725 SAMEA13465139 196153039
    ## 2 ERR9321403      OTHER        501 108988401 PRJEB51725 SAMEA13465151  82597893
    ## 3 ERR9321408      OTHER        501 258263353 PRJEB51725 SAMEA13465140 189149611
    ## 4 ERR9321417      OTHER        501 258339970 PRJEB51725 SAMEA13465141 185266909
    ## 5 ERR9321418      OTHER        501 344352440 PRJEB51725 SAMEA13465144 245957565
    ## 6 ERR9321419      OTHER        501 303120296 PRJEB51725 SAMEA13465147 211182206
    ##                          Center.Name common_name Consent DATASTORE.filetype
    ## 1 CONSIGLIO NAZIONALE DELLE RICERCHE       human  public   fastq,run.zq,sra
    ## 2 CONSIGLIO NAZIONALE DELLE RICERCHE       human  public   fastq,run.zq,sra
    ## 3 CONSIGLIO NAZIONALE DELLE RICERCHE       human  public   fastq,run.zq,sra
    ## 4 CONSIGLIO NAZIONALE DELLE RICERCHE       human  public   fastq,run.zq,sra
    ## 5 CONSIGLIO NAZIONALE DELLE RICERCHE       human  public   fastq,run.zq,sra
    ## 6 CONSIGLIO NAZIONALE DELLE RICERCHE       human  public   fastq,run.zq,sra
    ##   DATASTORE.provider                     DATASTORE.region
    ## 1         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   ENA.FIRST.PUBLIC..run. ena_first_public ENA.LAST.UPDATE..run. ena_last_update
    ## 1             2023-03-17       2023-03-17            2023-03-17      2023-03-17
    ## 2             2023-03-17       2023-03-17            2023-03-17      2023-03-17
    ## 3             2023-03-17       2023-03-17            2023-03-17      2023-03-17
    ## 4             2023-03-17       2023-03-17            2023-03-17      2023-03-17
    ## 5             2023-03-17       2023-03-17            2023-03-17      2023-03-17
    ## 6             2023-03-17       2023-03-17            2023-03-17      2023-03-17
    ##   Experiment   external_id                 INSDC_center_alias
    ## 1 ERX8863398 SAMEA13465139 Consiglio Nazionale delle Ricerche
    ## 2 ERX8863402 SAMEA13465151 Consiglio Nazionale delle Ricerche
    ## 3 ERX8863407 SAMEA13465140 Consiglio Nazionale delle Ricerche
    ## 4 ERX8863416 SAMEA13465141 Consiglio Nazionale delle Ricerche
    ## 5 ERX8863417 SAMEA13465144 Consiglio Nazionale delle Ricerche
    ## 6 ERX8863418 SAMEA13465147 Consiglio Nazionale delle Ricerche
    ##                    INSDC_center_name   INSDC_first_public    INSDC_last_update
    ## 1 Consiglio Nazionale delle Ricerche 2023-03-17T00:30:27Z 2023-03-17T00:30:27Z
    ## 2 Consiglio Nazionale delle Ricerche 2023-03-17T00:30:27Z 2023-03-17T00:30:27Z
    ## 3 Consiglio Nazionale delle Ricerche 2023-03-17T00:30:27Z 2023-03-17T00:30:27Z
    ## 4 Consiglio Nazionale delle Ricerche 2023-03-17T00:30:27Z 2023-03-17T00:30:27Z
    ## 5 Consiglio Nazionale delle Ricerche 2023-03-17T00:31:31Z 2023-03-17T00:31:31Z
    ## 6 Consiglio Nazionale delle Ricerche 2023-03-17T00:30:27Z 2023-03-17T00:30:27Z
    ##   INSDC_status     Instrument Library.Name LibraryLayout LibrarySelection
    ## 1       public Illumina MiSeq           A4        PAIRED            other
    ## 2       public Illumina MiSeq           A8        PAIRED            other
    ## 3       public Illumina MiSeq           B4        PAIRED            other
    ## 4       public Illumina MiSeq           C4        PAIRED            other
    ## 5       public Illumina MiSeq           C5        PAIRED            other
    ## 6       public Illumina MiSeq           C6        PAIRED            other
    ##   LibrarySource     Organism Platform          ReleaseDate          create_date
    ## 1   METAGENOMIC Homo sapiens ILLUMINA 2023-03-18T00:00:00Z 2023-07-18T10:41:00Z
    ## 2   METAGENOMIC Homo sapiens ILLUMINA 2023-03-18T00:00:00Z 2023-07-18T10:41:00Z
    ## 3   METAGENOMIC Homo sapiens ILLUMINA 2023-03-18T00:00:00Z 2023-07-18T10:41:00Z
    ## 4   METAGENOMIC Homo sapiens ILLUMINA 2023-03-18T00:00:00Z 2023-07-18T10:41:00Z
    ## 5   METAGENOMIC Homo sapiens ILLUMINA 2023-03-18T00:00:00Z 2023-07-18T10:42:00Z
    ## 6   METAGENOMIC Homo sapiens ILLUMINA 2023-03-18T00:00:00Z 2023-07-18T10:42:00Z
    ##   version   Sample.Name Sample_name Scientific_Name SRA.Study Submitter_Id
    ## 1       1 SAMEA13465139          A4    Homo sapiens ERP136380           A4
    ## 2       1 SAMEA13465151          A8    Homo sapiens ERP136380           A8
    ## 3       1 SAMEA13465140          B4    Homo sapiens ERP136380           B4
    ## 4       1 SAMEA13465141          C4    Homo sapiens ERP136380           C4
    ## 5       1 SAMEA13465144          C5    Homo sapiens ERP136380           C5
    ## 6       1 SAMEA13465147          C6    Homo sapiens ERP136380           C6
    ##   Collection_Date environment_.biome. environment_.feature.
    ## 1            2017                oral                 mouth
    ## 2            2017                oral                 mouth
    ## 3            2017                oral                 mouth
    ## 4            2017                oral                 mouth
    ## 5            2017                oral                 mouth
    ## 6            2017                oral                 mouth
    ##   environment_.material. geo_loc_name_country geo_loc_name_country_continent
    ## 1   supragingival plaque                Italy                         Europe
    ## 2   supragingival plaque                Italy                         Europe
    ## 3   supragingival plaque                Italy                         Europe
    ## 4   supragingival plaque                Italy                         Europe
    ## 5   supragingival plaque                Italy                         Europe
    ## 6   supragingival plaque                Italy                         Europe
    ##   geographic_location_.country_and.or_sea. geographic_location_.latitude.
    ## 1                                    Italy                   not provided
    ## 2                                    Italy                   not provided
    ## 3                                    Italy                   not provided
    ## 4                                    Italy                   not provided
    ## 5                                    Italy                   not provided
    ## 6                                    Italy                   not provided
    ##   geographic_location_.longitude. human_oral_environmental_package
    ## 1                    not provided                       human-oral
    ## 2                    not provided                       human-oral
    ## 3                    not provided                       human-oral
    ## 4                    not provided                       human-oral
    ## 5                    not provided                       human-oral
    ## 6                    not provided                       human-oral
    ##   investigation_type project_name Sequencing_method
    ## 1         metagenome   PRJEB51725          Illumina
    ## 2         metagenome   PRJEB51725          Illumina
    ## 3         metagenome   PRJEB51725          Illumina
    ## 4         metagenome   PRJEB51725          Illumina
    ## 5         metagenome   PRJEB51725          Illumina
    ## 6         metagenome   PRJEB51725          Illumina
    ##   geographic_location_.region_and_locality. isolate receipt_date
    ## 1                                                             NA
    ## 2                                                             NA
    ## 3                                                             NA
    ## 4                                                             NA
    ## 5                                                             NA
    ## 6                                                             NA

# Step 2: Create a Uniform Metadata Template

The uniform metadata should include the following columns:

1.  **BioProject**: Project ID (e.g., “PRJNA339212”)
2.  **Run**: Sample ID (e.g., “ERRXXXXXX” or “SRRXXXXXX”)
3.  **sample_name**: Original sample name submitted by the authors (not
    the database-assigned ID)
4.  **Age**: Age of the subject (if available)
5.  **BMI**: Body Mass Index (if available)
6.  **Sex**: Sex of the subject (“F” or “M”)
7.  **Region**: Targeted 16S rRNA region (e.g., “V3-V4”)
8.  **Platform**: Sequencing platform (e.g., “Illumina MiSeq”)
9.  **Country**: Country of sample origin (e.g., “China”)
10. **Site**: Sample site (e.g., “saliva”, “dental plaque”)
11. **public_time**: Data release date (format: “YYYY-MM-DD”)
12. **PCR_primers**: PCR primers used (e.g., “515F-806R”)
13. **title**: Title of the associated publication
14. **Disease_study**: Whether the study involves a disease (“YES” or
    “NO”)
15. **Status**: Disease status (e.g., “T2DM”, “HC”). Note: If the study
    includes both disease and healthy samples, this must be specified
    per sample.
16. **matched_samples**: matched samples (“YES” or “NO”)

``` r
# Create a uniform metadata dataframe
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Sample.Name,
  Age = NA,
  BMI = NA,
  Sex=NA,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$environment_.material.,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "341F-805R"
metadata$title <- "Oral ecological environment modifications by hard-cheese: from pH to microbiome: a prospective cohort study based on 16S rRNA metabarcoding approach"
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status<-"H"

keep_subjects <- unique(subset(metadata, Site == "supragingival plaque")$sample_name)
metadata <- subset(metadata, sample_name %in% keep_subjects)

metadata$matched_samples <- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 27"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 27"

# Step 5: Visualize the Metadata

``` r
knitr::include_graphics("metadata_summary.png")
```

![](metadata_summary.png)<!-- -->

# Step 6: Finalize and Export the Metadata

``` r
# Reorder columns to match the desired format
metadata <- metadata[, c(
  "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
  "Region", "Platform", "Country", "Site", "public_time", 
  "PCR_primers", "title", "Disease_study", "Status","matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = TRUE)
head(metadata)
```

    ##   BioProject        Run   sample_name Age BMI Sex Region       Platform Country
    ## 1 PRJEB51725 ERR9321399 SAMEA13465139  NA  NA  NA  V3-V4 Illumina MiSeq   Italy
    ## 2 PRJEB51725 ERR9321403 SAMEA13465151  NA  NA  NA  V3-V4 Illumina MiSeq   Italy
    ## 3 PRJEB51725 ERR9321408 SAMEA13465140  NA  NA  NA  V3-V4 Illumina MiSeq   Italy
    ## 4 PRJEB51725 ERR9321417 SAMEA13465141  NA  NA  NA  V3-V4 Illumina MiSeq   Italy
    ## 5 PRJEB51725 ERR9321418 SAMEA13465144  NA  NA  NA  V3-V4 Illumina MiSeq   Italy
    ## 6 PRJEB51725 ERR9321419 SAMEA13465147  NA  NA  NA  V3-V4 Illumina MiSeq   Italy
    ##                   Site          public_time PCR_primers
    ## 1 supragingival plaque 2023-03-18T00:00:00Z   341F-805R
    ## 2 supragingival plaque 2023-03-18T00:00:00Z   341F-805R
    ## 3 supragingival plaque 2023-03-18T00:00:00Z   341F-805R
    ## 4 supragingival plaque 2023-03-18T00:00:00Z   341F-805R
    ## 5 supragingival plaque 2023-03-18T00:00:00Z   341F-805R
    ## 6 supragingival plaque 2023-03-18T00:00:00Z   341F-805R
    ##                                                                                                                                                  title
    ## 1 Oral ecological environment modifications by hard-cheese: from pH to microbiome: a prospective cohort study based on 16S rRNA metabarcoding approach
    ## 2 Oral ecological environment modifications by hard-cheese: from pH to microbiome: a prospective cohort study based on 16S rRNA metabarcoding approach
    ## 3 Oral ecological environment modifications by hard-cheese: from pH to microbiome: a prospective cohort study based on 16S rRNA metabarcoding approach
    ## 4 Oral ecological environment modifications by hard-cheese: from pH to microbiome: a prospective cohort study based on 16S rRNA metabarcoding approach
    ## 5 Oral ecological environment modifications by hard-cheese: from pH to microbiome: a prospective cohort study based on 16S rRNA metabarcoding approach
    ## 6 Oral ecological environment modifications by hard-cheese: from pH to microbiome: a prospective cohort study based on 16S rRNA metabarcoding approach
    ##   Disease_study Status matched_samples
    ## 1            NO      H              NO
    ## 2            NO      H              NO
    ## 3            NO      H              NO
    ## 4            NO      H              NO
    ## 5            NO      H              NO
    ## 6            NO      H              NO

## Notes for Students:

1.  **Metadata Accuracy**: Always cross-check the metadata with the
    original publication to ensure accuracy.
2.  **Manual Entries**: Columns like `PCR_primers`, and `Status` often
    require manual entry based on the publication.
3.  **Disease Status**: Be careful when assigning `Status`—some studies
    include both disease and healthy samples.
4.  **Visualization**: The plots help quickly identify the distribution
    of values in each column, which is useful for quality control.
5.  **Reproducibility**: Update the working directory and file paths to
    match your local environment.
