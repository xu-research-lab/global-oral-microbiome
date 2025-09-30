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
setwd("D:/oral_project/PRJNA316469")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries where Organism is "gut metagenome" (non-oral samples)
#initial_meta <- subset(initial_meta,isolation_source!= "feces")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR3991812   AMPLICON        544 13423259 PRJNA316469 SAMN04580922
    ## 2 SRR3991813   AMPLICON        526 20439344 PRJNA316469 SAMN04580956
    ## 3 SRR3991814   AMPLICON        548  7224060 PRJNA316469 SAMN04580993
    ## 4 SRR3991815   AMPLICON        561  8650460 PRJNA316469 SAMN04581022
    ## 5 SRR3991816   AMPLICON        557  3466449 PRJNA316469 SAMN04581064
    ## 6 SRR3991817   AMPLICON        562 16283513 PRJNA316469 SAMN04581095
    ##                BioSampleModel    Bytes               Center.Name Consent
    ## 1 Metagenome or environmental  6979848 NATIONAL CANCER INSTITUTE  public
    ## 2 Metagenome or environmental 11863981 NATIONAL CANCER INSTITUTE  public
    ## 3 Metagenome or environmental  3737750 NATIONAL CANCER INSTITUTE  public
    ## 4 Metagenome or environmental  4474865 NATIONAL CANCER INSTITUTE  public
    ## 5 Metagenome or environmental  1833751 NATIONAL CANCER INSTITUTE  public
    ## 6 Metagenome or environmental  8261482 NATIONAL CANCER INSTITUTE  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##                        env_biome     env_feature                   env_material
    ## 1                   human saliva healthy subject                   human saliva
    ## 2       human tongue dorsum swab healthy subject       human tongue dorsum swab
    ## 3         human hard palate swab healthy subject         human hard palate swab
    ## 4       human buccal mucosa swab healthy subject       human buccal mucosa swab
    ## 5 human Keratinized gingiva swab healthy subject human Keratinized gingiva swab
    ## 6    human Palatine tonsils swab healthy subject    human Palatine tonsils swab
    ##   ETHNICITY Experiment geo_loc_name_country geo_loc_name_country_continent
    ## 1 caucasian SRX1993057                  USA                  North America
    ## 2 caucasian SRX1993058                  USA                  North America
    ## 3 caucasian SRX1993059                  USA                  North America
    ## 4 caucasian SRX1993060                  USA                  North America
    ## 5 caucasian SRX1993061                  USA                  North America
    ## 6 caucasian SRX1993062                  USA                  North America
    ##          geo_loc_name host_age         HOST host_sex host_subject_id
    ## 1 USA: Rochester\\,NY    45.33 Homo sapiens   female         GYU0007
    ## 2 USA: Rochester\\,NY    45.33 Homo sapiens   female         GYU0007
    ## 3 USA: Rochester\\,NY    45.33 Homo sapiens   female         GYU0007
    ## 4 USA: Rochester\\,NY    45.33 Homo sapiens   female         GYU0007
    ## 5 USA: Rochester\\,NY    45.33 Homo sapiens   female         GYU0007
    ## 6 USA: Rochester\\,NY    45.33 Homo sapiens   female         GYU0007
    ##       Instrument                              Library.Name LibraryLayout
    ## 1 Illumina MiSeq  TGGTCACACTAAGCCT: GYU0007.0102.H4.RAV411        PAIRED
    ## 2 Illumina MiSeq         GTCGTGATTAGACCTA: GYU0007.0201.H2        PAIRED
    ## 3 Illumina MiSeq  CATGCCTATGAACCTT: GYU0007.0301.A5.RAV409        PAIRED
    ## 4 Illumina MiSeq  AGGAGTCCCTAATCGA: GYU0007.0401.E2.RAV409        PAIRED
    ## 5 Illumina MiSeq CCACTCCTCTAAGCCT: GYU0007.0501.H10.RAV411        PAIRED
    ## 6 Illumina MiSeq  TTGACCCTTGTTCTCT: GYU0007.0601.C8.RAV411        PAIRED
    ##   LibrarySelection LibrarySource         Organism Platform          ReleaseDate
    ## 1              PCR   METAGENOMIC human metagenome ILLUMINA 2016-08-01T00:00:00Z
    ## 2              PCR   METAGENOMIC human metagenome ILLUMINA 2016-08-01T00:00:00Z
    ## 3              PCR   METAGENOMIC human metagenome ILLUMINA 2016-08-01T00:00:00Z
    ## 4              PCR   METAGENOMIC human metagenome ILLUMINA 2016-08-01T00:00:00Z
    ## 5              PCR   METAGENOMIC human metagenome ILLUMINA 2016-08-01T00:00:00Z
    ## 6              PCR   METAGENOMIC human metagenome ILLUMINA 2016-08-01T00:00:00Z
    ##               Sample.Name smoker      source_material_id SRA.Study
    ## 1  GYU0007.0102.H4.RAV411  never  GYU0007.0102.H4.RAV411 SRP080748
    ## 2         GYU0007.0201.H2  never         GYU0007.0201.H2 SRP080748
    ## 3  GYU0007.0301.A5.RAV409  never  GYU0007.0301.A5.RAV409 SRP080748
    ## 4  GYU0007.0401.E2.RAV409  never  GYU0007.0401.E2.RAV409 SRP080748
    ## 5 GYU0007.0501.H10.RAV411  never GYU0007.0501.H10.RAV411 SRP080748
    ## 6  GYU0007.0601.C8.RAV411  never  GYU0007.0601.C8.RAV411 SRP080748
    ##            create_date version
    ## 1 2016-08-01T14:20:00Z       1
    ## 2 2016-08-01T14:59:00Z       1
    ## 3 2016-08-01T14:20:00Z       1
    ## 4 2016-08-01T14:20:00Z       1
    ## 5 2016-08-01T14:20:00Z       1
    ## 6 2016-08-01T14:20:00Z       1

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

``` r
# Create a uniform metadata dataframe
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$host_subject_id,
  Age = initial_meta$host_age,
  BMI = NA,
  Sex=initial_meta$host_sex,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$env_material,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- NA
metadata$title <- "The effect of cigarette smoking on the oral and nasal microbiota"
#metadata<-subset(uniform_meta,Site!="feces")
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status <- "H"
#extract subject id by remove prefix of sample_name
#metadata$sample_name <- gsub("[FSP].*", "", metadata$sample_name)
#get subjects id with Saliva available
keep_subjects<-subset(metadata,Site %in% c("human hard palate swab","human Keratinized gingiva swab","human oral wash","human saliva","human Supragingival plaque scrape","human tongue dorsum swab"))
#remove subjects without Saliva samples
metadata<-subset(metadata,sample_name %in% keep_subjects$sample_name)
#get subjects sampling counts
sample_counts <- table(metadata$sample_name)
# Create matched_samples column, any subjects have >=2 samples will be label "YES", otherwise "NO"
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")

metadata$Site<-gsub("human ","",metadata$Site)
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 328"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 43"

# Step 5: Visualize the Metadata

``` r
knitr::include_graphics("Metadata_Summary.png")
```

<img src="Metadata_Summary.png" width="1200" />

# Step 6: Finalize and Export the Metadata

``` r
# Reorder columns to match the desired format
metadata <- metadata[, c(
  "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
  "Region", "Platform", "Country", "Site", "public_time", 
  "PCR_primers", "title", "Disease_study", "Status","matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##    BioProject        Run sample_name   Age BMI    Sex Region       Platform
    ## 1 PRJNA316469 SRR3991812     GYU0007 45.33  NA female  V3-V4 Illumina MiSeq
    ## 2 PRJNA316469 SRR3991813     GYU0007 45.33  NA female  V3-V4 Illumina MiSeq
    ## 3 PRJNA316469 SRR3991814     GYU0007 45.33  NA female  V3-V4 Illumina MiSeq
    ## 4 PRJNA316469 SRR3991815     GYU0007 45.33  NA female  V3-V4 Illumina MiSeq
    ## 5 PRJNA316469 SRR3991816     GYU0007 45.33  NA female  V3-V4 Illumina MiSeq
    ## 6 PRJNA316469 SRR3991817     GYU0007 45.33  NA female  V3-V4 Illumina MiSeq
    ##   Country                     Site          public_time PCR_primers
    ## 1     USA                   saliva 2016-08-01T00:00:00Z          NA
    ## 2     USA       tongue dorsum swab 2016-08-01T00:00:00Z          NA
    ## 3     USA         hard palate swab 2016-08-01T00:00:00Z          NA
    ## 4     USA       buccal mucosa swab 2016-08-01T00:00:00Z          NA
    ## 5     USA Keratinized gingiva swab 2016-08-01T00:00:00Z          NA
    ## 6     USA    Palatine tonsils swab 2016-08-01T00:00:00Z          NA
    ##                                                              title
    ## 1 The effect of cigarette smoking on the oral and nasal microbiota
    ## 2 The effect of cigarette smoking on the oral and nasal microbiota
    ## 3 The effect of cigarette smoking on the oral and nasal microbiota
    ## 4 The effect of cigarette smoking on the oral and nasal microbiota
    ## 5 The effect of cigarette smoking on the oral and nasal microbiota
    ## 6 The effect of cigarette smoking on the oral and nasal microbiota
    ##   Disease_study Status matched_samples
    ## 1            NO      H             YES
    ## 2            NO      H             YES
    ## 3            NO      H             YES
    ## 4            NO      H             YES
    ## 5            NO      H             YES
    ## 6            NO      H             YES

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
