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
setwd("D:/metadata/2/PRJEB20889")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries non-oral samples
#initial_meta <- subset(initial_meta,BioSampleModel != "Metagenome or environmental")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen   Bases BioProject      BioSample   Bytes
    ## 1 ERR1983891   AMPLICON        302 5527204 PRJEB20889 SAMEA104079181 3262643
    ## 2 ERR1983892   AMPLICON        302 3735136 PRJEB20889 SAMEA104079182 2153188
    ## 3 ERR1983893   AMPLICON        302 6441358 PRJEB20889 SAMEA104079183 3877544
    ## 4 ERR1983894   AMPLICON        302 4417052 PRJEB20889 SAMEA104079184 2580728
    ## 5 ERR1983895   AMPLICON        302 2140576 PRJEB20889 SAMEA104079185 1260573
    ## 6 ERR1983896   AMPLICON        302 5215540 PRJEB20889 SAMEA104079186 3275267
    ##                        Center.Name Collection_Date Consent DATASTORE.filetype
    ## 1 UNIVERSITY OF COLORADO - BOULDER            2003  public   fastq,run.zq,sra
    ## 2 UNIVERSITY OF COLORADO - BOULDER            2003  public   fastq,run.zq,sra
    ## 3 UNIVERSITY OF COLORADO - BOULDER            2003  public   fastq,run.zq,sra
    ## 4 UNIVERSITY OF COLORADO - BOULDER            2003  public   fastq,run.zq,sra
    ## 5 UNIVERSITY OF COLORADO - BOULDER            2003  public   fastq,run.zq,sra
    ## 6 UNIVERSITY OF COLORADO - BOULDER            2003  public   fastq,run.zq,sra
    ##   DATASTORE.provider                         DATASTORE.region
    ## 1     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1
    ## 2     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1
    ## 3     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1
    ## 4     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1
    ## 5     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1
    ## 6     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1
    ##   ENA.FIRST.PUBLIC..run. ENA.LAST.UPDATE..run. environment_.biome.
    ## 1             2017-06-29            2018-11-16                oral
    ## 2             2017-06-29            2018-11-16                oral
    ## 3             2017-06-29            2018-11-16                oral
    ## 4             2017-06-29            2018-11-16                oral
    ## 5             2017-06-29            2018-11-16                oral
    ## 6             2017-06-29            2018-11-16                oral
    ##   environment_.feature. environment_.material. Experiment    external_id
    ## 1                  oral                 saliva ERX2044267 SAMEA104079181
    ## 2                  oral                 saliva ERX2044268 SAMEA104079182
    ## 3                  oral                 saliva ERX2044269 SAMEA104079183
    ## 4                  oral                 saliva ERX2044270 SAMEA104079184
    ## 5                  oral                 saliva ERX2044271 SAMEA104079185
    ## 6                  oral                 saliva ERX2044272 SAMEA104079186
    ##   geo_loc_name_country geo_loc_name_country_continent
    ## 1                  USA                  North America
    ## 2                  USA                  North America
    ## 3                  USA                  North America
    ## 4                  USA                  North America
    ## 5                  USA                  North America
    ## 6                  USA                  North America
    ##   geographic_location_.country_and.or_sea. geographic_location_.latitude.
    ## 1                                      USA                       40.74189
    ## 2                                      USA                       40.74189
    ## 3                                      USA                       40.74189
    ## 4                                      USA                       40.74189
    ## 5                                      USA                       40.74189
    ## 6                                      USA                       40.74189
    ##   geographic_location_.longitude. humam_oral_environmental_package
    ## 1                       -73.98931                       human-oral
    ## 2                       -73.98931                       human-oral
    ## 3                       -73.98931                       human-oral
    ## 4                       -73.98931                       human-oral
    ## 5                       -73.98931                       human-oral
    ## 6                       -73.98931                       human-oral
    ##                  INSDC_center_name   INSDC_first_public    INSDC_last_update
    ## 1 UNIVERSITY OF COLORADO - BOULDER 2017-06-29T17:02:36Z 2017-05-16T19:48:02Z
    ## 2 UNIVERSITY OF COLORADO - BOULDER 2017-06-29T17:02:36Z 2017-05-16T19:48:02Z
    ## 3 UNIVERSITY OF COLORADO - BOULDER 2017-06-29T17:02:36Z 2017-05-16T19:48:02Z
    ## 4 UNIVERSITY OF COLORADO - BOULDER 2017-06-29T17:02:36Z 2017-05-16T19:48:02Z
    ## 5 UNIVERSITY OF COLORADO - BOULDER 2017-06-29T17:02:36Z 2017-05-16T19:48:02Z
    ## 6 UNIVERSITY OF COLORADO - BOULDER 2017-06-29T17:02:36Z 2017-05-16T19:48:02Z
    ##   INSDC_status     Instrument investigation_type Library.Name LibraryLayout
    ## 1       public Illumina MiSeq         metagenome  unspecified        PAIRED
    ## 2       public Illumina MiSeq         metagenome  unspecified        PAIRED
    ## 3       public Illumina MiSeq         metagenome  unspecified        PAIRED
    ## 4       public Illumina MiSeq         metagenome  unspecified        PAIRED
    ## 5       public Illumina MiSeq         metagenome  unspecified        PAIRED
    ## 6       public Illumina MiSeq         metagenome  unspecified        PAIRED
    ##   LibrarySelection LibrarySource        Organism Platform
    ## 1              PCR   METAGENOMIC oral metagenome ILLUMINA
    ## 2              PCR   METAGENOMIC oral metagenome ILLUMINA
    ## 3              PCR   METAGENOMIC oral metagenome ILLUMINA
    ## 4              PCR   METAGENOMIC oral metagenome ILLUMINA
    ## 5              PCR   METAGENOMIC oral metagenome ILLUMINA
    ## 6              PCR   METAGENOMIC oral metagenome ILLUMINA
    ##                                     project_name          ReleaseDate
    ## 1 krauter_salivary_microbiome_part_1_of_6_study2 2017-07-01T00:00:00Z
    ## 2 krauter_salivary_microbiome_part_1_of_6_study2 2017-07-01T00:00:00Z
    ## 3 krauter_salivary_microbiome_part_1_of_6_study2 2017-07-01T00:00:00Z
    ## 4 krauter_salivary_microbiome_part_1_of_6_study2 2017-07-01T00:00:00Z
    ## 5 krauter_salivary_microbiome_part_1_of_6_study2 2017-07-01T00:00:00Z
    ## 6 krauter_salivary_microbiome_part_1_of_6_study2 2017-07-01T00:00:00Z
    ##      Sample.Name Sample_name Sequencing_method SRA.Study Submitter_Id
    ## 1 SAMEA104079181  bg.s2.1007    Illumina MiSeq ERP023086   bg.s2.1007
    ## 2 SAMEA104079182  bg.s2.1022    Illumina MiSeq ERP023086   bg.s2.1022
    ## 3 SAMEA104079183  bg.s2.1024    Illumina MiSeq ERP023086   bg.s2.1024
    ## 4 SAMEA104079184  bg.s2.1027    Illumina MiSeq ERP023086   bg.s2.1027
    ## 5 SAMEA104079185  bg.s2.1039    Illumina MiSeq ERP023086   bg.s2.1039
    ## 6 SAMEA104079186  bg.s2.1040    Illumina MiSeq ERP023086   bg.s2.1040
    ##   ENA_first_public ena_last_update               INSDC_center_alias
    ## 1       2017-06-29      2017-05-16 UNIVERSITY OF COLORADO - BOULDER
    ## 2       2017-06-29      2017-05-16 UNIVERSITY OF COLORADO - BOULDER
    ## 3       2017-06-29      2017-05-16 UNIVERSITY OF COLORADO - BOULDER
    ## 4       2017-06-29      2017-05-16 UNIVERSITY OF COLORADO - BOULDER
    ## 5       2017-06-29      2017-05-16 UNIVERSITY OF COLORADO - BOULDER
    ## 6       2017-06-29      2017-05-16 UNIVERSITY OF COLORADO - BOULDER
    ##            create_date version ena_first_public ena_last_update.1
    ## 1 2017-07-01T16:55:00Z       1                                   
    ## 2 2017-07-01T16:55:00Z       1                                   
    ## 3 2017-07-01T16:55:00Z       1                                   
    ## 4 2017-07-01T16:55:00Z       1                                   
    ## 5 2017-07-01T16:56:00Z       1                                   
    ## 6 2017-07-01T16:56:00Z       1                                   
    ##   Scientific_Name
    ## 1                
    ## 2                
    ## 3                
    ## 4                
    ## 5                
    ## 6

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
  Region = "V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$environment_.material.,
  public_time = initial_meta$ReleaseDate,
  Status = NA,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "515F-806R"
metadata$title <- "Genetic influences on the human oral microbiome"
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
#metadata$Status<-"H"

keep_subjects <- unique(subset(metadata, Site == "saliva")$sample_name)
metadata <- subset(metadata, sample_name %in% keep_subjects)

metadata$matched_samples <- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 397"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 397"

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

    ##   BioProject        Run    sample_name Age BMI Sex Region       Platform
    ## 1 PRJEB20889 ERR1983891 SAMEA104079181  NA  NA  NA     V4 Illumina MiSeq
    ## 2 PRJEB20889 ERR1983892 SAMEA104079182  NA  NA  NA     V4 Illumina MiSeq
    ## 3 PRJEB20889 ERR1983893 SAMEA104079183  NA  NA  NA     V4 Illumina MiSeq
    ## 4 PRJEB20889 ERR1983894 SAMEA104079184  NA  NA  NA     V4 Illumina MiSeq
    ## 5 PRJEB20889 ERR1983895 SAMEA104079185  NA  NA  NA     V4 Illumina MiSeq
    ## 6 PRJEB20889 ERR1983896 SAMEA104079186  NA  NA  NA     V4 Illumina MiSeq
    ##   Country   Site          public_time PCR_primers
    ## 1     USA saliva 2017-07-01T00:00:00Z   515F-806R
    ## 2     USA saliva 2017-07-01T00:00:00Z   515F-806R
    ## 3     USA saliva 2017-07-01T00:00:00Z   515F-806R
    ## 4     USA saliva 2017-07-01T00:00:00Z   515F-806R
    ## 5     USA saliva 2017-07-01T00:00:00Z   515F-806R
    ## 6     USA saliva 2017-07-01T00:00:00Z   515F-806R
    ##                                             title Disease_study Status
    ## 1 Genetic influences on the human oral microbiome            NO     NA
    ## 2 Genetic influences on the human oral microbiome            NO     NA
    ## 3 Genetic influences on the human oral microbiome            NO     NA
    ## 4 Genetic influences on the human oral microbiome            NO     NA
    ## 5 Genetic influences on the human oral microbiome            NO     NA
    ## 6 Genetic influences on the human oral microbiome            NO     NA
    ##   matched_samples
    ## 1              NO
    ## 2              NO
    ## 3              NO
    ## 4              NO
    ## 5              NO
    ## 6              NO

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
