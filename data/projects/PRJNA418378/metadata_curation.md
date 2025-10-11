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
library(readxl)

# Set working directory (update this path to your local directory)
setwd("/home/dongbiao/word_embedding_microbiome/all_data/oral/global-oral-microbiome/data/projects/PRJNA418378")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv")

head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR6318146   AMPLICON        602 34167714 PRJNA418378 SAMN08042376
    ## 2 SRR6318147   AMPLICON        602 42858186 PRJNA418378 SAMN08042377
    ## 3 SRR6318148   AMPLICON        602 32042654 PRJNA418378 SAMN08042372
    ## 4 SRR6318149   AMPLICON        602 37423330 PRJNA418378 SAMN08042373
    ## 5 SRR6318151   AMPLICON        602 31651354 PRJNA418378 SAMN08042375
    ## 6 SRR6318152   AMPLICON        602 38912076 PRJNA418378 SAMN08042368
    ##                         BioSampleModel    Bytes
    ## 1 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 16278075
    ## 2 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 20629558
    ## 3 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 15404747
    ## 4 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 17649759
    ## 5 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 15199127
    ## 6 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 18750111
    ##                                                                   Center.Name
    ## 1 AFFILIATED HOSPITAL OF STOMATOLOGY,COLLEGE OF MEDICINE, ZHEJIANG UNIVERSITY
    ## 2 AFFILIATED HOSPITAL OF STOMATOLOGY,COLLEGE OF MEDICINE, ZHEJIANG UNIVERSITY
    ## 3 AFFILIATED HOSPITAL OF STOMATOLOGY,COLLEGE OF MEDICINE, ZHEJIANG UNIVERSITY
    ## 4 AFFILIATED HOSPITAL OF STOMATOLOGY,COLLEGE OF MEDICINE, ZHEJIANG UNIVERSITY
    ## 5 AFFILIATED HOSPITAL OF STOMATOLOGY,COLLEGE OF MEDICINE, ZHEJIANG UNIVERSITY
    ## 6 AFFILIATED HOSPITAL OF STOMATOLOGY,COLLEGE OF MEDICINE, ZHEJIANG UNIVERSITY
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1      2013-04-12  public         run.zq,sra         gs,ncbi,s3
    ## 2      2013-04-12  public         run.zq,sra         gs,ncbi,s3
    ## 3      2013-04-12  public         run.zq,sra         gs,ncbi,s3
    ## 4      2013-04-12  public         run.zq,sra         gs,ncbi,s3
    ## 5      2013-04-12  public         run.zq,sra         gs,ncbi,s3
    ## 6      2013-04-12  public         run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region      env_biome    env_feature
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 not applicable not applicable
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 not applicable not applicable
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 not applicable not applicable
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 not applicable not applicable
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 not applicable not applicable
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 not applicable not applicable
    ##     env_material Experiment geo_loc_name_country geo_loc_name_country_continent
    ## 1 not applicable SRX3418091                China                           Asia
    ## 2 not applicable SRX3418090                China                           Asia
    ## 3 not applicable SRX3418089                China                           Asia
    ## 4 not applicable SRX3418088                China                           Asia
    ## 5 not applicable SRX3418086                China                           Asia
    ## 6 not applicable SRX3418085                China                           Asia
    ##     geo_loc_name         HOST     Instrument          lat_lon Library.Name
    ## 1 China:Zhejiang Homo sapiens Illumina MiSeq 30.23 N 119.90 E         bt23
    ## 2 China:Zhejiang Homo sapiens Illumina MiSeq 30.23 N 119.91 E         bt25
    ## 3 China:Zhejiang Homo sapiens Illumina MiSeq 30.23 N 119.86 E         bt11
    ## 4 China:Zhejiang Homo sapiens Illumina MiSeq 30.23 N 119.87 E         bt13
    ## 5 China:Zhejiang Homo sapiens Illumina MiSeq 30.23 N 119.89 E         bt19
    ## 6 China:Zhejiang Homo sapiens Illumina MiSeq 30.23 N 119.82 E         bk01
    ##   LibraryLayout LibrarySelection LibrarySource              Organism Platform
    ## 1        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 2        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 3        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 4        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 5        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 6        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ##            ReleaseDate          create_date version Sample.Name SRA.Study
    ## 1 2021-11-20T00:00:00Z 2017-11-22T21:51:00Z       1        bt23 SRP125517
    ## 2 2021-11-20T00:00:00Z 2017-11-22T21:51:00Z       1        bt25 SRP125517
    ## 3 2021-11-20T00:00:00Z 2017-11-22T21:51:00Z       1        bt11 SRP125517
    ## 4 2021-11-20T00:00:00Z 2017-11-22T21:51:00Z       1        bt13 SRP125517
    ## 5 2021-11-20T00:00:00Z 2017-11-22T21:51:00Z       1        bt19 SRP125517
    ## 6 2021-11-20T00:00:00Z 2017-11-22T21:51:00Z       1        bk01 SRP125517

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
  sample_name = initial_meta$Sample.Name,
  Age = NA,
  BMI = NA,
  Sex=NA,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = "saliva",
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
metadata<-unique(metadata)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "338F-806R"
metadata$title <- "Dynamic Alterations in Salivary Microbiota Related to Dental Caries and Age in Preschool Children With Deciduous Dentition: A 2-Year Follow-Up Study"
metadata$Disease_study <- "YES"

# Assign disease status based on sample_name patterns
metadata$Status <- NA
metadata$matched_samples <- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 100"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 100"

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
metadata<-unique(metadata)
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##    BioProject        Run sample_name Age BMI Sex Region       Platform Country
    ## 1 PRJNA418378 SRR6318146        bt23  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 2 PRJNA418378 SRR6318147        bt25  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 3 PRJNA418378 SRR6318148        bt11  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 4 PRJNA418378 SRR6318149        bt13  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 5 PRJNA418378 SRR6318151        bt19  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 6 PRJNA418378 SRR6318152        bk01  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ##     Site          public_time PCR_primers
    ## 1 saliva 2021-11-20T00:00:00Z   338F-806R
    ## 2 saliva 2021-11-20T00:00:00Z   338F-806R
    ## 3 saliva 2021-11-20T00:00:00Z   338F-806R
    ## 4 saliva 2021-11-20T00:00:00Z   338F-806R
    ## 5 saliva 2021-11-20T00:00:00Z   338F-806R
    ## 6 saliva 2021-11-20T00:00:00Z   338F-806R
    ##                                                                                                                                                  title
    ## 1 Dynamic Alterations in Salivary Microbiota Related to Dental Caries and Age in Preschool Children With Deciduous Dentition: A 2-Year Follow-Up Study
    ## 2 Dynamic Alterations in Salivary Microbiota Related to Dental Caries and Age in Preschool Children With Deciduous Dentition: A 2-Year Follow-Up Study
    ## 3 Dynamic Alterations in Salivary Microbiota Related to Dental Caries and Age in Preschool Children With Deciduous Dentition: A 2-Year Follow-Up Study
    ## 4 Dynamic Alterations in Salivary Microbiota Related to Dental Caries and Age in Preschool Children With Deciduous Dentition: A 2-Year Follow-Up Study
    ## 5 Dynamic Alterations in Salivary Microbiota Related to Dental Caries and Age in Preschool Children With Deciduous Dentition: A 2-Year Follow-Up Study
    ## 6 Dynamic Alterations in Salivary Microbiota Related to Dental Caries and Age in Preschool Children With Deciduous Dentition: A 2-Year Follow-Up Study
    ##   Disease_study Status matched_samples
    ## 1           YES     NA              NO
    ## 2           YES     NA              NO
    ## 3           YES     NA              NO
    ## 4           YES     NA              NO
    ## 5           YES     NA              NO
    ## 6           YES     NA              NO

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
