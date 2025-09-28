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

setwd("D:/PRJNA714139")
initial_meta<-read.csv("SraRunTable .csv",header = T)
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen     Bases  BioProject    BioSample
    ## 1 SRR13964430   AMPLICON        500 127206000 PRJNA714139 SAMN18290700
    ## 2 SRR13964431   AMPLICON        500 120534000 PRJNA714139 SAMN18290699
    ## 3 SRR13964432   AMPLICON        500 121112000 PRJNA714139 SAMN18290770
    ## 4 SRR13964433   AMPLICON        500 136951000 PRJNA714139 SAMN18290769
    ## 5 SRR13964434   AMPLICON        500 145208500 PRJNA714139 SAMN18290768
    ## 6 SRR13964435   AMPLICON        500 106234000 PRJNA714139 SAMN18290767
    ##                                BioSampleModel    Bytes            Center.Name
    ## 1 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 85578999 UNIVERSITY OF MANITOBA
    ## 2 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 78055890 UNIVERSITY OF MANITOBA
    ## 3 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 79941846 UNIVERSITY OF MANITOBA
    ## 4 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 87591392 UNIVERSITY OF MANITOBA
    ## 5 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 95413506 UNIVERSITY OF MANITOBA
    ## 6 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 67416666 UNIVERSITY OF MANITOBA
    ##         Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1 2017-12-13/2018-08-14  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2 2017-12-13/2018-08-14  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3 2017-12-13/2018-08-14  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4 2017-12-13/2018-08-14  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5 2017-12-13/2018-08-14  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6 2017-12-13/2018-08-14  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region        env_broad_scale
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 mouth [UBERON:0000165]
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 mouth [UBERON:0000165]
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 mouth [UBERON:0000165]
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 mouth [UBERON:0000165]
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 mouth [UBERON:0000165]
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 mouth [UBERON:0000165]
    ##              env_local_scale             env_medium  Experiment
    ## 1 carious teeth [MP:0030101] Saliva [ENVO:02000036] SRX10342501
    ## 2 carious teeth [MP:0030101] Saliva [ENVO:02000036] SRX10342500
    ## 3 carious teeth [MP:0030101] Saliva [ENVO:02000036] SRX10342499
    ## 4 carious teeth [MP:0030101] Saliva [ENVO:02000036] SRX10342498
    ## 5 carious teeth [MP:0030101] Saliva [ENVO:02000036] SRX10342497
    ## 6 carious teeth [MP:0030101] Saliva [ENVO:02000036] SRX10342496
    ##   geo_loc_name_country geo_loc_name_country_continent     geo_loc_name host_age
    ## 1               Canada                  North America Canada: Winnipeg 2.583333
    ## 2               Canada                  North America Canada: Winnipeg 3.083333
    ## 3               Canada                  North America Canada: Winnipeg 3.500000
    ## 4               Canada                  North America Canada: Winnipeg 4.166667
    ## 5               Canada                  North America Canada: Winnipeg 5.833333
    ## 6               Canada                  North America Canada: Winnipeg 3.333333
    ##   Host_disease         HOST host_sex host_subject_id     Instrument
    ## 1        S-ECC Homo sapiens   female             S18 Illumina MiSeq
    ## 2        S-ECC Homo sapiens     male             S17 Illumina MiSeq
    ## 3  Caries-free Homo sapiens   female              S9 Illumina MiSeq
    ## 4  Caries-free Homo sapiens     male             S80 Illumina MiSeq
    ## 5        S-ECC Homo sapiens     male              S8 Illumina MiSeq
    ## 6  Caries-free Homo sapiens     male             S79 Illumina MiSeq
    ##   isolation_source         lat_lon Library.Name LibraryLayout LibrarySelection
    ## 1       Oral swabs 49.89 N 97.14 W          S18        PAIRED              PCR
    ## 2       Oral swabs 49.89 N 97.14 W          S17        PAIRED              PCR
    ## 3       Oral swabs 49.89 N 97.14 W           S9        PAIRED              PCR
    ## 4       Oral swabs 49.89 N 97.14 W          S80        PAIRED              PCR
    ## 5       Oral swabs 49.89 N 97.14 W           S8        PAIRED              PCR
    ## 6       Oral swabs 49.89 N 97.14 W          S79        PAIRED              PCR
    ##   LibrarySource              Organism Platform          ReleaseDate
    ## 1   METAGENOMIC human oral metagenome ILLUMINA 2021-03-16T00:00:00Z
    ## 2   METAGENOMIC human oral metagenome ILLUMINA 2021-03-16T00:00:00Z
    ## 3   METAGENOMIC human oral metagenome ILLUMINA 2021-03-16T00:00:00Z
    ## 4   METAGENOMIC human oral metagenome ILLUMINA 2021-03-16T00:00:00Z
    ## 5   METAGENOMIC human oral metagenome ILLUMINA 2021-03-16T00:00:00Z
    ## 6   METAGENOMIC human oral metagenome ILLUMINA 2021-03-16T00:00:00Z
    ##            create_date version samp_collect_device Sample.Name SRA.Study
    ## 1 2021-03-15T15:47:00Z       1           Oral swab         S18 SRP310674
    ## 2 2021-03-15T15:47:00Z       1           Oral swab         S17 SRP310674
    ## 3 2021-03-15T15:47:00Z       1           Oral swab          S9 SRP310674
    ## 4 2021-03-15T15:47:00Z       1           Oral swab         S80 SRP310674
    ## 5 2021-03-15T15:47:00Z       1           Oral swab          S8 SRP310674
    ## 6 2021-03-15T15:47:00Z       1           Oral swab         S79 SRP310674

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
  Age = initial_meta$host_age,
  Sex = initial_meta$host_sex,
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Host_disease= initial_meta$Host_disease,
  Site = initial_meta$env_medium,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
metadata <- metadata %>%
  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
metadata$BMI <- NA_character_
metadata$PCR_primers <- "515F-806R"
metadata$Region <- "V4"
metadata$title <-"Characterization of Supragingival Plaque and Oral Swab Microbiomes in Children With Severe Early Childhood Caries"
metadata$Disease_study <- "YES"
metadata<-metadata%>%
  mutate(Status=case_when(grepl("S-ECC",Host_disease)~"DC",grepl("Caries-free",Host_disease)~"HC"))
metadata$matched_samples <-"NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 80"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 80"

# Step 5: Visualize the Metadata

``` r
knitr::include_graphics("Metadata_Summary.png")
```

![](Metadata_Summary.png)<!-- -->

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

    ##    BioProject         Run sample_name      Age  BMI Sex Region       Platform
    ## 1 PRJNA714139 SRR13964430         S18 2.583333 <NA>   F     V4 Illumina MiSeq
    ## 2 PRJNA714139 SRR13964431         S17 3.083333 <NA>   M     V4 Illumina MiSeq
    ## 3 PRJNA714139 SRR13964432          S9 3.500000 <NA>   F     V4 Illumina MiSeq
    ## 4 PRJNA714139 SRR13964433         S80 4.166667 <NA>   M     V4 Illumina MiSeq
    ## 5 PRJNA714139 SRR13964434          S8 5.833333 <NA>   M     V4 Illumina MiSeq
    ## 6 PRJNA714139 SRR13964435         S79 3.333333 <NA>   M     V4 Illumina MiSeq
    ##   Country                   Site          public_time PCR_primers
    ## 1  Canada Saliva [ENVO:02000036] 2021-03-16T00:00:00Z   515F-806R
    ## 2  Canada Saliva [ENVO:02000036] 2021-03-16T00:00:00Z   515F-806R
    ## 3  Canada Saliva [ENVO:02000036] 2021-03-16T00:00:00Z   515F-806R
    ## 4  Canada Saliva [ENVO:02000036] 2021-03-16T00:00:00Z   515F-806R
    ## 5  Canada Saliva [ENVO:02000036] 2021-03-16T00:00:00Z   515F-806R
    ## 6  Canada Saliva [ENVO:02000036] 2021-03-16T00:00:00Z   515F-806R
    ##                                                                                                               title
    ## 1 Characterization of Supragingival Plaque and Oral Swab Microbiomes in Children With Severe Early Childhood Caries
    ## 2 Characterization of Supragingival Plaque and Oral Swab Microbiomes in Children With Severe Early Childhood Caries
    ## 3 Characterization of Supragingival Plaque and Oral Swab Microbiomes in Children With Severe Early Childhood Caries
    ## 4 Characterization of Supragingival Plaque and Oral Swab Microbiomes in Children With Severe Early Childhood Caries
    ## 5 Characterization of Supragingival Plaque and Oral Swab Microbiomes in Children With Severe Early Childhood Caries
    ## 6 Characterization of Supragingival Plaque and Oral Swab Microbiomes in Children With Severe Early Childhood Caries
    ##   Disease_study Status matched_samples
    ## 1           YES     DC              NO
    ## 2           YES     DC              NO
    ## 3           YES     HC              NO
    ## 4           YES     HC              NO
    ## 5           YES     DC              NO
    ## 6           YES     HC              NO
