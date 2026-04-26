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
library(stringr)

# Set working directory (update this path to your local directory)
setwd("/home/hby/oral_project/PRJEB39064")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = T)

# Combine the two metadata files
# meta <- rbind(initial_meta, initial_meta1)

# Remove entries non-oral samples
# initial_meta <- subset(initial_meta,isolate %in% c("saliva","tongue coat"))

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
# nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen     Bases BioProject    BioSample    Bytes Center.Name Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1 ERR4296609   AMPLICON        602  99615950 PRJEB39064 SAMEA7016816 63294440        NCKU        2013/6/4  public   run.zq,sra,fastq         gs,ncbi,s3
    ## 2 ERR4296610   AMPLICON        602 127905134 PRJEB39064 SAMEA7016817 81852204        NCKU        2013/6/5  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3 ERR4296611   AMPLICON        602 107340814 PRJEB39064 SAMEA7016818 67172464        NCKU       2013/8/23  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4 ERR4296612   AMPLICON        602  86025198 PRJEB39064 SAMEA7016819 54966027        NCKU       2013/9/12  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5 ERR4296613   AMPLICON        602 102600666 PRJEB39064 SAMEA7016820 67349647        NCKU      2013/10/15  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6 ERR4296614   AMPLICON        602 105632940 PRJEB39064 SAMEA7016821 66770558        NCKU      2013/10/23  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region ENA.FIRST.PUBLIC..run. ENA_first_public ENA.LAST.UPDATE..run. ENA_last_update environment_.biome. environment_.feature. environment_material
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1              2021/9/10        2021/9/10             2021/9/10       2021/9/10               Human       Oral Metagenome               Saliva
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1              2021/9/10        2021/9/10             2021/9/10       2021/9/10               Human       Oral Metagenome               Saliva
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1              2021/9/10        2021/9/10             2021/9/10       2021/9/10               Human       Oral Metagenome               Saliva
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1              2021/9/10        2021/9/10             2021/9/10       2021/9/10               Human       Oral Metagenome               Saliva
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1              2021/9/10        2021/9/10             2021/9/10       2021/9/10               Human       Oral Metagenome               Saliva
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1              2021/9/10        2021/9/10             2021/9/10       2021/9/10               Human       Oral Metagenome               Saliva
    ##   Experiment  external_id geo_loc_name_country geo_loc_name_country_continent geographic_location_.country_and.or_sea. geographic_location_.latitude. geographic_location_.longitude.
    ## 1 ERX4245364 SAMEA7016816                China                           Asia                                    China                       23.28943                        120.3254
    ## 2 ERX4245365 SAMEA7016817                China                           Asia                                    China                       23.28943                        120.3254
    ## 3 ERX4245366 SAMEA7016818                China                           Asia                                    China                       23.28943                        120.3254
    ## 4 ERX4245367 SAMEA7016819                China                           Asia                                    China                       23.28943                        120.3254
    ## 5 ERX4245368 SAMEA7016820                China                           Asia                                    China                       23.28943                        120.3254
    ## 6 ERX4245369 SAMEA7016821                China                           Asia                                    China                       23.28943                        120.3254
    ##   Host_Age host_disease_status host_sex human_oral_environmental_package             INSDC_center_alias              INSDC_center_name   INSDC_first_public    INSDC_last_update
    ## 1       64                   H     male                       human-oral NATIONAL CHENG KUNG UNIVERSITY NATIONAL CHENG KUNG UNIVERSITY 2021-09-10T04:15:15Z 2021-09-10T04:15:15Z
    ## 2       50                   H     male                       human-oral NATIONAL CHENG KUNG UNIVERSITY NATIONAL CHENG KUNG UNIVERSITY 2021-09-10T04:15:15Z 2021-09-10T04:15:15Z
    ## 3       55                   H     male                       human-oral NATIONAL CHENG KUNG UNIVERSITY NATIONAL CHENG KUNG UNIVERSITY 2021-09-10T04:15:15Z 2021-09-10T04:15:15Z
    ## 4       48                   H     male                       human-oral NATIONAL CHENG KUNG UNIVERSITY NATIONAL CHENG KUNG UNIVERSITY 2021-09-10T04:15:15Z 2021-09-10T04:15:15Z
    ## 5       49                   H     male                       human-oral NATIONAL CHENG KUNG UNIVERSITY NATIONAL CHENG KUNG UNIVERSITY 2021-09-10T04:15:15Z 2021-09-10T04:15:15Z
    ## 6       53                   H     male                       human-oral NATIONAL CHENG KUNG UNIVERSITY NATIONAL CHENG KUNG UNIVERSITY 2021-09-10T04:15:15Z 2021-09-10T04:15:15Z
    ##   INSDC_status     Instrument investigation_type Library.Name LibraryLayout LibrarySelection LibrarySource        Organism Platform           project_name          ReleaseDate
    ## 1       public Illumina MiSeq         metagenome  unspecified        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA Oral Cancer Microbiome 2021-09-30T00:00:00Z
    ## 2       public Illumina MiSeq         metagenome  unspecified        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA Oral Cancer Microbiome 2021-09-30T00:00:00Z
    ## 3       public Illumina MiSeq         metagenome  unspecified        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA Oral Cancer Microbiome 2021-09-30T00:00:00Z
    ## 4       public Illumina MiSeq         metagenome  unspecified        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA Oral Cancer Microbiome 2021-09-30T00:00:00Z
    ## 5       public Illumina MiSeq         metagenome  unspecified        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA Oral Cancer Microbiome 2021-09-30T00:00:00Z
    ## 6       public Illumina MiSeq         metagenome  unspecified        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA Oral Cancer Microbiome 2021-09-30T00:00:00Z
    ##            create_date version  Sample.Name Sample_name Sequencing_method SRA.Study Submitter_Id
    ## 1 2021-10-10T17:07:00Z       1 SAMEA7016816  Normal-132    Illumina MiSeq ERP122533   Normal-132
    ## 2 2021-10-10T17:11:00Z       1 SAMEA7016817  Normal-136    Illumina MiSeq ERP122533   Normal-136
    ## 3 2021-10-10T17:07:00Z       1 SAMEA7016818  Normal-193    Illumina MiSeq ERP122533   Normal-193
    ## 4 2021-10-10T17:07:00Z       1 SAMEA7016819  Normal-219    Illumina MiSeq ERP122533   Normal-219
    ## 5 2021-10-10T17:10:00Z       1 SAMEA7016820  Normal-251    Illumina MiSeq ERP122533   Normal-251
    ## 6 2021-10-10T17:09:00Z       1 SAMEA7016821  Normal-252    Illumina MiSeq ERP122533   Normal-252

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
  sample_name = initial_meta$Sample_name,
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$environment_material,
  public_time = initial_meta$ReleaseDate,
  Status = initial_meta$host_disease_status,
  Sex = initial_meta$host_sex,
  Age = initial_meta$Host_Age,
  BMI = NA,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Site = case_when(
      Site == "Saliva" ~ "saliva",
      TRUE ~ as.character(Status)
    )
  )
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Region <- "V3-V4"
metadata$PCR_primers <- "341F-806R"
metadata$title <- "Taxonomic and Functional Dysregulation in Salivary Microbiomes During Oral Carcinogenesis"
metadata$Disease_study <- "YES"
metadata$matched_samples <- "No"

cat("\nFinal sample distribution:\n")
```

    ## 
    ## Final sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ## Normal-132 Normal-136 Normal-193 Normal-219 Normal-251 Normal-252 Normal-284 Normal-307 Normal-315 Normal-330 Normal-335 Normal-378 Normal-381 Normal-412 Normal-427 Normal-485 
    ##          1          1          1          1          1          1          1          1          1          1          1          1          1          1          1          1 
    ## Normal-487 Normal-550 Normal-682 Normal-724 Normal-728 Normal-729 Normal-739 Normal-770 Normal-788 Normal-813 Normal-818    NRVH-12   NRVH-121   NRVH-139   NRVH-208   NRVH-226 
    ##          1          1          1          1          1          1          1          1          1          1          1          1          1          1          1          1 
    ##   NRVH-282   NRVH-364   NRVH-451   NRVH-481   NRVH-544   NRVH-631   NRVH-653     NRVH-7   NRVH-700   NRVH-714   NRVH-742   NRVH-803   NRVH-895   NRVH-913   NRVH-932   NRVH-986 
    ##          1          1          1          1          1          1          1          1          1          1          1          1          1          1          1          1 
    ##    OSCC-10  OSCC-1003  OSCC-1017  OSCC-1020  OSCC-1056   OSCC-106   OSCC-117   OSCC-146   OSCC-151   OSCC-158   OSCC-159   OSCC-160    OSCC-18   OSCC-204   OSCC-230    OSCC-38 
    ##          1          1          1          1          1          1          1          1          1          1          1          1          1          1          1          1 
    ##   OSCC-433   OSCC-489   OSCC-519   OSCC-579   OSCC-587    OSCC-62   OSCC-673    OSCC-87   OSCC-887   OSCC-982   OSCC-997 
    ##          1          1          1          1          1          1          1          1          1          1          1

``` r
table(metadata$Site)
```

    ## 
    ## saliva 
    ##     75

``` r
table(metadata$matched_samples)
```

    ## 
    ## No 
    ## 75

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 75"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 75"

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
  "PCR_primers", "title", "Disease_study", "Status", "matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##   BioProject        Run sample_name Age BMI  Sex Region       Platform Country   Site          public_time PCR_primers
    ## 1 PRJEB39064 ERR4296609  Normal-132  64  NA male  V3-V4 Illumina MiSeq   China saliva 2021-09-30T00:00:00Z   341F-806R
    ## 2 PRJEB39064 ERR4296610  Normal-136  50  NA male  V3-V4 Illumina MiSeq   China saliva 2021-09-30T00:00:00Z   341F-806R
    ## 3 PRJEB39064 ERR4296611  Normal-193  55  NA male  V3-V4 Illumina MiSeq   China saliva 2021-09-30T00:00:00Z   341F-806R
    ## 4 PRJEB39064 ERR4296612  Normal-219  48  NA male  V3-V4 Illumina MiSeq   China saliva 2021-09-30T00:00:00Z   341F-806R
    ## 5 PRJEB39064 ERR4296613  Normal-251  49  NA male  V3-V4 Illumina MiSeq   China saliva 2021-09-30T00:00:00Z   341F-806R
    ## 6 PRJEB39064 ERR4296614  Normal-252  53  NA male  V3-V4 Illumina MiSeq   China saliva 2021-09-30T00:00:00Z   341F-806R
    ##                                                                                       title Disease_study Status matched_samples
    ## 1 Taxonomic and Functional Dysregulation in Salivary Microbiomes During Oral Carcinogenesis           YES      H              No
    ## 2 Taxonomic and Functional Dysregulation in Salivary Microbiomes During Oral Carcinogenesis           YES      H              No
    ## 3 Taxonomic and Functional Dysregulation in Salivary Microbiomes During Oral Carcinogenesis           YES      H              No
    ## 4 Taxonomic and Functional Dysregulation in Salivary Microbiomes During Oral Carcinogenesis           YES      H              No
    ## 5 Taxonomic and Functional Dysregulation in Salivary Microbiomes During Oral Carcinogenesis           YES      H              No
    ## 6 Taxonomic and Functional Dysregulation in Salivary Microbiomes During Oral Carcinogenesis           YES      H              No

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
