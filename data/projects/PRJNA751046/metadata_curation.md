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
setwd("D:/studylearning/董彪师兄项目/任务/张辰飞/PRJNA751046")

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

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR21474442   AMPLICON        602 17958262 PRJNA751046 SAMN30722499
    ## 2 SRR21474443   AMPLICON        602 30789290 PRJNA751046 SAMN30722498
    ## 3 SRR21474450   AMPLICON        602 36920058 PRJNA751046 SAMN30722509
    ## 4 SRR21474451   AMPLICON        602 28743694 PRJNA751046 SAMN30722508
    ## 5 SRR21484724   AMPLICON        602 56516362 PRJNA751046 SAMN30725542
    ## 6 SRR21484725   AMPLICON        602 26359774 PRJNA751046 SAMN30725541
    ##                                BioSampleModel    Bytes
    ## 1 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 10385393
    ## 2 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 17344813
    ## 3 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 20145434
    ## 4 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 16293854
    ## 5 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 32135564
    ## 6 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 15129962
    ##                                             Center.Name Collection_Date Consent
    ## 1 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-02-10  public
    ## 2 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-02-05  public
    ## 3 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-04-10  public
    ## 4 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-03-27  public
    ## 5 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-08-15  public
    ## 6 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-08-08  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##         env_broad_scale env_local_scale env_medium  Experiment
    ## 1           oral cavity  not applicable       <NA> SRX17478091
    ## 2           oral cavity  not applicable       <NA> SRX17478090
    ## 3           oral cavity  not applicable       <NA> SRX17478083
    ## 4           oral cavity  not applicable       <NA> SRX17478082
    ## 5 Human Oral microbiome            Oral    missing SRX17488295
    ## 6 Human Oral microbiome            Oral    missing SRX17488294
    ##   geo_loc_name_country geo_loc_name_country_continent geo_loc_name         HOST
    ## 1                India                           Asia        India Homo sapiens
    ## 2                India                           Asia        India Homo sapiens
    ## 3                India                           Asia        India Homo sapiens
    ## 4                India                           Asia        India Homo sapiens
    ## 5                India                           Asia        India Homo sapiens
    ## 6                India                           Asia        India Homo sapiens
    ##       Instrument lat_lon Library.Name LibraryLayout LibrarySelection
    ## 1 Illumina MiSeq    <NA>         S207        PAIRED              PCR
    ## 2 Illumina MiSeq    <NA>         S195        PAIRED              PCR
    ## 3 Illumina MiSeq    <NA>         S244        PAIRED              PCR
    ## 4 Illumina MiSeq    <NA>         S232        PAIRED              PCR
    ## 5 Illumina MiSeq missing         S214        PAIRED              PCR
    ## 6 Illumina MiSeq missing         S202        PAIRED              PCR
    ##   LibrarySource              Organism Platform          ReleaseDate
    ## 1   METAGENOMIC human oral metagenome ILLUMINA 2022-09-08T00:00:00Z
    ## 2   METAGENOMIC human oral metagenome ILLUMINA 2022-09-08T00:00:00Z
    ## 3   METAGENOMIC human oral metagenome ILLUMINA 2022-09-08T00:00:00Z
    ## 4   METAGENOMIC human oral metagenome ILLUMINA 2022-09-08T00:00:00Z
    ## 5   METAGENOMIC human oral metagenome ILLUMINA 2022-09-09T00:00:00Z
    ## 6   METAGENOMIC human oral metagenome ILLUMINA 2022-09-09T00:00:00Z
    ##            create_date version Sample.Name SRA.Study Host_Age
    ## 1 2022-09-08T17:38:00Z       1         C21 SRP330565       NA
    ## 2 2022-09-08T18:18:00Z       1         C20 SRP330565       NA
    ## 3 2022-09-08T17:22:00Z       1         C36 SRP330565       NA
    ## 4 2022-09-08T17:09:00Z       1         C34 SRP330565       NA
    ## 5 2022-09-09T02:43:00Z       1          O7 SRP330565       NA
    ## 6 2022-09-09T02:15:00Z       1          O1 SRP330565       NA

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
  Sex= NA,
  Region = "V6-V8",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = "oral rinse",
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "969F-1406R"
metadata$title <- "Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata <- metadata %>% 
  mutate(Status = case_when(
    grepl("C",sample_name) ~ "HC",
    grepl("T",sample_name) ~ "Tobacco",
    TRUE ~ "OSCC"
  ))


#extract subject id by remove prefix of sample_name
#metadata$sample_name <- gsub("T|CR|CS|C|N|S", "", metadata$sample_name)

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##   C1  C10  C13  C14  C15  C17  C18  C19   C2  C20  C21  C22  C23  C24  C27  C29 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  C31  C32  C33  C34  C36  C37  C39   C4  C43  C44  C45  C46  C47  C48  C49   C5 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  C50  C51  C52  C53   C6   O1  O23  O25  O28  O31  O32  O33  O34  O36  O37  O38 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  O39  O40  O47  O49  O50  O51  O53  O54  O55  O57  O60  O61  O62  O65   O7   O9 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## OC13 OC15 OC16 OC17  OC2 OC21 OC22  OC4  OC6   T1  T10  T11  T12  T13  T14  T15 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  T16  T18  T19   T2  T20  T21  T22  T23  T24  T25  T26  T27  T28  T29   T3  T30 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  T31  T32  T33  T34  T35  T36  T37  T39   T4  T41  T43  T44  T45   T6   T7   T9 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1

``` r
table(metadata$Site)
```

    ## 
    ## oral rinse 
    ##        112

``` r
# Calculate sampling frequency per subject
sample_counts <- table(metadata$sample_name)
cat("\nNumber of samples per subject:\n")
```

    ## 
    ## Number of samples per subject:

``` r
print(sample_counts)
```

    ## 
    ##   C1  C10  C13  C14  C15  C17  C18  C19   C2  C20  C21  C22  C23  C24  C27  C29 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  C31  C32  C33  C34  C36  C37  C39   C4  C43  C44  C45  C46  C47  C48  C49   C5 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  C50  C51  C52  C53   C6   O1  O23  O25  O28  O31  O32  O33  O34  O36  O37  O38 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  O39  O40  O47  O49  O50  O51  O53  O54  O55  O57  O60  O61  O62  O65   O7   O9 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## OC13 OC15 OC16 OC17  OC2 OC21 OC22  OC4  OC6   T1  T10  T11  T12  T13  T14  T15 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  T16  T18  T19   T2  T20  T21  T22  T23  T24  T25  T26  T27  T28  T29   T3  T30 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  T31  T32  T33  T34  T35  T36  T37  T39   T4  T41  T43  T44  T45   T6   T7   T9 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1

``` r
# Create matched_samples column
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")

cat("\nMatched samples distribution:\n")
```

    ## 
    ## Matched samples distribution:

``` r
table(metadata$matched_samples)
```

    ## 
    ##  NO 
    ## 112

``` r
# Identify samples for removal: single samples excluding oral rinse
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & !(metadata$Site %in% c("oral rinse"))

cat("\nNumber of samples to be removed:", sum(samples_to_remove), "\n")
```

    ## 
    ## Number of samples to be removed: 0

``` r
# Remove identified samples
metadata_clean <- metadata[!samples_to_remove, ]

cat("\nFinal sample distribution after cleaning:\n")
```

    ## 
    ## Final sample distribution after cleaning:

``` r
table(metadata_clean$sample_name)
```

    ## 
    ##   C1  C10  C13  C14  C15  C17  C18  C19   C2  C20  C21  C22  C23  C24  C27  C29 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  C31  C32  C33  C34  C36  C37  C39   C4  C43  C44  C45  C46  C47  C48  C49   C5 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  C50  C51  C52  C53   C6   O1  O23  O25  O28  O31  O32  O33  O34  O36  O37  O38 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  O39  O40  O47  O49  O50  O51  O53  O54  O55  O57  O60  O61  O62  O65   O7   O9 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## OC13 OC15 OC16 OC17  OC2 OC21 OC22  OC4  OC6   T1  T10  T11  T12  T13  T14  T15 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  T16  T18  T19   T2  T20  T21  T22  T23  T24  T25  T26  T27  T28  T29   T3  T30 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  T31  T32  T33  T34  T35  T36  T37  T39   T4  T41  T43  T44  T45   T6   T7   T9 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1

``` r
table(metadata_clean$Site)
```

    ## 
    ## oral rinse 
    ##        112

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ##  NO 
    ## 112

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 112"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 112"

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

    ##    BioProject         Run sample_name Age BMI Sex Region       Platform Country
    ## 1 PRJNA751046 SRR21474442         C21  NA  NA  NA  V6-V8 Illumina MiSeq   India
    ## 2 PRJNA751046 SRR21474443         C20  NA  NA  NA  V6-V8 Illumina MiSeq   India
    ## 3 PRJNA751046 SRR21474450         C36  NA  NA  NA  V6-V8 Illumina MiSeq   India
    ## 4 PRJNA751046 SRR21474451         C34  NA  NA  NA  V6-V8 Illumina MiSeq   India
    ## 5 PRJNA751046 SRR21484724          O7  NA  NA  NA  V6-V8 Illumina MiSeq   India
    ## 6 PRJNA751046 SRR21484725          O1  NA  NA  NA  V6-V8 Illumina MiSeq   India
    ##         Site          public_time PCR_primers
    ## 1 oral rinse 2022-09-08T00:00:00Z  969F-1406R
    ## 2 oral rinse 2022-09-08T00:00:00Z  969F-1406R
    ## 3 oral rinse 2022-09-08T00:00:00Z  969F-1406R
    ## 4 oral rinse 2022-09-08T00:00:00Z  969F-1406R
    ## 5 oral rinse 2022-09-09T00:00:00Z  969F-1406R
    ## 6 oral rinse 2022-09-09T00:00:00Z  969F-1406R
    ##                                                                                                                            title
    ## 1 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study
    ## 2 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study
    ## 3 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study
    ## 4 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study
    ## 5 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study
    ## 6 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study
    ##   Disease_study Status matched_samples
    ## 1           YES     HC              NO
    ## 2           YES     HC              NO
    ## 3           YES     HC              NO
    ## 4           YES     HC              NO
    ## 5           YES   OSCC              NO
    ## 6           YES   OSCC              NO

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
