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
setwd("D:/studylearning/董彪师兄项目/任务/修改20250922/PRJNA575550")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen     Bases  BioProject    BioSample
    ## 1 SRR10226385   AMPLICON        602 175814702 PRJNA575550 SAMN12900412
    ## 2 SRR10226386   AMPLICON        602 153531070 PRJNA575550 SAMN12900313
    ## 3 SRR10226387   AMPLICON        602  96073782 PRJNA575550 SAMN12900329
    ## 4 SRR10226388   AMPLICON        602  81838288 PRJNA575550 SAMN12900349
    ## 5 SRR10226389   AMPLICON        602 102554312 PRJNA575550 SAMN12900432
    ## 6 SRR10226390   AMPLICON        602  95747498 PRJNA575550 SAMN12900428
    ##                         BioSampleModel     Bytes        Center.Name
    ## 1 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 124417547 UPPSALA UNIVERSITY
    ## 2 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 108726972 UPPSALA UNIVERSITY
    ## 3 MIMS.me,MIGS/MIMS/MIMARKS.human-oral  67939775 UPPSALA UNIVERSITY
    ## 4 MIMS.me,MIGS/MIMS/MIMARKS.human-oral  58380009 UPPSALA UNIVERSITY
    ## 5 MIMS.me,MIGS/MIMS/MIMARKS.human-oral  72657066 UPPSALA UNIVERSITY
    ## 6 MIMS.me,MIGS/MIMS/MIMARKS.human-oral  68227568 UPPSALA UNIVERSITY
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1      2012-06-19  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2      2012-01-23  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3      2012-09-13  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4      2013-02-12  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5      2013-01-09  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6      2014-06-17  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region env_broad_scale env_local_scale
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1      human_oral          saliva
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1      human_oral          saliva
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1      human_oral          saliva
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1      human_oral          saliva
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1      human_oral          saliva
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1      human_oral          saliva
    ##   env_medium Experiment geo_loc_name_country geo_loc_name_country_continent
    ## 1     saliva SRX6946109               Sweden                         Europe
    ## 2     saliva SRX6946108               Sweden                         Europe
    ## 3     saliva SRX6946107               Sweden                         Europe
    ## 4     saliva SRX6946106               Sweden                         Europe
    ## 5     saliva SRX6946105               Sweden                         Europe
    ## 6     saliva SRX6946104               Sweden                         Europe
    ##   geo_loc_name         HOST     Instrument isolate isolation_source
    ## 1       Sweden Homo_sapiens Illumina MiSeq     129      oral_cavity
    ## 2       Sweden Homo_sapiens Illumina MiSeq      28      oral_cavity
    ## 3       Sweden Homo_sapiens Illumina MiSeq      45      oral_cavity
    ## 4       Sweden Homo_sapiens Illumina MiSeq      66      oral_cavity
    ## 5       Sweden Homo_sapiens Illumina MiSeq     149      oral_cavity
    ## 6       Sweden Homo_sapiens Illumina MiSeq     145      oral_cavity
    ##           lat_lon                  Library.Name LibraryLayout LibrarySelection
    ## 1 59.30 N 18.09 E SAL210_S103_L001_R1_001.fastq        PAIRED              PCR
    ## 2 59.30 N 18.09 E    SAL11_S5_L001_R1_001.fastq        PAIRED              PCR
    ## 3 59.30 N 18.09 E   SAL21_S39_L001_R1_001.fastq        PAIRED              PCR
    ## 4 59.30 N 18.09 E SAL208_S102_L001_R1_001.fastq        PAIRED              PCR
    ## 5 59.30 N 18.09 E SAL206_S101_L001_R1_001.fastq        PAIRED              PCR
    ## 6 59.30 N 18.09 E SAL205_S100_L001_R1_001.fastq        PAIRED              PCR
    ##   LibrarySource              Organism Platform          ReleaseDate
    ## 1   METAGENOMIC human oral metagenome ILLUMINA 2019-10-04T00:00:00Z
    ## 2   METAGENOMIC human oral metagenome ILLUMINA 2019-10-04T00:00:00Z
    ## 3   METAGENOMIC human oral metagenome ILLUMINA 2019-10-04T00:00:00Z
    ## 4   METAGENOMIC human oral metagenome ILLUMINA 2019-10-04T00:00:00Z
    ## 5   METAGENOMIC human oral metagenome ILLUMINA 2019-10-04T00:00:00Z
    ## 6   METAGENOMIC human oral metagenome ILLUMINA 2019-10-04T00:00:00Z
    ##            create_date version Sample.Name SRA.Study
    ## 1 2019-10-03T13:22:00Z       1      SAL210 SRP224096
    ## 2 2019-10-03T13:22:00Z       1       SAL11 SRP224096
    ## 3 2019-10-03T13:22:00Z       1       SAL21 SRP224096
    ## 4 2019-10-03T13:22:00Z       1      SAL208 SRP224096
    ## 5 2019-10-03T13:22:00Z       1      SAL206 SRP224096
    ## 6 2019-10-03T13:22:00Z       1      SAL205 SRP224096

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
  Sex= NA,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$env_medium,
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
metadata$title <- "Oral Microbiota Development in Early Childhood"
metadata$Disease_study <- "NO"
metadata$Status <- NA_character_

#extract subject id by remove prefix of sample_name
metadata$sample_name <- gsub("SAL", "", metadata$sample_name)

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##   1  10 100 102 105 106 107 108 109  11 110 113 114 115 116 117 118 119  12 120 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 125 126 128 129  13 131 132 133 134 135 136 137 139  14 140 142 143 144 146 148 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 149  15 150 151 152 153 155 156 157 158 159 160 161 163 164 166 167  17 170 172 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 173 174 175 176 178 179 180 181 183 184 185 186 187 189 190 191 192 193 195 196 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 197 199  20 200 202 203 205 206 208  21 210 211 212 214 217  22 220 223 224 225 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 226 227 229  23 231 232 233 235 236 237 238  24 240 242 243 244 245 246  25  26 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  29  31  33  34  36  37  38  39   4  40  42  43  45  47  48   5  50  51  52  53 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  54  55  57  58  59   6  60  61  64  66  68  69  70  72  73  74  76  78  81  83 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  84  85  86  87  88  89   9  90  91  93  95  97  98 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1

``` r
table(metadata$Site)
```

    ## 
    ## saliva 
    ##    173

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
    ##   1  10 100 102 105 106 107 108 109  11 110 113 114 115 116 117 118 119  12 120 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 125 126 128 129  13 131 132 133 134 135 136 137 139  14 140 142 143 144 146 148 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 149  15 150 151 152 153 155 156 157 158 159 160 161 163 164 166 167  17 170 172 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 173 174 175 176 178 179 180 181 183 184 185 186 187 189 190 191 192 193 195 196 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 197 199  20 200 202 203 205 206 208  21 210 211 212 214 217  22 220 223 224 225 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 226 227 229  23 231 232 233 235 236 237 238  24 240 242 243 244 245 246  25  26 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  29  31  33  34  36  37  38  39   4  40  42  43  45  47  48   5  50  51  52  53 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  54  55  57  58  59   6  60  61  64  66  68  69  70  72  73  74  76  78  81  83 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  84  85  86  87  88  89   9  90  91  93  95  97  98 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1

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
    ## 173

``` r
# Identify samples for removal: single samples excluding saliva 
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & !(metadata$Site %in% c("saliva"))

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
    ##   1  10 100 102 105 106 107 108 109  11 110 113 114 115 116 117 118 119  12 120 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 125 126 128 129  13 131 132 133 134 135 136 137 139  14 140 142 143 144 146 148 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 149  15 150 151 152 153 155 156 157 158 159 160 161 163 164 166 167  17 170 172 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 173 174 175 176 178 179 180 181 183 184 185 186 187 189 190 191 192 193 195 196 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 197 199  20 200 202 203 205 206 208  21 210 211 212 214 217  22 220 223 224 225 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 226 227 229  23 231 232 233 235 236 237 238  24 240 242 243 244 245 246  25  26 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  29  31  33  34  36  37  38  39   4  40  42  43  45  47  48   5  50  51  52  53 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  54  55  57  58  59   6  60  61  64  66  68  69  70  72  73  74  76  78  81  83 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  84  85  86  87  88  89   9  90  91  93  95  97  98 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1

``` r
table(metadata_clean$Site)
```

    ## 
    ## saliva 
    ##    173

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ##  NO 
    ## 173

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 173"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 173"

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
    ## 1 PRJNA575550 SRR10226385         210  NA  NA  NA  V3-V4 Illumina MiSeq  Sweden
    ## 2 PRJNA575550 SRR10226386          11  NA  NA  NA  V3-V4 Illumina MiSeq  Sweden
    ## 3 PRJNA575550 SRR10226387          21  NA  NA  NA  V3-V4 Illumina MiSeq  Sweden
    ## 4 PRJNA575550 SRR10226388         208  NA  NA  NA  V3-V4 Illumina MiSeq  Sweden
    ## 5 PRJNA575550 SRR10226389         206  NA  NA  NA  V3-V4 Illumina MiSeq  Sweden
    ## 6 PRJNA575550 SRR10226390         205  NA  NA  NA  V3-V4 Illumina MiSeq  Sweden
    ##     Site          public_time PCR_primers
    ## 1 saliva 2019-10-04T00:00:00Z   341F-805R
    ## 2 saliva 2019-10-04T00:00:00Z   341F-805R
    ## 3 saliva 2019-10-04T00:00:00Z   341F-805R
    ## 4 saliva 2019-10-04T00:00:00Z   341F-805R
    ## 5 saliva 2019-10-04T00:00:00Z   341F-805R
    ## 6 saliva 2019-10-04T00:00:00Z   341F-805R
    ##                                            title Disease_study Status
    ## 1 Oral Microbiota Development in Early Childhood            NO   <NA>
    ## 2 Oral Microbiota Development in Early Childhood            NO   <NA>
    ## 3 Oral Microbiota Development in Early Childhood            NO   <NA>
    ## 4 Oral Microbiota Development in Early Childhood            NO   <NA>
    ## 5 Oral Microbiota Development in Early Childhood            NO   <NA>
    ## 6 Oral Microbiota Development in Early Childhood            NO   <NA>
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
