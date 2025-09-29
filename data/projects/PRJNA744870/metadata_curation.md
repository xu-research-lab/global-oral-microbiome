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
setwd("D:/studylearning/董彪师兄项目/任务/张辰飞/PRJNA744870")

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
    ## 1 SRR15719926   AMPLICON        602 23042754 PRJNA744870 SAMN21215907
    ## 2 SRR15719937   AMPLICON        602 23526160 PRJNA744870 SAMN21215897
    ## 3 SRR15719951   AMPLICON        602 28818944 PRJNA744870 SAMN21215884
    ## 4 SRR15719961   AMPLICON        602 22373932 PRJNA744870 SAMN21215816
    ## 5 SRR15719966   AMPLICON        602  4212796 PRJNA744870 SAMN21215815
    ## 6 SRR15719973   AMPLICON        602 36488424 PRJNA744870 SAMN21215951
    ##                BioSampleModel    Bytes                         Center.Name
    ## 1 Metagenome or environmental 15416009 THE CHINESE UNIVERSITY OF HONG KONG
    ## 2 Metagenome or environmental 15389328 THE CHINESE UNIVERSITY OF HONG KONG
    ## 3 Metagenome or environmental 17930945 THE CHINESE UNIVERSITY OF HONG KONG
    ## 4 Metagenome or environmental 13334903 THE CHINESE UNIVERSITY OF HONG KONG
    ## 5 Metagenome or environmental  2130157 THE CHINESE UNIVERSITY OF HONG KONG
    ## 6 Metagenome or environmental 22725477 THE CHINESE UNIVERSITY OF HONG KONG
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region  Experiment geo_loc_name_country
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX12015501            Hong Kong
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX12015490            Hong Kong
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX12015476            Hong Kong
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX12015466            Hong Kong
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX12015461            Hong Kong
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX12015834            Hong Kong
    ##   geo_loc_name_country_continent geo_loc_name         HOST     Instrument
    ## 1                           Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 2                           Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 3                           Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 4                           Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 5                           Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 6                           Asia    Hong Kong Homo sapiens Illumina MiSeq
    ##   lat_lon             Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1 missing    Health26-Oral-16SV3V4        PAIRED              PCR       GENOMIC
    ## 2 missing    Health16-Oral-16SV3V4        PAIRED              PCR       GENOMIC
    ## 3 missing    Health03-Oral-16SV3V4        PAIRED              PCR       GENOMIC
    ## 4 missing Patient11-PostM3-16SV3V4        PAIRED              PCR       GENOMIC
    ## 5 missing Patient10-PostM3-16SV3V4        PAIRED              PCR       GENOMIC
    ## 6 missing    Health70-Oral-16SV3V4        PAIRED              PCR       GENOMIC
    ##          Organism Platform          ReleaseDate          create_date version
    ## 1 oral metagenome ILLUMINA 2021-09-28T00:00:00Z 2021-09-04T12:36:00Z       1
    ## 2 oral metagenome ILLUMINA 2021-09-28T00:00:00Z 2021-09-04T12:36:00Z       1
    ## 3 oral metagenome ILLUMINA 2021-09-28T00:00:00Z 2021-09-04T12:37:00Z       1
    ## 4 oral metagenome ILLUMINA 2021-09-28T00:00:00Z 2021-09-04T12:36:00Z       1
    ## 5 oral metagenome ILLUMINA 2021-09-28T00:00:00Z 2021-09-04T12:36:00Z       1
    ## 6 oral metagenome ILLUMINA 2021-09-28T00:00:00Z 2021-09-04T12:36:00Z       1
    ##   Sample.Name SRA.Study   isolate sample_type   strain isolation_source
    ## 1    Health26 SRP328065  Health26  oral rinse Bacteria               NA
    ## 2    Health16 SRP328065  Health16  oral rinse Bacteria               NA
    ## 3    Health03 SRP328065  Health03  oral rinse Bacteria               NA
    ## 4   Patient11 SRP328065 Patient11  oral rinse Bacteria               NA
    ## 5   Patient10 SRP328065 Patient10  oral rinse Bacteria               NA
    ## 6    Health70 SRP328065  Health70  oral rinse Bacteria               NA

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
  sample_name = initial_meta$Library.Name,
  Age = NA,
  BMI = NA,
  Sex= NA,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = "China",
  Site = initial_meta$sample_type,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "341F-806R"
metadata$title <- "Restoration of the Oral Microbiota After Surgery for Head and Neck Squamous Cell Carcinoma Is Associated With Patient Outcomes"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata <- metadata %>% 
  mutate(Status = case_when(
    grepl("Patient",sample_name) ~ "HNSCC",
    grepl("Health",sample_name) ~ "HC"
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
    ##     Health01-Oral-16SV3V4     Health02-Oral-16SV3V4     Health03-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health04-Oral-16SV3V4     Health05-Oral-16SV3V4     Health06-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health07-Oral-16SV3V4     Health08-Oral-16SV3V4     Health09-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health10-Oral-16SV3V4     Health11-Oral-16SV3V4     Health12-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health13-Oral-16SV3V4     Health14-Oral-16SV3V4     Health15-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health16-Oral-16SV3V4     Health17-Oral-16SV3V4     Health18-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health19-Oral-16SV3V4     Health20-Oral-16SV3V4     Health21-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health22-Oral-16SV3V4     Health23-Oral-16SV3V4     Health24-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health25-Oral-16SV3V4     Health26-Oral-16SV3V4     Health27-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health28-Oral-16SV3V4     Health29-Oral-16SV3V4     Health30-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health31-Oral-16SV3V4     Health32-Oral-16SV3V4     Health33-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health34-Oral-16SV3V4     Health35-Oral-16SV3V4     Health36-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health37-Oral-16SV3V4     Health38-Oral-16SV3V4     Health39-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health40-Oral-16SV3V4     Health41-Oral-16SV3V4     Health42-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health43-Oral-16SV3V4     Health44-Oral-16SV3V4     Health45-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health46-Oral-16SV3V4     Health47-Oral-16SV3V4     Health48-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health49-Oral-16SV3V4     Health50-Oral-16SV3V4     Health51-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health52-Oral-16SV3V4     Health53-Oral-16SV3V4     Health54-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health55-Oral-16SV3V4     Health56-Oral-16SV3V4     Health57-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health58-Oral-16SV3V4     Health59-Oral-16SV3V4     Health60-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health61-Oral-16SV3V4     Health62-Oral-16SV3V4     Health63-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health64-Oral-16SV3V4     Health65-Oral-16SV3V4     Health66-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health67-Oral-16SV3V4     Health68-Oral-16SV3V4     Health69-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health70-Oral-16SV3V4     Health71-Oral-16SV3V4     Health72-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health73-Oral-16SV3V4     Health74-Oral-16SV3V4     Health75-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health76-Oral-16SV3V4  Patient01-PostM1-16SV3V4  Patient01-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient01-PostM6-16SV3V4     Patient01-Pre-16SV3V4 Patient01R-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ## Patient01R-PostM3-16SV3V4 Patient01R-PostM6-16SV3V4    Patient01R-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient02-PostM1-16SV3V4  Patient02-PostM3-16SV3V4  Patient02-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient02-Pre-16SV3V4  Patient03-PostM1-16SV3V4  Patient03-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient03-PostM6-16SV3V4     Patient03-Pre-16SV3V4  Patient04-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient04-PostM3-16SV3V4  Patient04-PostM6-16SV3V4     Patient04-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient05-PostM1-16SV3V4  Patient05-PostM3-16SV3V4  Patient05-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient05-Pre-16SV3V4  Patient06-PostM1-16SV3V4  Patient06-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient06-PostM6-16SV3V4     Patient06-Pre-16SV3V4  Patient07-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient07-PostM3-16SV3V4  Patient07-PostM6-16SV3V4     Patient07-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient08-PostM1-16SV3V4  Patient08-PostM3-16SV3V4  Patient08-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient08-Pre-16SV3V4  Patient09-PostM1-16SV3V4  Patient09-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient09-PostM6-16SV3V4     Patient09-Pre-16SV3V4  Patient10-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient10-PostM3-16SV3V4  Patient10-PostM6-16SV3V4     Patient10-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient11-PostM1-16SV3V4  Patient11-PostM3-16SV3V4  Patient11-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient11-Pre-16SV3V4  Patient12-PostM1-16SV3V4  Patient12-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient12-PostM6-16SV3V4     Patient12-Pre-16SV3V4  Patient13-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient13-PostM3-16SV3V4  Patient13-PostM6-16SV3V4     Patient13-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient14-PostM1-16SV3V4  Patient14-PostM3-16SV3V4  Patient14-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient14-Pre-16SV3V4  Patient15-PostM1-16SV3V4  Patient15-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient15-PostM6-16SV3V4     Patient15-Pre-16SV3V4  Patient16-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient16-PostM3-16SV3V4  Patient16-PostM6-16SV3V4     Patient16-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient17-PostM1-16SV3V4  Patient17-PostM3-16SV3V4  Patient17-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient17-Pre-16SV3V4  Patient18-PostM1-16SV3V4  Patient18-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient18-PostM6-16SV3V4     Patient18-Pre-16SV3V4  Patient19-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient19-PostM3-16SV3V4  Patient19-PostM6-16SV3V4     Patient19-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient20-PostM1-16SV3V4  Patient20-PostM3-16SV3V4  Patient20-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient20-Pre-16SV3V4  Patient21-PostM1-16SV3V4  Patient21-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient21-PostM6-16SV3V4     Patient21-Pre-16SV3V4  Patient22-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient22-PostM3-16SV3V4  Patient22-PostM6-16SV3V4     Patient22-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient23-PostM1-16SV3V4  Patient23-PostM3-16SV3V4  Patient23-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient23-Pre-16SV3V4  Patient24-PostM1-16SV3V4  Patient24-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient24-PostM6-16SV3V4     Patient24-Pre-16SV3V4  Patient25-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient25-PostM3-16SV3V4  Patient25-PostM6-16SV3V4     Patient25-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient26-PostM1-16SV3V4  Patient26-PostM3-16SV3V4  Patient26-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient26-Pre-16SV3V4  Patient27-PostM1-16SV3V4  Patient27-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient27-PostM6-16SV3V4     Patient27-Pre-16SV3V4  Patient28-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient28-PostM3-16SV3V4  Patient28-PostM6-16SV3V4     Patient28-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient29-PostM1-16SV3V4  Patient29-PostM3-16SV3V4  Patient29-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient29-Pre-16SV3V4  Patient30-PostM1-16SV3V4  Patient30-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient30-PostM6-16SV3V4     Patient30-Pre-16SV3V4  Patient31-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient31-PostM3-16SV3V4  Patient31-PostM6-16SV3V4     Patient31-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient32-PostM1-16SV3V4  Patient32-PostM3-16SV3V4  Patient32-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient32-Pre-16SV3V4  Patient33-PostM1-16SV3V4  Patient33-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient33-PostM6-16SV3V4     Patient33-Pre-16SV3V4  Patient34-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient34-PostM3-16SV3V4  Patient34-PostM6-16SV3V4     Patient34-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient35-PostM1-16SV3V4  Patient35-PostM3-16SV3V4  Patient35-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient35-Pre-16SV3V4  Patient36-PostM1-16SV3V4  Patient36-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient36-PostM6-16SV3V4     Patient36-Pre-16SV3V4  Patient37-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient37-PostM3-16SV3V4  Patient37-PostM6-16SV3V4     Patient37-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient38-PostM1-16SV3V4  Patient38-PostM3-16SV3V4  Patient38-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient38-Pre-16SV3V4  Patient39-PostM1-16SV3V4  Patient39-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient39-PostM6-16SV3V4     Patient39-Pre-16SV3V4  Patient40-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient40-PostM3-16SV3V4  Patient40-PostM6-16SV3V4     Patient40-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient41-PostM1-16SV3V4  Patient41-PostM3-16SV3V4  Patient41-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient41-Pre-16SV3V4  Patient42-PostM1-16SV3V4  Patient42-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient42-PostM6-16SV3V4     Patient42-Pre-16SV3V4  Patient43-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient43-PostM3-16SV3V4  Patient43-PostM6-16SV3V4     Patient43-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient44-PostM1-16SV3V4  Patient44-PostM3-16SV3V4  Patient44-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient44-Pre-16SV3V4  Patient45-PostM1-16SV3V4  Patient45-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient45-PostM6-16SV3V4     Patient45-Pre-16SV3V4  Patient46-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient46-PostM3-16SV3V4  Patient46-PostM6-16SV3V4     Patient46-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient47-PostM1-16SV3V4  Patient47-PostM3-16SV3V4  Patient47-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient47-Pre-16SV3V4  Patient48-PostM1-16SV3V4  Patient48-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient48-PostM6-16SV3V4     Patient48-Pre-16SV3V4  Patient49-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient49-PostM3-16SV3V4  Patient49-PostM6-16SV3V4     Patient49-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient50-PostM1-16SV3V4  Patient50-PostM3-16SV3V4  Patient50-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient50-Pre-16SV3V4  Patient51-PostM1-16SV3V4  Patient51-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient51-PostM6-16SV3V4     Patient51-Pre-16SV3V4  Patient52-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient52-PostM3-16SV3V4  Patient52-PostM6-16SV3V4     Patient52-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient53-PostM1-16SV3V4  Patient53-PostM3-16SV3V4  Patient53-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient53-Pre-16SV3V4  Patient54-PostM1-16SV3V4  Patient54-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient54-PostM6-16SV3V4     Patient54-Pre-16SV3V4  Patient55-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient55-PostM3-16SV3V4  Patient55-PostM6-16SV3V4     Patient55-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient56-PostM1-16SV3V4  Patient56-PostM3-16SV3V4  Patient56-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient56-Pre-16SV3V4  Patient57-PostM1-16SV3V4  Patient57-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient57-PostM6-16SV3V4     Patient57-Pre-16SV3V4  Patient58-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient58-PostM3-16SV3V4  Patient58-PostM6-16SV3V4     Patient58-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient59-PostM1-16SV3V4  Patient59-PostM3-16SV3V4  Patient59-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient59-Pre-16SV3V4  Patient60-PostM1-16SV3V4  Patient60-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient60-PostM6-16SV3V4     Patient60-Pre-16SV3V4  Patient61-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient61-PostM3-16SV3V4  Patient61-PostM6-16SV3V4     Patient61-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient62-PostM1-16SV3V4  Patient62-PostM3-16SV3V4  Patient62-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient62-Pre-16SV3V4  Patient63-PostM1-16SV3V4  Patient63-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient63-PostM6-16SV3V4     Patient63-Pre-16SV3V4  Patient64-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient64-PostM3-16SV3V4  Patient64-PostM6-16SV3V4     Patient64-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient65-PostM1-16SV3V4  Patient65-PostM3-16SV3V4  Patient65-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient65-Pre-16SV3V4  Patient66-PostM1-16SV3V4  Patient66-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient66-PostM6-16SV3V4     Patient66-Pre-16SV3V4  Patient67-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient67-PostM3-16SV3V4  Patient67-PostM6-16SV3V4     Patient67-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient68-PostM1-16SV3V4  Patient68-PostM3-16SV3V4  Patient68-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient68-Pre-16SV3V4  Patient69-PostM1-16SV3V4  Patient69-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient69-PostM6-16SV3V4     Patient69-Pre-16SV3V4  Patient70-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient70-PostM3-16SV3V4  Patient70-PostM6-16SV3V4     Patient70-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient71-PostM1-16SV3V4  Patient71-PostM3-16SV3V4  Patient71-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient71-Pre-16SV3V4  Patient72-PostM1-16SV3V4  Patient72-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient72-PostM6-16SV3V4     Patient72-Pre-16SV3V4  Patient73-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient73-PostM3-16SV3V4  Patient73-PostM6-16SV3V4     Patient73-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient74-PostM1-16SV3V4  Patient74-PostM3-16SV3V4  Patient74-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient74-Pre-16SV3V4  Patient75-PostM1-16SV3V4  Patient75-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient75-PostM6-16SV3V4     Patient75-Pre-16SV3V4  Patient76-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient76-PostM3-16SV3V4  Patient76-PostM6-16SV3V4     Patient76-Pre-16SV3V4 
    ##                         1                         1                         1

``` r
table(metadata$Site)
```

    ## 
    ## oral rinse 
    ##        384

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
    ##     Health01-Oral-16SV3V4     Health02-Oral-16SV3V4     Health03-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health04-Oral-16SV3V4     Health05-Oral-16SV3V4     Health06-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health07-Oral-16SV3V4     Health08-Oral-16SV3V4     Health09-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health10-Oral-16SV3V4     Health11-Oral-16SV3V4     Health12-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health13-Oral-16SV3V4     Health14-Oral-16SV3V4     Health15-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health16-Oral-16SV3V4     Health17-Oral-16SV3V4     Health18-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health19-Oral-16SV3V4     Health20-Oral-16SV3V4     Health21-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health22-Oral-16SV3V4     Health23-Oral-16SV3V4     Health24-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health25-Oral-16SV3V4     Health26-Oral-16SV3V4     Health27-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health28-Oral-16SV3V4     Health29-Oral-16SV3V4     Health30-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health31-Oral-16SV3V4     Health32-Oral-16SV3V4     Health33-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health34-Oral-16SV3V4     Health35-Oral-16SV3V4     Health36-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health37-Oral-16SV3V4     Health38-Oral-16SV3V4     Health39-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health40-Oral-16SV3V4     Health41-Oral-16SV3V4     Health42-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health43-Oral-16SV3V4     Health44-Oral-16SV3V4     Health45-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health46-Oral-16SV3V4     Health47-Oral-16SV3V4     Health48-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health49-Oral-16SV3V4     Health50-Oral-16SV3V4     Health51-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health52-Oral-16SV3V4     Health53-Oral-16SV3V4     Health54-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health55-Oral-16SV3V4     Health56-Oral-16SV3V4     Health57-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health58-Oral-16SV3V4     Health59-Oral-16SV3V4     Health60-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health61-Oral-16SV3V4     Health62-Oral-16SV3V4     Health63-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health64-Oral-16SV3V4     Health65-Oral-16SV3V4     Health66-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health67-Oral-16SV3V4     Health68-Oral-16SV3V4     Health69-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health70-Oral-16SV3V4     Health71-Oral-16SV3V4     Health72-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health73-Oral-16SV3V4     Health74-Oral-16SV3V4     Health75-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health76-Oral-16SV3V4  Patient01-PostM1-16SV3V4  Patient01-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient01-PostM6-16SV3V4     Patient01-Pre-16SV3V4 Patient01R-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ## Patient01R-PostM3-16SV3V4 Patient01R-PostM6-16SV3V4    Patient01R-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient02-PostM1-16SV3V4  Patient02-PostM3-16SV3V4  Patient02-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient02-Pre-16SV3V4  Patient03-PostM1-16SV3V4  Patient03-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient03-PostM6-16SV3V4     Patient03-Pre-16SV3V4  Patient04-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient04-PostM3-16SV3V4  Patient04-PostM6-16SV3V4     Patient04-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient05-PostM1-16SV3V4  Patient05-PostM3-16SV3V4  Patient05-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient05-Pre-16SV3V4  Patient06-PostM1-16SV3V4  Patient06-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient06-PostM6-16SV3V4     Patient06-Pre-16SV3V4  Patient07-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient07-PostM3-16SV3V4  Patient07-PostM6-16SV3V4     Patient07-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient08-PostM1-16SV3V4  Patient08-PostM3-16SV3V4  Patient08-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient08-Pre-16SV3V4  Patient09-PostM1-16SV3V4  Patient09-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient09-PostM6-16SV3V4     Patient09-Pre-16SV3V4  Patient10-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient10-PostM3-16SV3V4  Patient10-PostM6-16SV3V4     Patient10-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient11-PostM1-16SV3V4  Patient11-PostM3-16SV3V4  Patient11-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient11-Pre-16SV3V4  Patient12-PostM1-16SV3V4  Patient12-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient12-PostM6-16SV3V4     Patient12-Pre-16SV3V4  Patient13-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient13-PostM3-16SV3V4  Patient13-PostM6-16SV3V4     Patient13-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient14-PostM1-16SV3V4  Patient14-PostM3-16SV3V4  Patient14-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient14-Pre-16SV3V4  Patient15-PostM1-16SV3V4  Patient15-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient15-PostM6-16SV3V4     Patient15-Pre-16SV3V4  Patient16-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient16-PostM3-16SV3V4  Patient16-PostM6-16SV3V4     Patient16-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient17-PostM1-16SV3V4  Patient17-PostM3-16SV3V4  Patient17-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient17-Pre-16SV3V4  Patient18-PostM1-16SV3V4  Patient18-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient18-PostM6-16SV3V4     Patient18-Pre-16SV3V4  Patient19-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient19-PostM3-16SV3V4  Patient19-PostM6-16SV3V4     Patient19-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient20-PostM1-16SV3V4  Patient20-PostM3-16SV3V4  Patient20-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient20-Pre-16SV3V4  Patient21-PostM1-16SV3V4  Patient21-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient21-PostM6-16SV3V4     Patient21-Pre-16SV3V4  Patient22-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient22-PostM3-16SV3V4  Patient22-PostM6-16SV3V4     Patient22-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient23-PostM1-16SV3V4  Patient23-PostM3-16SV3V4  Patient23-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient23-Pre-16SV3V4  Patient24-PostM1-16SV3V4  Patient24-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient24-PostM6-16SV3V4     Patient24-Pre-16SV3V4  Patient25-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient25-PostM3-16SV3V4  Patient25-PostM6-16SV3V4     Patient25-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient26-PostM1-16SV3V4  Patient26-PostM3-16SV3V4  Patient26-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient26-Pre-16SV3V4  Patient27-PostM1-16SV3V4  Patient27-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient27-PostM6-16SV3V4     Patient27-Pre-16SV3V4  Patient28-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient28-PostM3-16SV3V4  Patient28-PostM6-16SV3V4     Patient28-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient29-PostM1-16SV3V4  Patient29-PostM3-16SV3V4  Patient29-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient29-Pre-16SV3V4  Patient30-PostM1-16SV3V4  Patient30-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient30-PostM6-16SV3V4     Patient30-Pre-16SV3V4  Patient31-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient31-PostM3-16SV3V4  Patient31-PostM6-16SV3V4     Patient31-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient32-PostM1-16SV3V4  Patient32-PostM3-16SV3V4  Patient32-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient32-Pre-16SV3V4  Patient33-PostM1-16SV3V4  Patient33-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient33-PostM6-16SV3V4     Patient33-Pre-16SV3V4  Patient34-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient34-PostM3-16SV3V4  Patient34-PostM6-16SV3V4     Patient34-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient35-PostM1-16SV3V4  Patient35-PostM3-16SV3V4  Patient35-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient35-Pre-16SV3V4  Patient36-PostM1-16SV3V4  Patient36-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient36-PostM6-16SV3V4     Patient36-Pre-16SV3V4  Patient37-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient37-PostM3-16SV3V4  Patient37-PostM6-16SV3V4     Patient37-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient38-PostM1-16SV3V4  Patient38-PostM3-16SV3V4  Patient38-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient38-Pre-16SV3V4  Patient39-PostM1-16SV3V4  Patient39-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient39-PostM6-16SV3V4     Patient39-Pre-16SV3V4  Patient40-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient40-PostM3-16SV3V4  Patient40-PostM6-16SV3V4     Patient40-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient41-PostM1-16SV3V4  Patient41-PostM3-16SV3V4  Patient41-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient41-Pre-16SV3V4  Patient42-PostM1-16SV3V4  Patient42-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient42-PostM6-16SV3V4     Patient42-Pre-16SV3V4  Patient43-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient43-PostM3-16SV3V4  Patient43-PostM6-16SV3V4     Patient43-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient44-PostM1-16SV3V4  Patient44-PostM3-16SV3V4  Patient44-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient44-Pre-16SV3V4  Patient45-PostM1-16SV3V4  Patient45-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient45-PostM6-16SV3V4     Patient45-Pre-16SV3V4  Patient46-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient46-PostM3-16SV3V4  Patient46-PostM6-16SV3V4     Patient46-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient47-PostM1-16SV3V4  Patient47-PostM3-16SV3V4  Patient47-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient47-Pre-16SV3V4  Patient48-PostM1-16SV3V4  Patient48-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient48-PostM6-16SV3V4     Patient48-Pre-16SV3V4  Patient49-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient49-PostM3-16SV3V4  Patient49-PostM6-16SV3V4     Patient49-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient50-PostM1-16SV3V4  Patient50-PostM3-16SV3V4  Patient50-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient50-Pre-16SV3V4  Patient51-PostM1-16SV3V4  Patient51-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient51-PostM6-16SV3V4     Patient51-Pre-16SV3V4  Patient52-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient52-PostM3-16SV3V4  Patient52-PostM6-16SV3V4     Patient52-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient53-PostM1-16SV3V4  Patient53-PostM3-16SV3V4  Patient53-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient53-Pre-16SV3V4  Patient54-PostM1-16SV3V4  Patient54-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient54-PostM6-16SV3V4     Patient54-Pre-16SV3V4  Patient55-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient55-PostM3-16SV3V4  Patient55-PostM6-16SV3V4     Patient55-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient56-PostM1-16SV3V4  Patient56-PostM3-16SV3V4  Patient56-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient56-Pre-16SV3V4  Patient57-PostM1-16SV3V4  Patient57-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient57-PostM6-16SV3V4     Patient57-Pre-16SV3V4  Patient58-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient58-PostM3-16SV3V4  Patient58-PostM6-16SV3V4     Patient58-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient59-PostM1-16SV3V4  Patient59-PostM3-16SV3V4  Patient59-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient59-Pre-16SV3V4  Patient60-PostM1-16SV3V4  Patient60-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient60-PostM6-16SV3V4     Patient60-Pre-16SV3V4  Patient61-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient61-PostM3-16SV3V4  Patient61-PostM6-16SV3V4     Patient61-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient62-PostM1-16SV3V4  Patient62-PostM3-16SV3V4  Patient62-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient62-Pre-16SV3V4  Patient63-PostM1-16SV3V4  Patient63-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient63-PostM6-16SV3V4     Patient63-Pre-16SV3V4  Patient64-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient64-PostM3-16SV3V4  Patient64-PostM6-16SV3V4     Patient64-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient65-PostM1-16SV3V4  Patient65-PostM3-16SV3V4  Patient65-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient65-Pre-16SV3V4  Patient66-PostM1-16SV3V4  Patient66-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient66-PostM6-16SV3V4     Patient66-Pre-16SV3V4  Patient67-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient67-PostM3-16SV3V4  Patient67-PostM6-16SV3V4     Patient67-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient68-PostM1-16SV3V4  Patient68-PostM3-16SV3V4  Patient68-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient68-Pre-16SV3V4  Patient69-PostM1-16SV3V4  Patient69-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient69-PostM6-16SV3V4     Patient69-Pre-16SV3V4  Patient70-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient70-PostM3-16SV3V4  Patient70-PostM6-16SV3V4     Patient70-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient71-PostM1-16SV3V4  Patient71-PostM3-16SV3V4  Patient71-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient71-Pre-16SV3V4  Patient72-PostM1-16SV3V4  Patient72-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient72-PostM6-16SV3V4     Patient72-Pre-16SV3V4  Patient73-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient73-PostM3-16SV3V4  Patient73-PostM6-16SV3V4     Patient73-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient74-PostM1-16SV3V4  Patient74-PostM3-16SV3V4  Patient74-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient74-Pre-16SV3V4  Patient75-PostM1-16SV3V4  Patient75-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient75-PostM6-16SV3V4     Patient75-Pre-16SV3V4  Patient76-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient76-PostM3-16SV3V4  Patient76-PostM6-16SV3V4     Patient76-Pre-16SV3V4 
    ##                         1                         1                         1

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
    ## 384

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
    ##     Health01-Oral-16SV3V4     Health02-Oral-16SV3V4     Health03-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health04-Oral-16SV3V4     Health05-Oral-16SV3V4     Health06-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health07-Oral-16SV3V4     Health08-Oral-16SV3V4     Health09-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health10-Oral-16SV3V4     Health11-Oral-16SV3V4     Health12-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health13-Oral-16SV3V4     Health14-Oral-16SV3V4     Health15-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health16-Oral-16SV3V4     Health17-Oral-16SV3V4     Health18-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health19-Oral-16SV3V4     Health20-Oral-16SV3V4     Health21-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health22-Oral-16SV3V4     Health23-Oral-16SV3V4     Health24-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health25-Oral-16SV3V4     Health26-Oral-16SV3V4     Health27-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health28-Oral-16SV3V4     Health29-Oral-16SV3V4     Health30-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health31-Oral-16SV3V4     Health32-Oral-16SV3V4     Health33-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health34-Oral-16SV3V4     Health35-Oral-16SV3V4     Health36-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health37-Oral-16SV3V4     Health38-Oral-16SV3V4     Health39-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health40-Oral-16SV3V4     Health41-Oral-16SV3V4     Health42-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health43-Oral-16SV3V4     Health44-Oral-16SV3V4     Health45-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health46-Oral-16SV3V4     Health47-Oral-16SV3V4     Health48-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health49-Oral-16SV3V4     Health50-Oral-16SV3V4     Health51-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health52-Oral-16SV3V4     Health53-Oral-16SV3V4     Health54-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health55-Oral-16SV3V4     Health56-Oral-16SV3V4     Health57-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health58-Oral-16SV3V4     Health59-Oral-16SV3V4     Health60-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health61-Oral-16SV3V4     Health62-Oral-16SV3V4     Health63-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health64-Oral-16SV3V4     Health65-Oral-16SV3V4     Health66-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health67-Oral-16SV3V4     Health68-Oral-16SV3V4     Health69-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health70-Oral-16SV3V4     Health71-Oral-16SV3V4     Health72-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health73-Oral-16SV3V4     Health74-Oral-16SV3V4     Health75-Oral-16SV3V4 
    ##                         1                         1                         1 
    ##     Health76-Oral-16SV3V4  Patient01-PostM1-16SV3V4  Patient01-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient01-PostM6-16SV3V4     Patient01-Pre-16SV3V4 Patient01R-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ## Patient01R-PostM3-16SV3V4 Patient01R-PostM6-16SV3V4    Patient01R-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient02-PostM1-16SV3V4  Patient02-PostM3-16SV3V4  Patient02-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient02-Pre-16SV3V4  Patient03-PostM1-16SV3V4  Patient03-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient03-PostM6-16SV3V4     Patient03-Pre-16SV3V4  Patient04-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient04-PostM3-16SV3V4  Patient04-PostM6-16SV3V4     Patient04-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient05-PostM1-16SV3V4  Patient05-PostM3-16SV3V4  Patient05-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient05-Pre-16SV3V4  Patient06-PostM1-16SV3V4  Patient06-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient06-PostM6-16SV3V4     Patient06-Pre-16SV3V4  Patient07-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient07-PostM3-16SV3V4  Patient07-PostM6-16SV3V4     Patient07-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient08-PostM1-16SV3V4  Patient08-PostM3-16SV3V4  Patient08-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient08-Pre-16SV3V4  Patient09-PostM1-16SV3V4  Patient09-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient09-PostM6-16SV3V4     Patient09-Pre-16SV3V4  Patient10-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient10-PostM3-16SV3V4  Patient10-PostM6-16SV3V4     Patient10-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient11-PostM1-16SV3V4  Patient11-PostM3-16SV3V4  Patient11-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient11-Pre-16SV3V4  Patient12-PostM1-16SV3V4  Patient12-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient12-PostM6-16SV3V4     Patient12-Pre-16SV3V4  Patient13-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient13-PostM3-16SV3V4  Patient13-PostM6-16SV3V4     Patient13-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient14-PostM1-16SV3V4  Patient14-PostM3-16SV3V4  Patient14-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient14-Pre-16SV3V4  Patient15-PostM1-16SV3V4  Patient15-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient15-PostM6-16SV3V4     Patient15-Pre-16SV3V4  Patient16-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient16-PostM3-16SV3V4  Patient16-PostM6-16SV3V4     Patient16-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient17-PostM1-16SV3V4  Patient17-PostM3-16SV3V4  Patient17-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient17-Pre-16SV3V4  Patient18-PostM1-16SV3V4  Patient18-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient18-PostM6-16SV3V4     Patient18-Pre-16SV3V4  Patient19-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient19-PostM3-16SV3V4  Patient19-PostM6-16SV3V4     Patient19-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient20-PostM1-16SV3V4  Patient20-PostM3-16SV3V4  Patient20-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient20-Pre-16SV3V4  Patient21-PostM1-16SV3V4  Patient21-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient21-PostM6-16SV3V4     Patient21-Pre-16SV3V4  Patient22-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient22-PostM3-16SV3V4  Patient22-PostM6-16SV3V4     Patient22-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient23-PostM1-16SV3V4  Patient23-PostM3-16SV3V4  Patient23-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient23-Pre-16SV3V4  Patient24-PostM1-16SV3V4  Patient24-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient24-PostM6-16SV3V4     Patient24-Pre-16SV3V4  Patient25-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient25-PostM3-16SV3V4  Patient25-PostM6-16SV3V4     Patient25-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient26-PostM1-16SV3V4  Patient26-PostM3-16SV3V4  Patient26-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient26-Pre-16SV3V4  Patient27-PostM1-16SV3V4  Patient27-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient27-PostM6-16SV3V4     Patient27-Pre-16SV3V4  Patient28-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient28-PostM3-16SV3V4  Patient28-PostM6-16SV3V4     Patient28-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient29-PostM1-16SV3V4  Patient29-PostM3-16SV3V4  Patient29-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient29-Pre-16SV3V4  Patient30-PostM1-16SV3V4  Patient30-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient30-PostM6-16SV3V4     Patient30-Pre-16SV3V4  Patient31-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient31-PostM3-16SV3V4  Patient31-PostM6-16SV3V4     Patient31-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient32-PostM1-16SV3V4  Patient32-PostM3-16SV3V4  Patient32-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient32-Pre-16SV3V4  Patient33-PostM1-16SV3V4  Patient33-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient33-PostM6-16SV3V4     Patient33-Pre-16SV3V4  Patient34-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient34-PostM3-16SV3V4  Patient34-PostM6-16SV3V4     Patient34-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient35-PostM1-16SV3V4  Patient35-PostM3-16SV3V4  Patient35-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient35-Pre-16SV3V4  Patient36-PostM1-16SV3V4  Patient36-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient36-PostM6-16SV3V4     Patient36-Pre-16SV3V4  Patient37-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient37-PostM3-16SV3V4  Patient37-PostM6-16SV3V4     Patient37-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient38-PostM1-16SV3V4  Patient38-PostM3-16SV3V4  Patient38-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient38-Pre-16SV3V4  Patient39-PostM1-16SV3V4  Patient39-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient39-PostM6-16SV3V4     Patient39-Pre-16SV3V4  Patient40-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient40-PostM3-16SV3V4  Patient40-PostM6-16SV3V4     Patient40-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient41-PostM1-16SV3V4  Patient41-PostM3-16SV3V4  Patient41-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient41-Pre-16SV3V4  Patient42-PostM1-16SV3V4  Patient42-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient42-PostM6-16SV3V4     Patient42-Pre-16SV3V4  Patient43-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient43-PostM3-16SV3V4  Patient43-PostM6-16SV3V4     Patient43-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient44-PostM1-16SV3V4  Patient44-PostM3-16SV3V4  Patient44-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient44-Pre-16SV3V4  Patient45-PostM1-16SV3V4  Patient45-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient45-PostM6-16SV3V4     Patient45-Pre-16SV3V4  Patient46-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient46-PostM3-16SV3V4  Patient46-PostM6-16SV3V4     Patient46-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient47-PostM1-16SV3V4  Patient47-PostM3-16SV3V4  Patient47-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient47-Pre-16SV3V4  Patient48-PostM1-16SV3V4  Patient48-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient48-PostM6-16SV3V4     Patient48-Pre-16SV3V4  Patient49-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient49-PostM3-16SV3V4  Patient49-PostM6-16SV3V4     Patient49-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient50-PostM1-16SV3V4  Patient50-PostM3-16SV3V4  Patient50-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient50-Pre-16SV3V4  Patient51-PostM1-16SV3V4  Patient51-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient51-PostM6-16SV3V4     Patient51-Pre-16SV3V4  Patient52-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient52-PostM3-16SV3V4  Patient52-PostM6-16SV3V4     Patient52-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient53-PostM1-16SV3V4  Patient53-PostM3-16SV3V4  Patient53-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient53-Pre-16SV3V4  Patient54-PostM1-16SV3V4  Patient54-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient54-PostM6-16SV3V4     Patient54-Pre-16SV3V4  Patient55-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient55-PostM3-16SV3V4  Patient55-PostM6-16SV3V4     Patient55-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient56-PostM1-16SV3V4  Patient56-PostM3-16SV3V4  Patient56-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient56-Pre-16SV3V4  Patient57-PostM1-16SV3V4  Patient57-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient57-PostM6-16SV3V4     Patient57-Pre-16SV3V4  Patient58-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient58-PostM3-16SV3V4  Patient58-PostM6-16SV3V4     Patient58-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient59-PostM1-16SV3V4  Patient59-PostM3-16SV3V4  Patient59-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient59-Pre-16SV3V4  Patient60-PostM1-16SV3V4  Patient60-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient60-PostM6-16SV3V4     Patient60-Pre-16SV3V4  Patient61-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient61-PostM3-16SV3V4  Patient61-PostM6-16SV3V4     Patient61-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient62-PostM1-16SV3V4  Patient62-PostM3-16SV3V4  Patient62-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient62-Pre-16SV3V4  Patient63-PostM1-16SV3V4  Patient63-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient63-PostM6-16SV3V4     Patient63-Pre-16SV3V4  Patient64-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient64-PostM3-16SV3V4  Patient64-PostM6-16SV3V4     Patient64-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient65-PostM1-16SV3V4  Patient65-PostM3-16SV3V4  Patient65-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient65-Pre-16SV3V4  Patient66-PostM1-16SV3V4  Patient66-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient66-PostM6-16SV3V4     Patient66-Pre-16SV3V4  Patient67-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient67-PostM3-16SV3V4  Patient67-PostM6-16SV3V4     Patient67-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient68-PostM1-16SV3V4  Patient68-PostM3-16SV3V4  Patient68-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient68-Pre-16SV3V4  Patient69-PostM1-16SV3V4  Patient69-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient69-PostM6-16SV3V4     Patient69-Pre-16SV3V4  Patient70-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient70-PostM3-16SV3V4  Patient70-PostM6-16SV3V4     Patient70-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient71-PostM1-16SV3V4  Patient71-PostM3-16SV3V4  Patient71-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient71-Pre-16SV3V4  Patient72-PostM1-16SV3V4  Patient72-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient72-PostM6-16SV3V4     Patient72-Pre-16SV3V4  Patient73-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient73-PostM3-16SV3V4  Patient73-PostM6-16SV3V4     Patient73-Pre-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient74-PostM1-16SV3V4  Patient74-PostM3-16SV3V4  Patient74-PostM6-16SV3V4 
    ##                         1                         1                         1 
    ##     Patient74-Pre-16SV3V4  Patient75-PostM1-16SV3V4  Patient75-PostM3-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient75-PostM6-16SV3V4     Patient75-Pre-16SV3V4  Patient76-PostM1-16SV3V4 
    ##                         1                         1                         1 
    ##  Patient76-PostM3-16SV3V4  Patient76-PostM6-16SV3V4     Patient76-Pre-16SV3V4 
    ##                         1                         1                         1

``` r
table(metadata_clean$Site)
```

    ## 
    ## oral rinse 
    ##        384

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ##  NO 
    ## 384

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 384"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 384"

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

    ##    BioProject         Run              sample_name Age BMI Sex Region
    ## 1 PRJNA744870 SRR15719926    Health26-Oral-16SV3V4  NA  NA  NA  V3-V4
    ## 2 PRJNA744870 SRR15719937    Health16-Oral-16SV3V4  NA  NA  NA  V3-V4
    ## 3 PRJNA744870 SRR15719951    Health03-Oral-16SV3V4  NA  NA  NA  V3-V4
    ## 4 PRJNA744870 SRR15719961 Patient11-PostM3-16SV3V4  NA  NA  NA  V3-V4
    ## 5 PRJNA744870 SRR15719966 Patient10-PostM3-16SV3V4  NA  NA  NA  V3-V4
    ## 6 PRJNA744870 SRR15719973    Health70-Oral-16SV3V4  NA  NA  NA  V3-V4
    ##         Platform Country       Site          public_time PCR_primers
    ## 1 Illumina MiSeq   China oral rinse 2021-09-28T00:00:00Z   341F-806R
    ## 2 Illumina MiSeq   China oral rinse 2021-09-28T00:00:00Z   341F-806R
    ## 3 Illumina MiSeq   China oral rinse 2021-09-28T00:00:00Z   341F-806R
    ## 4 Illumina MiSeq   China oral rinse 2021-09-28T00:00:00Z   341F-806R
    ## 5 Illumina MiSeq   China oral rinse 2021-09-28T00:00:00Z   341F-806R
    ## 6 Illumina MiSeq   China oral rinse 2021-09-28T00:00:00Z   341F-806R
    ##                                                                                                                            title
    ## 1 Restoration of the Oral Microbiota After Surgery for Head and Neck Squamous Cell Carcinoma Is Associated With Patient Outcomes
    ## 2 Restoration of the Oral Microbiota After Surgery for Head and Neck Squamous Cell Carcinoma Is Associated With Patient Outcomes
    ## 3 Restoration of the Oral Microbiota After Surgery for Head and Neck Squamous Cell Carcinoma Is Associated With Patient Outcomes
    ## 4 Restoration of the Oral Microbiota After Surgery for Head and Neck Squamous Cell Carcinoma Is Associated With Patient Outcomes
    ## 5 Restoration of the Oral Microbiota After Surgery for Head and Neck Squamous Cell Carcinoma Is Associated With Patient Outcomes
    ## 6 Restoration of the Oral Microbiota After Surgery for Head and Neck Squamous Cell Carcinoma Is Associated With Patient Outcomes
    ##   Disease_study Status matched_samples
    ## 1           YES     HC              NO
    ## 2           YES     HC              NO
    ## 3           YES     HC              NO
    ## 4           YES  HNSCC              NO
    ## 5           YES  HNSCC              NO
    ## 6           YES     HC              NO

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
