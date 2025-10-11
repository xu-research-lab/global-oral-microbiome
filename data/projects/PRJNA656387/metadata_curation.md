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
setwd("D:/oral_project/PRJNA656387")

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
    ## 1 SRR12488523   AMPLICON        298 27838141 PRJNA656387 SAMN15852592
    ## 2 SRR12488524   AMPLICON        297 19372419 PRJNA656387 SAMN15852591
    ## 3 SRR12488525   AMPLICON        297 29310790 PRJNA656387 SAMN15852590
    ## 4 SRR12488526   AMPLICON        296 33249715 PRJNA656387 SAMN15852589
    ## 5 SRR12488527   AMPLICON        297 24064616 PRJNA656387 SAMN15852588
    ## 6 SRR12488528   AMPLICON        296 20628159 PRJNA656387 SAMN15852587
    ##                                BioSampleModel    Bytes              Center.Name
    ## 1 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 16967411 UNIVERSITY OF QUEENSLAND
    ## 2 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 11757366 UNIVERSITY OF QUEENSLAND
    ## 3 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 17847286 UNIVERSITY OF QUEENSLAND
    ## 4 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 20516979 UNIVERSITY OF QUEENSLAND
    ## 5 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 14366077 UNIVERSITY OF QUEENSLAND
    ## 6 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 12378219 UNIVERSITY OF QUEENSLAND
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region  env_broad_scale       env_local_scale
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome Human oral microbiome
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome Human oral microbiome
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome Human oral microbiome
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome Human oral microbiome
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome Human oral microbiome
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome Human oral microbiome
    ##                env_medium Experiment geo_loc_name_country
    ## 1 Human tongue microbiome SRX8980475            Australia
    ## 2 Human tongue microbiome SRX8980474            Australia
    ## 3 Human tongue microbiome SRX8980473            Australia
    ## 4 Human tongue microbiome SRX8980472            Australia
    ## 5 Human tongue microbiome SRX8980471            Australia
    ## 6 Human tongue microbiome SRX8980470            Australia
    ##   geo_loc_name_country_continent        geo_loc_name host_disease
    ## 1                        Oceania Australia: Brisbane      Healthy
    ## 2                        Oceania Australia: Brisbane      Healthy
    ## 3                        Oceania Australia: Brisbane      Healthy
    ## 4                        Oceania Australia: Brisbane      Healthy
    ## 5                        Oceania Australia: Brisbane      Healthy
    ## 6                        Oceania Australia: Brisbane      Healthy
    ##              host_family_relationship         HOST host_tissue_sampled
    ## 1 First-degree relative of RA patient Homo sapiens     Tongue scraping
    ## 2 First-degree relative of RA patient Homo sapiens     Tongue scraping
    ## 3 First-degree relative of RA patient Homo sapiens     Tongue scraping
    ## 4 First-degree relative of RA patient Homo sapiens     Tongue scraping
    ## 5 First-degree relative of RA patient Homo sapiens     Tongue scraping
    ## 6 First-degree relative of RA patient Homo sapiens     Tongue scraping
    ##       Instrument         lat_lon Library.Name LibraryLayout LibrarySelection
    ## 1 Illumina MiSeq 27.28 S 153.1 E        S9074        SINGLE              PCR
    ## 2 Illumina MiSeq 27.28 S 153.1 E        S9073        SINGLE              PCR
    ## 3 Illumina MiSeq 27.28 S 153.1 E        S9072        SINGLE              PCR
    ## 4 Illumina MiSeq 27.28 S 153.1 E        S9071        SINGLE              PCR
    ## 5 Illumina MiSeq 27.28 S 153.1 E        S9070        SINGLE              PCR
    ## 6 Illumina MiSeq 27.28 S 153.1 E        S9069        SINGLE              PCR
    ##   LibrarySource              Organism Platform          ReleaseDate
    ## 1   METAGENOMIC human oral metagenome ILLUMINA 2022-09-29T00:00:00Z
    ## 2   METAGENOMIC human oral metagenome ILLUMINA 2022-09-29T00:00:00Z
    ## 3   METAGENOMIC human oral metagenome ILLUMINA 2022-09-29T00:00:00Z
    ## 4   METAGENOMIC human oral metagenome ILLUMINA 2022-09-29T00:00:00Z
    ## 5   METAGENOMIC human oral metagenome ILLUMINA 2022-09-29T00:00:00Z
    ## 6   METAGENOMIC human oral metagenome ILLUMINA 2022-09-29T00:00:00Z
    ##            create_date version sample_description Sample.Name SRA.Study
    ## 1 2020-08-21T13:41:00Z       1            T_F3205       S9074 SRP278291
    ## 2 2020-08-21T13:42:00Z       1            T_F1104       S9073 SRP278291
    ## 3 2020-08-21T13:42:00Z       1            T_F2408       S9072 SRP278291
    ## 4 2020-08-21T13:42:00Z       1            T_F2404       S9071 SRP278291
    ## 5 2020-08-21T13:41:00Z       1            T_F2212       S9070 SRP278291
    ## 6 2020-08-21T13:42:00Z       1            T_F4102       S9069 SRP278291

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
  Sex = NA,
  Region = NA,
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = replace(initial_meta$Organism, initial_meta$Organism == "human oral metagenome", "tongue"),
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
metadata$title <- "Streptococcus species enriched in the oral cavity of patients with RA are a source of peptidoglycan-polysaccharide polymers that can induce arthritis in mice"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata <- metadata %>%
  mutate(Status = case_when(
    grepl("Healthy", initial_meta$host_disease) ~ "HC",
    grepl("Rheumatoid arthritis", initial_meta$host_disease) ~ "RA"
  ))

#get subjects id according to sample id
metadata$sample_name<-gsub("S9|SA|S8","",metadata$sample_name)
#get subjects id with Saliva available
keep_subjects<-subset(metadata,Site %in% "tongue")
#remove subjects without Saliva samples
metadata<-subset(metadata,sample_name %in% keep_subjects$sample_name)
#get subjects sampling counts
sample_counts <- table(metadata$sample_name)
# Create matched_samples column, any subjects have >=2 samples will be label "YES", otherwise "NO"
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 199"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 199"

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
write.csv(metadata, file = "metadata.csv", quote = TRUE)
head(metadata)
```

    ##    BioProject         Run sample_name Age BMI Sex Region       Platform
    ## 1 PRJNA656387 SRR12488523         074  NA  NA  NA     NA Illumina MiSeq
    ## 2 PRJNA656387 SRR12488524         073  NA  NA  NA     NA Illumina MiSeq
    ## 3 PRJNA656387 SRR12488525         072  NA  NA  NA     NA Illumina MiSeq
    ## 4 PRJNA656387 SRR12488526         071  NA  NA  NA     NA Illumina MiSeq
    ## 5 PRJNA656387 SRR12488527         070  NA  NA  NA     NA Illumina MiSeq
    ## 6 PRJNA656387 SRR12488528         069  NA  NA  NA     NA Illumina MiSeq
    ##     Country   Site          public_time PCR_primers
    ## 1 Australia tongue 2022-09-29T00:00:00Z          NA
    ## 2 Australia tongue 2022-09-29T00:00:00Z          NA
    ## 3 Australia tongue 2022-09-29T00:00:00Z          NA
    ## 4 Australia tongue 2022-09-29T00:00:00Z          NA
    ## 5 Australia tongue 2022-09-29T00:00:00Z          NA
    ## 6 Australia tongue 2022-09-29T00:00:00Z          NA
    ##                                                                                                                                                           title
    ## 1 Streptococcus species enriched in the oral cavity of patients with RA are a source of peptidoglycan-polysaccharide polymers that can induce arthritis in mice
    ## 2 Streptococcus species enriched in the oral cavity of patients with RA are a source of peptidoglycan-polysaccharide polymers that can induce arthritis in mice
    ## 3 Streptococcus species enriched in the oral cavity of patients with RA are a source of peptidoglycan-polysaccharide polymers that can induce arthritis in mice
    ## 4 Streptococcus species enriched in the oral cavity of patients with RA are a source of peptidoglycan-polysaccharide polymers that can induce arthritis in mice
    ## 5 Streptococcus species enriched in the oral cavity of patients with RA are a source of peptidoglycan-polysaccharide polymers that can induce arthritis in mice
    ## 6 Streptococcus species enriched in the oral cavity of patients with RA are a source of peptidoglycan-polysaccharide polymers that can induce arthritis in mice
    ##   Disease_study Status matched_samples
    ## 1           YES     HC              NO
    ## 2           YES     HC              NO
    ## 3           YES     HC              NO
    ## 4           YES     HC              NO
    ## 5           YES     HC              NO
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
