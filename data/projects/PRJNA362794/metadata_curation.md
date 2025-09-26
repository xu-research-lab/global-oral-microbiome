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
library(tidyverse)
library(ggplot2)
library(cowplot)

setwd("C:/Users/zyc/Desktop/oral_project/data/projects/PRJNA362794")
initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample                       BioSampleModel    Bytes
    ## 1 SRR5196235   AMPLICON        494 27033131 PRJNA362794 SAMN06240299 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 17324281
    ## 2 SRR5196236   AMPLICON        494 25527276 PRJNA362794 SAMN06241860 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 16271666
    ## 3 SRR5196237   AMPLICON        494 24645185 PRJNA362794 SAMN06240331 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 15571251
    ## 4 SRR5196238   AMPLICON        494 29796878 PRJNA362794 SAMN06240332 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 18771287
    ## 5 SRR5196239   AMPLICON        494 24372247 PRJNA362794 SAMN06240339 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 15726825
    ## 6 SRR5196240   AMPLICON        494 26370965 PRJNA362794 SAMN06240348 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 16891834
    ##                                                                  Center.Name Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1 NINTH PEOPLE'S HOSPITAL, SCHOOL OF MEDICINE, SHANGHAI JIAO TONG UNIVERSITY      2016-03-18  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2 NINTH PEOPLE'S HOSPITAL, SCHOOL OF MEDICINE, SHANGHAI JIAO TONG UNIVERSITY      2016-03-22  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3 NINTH PEOPLE'S HOSPITAL, SCHOOL OF MEDICINE, SHANGHAI JIAO TONG UNIVERSITY      2016-03-30  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4 NINTH PEOPLE'S HOSPITAL, SCHOOL OF MEDICINE, SHANGHAI JIAO TONG UNIVERSITY      2016-04-01  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5 NINTH PEOPLE'S HOSPITAL, SCHOOL OF MEDICINE, SHANGHAI JIAO TONG UNIVERSITY      2016-04-02  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6 NINTH PEOPLE'S HOSPITAL, SCHOOL OF MEDICINE, SHANGHAI JIAO TONG UNIVERSITY      2016-04-03  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region env_biome env_feature env_material Experiment geo_loc_name_country geo_loc_name_country_continent
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1      oral      mucosa suface swabs SRX2511896                China                           Asia
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1      oral      mucosa suface swabs SRX2511897                China                           Asia
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1      oral      mucosa suface swabs SRX2511898                China                           Asia
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1      oral      mucosa suface swabs SRX2511899                China                           Asia
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1      oral      mucosa suface swabs SRX2511900                China                           Asia
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1      oral      mucosa suface swabs SRX2511901                China                           Asia
    ##     geo_loc_name         HOST     Instrument       lat_lon Library.Name LibraryLayout LibrarySelection LibrarySource              Organism
    ## 1 China:Shanghai Homo sapiens Illumina MiSeq not collected          1_N        PAIRED              PCR   METAGENOMIC human oral metagenome
    ## 2 China:Shanghai Homo sapiens Illumina MiSeq not collected          3_N        PAIRED              PCR   METAGENOMIC human oral metagenome
    ## 3 China:Shanghai Homo sapiens Illumina MiSeq not collected          4_N        PAIRED              PCR   METAGENOMIC human oral metagenome
    ## 4 China:Shanghai Homo sapiens Illumina MiSeq not collected          5_N        PAIRED              PCR   METAGENOMIC human oral metagenome
    ## 5 China:Shanghai Homo sapiens Illumina MiSeq not collected          6_N        PAIRED              PCR   METAGENOMIC human oral metagenome
    ## 6 China:Shanghai Homo sapiens Illumina MiSeq not collected          7_N        PAIRED              PCR   METAGENOMIC human oral metagenome
    ##   Platform          ReleaseDate          create_date version Sample.Name SRA.Study
    ## 1 ILLUMINA 2017-01-28T00:00:00Z 2017-01-23T21:18:00Z       1         1_N SRP097643
    ## 2 ILLUMINA 2017-01-28T00:00:00Z 2017-01-23T21:18:00Z       1         3_N SRP097643
    ## 3 ILLUMINA 2017-01-28T00:00:00Z 2017-01-23T21:18:00Z       1         4_N SRP097643
    ## 4 ILLUMINA 2017-01-28T00:00:00Z 2017-01-23T21:18:00Z       1         5_N SRP097643
    ## 5 ILLUMINA 2017-01-28T00:00:00Z 2017-01-23T21:18:00Z       1         6_N SRP097643
    ## 6 ILLUMINA 2017-01-28T00:00:00Z 2017-01-23T21:18:00Z       1         7_N SRP097643

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
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Sample.Name,
  Age = NA, 
  BMI = NA,  
  Sex = NA,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$env_feature,
  public_time = initial_meta$ReleaseDate,
  matched_samples = "NO",
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# metadata <- metadata %>%
#  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
```

``` r
metadata$title <- "Signatures of Mucosal Microbiome in Oral Squamous Cell Carcinoma Identified Using a Random Forest Model"
metadata$Disease_study <- "YES"
metadata$PCR_primers <- "341F-806R"
metadata <- metadata%>%mutate(Status = sapply(str_split(sample_name, "_"), function(x) tail(x, 1)))%>% mutate(Status = case_when(
    Status == "N" ~ "OSCC_normal",
    TRUE ~ "OSCC_tumor" # 其他所有情况
  ))
head(metadata)
```

    ##    BioProject        Run sample_name Age BMI Sex Region       Platform Country   Site          public_time matched_samples
    ## 1 PRJNA362794 SRR5196235         1_N  NA  NA  NA  V3-V4 Illumina MiSeq   China mucosa 2017-01-28T00:00:00Z              NO
    ## 2 PRJNA362794 SRR5196236         3_N  NA  NA  NA  V3-V4 Illumina MiSeq   China mucosa 2017-01-28T00:00:00Z              NO
    ## 3 PRJNA362794 SRR5196237         4_N  NA  NA  NA  V3-V4 Illumina MiSeq   China mucosa 2017-01-28T00:00:00Z              NO
    ## 4 PRJNA362794 SRR5196238         5_N  NA  NA  NA  V3-V4 Illumina MiSeq   China mucosa 2017-01-28T00:00:00Z              NO
    ## 5 PRJNA362794 SRR5196239         6_N  NA  NA  NA  V3-V4 Illumina MiSeq   China mucosa 2017-01-28T00:00:00Z              NO
    ## 6 PRJNA362794 SRR5196240         7_N  NA  NA  NA  V3-V4 Illumina MiSeq   China mucosa 2017-01-28T00:00:00Z              NO
    ##                                                                                                     title Disease_study PCR_primers
    ## 1 Signatures of Mucosal Microbiome in Oral Squamous Cell Carcinoma Identified Using a Random Forest Model           YES   341F-806R
    ## 2 Signatures of Mucosal Microbiome in Oral Squamous Cell Carcinoma Identified Using a Random Forest Model           YES   341F-806R
    ## 3 Signatures of Mucosal Microbiome in Oral Squamous Cell Carcinoma Identified Using a Random Forest Model           YES   341F-806R
    ## 4 Signatures of Mucosal Microbiome in Oral Squamous Cell Carcinoma Identified Using a Random Forest Model           YES   341F-806R
    ## 5 Signatures of Mucosal Microbiome in Oral Squamous Cell Carcinoma Identified Using a Random Forest Model           YES   341F-806R
    ## 6 Signatures of Mucosal Microbiome in Oral Squamous Cell Carcinoma Identified Using a Random Forest Model           YES   341F-806R
    ##        Status
    ## 1 OSCC_normal
    ## 2 OSCC_normal
    ## 3 OSCC_normal
    ## 4 OSCC_normal
    ## 5 OSCC_normal
    ## 6 OSCC_normal

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
plot_data <- metadata %>%
  select(-BioProject, -Run, -sample_name, -title)

# Function to create simple plots for each column
create_simple_plots <- function(data) {
  plots <- list()
  
  for (col in colnames(data)) {
    # Check if the column is all NA
    if (all(is.na(data[[col]]))) {
      # Create an empty plot with text indicating all NA
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "this column is NA", size = 6) +
        labs(title = col) +
        theme_void()
    } else if (is.numeric(data[[col]])) {
      # Numeric variable: histogram
      p <- ggplot(data, aes(x = .data[[col]])) +
        geom_histogram(fill = "blue", alpha = 0.6, bins = 20, na.rm = TRUE) +
        labs(title = col) +
        theme_minimal()
    } else {
      # Categorical variable: bar plot (remove NA values)
      count_data <- as.data.frame(table(data[[col]], useNA = "ifany"))
      p <- ggplot(count_data, aes(x = reorder(Var1, -Freq), y = Freq)) +
        geom_col(fill = "red", alpha = 0.6) +
        labs(title = col, x = NULL, y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 15, hjust = 1))
    }
    plots[[col]] <- p
  }
  
  return(plots)
}


plots <- create_simple_plots(plot_data)
p<-plot_grid(plots$Age, 
          plots$BMI, 
          plots$Sex, 
          plots$Region, 
          plots$Platform, 
          plots$Country, 
          plots$Site,
          plots$public_time, 
          plots$PCR_primers,
          plots$Disease_study,
          plots$Status,
          plots$matched_samples,# add 
          ncol = 3, align = "hv", axis = "lb")
ggsave("Metadata_Summary.png", plot = p, width = 8, height = 6, dpi = 150)
```

``` r
knitr::include_graphics("Metadata_Summary.png")
```

![](Metadata_Summary.png)<!-- -->

# Step 6: Finalize and Export the Metadata

``` r
metadata <- metadata[, c(
  "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
  "Region", "Platform", "Country", "Site", "public_time", 
  "PCR_primers", "title", "Disease_study", "Status","matched_samples"
)]

write.csv(metadata, file = "metadata.csv", quote = FALSE)
```
