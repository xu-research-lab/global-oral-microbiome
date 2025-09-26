# Step 1: Obtain Initial Metadata from Each Project

To obtain the initial metadata: 1. Search for the BioProject ID (e.g.,
“PRJEBXXXX”) on the [NCBI SRA
website](https://www.ncbi.nlm.nih.gov/Traces/study/). 2. In the “Select”
table, click on “Metadata” and download the file. 3. The downloaded file
is typically named “SraRunTable.csv” and contains metadata for the
project.

    library(tidyverse)
    library(ggplot2)
    library(cowplot)

    setwd("C:/Users/zyc/Desktop/oral_project/data/projects/PRJNA489752")
    initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
    head(initial_meta)

    ##          Run Assay.Type AvgSpotLen     Bases  BioProject    BioSample              BioSampleModel    Bytes    Center.Name Collection_Date
    ## 1 SRR7820324   AMPLICON        599  96788761 PRJNA489752 SAMN09988025 Metagenome or environmental 50475102 MACROGEN, INC.      2017-12-21
    ## 2 SRR7820325   AMPLICON        599  86571955 PRJNA489752 SAMN09988026 Metagenome or environmental 45394452 MACROGEN, INC.      2017-12-21
    ## 3 SRR7820326   AMPLICON        599 100520236 PRJNA489752 SAMN09988023 Metagenome or environmental 52267415 MACROGEN, INC.      2017-12-21
    ## 4 SRR7820327   AMPLICON        599 102977490 PRJNA489752 SAMN09988024 Metagenome or environmental 54010454 MACROGEN, INC.      2017-12-21
    ## 5 SRR7820328   AMPLICON        599 100634029 PRJNA489752 SAMN09988077 Metagenome or environmental 52596721 MACROGEN, INC.      2017-12-21
    ## 6 SRR7820329   AMPLICON        599 104677345 PRJNA489752 SAMN09988078 Metagenome or environmental 54028494 MACROGEN, INC.      2017-12-21
    ##   Consent DATASTORE.filetype DATASTORE.provider                     DATASTORE.region Experiment geo_loc_name_country
    ## 1  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX4671707         uncalculated
    ## 2  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX4671706         uncalculated
    ## 3  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX4671705         uncalculated
    ## 4  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX4671704         uncalculated
    ## 5  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX4671703         uncalculated
    ## 6  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX4671702         uncalculated
    ##   geo_loc_name_country_continent       geo_loc_name         HOST     Instrument isolation_source          lat_lon Library.Name LibraryLayout
    ## 1                   uncalculated Korea: Republic of Homo sapiens Illumina MiSeq           saliva 37.48 N 126.89 E  S-03-Saliva        PAIRED
    ## 2                   uncalculated Korea: Republic of Homo sapiens Illumina MiSeq           saliva 37.48 N 126.89 E  S-04-Saliva        PAIRED
    ## 3                   uncalculated Korea: Republic of Homo sapiens Illumina MiSeq           saliva 37.48 N 126.89 E  S-01-Saliva        PAIRED
    ## 4                   uncalculated Korea: Republic of Homo sapiens Illumina MiSeq           saliva 37.48 N 126.89 E  S-02-Saliva        PAIRED
    ## 5                   uncalculated Korea: Republic of Homo sapiens Illumina MiSeq           saliva 37.48 N 126.89 E  S-07-Saliva        PAIRED
    ## 6                   uncalculated Korea: Republic of Homo sapiens Illumina MiSeq           saliva 37.48 N 126.89 E  S-08-Saliva        PAIRED
    ##   LibrarySelection LibrarySource         Organism Platform          ReleaseDate          create_date version Sample.Name SRA.Study
    ## 1              PCR   METAGENOMIC human metagenome ILLUMINA 2020-09-13T00:00:00Z 2018-09-11T21:14:00Z       1 S-03-Saliva SRP161532
    ## 2              PCR   METAGENOMIC human metagenome ILLUMINA 2020-09-13T00:00:00Z 2018-09-11T21:14:00Z       1 S-04-Saliva SRP161532
    ## 3              PCR   METAGENOMIC human metagenome ILLUMINA 2020-09-13T00:00:00Z 2018-09-11T21:14:00Z       1 S-01-Saliva SRP161532
    ## 4              PCR   METAGENOMIC human metagenome ILLUMINA 2020-09-13T00:00:00Z 2018-09-11T21:14:00Z       1 S-02-Saliva SRP161532
    ## 5              PCR   METAGENOMIC human metagenome ILLUMINA 2020-09-13T00:00:00Z 2018-09-11T21:14:00Z       1 S-07-Saliva SRP161532
    ## 6              PCR   METAGENOMIC human metagenome ILLUMINA 2020-09-13T00:00:00Z 2018-09-11T21:14:00Z       1 S-08-Saliva SRP161532

# Step 2: Create a Uniform Metadata Template

The uniform metadata should include the following columns:

1.  **BioProject**: Project ID (e.g., “PRJNA339212”)
2.  **Run**: Sample ID (e.g., “ERRXXXXXX” or “SRRXXXXXX”)
3.  **sample\_name**: Original sample name submitted by the authors (not
    the database-assigned ID)
4.  **Age**: Age of the subject (if available)
5.  **BMI**: Body Mass Index (if available)
6.  **Sex**: Sex of the subject (“F” or “M”)
7.  **Region**: Targeted 16S rRNA region (e.g., “V3-V4”)
8.  **Platform**: Sequencing platform (e.g., “Illumina MiSeq”)
9.  **Country**: Country of sample origin (e.g., “China”)
10. **Site**: Sample site (e.g., “saliva”, “dental plaque”)
11. **public\_time**: Data release date (format: “YYYY-MM-DD”)
12. **PCR\_primers**: PCR primers used (e.g., “515F-806R”)
13. **title**: Title of the associated publication
14. **Disease\_study**: Whether the study involves a disease (“YES” or
    “NO”)
15. **Status**: Disease status (e.g., “T2DM”, “HC”). Note: If the study
    includes both disease and healthy samples, this must be specified
    per sample.
16. **matched\_samples**: matched samples (“YES” or “NO”)

<!-- -->

    metadata <- data.frame(
      BioProject = initial_meta$BioProject,
      Run = initial_meta$Run,
      sample_name = initial_meta$Sample.Name,
      Age = NA, 
      BMI = NA,  
      Sex = NA,
      Platform = initial_meta$Instrument,
      Country = initial_meta$geo_loc_name,
      Site = initial_meta$isolation_source,
      public_time = initial_meta$ReleaseDate,
      stringsAsFactors = FALSE
    )

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

     # metadata <- metadata %>%
     #  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
    metadata <- metadata %>% mutate(Country = sapply(str_split(Country, ":"), function(x) head(x, 1)))

    metadata$matched_samples = "YES"
    metadata$Region = "V3-V4"
    metadata$Status = "periodontal disease"
    metadata$title <- "Metagenomic characterization of the microbiomes in five different body habitats of otherwise healthy individuals with periodontal disease"
    metadata$Disease_study <- "YES"
    metadata$PCR_primers <- "519F-806R"
    head(metadata)

    ##    BioProject        Run sample_name Age BMI Sex       Platform Country   Site          public_time matched_samples Region
    ## 1 PRJNA489752 SRR7820324 S-03-Saliva  NA  NA  NA Illumina MiSeq   Korea saliva 2020-09-13T00:00:00Z             YES  V3-V4
    ## 2 PRJNA489752 SRR7820325 S-04-Saliva  NA  NA  NA Illumina MiSeq   Korea saliva 2020-09-13T00:00:00Z             YES  V3-V4
    ## 3 PRJNA489752 SRR7820326 S-01-Saliva  NA  NA  NA Illumina MiSeq   Korea saliva 2020-09-13T00:00:00Z             YES  V3-V4
    ## 4 PRJNA489752 SRR7820327 S-02-Saliva  NA  NA  NA Illumina MiSeq   Korea saliva 2020-09-13T00:00:00Z             YES  V3-V4
    ## 5 PRJNA489752 SRR7820328 S-07-Saliva  NA  NA  NA Illumina MiSeq   Korea saliva 2020-09-13T00:00:00Z             YES  V3-V4
    ## 6 PRJNA489752 SRR7820329 S-08-Saliva  NA  NA  NA Illumina MiSeq   Korea saliva 2020-09-13T00:00:00Z             YES  V3-V4
    ##                Status
    ## 1 periodontal disease
    ## 2 periodontal disease
    ## 3 periodontal disease
    ## 4 periodontal disease
    ## 5 periodontal disease
    ## 6 periodontal disease
    ##                                                                                                                                       title
    ## 1 Metagenomic characterization of the microbiomes in five different body habitats of otherwise healthy individuals with periodontal disease
    ## 2 Metagenomic characterization of the microbiomes in five different body habitats of otherwise healthy individuals with periodontal disease
    ## 3 Metagenomic characterization of the microbiomes in five different body habitats of otherwise healthy individuals with periodontal disease
    ## 4 Metagenomic characterization of the microbiomes in five different body habitats of otherwise healthy individuals with periodontal disease
    ## 5 Metagenomic characterization of the microbiomes in five different body habitats of otherwise healthy individuals with periodontal disease
    ## 6 Metagenomic characterization of the microbiomes in five different body habitats of otherwise healthy individuals with periodontal disease
    ##   Disease_study PCR_primers
    ## 1           YES   519F-806R
    ## 2           YES   519F-806R
    ## 3           YES   519F-806R
    ## 4           YES   519F-806R
    ## 5           YES   519F-806R
    ## 6           YES   519F-806R

# Step 4: Summarize and Validate the Metadata

    print(paste("Number of unique runs:", length(unique(metadata$Run))))

    ## [1] "Number of unique runs: 50"

    print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))

    ## [1] "Number of unique sample names: 50"

# Step 5: Visualize the Metadata

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
            geom_histogram(fill = "blue", alpha = 0.6, bins = 20, na.rm = T) +
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
    ggsave("metadata_summary.png", plot = p, width = 8, height = 6, dpi = 150)

    knitr::include_graphics("metadata_summary.png")

![](metadata_summary.png) \# Step 6: Finalize and Export the Metadata

    metadata <- metadata[, c(
      "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
      "Region", "Platform", "Country", "Site", "public_time", 
      "PCR_primers", "title", "Disease_study", "Status","matched_samples"
    )]

    write.csv(metadata, file = "metadata.csv", quote = FALSE)
