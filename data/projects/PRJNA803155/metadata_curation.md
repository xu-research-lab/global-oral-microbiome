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

    setwd("C:/Users/zyc/Desktop/oral_project/data/projects/PRJNA803155")
    initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
    head(initial_meta)

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample              BioSampleModel    Bytes
    ## 1 SRR17871271   AMPLICON        495 81270090 PRJNA803155 SAMN25635642 Metagenome or environmental 26913391
    ## 2 SRR17871272   AMPLICON        492 86192988 PRJNA803155 SAMN25635641 Metagenome or environmental 29056783
    ## 3 SRR17871273   AMPLICON        493 80608951 PRJNA803155 SAMN25635640 Metagenome or environmental 27073955
    ## 4 SRR17871274   AMPLICON        494 77877124 PRJNA803155 SAMN25635639 Metagenome or environmental 25358587
    ## 5 SRR17871275   AMPLICON        495 71995275 PRJNA803155 SAMN25635638 Metagenome or environmental 23678018
    ## 6 SRR17871276   AMPLICON        492 77198244 PRJNA803155 SAMN25635637 Metagenome or environmental 25529223
    ##                                                            Center.Name Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1 SHANGHAI JIAOTONG UNIVERSITY SCHOOL OF MEDCINE 9TH PEOPLE'S HOSPITAL      2020-10-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2 SHANGHAI JIAOTONG UNIVERSITY SCHOOL OF MEDCINE 9TH PEOPLE'S HOSPITAL      2020-10-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3 SHANGHAI JIAOTONG UNIVERSITY SCHOOL OF MEDCINE 9TH PEOPLE'S HOSPITAL      2020-10-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4 SHANGHAI JIAOTONG UNIVERSITY SCHOOL OF MEDCINE 9TH PEOPLE'S HOSPITAL      2020-10-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5 SHANGHAI JIAOTONG UNIVERSITY SCHOOL OF MEDCINE 9TH PEOPLE'S HOSPITAL      2020-10-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6 SHANGHAI JIAOTONG UNIVERSITY SCHOOL OF MEDCINE 9TH PEOPLE'S HOSPITAL      2020-10-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region  Experiment geo_loc_name_country geo_loc_name_country_continent   geo_loc_name         HOST
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX14031031                China                           Asia China:Shanghai Homo sapiens
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX14031030                China                           Asia China:Shanghai Homo sapiens
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX14031029                China                           Asia China:Shanghai Homo sapiens
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX14031028                China                           Asia China:Shanghai Homo sapiens
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX14031027                China                           Asia China:Shanghai Homo sapiens
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX14031026                China                           Asia China:Shanghai Homo sapiens
    ##       Instrument isolation_source      lat_lon Library.Name LibraryLayout LibrarySelection LibrarySource   Organism Platform
    ## 1 Illumina MiSeq          missing 31 N 121.5 E         YS17        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ## 2 Illumina MiSeq          missing 31 N 121.5 E         YS16        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ## 3 Illumina MiSeq          missing 31 N 121.5 E         YS15        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ## 4 Illumina MiSeq          missing 31 N 121.5 E         YS14        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ## 5 Illumina MiSeq          missing 31 N 121.5 E         YS13        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ## 6 Illumina MiSeq          missing 31 N 121.5 E         YS12        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ##            ReleaseDate          create_date version Sample.Name SRA.Study
    ## 1 2022-02-03T00:00:00Z 2022-02-03T11:35:00Z       1        YS17 SRP358074
    ## 2 2022-02-03T00:00:00Z 2022-02-03T11:35:00Z       1        YS16 SRP358074
    ## 3 2022-02-03T00:00:00Z 2022-02-03T11:35:00Z       1        YS15 SRP358074
    ## 4 2022-02-03T00:00:00Z 2022-02-03T11:35:00Z       1        YS14 SRP358074
    ## 5 2022-02-03T00:00:00Z 2022-02-03T11:35:00Z       1        YS13 SRP358074
    ## 6 2022-02-03T00:00:00Z 2022-02-03T11:35:00Z       1        YS12 SRP358074

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
      Country = initial_meta$geo_loc_name_country,
      # Site = initial_meta$env_medium,
      public_time = initial_meta$ReleaseDate,
      # Status = initial_meta$Sample.Name,
      stringsAsFactors = FALSE
    )

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

     # metadata <- metadata %>%
     #  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
    # metadata <-metadata%>%
    #   mutate(Status = case_when(Site == "mucous" ~ "OSCC",
    #                             Site == "dental plaque" ~ "OSCC",
    #                             Site == "tumor" ~ "OSCC",
    #                             TRUE ~ NA))

    metadata$Site <- "tumor"
    metadata$Status <- "OSCC"
    metadata$matched_samples <- "NO"
    metadata$Region <- "V4-V5"
    metadata$title <- "Age-Related Cancer-Associated Microbiota Potentially Promotes Oral Squamous Cell Cancer Tumorigenesis by Distinct Mechanisms"
    metadata$Disease_study <- "YES"
    metadata$PCR_primers <- "515F-907R"
    head(metadata)

    ##    BioProject         Run sample_name Age BMI Sex       Platform Country          public_time  Site Status matched_samples Region
    ## 1 PRJNA803155 SRR17871271        YS17  NA  NA  NA Illumina MiSeq   China 2022-02-03T00:00:00Z tumor   OSCC              NO  V4-V5
    ## 2 PRJNA803155 SRR17871272        YS16  NA  NA  NA Illumina MiSeq   China 2022-02-03T00:00:00Z tumor   OSCC              NO  V4-V5
    ## 3 PRJNA803155 SRR17871273        YS15  NA  NA  NA Illumina MiSeq   China 2022-02-03T00:00:00Z tumor   OSCC              NO  V4-V5
    ## 4 PRJNA803155 SRR17871274        YS14  NA  NA  NA Illumina MiSeq   China 2022-02-03T00:00:00Z tumor   OSCC              NO  V4-V5
    ## 5 PRJNA803155 SRR17871275        YS13  NA  NA  NA Illumina MiSeq   China 2022-02-03T00:00:00Z tumor   OSCC              NO  V4-V5
    ## 6 PRJNA803155 SRR17871276        YS12  NA  NA  NA Illumina MiSeq   China 2022-02-03T00:00:00Z tumor   OSCC              NO  V4-V5
    ##                                                                                                                          title Disease_study
    ## 1 Age-Related Cancer-Associated Microbiota Potentially Promotes Oral Squamous Cell Cancer Tumorigenesis by Distinct Mechanisms           YES
    ## 2 Age-Related Cancer-Associated Microbiota Potentially Promotes Oral Squamous Cell Cancer Tumorigenesis by Distinct Mechanisms           YES
    ## 3 Age-Related Cancer-Associated Microbiota Potentially Promotes Oral Squamous Cell Cancer Tumorigenesis by Distinct Mechanisms           YES
    ## 4 Age-Related Cancer-Associated Microbiota Potentially Promotes Oral Squamous Cell Cancer Tumorigenesis by Distinct Mechanisms           YES
    ## 5 Age-Related Cancer-Associated Microbiota Potentially Promotes Oral Squamous Cell Cancer Tumorigenesis by Distinct Mechanisms           YES
    ## 6 Age-Related Cancer-Associated Microbiota Potentially Promotes Oral Squamous Cell Cancer Tumorigenesis by Distinct Mechanisms           YES
    ##   PCR_primers
    ## 1   515F-907R
    ## 2   515F-907R
    ## 3   515F-907R
    ## 4   515F-907R
    ## 5   515F-907R
    ## 6   515F-907R

# Step 4: Summarize and Validate the Metadata

    print(paste("Number of unique runs:", length(unique(metadata$Run))))

    ## [1] "Number of unique runs: 40"

    print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))

    ## [1] "Number of unique sample names: 40"

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
    ggsave("Metadata_Summary.png", plot = p, width = 8, height = 6, dpi = 150)

    knitr::include_graphics("Metadata_Summary.png")

![](Metadata_Summary.png) \# Step 6: Finalize and Export the Metadata

    metadata <- metadata[, c(
      "BioProject", "Run", "sample_name", "Age", "BMI", "Sex",
      "Region", "Platform", "Country", "Site", "public_time",
      "PCR_primers", "title", "Disease_study", "Status","matched_samples"
    )]



    write.csv(metadata, file = "metadata.csv", quote = FALSE)
