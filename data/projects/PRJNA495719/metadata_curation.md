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

    setwd("C:/Users/zyc/Desktop/oral_project/data/projects/PRJNA495719")
    initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
    head(initial_meta)

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample              BioSampleModel   Bytes
    ## 1 SRR8003020   AMPLICON        602 10425436 PRJNA495719 SAMN10230350 Metagenome or environmental 4798310
    ## 2 SRR8003021   AMPLICON        602 12638990 PRJNA495719 SAMN10230349 Metagenome or environmental 5784894
    ## 3 SRR8003022   AMPLICON        602 13898976 PRJNA495719 SAMN10230344 Metagenome or environmental 6421708
    ## 4 SRR8003023   AMPLICON        602 12446952 PRJNA495719 SAMN10230343 Metagenome or environmental 5708753
    ## 5 SRR8003024   AMPLICON        602 13874294 PRJNA495719 SAMN10230342 Metagenome or environmental 6374965
    ## 6 SRR8003025   AMPLICON        602 15226988 PRJNA495719 SAMN10230341 Metagenome or environmental 6993320
    ##                                                              Center.Name Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1 COLLEGE OF STOMATOLOGY, CHONGQING MEDICAL UNIVERSITY, CHONGQING, CHINA      2017-04-25  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2 COLLEGE OF STOMATOLOGY, CHONGQING MEDICAL UNIVERSITY, CHONGQING, CHINA      2017-04-25  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3 COLLEGE OF STOMATOLOGY, CHONGQING MEDICAL UNIVERSITY, CHONGQING, CHINA      2017-04-09  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4 COLLEGE OF STOMATOLOGY, CHONGQING MEDICAL UNIVERSITY, CHONGQING, CHINA      2017-04-09  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5 COLLEGE OF STOMATOLOGY, CHONGQING MEDICAL UNIVERSITY, CHONGQING, CHINA      2017-04-09  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6 COLLEGE OF STOMATOLOGY, CHONGQING MEDICAL UNIVERSITY, CHONGQING, CHINA      2017-04-09  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region Experiment geo_loc_name_country geo_loc_name_country_continent geo_loc_name         HOST
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX4833941         uncalculated                   uncalculated      missing Homo sapiens
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX4833940         uncalculated                   uncalculated      missing Homo sapiens
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX4833939         uncalculated                   uncalculated      missing Homo sapiens
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX4833938         uncalculated                   uncalculated      missing Homo sapiens
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX4833937         uncalculated                   uncalculated      missing Homo sapiens
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX4833936         uncalculated                   uncalculated      missing Homo sapiens
    ##       Instrument isolation_source lat_lon Library.Name LibraryLayout LibrarySelection LibrarySource         Organism Platform
    ## 1 Illumina MiSeq          missing missing         TP18        PAIRED              PCR   METAGENOMIC human metagenome ILLUMINA
    ## 2 Illumina MiSeq          missing missing         TP17        PAIRED              PCR   METAGENOMIC human metagenome ILLUMINA
    ## 3 Illumina MiSeq          missing missing         TP12        PAIRED              PCR   METAGENOMIC human metagenome ILLUMINA
    ## 4 Illumina MiSeq          missing missing         TP11        PAIRED              PCR   METAGENOMIC human metagenome ILLUMINA
    ## 5 Illumina MiSeq          missing missing         TP10        PAIRED              PCR   METAGENOMIC human metagenome ILLUMINA
    ## 6 Illumina MiSeq          missing missing          TP9        PAIRED              PCR   METAGENOMIC human metagenome ILLUMINA
    ##            ReleaseDate Sample.Name SRA.Study tmp          create_date version
    ## 1 2018-10-11T00:00:00Z        TP18 SRP165148  40 2018-10-11T04:56:00Z       1
    ## 2 2018-10-11T00:00:00Z        TP17 SRP165148  39 2018-10-11T04:57:00Z       1
    ## 3 2018-10-11T00:00:00Z        TP12 SRP165148  34 2018-10-11T04:57:00Z       1
    ## 4 2018-10-11T00:00:00Z        TP11 SRP165148  33 2018-10-11T04:57:00Z       1
    ## 5 2018-10-11T00:00:00Z        TP10 SRP165148  32 2018-10-11T04:57:00Z       1
    ## 6 2018-10-11T00:00:00Z         TP9 SRP165148  31 2018-10-11T04:57:00Z       1

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
      # Country = initial_meta$geo_loc_name,
      # Status = initial_meta$isolate,
      # Site = initial_meta$tissue,
      public_time = initial_meta$ReleaseDate,
      stringsAsFactors = FALSE
    )

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

     # metadata <- metadata %>%
     #  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
    metadata <- metadata %>% mutate(Status = case_when(
      grepl("CP|CS",sample_name) ~ "caries-free",
      grepl("TP|TS",sample_name) ~ "caries-active"
    ))
    metadata <- metadata %>% mutate(Site = case_when(
      grepl("CP|TP",sample_name) ~ "dental plaque",
      grepl("CS|TS",sample_name) ~ " saliva"
    ))

    metadata$Country <- "China"
    metadata$matched_samples = "NO"
    metadata$Region = "V3-V4"

    metadata$title <- "The Oral Microbiome in the Elderly With Dental Caries and Health"
    metadata$Disease_study <- "YES"
    metadata$PCR_primers <- "338F-806R"
    head(metadata)

    ##    BioProject        Run sample_name Age BMI Sex       Platform          public_time        Status          Site Country matched_samples
    ## 1 PRJNA495719 SRR8003020        TP18  NA  NA  NA Illumina MiSeq 2018-10-11T00:00:00Z caries-active dental plaque   China              NO
    ## 2 PRJNA495719 SRR8003021        TP17  NA  NA  NA Illumina MiSeq 2018-10-11T00:00:00Z caries-active dental plaque   China              NO
    ## 3 PRJNA495719 SRR8003022        TP12  NA  NA  NA Illumina MiSeq 2018-10-11T00:00:00Z caries-active dental plaque   China              NO
    ## 4 PRJNA495719 SRR8003023        TP11  NA  NA  NA Illumina MiSeq 2018-10-11T00:00:00Z caries-active dental plaque   China              NO
    ## 5 PRJNA495719 SRR8003024        TP10  NA  NA  NA Illumina MiSeq 2018-10-11T00:00:00Z caries-active dental plaque   China              NO
    ## 6 PRJNA495719 SRR8003025         TP9  NA  NA  NA Illumina MiSeq 2018-10-11T00:00:00Z caries-active dental plaque   China              NO
    ##   Region                                                            title Disease_study PCR_primers
    ## 1  V3-V4 The Oral Microbiome in the Elderly With Dental Caries and Health           YES   338F-806R
    ## 2  V3-V4 The Oral Microbiome in the Elderly With Dental Caries and Health           YES   338F-806R
    ## 3  V3-V4 The Oral Microbiome in the Elderly With Dental Caries and Health           YES   338F-806R
    ## 4  V3-V4 The Oral Microbiome in the Elderly With Dental Caries and Health           YES   338F-806R
    ## 5  V3-V4 The Oral Microbiome in the Elderly With Dental Caries and Health           YES   338F-806R
    ## 6  V3-V4 The Oral Microbiome in the Elderly With Dental Caries and Health           YES   338F-806R

# Step 4: Summarize and Validate the Metadata

    print(paste("Number of unique runs:", length(unique(metadata$Run))))

    ## [1] "Number of unique runs: 90"

    print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))

    ## [1] "Number of unique sample names: 90"

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
