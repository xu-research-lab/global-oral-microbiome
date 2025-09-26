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

    setwd("C:/Users/zyc/Desktop/oral_project/PRJNA421234")
    initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
    head(initial_meta)

    ##          Run Assay.Type AvgSpotLen     Bases  BioProject    BioSample              BioSampleModel     Bytes
    ## 1 SRR6364953   AMPLICON        584 192482296 PRJNA421234 SAMN08135458 Metagenome or environmental 125225726
    ## 2 SRR6364954   AMPLICON        599 154730497 PRJNA421234 SAMN08135457 Metagenome or environmental 101770082
    ## 3 SRR6364955   AMPLICON        596 107410059 PRJNA421234 SAMN08135460 Metagenome or environmental  70506452
    ## 4 SRR6364956   AMPLICON        579 189357173 PRJNA421234 SAMN08135459 Metagenome or environmental 120944546
    ## 5 SRR6364957   AMPLICON        599  88498101 PRJNA421234 SAMN08135454 Metagenome or environmental  57851331
    ## 6 SRR6364958   AMPLICON        582  84488425 PRJNA421234 SAMN08135453 Metagenome or environmental  54812380
    ##                  Center.Name Collection_Date Consent DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1 THE UNIVERSITY OF AUCKLAND         2016-02  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2 THE UNIVERSITY OF AUCKLAND         2016-02  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3 THE UNIVERSITY OF AUCKLAND         2017-03  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4 THE UNIVERSITY OF AUCKLAND         2017-03  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5 THE UNIVERSITY OF AUCKLAND         2016-02  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6 THE UNIVERSITY OF AUCKLAND         2016-03  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##                disease Experiment geo_loc_name_country geo_loc_name_country_continent          geo_loc_name         HOST
    ## 1 Dentally compromised SRX3460691          New Zealand                        Oceania New Zealand: Auckland Homo sapiens
    ## 2 Dentally compromised SRX3460690          New Zealand                        Oceania New Zealand: Auckland Homo sapiens
    ## 3 Dentally compromised SRX3460689          New Zealand                        Oceania New Zealand: Auckland Homo sapiens
    ## 4 Dentally compromised SRX3460688          New Zealand                        Oceania New Zealand: Auckland Homo sapiens
    ## 5 Dentally compromised SRX3460687          New Zealand                        Oceania New Zealand: Auckland Homo sapiens
    ## 6 Dentally compromised SRX3460686          New Zealand                        Oceania New Zealand: Auckland Homo sapiens
    ##       Instrument isolation_source       lat_lon             Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1 Illumina MiSeq           saliva not collected  CONTROL5_SALIVA_16S_S59        PAIRED              PCR   METAGENOMIC
    ## 2 Illumina MiSeq           saliva not collected  CONTROL4_SALIVA_16S_S58        PAIRED              PCR   METAGENOMIC
    ## 3 Illumina MiSeq           saliva not collected  CONTROL7_SALIVA_16S_S61        PAIRED              PCR   METAGENOMIC
    ## 4 Illumina MiSeq           saliva not collected  CONTROL6_SALIVA_16S_S60        PAIRED              PCR   METAGENOMIC
    ## 5 Illumina MiSeq           saliva not collected  CONTROL1_SALIVA_16S_S55        PAIRED              PCR   METAGENOMIC
    ## 6 Illumina MiSeq           saliva not collected CONTROL10_SALIVA_16S_S64        PAIRED              PCR   METAGENOMIC
    ##                Organism Platform          ReleaseDate          create_date version samp_collect_device samp_mat_process
    ## 1 human oral metagenome ILLUMINA 2019-01-01T00:00:00Z 2017-12-11T02:20:00Z       1 unstimulated saliva   Illumina MiSeq
    ## 2 human oral metagenome ILLUMINA 2019-01-01T00:00:00Z 2017-12-11T02:19:00Z       1 unstimulated saliva   Illumina MiSeq
    ## 3 human oral metagenome ILLUMINA 2019-01-01T00:00:00Z 2017-12-11T02:19:00Z       1 unstimulated saliva   Illumina MiSeq
    ## 4 human oral metagenome ILLUMINA 2019-01-01T00:00:00Z 2017-12-11T02:20:00Z       1 unstimulated saliva   Illumina MiSeq
    ## 5 human oral metagenome ILLUMINA 2019-01-01T00:00:00Z 2017-12-11T02:19:00Z       1 unstimulated saliva   Illumina MiSeq
    ## 6 human oral metagenome ILLUMINA 2019-01-01T00:00:00Z 2017-12-11T02:20:00Z       1 unstimulated saliva   Illumina MiSeq
    ##                Sample.Name SAMPLE SRA.Study
    ## 1  CONTROL5_SALIVA_16S_S59    S59 SRP126472
    ## 2  CONTROL4_SALIVA_16S_S58    S58 SRP126472
    ## 3  CONTROL7_SALIVA_16S_S61    S61 SRP126472
    ## 4  CONTROL6_SALIVA_16S_S60    S60 SRP126472
    ## 5  CONTROL1_SALIVA_16S_S55    S55 SRP126472
    ## 6 CONTROL10_SALIVA_16S_S64    S64 SRP126472

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
      Site = initial_meta$isolation_source,
      public_time = initial_meta$ReleaseDate,
      stringsAsFactors = FALSE
    )

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

     # metadata <- metadata %>%
     #  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
    metadata <-metadata%>%
      mutate(Status = sapply(str_split(sample_name, "_"), function(x) head(x, 1)))%>%
      mutate(Status = str_remove_all(Status,"[^A-Za-z]"))%>%
      mutate(Status = case_when(
        Status == "SALIVA" ~ "HC",
        Status == "OM" ~ "HNSCC",
        TRUE ~ "periodontal disease"
      ))

    metadata$matched_samples = "NO"
    metadata$Region = "V3-V4"
    metadata$title <- "The Clinical Potential of Oral Microbiota as a Screening Tool for Oral Squamous Cell Carcinomas"
    metadata$Disease_study <- "YES"
    metadata$PCR_primers <- "341F-806R"
    head(metadata)

    ##    BioProject        Run              sample_name Age BMI Sex       Platform     Country   Site          public_time
    ## 1 PRJNA421234 SRR6364953  CONTROL5_SALIVA_16S_S59  NA  NA  NA Illumina MiSeq New Zealand saliva 2019-01-01T00:00:00Z
    ## 2 PRJNA421234 SRR6364954  CONTROL4_SALIVA_16S_S58  NA  NA  NA Illumina MiSeq New Zealand saliva 2019-01-01T00:00:00Z
    ## 3 PRJNA421234 SRR6364955  CONTROL7_SALIVA_16S_S61  NA  NA  NA Illumina MiSeq New Zealand saliva 2019-01-01T00:00:00Z
    ## 4 PRJNA421234 SRR6364956  CONTROL6_SALIVA_16S_S60  NA  NA  NA Illumina MiSeq New Zealand saliva 2019-01-01T00:00:00Z
    ## 5 PRJNA421234 SRR6364957  CONTROL1_SALIVA_16S_S55  NA  NA  NA Illumina MiSeq New Zealand saliva 2019-01-01T00:00:00Z
    ## 6 PRJNA421234 SRR6364958 CONTROL10_SALIVA_16S_S64  NA  NA  NA Illumina MiSeq New Zealand saliva 2019-01-01T00:00:00Z
    ##                Status matched_samples Region
    ## 1 periodontal disease              NO  V3-V4
    ## 2 periodontal disease              NO  V3-V4
    ## 3 periodontal disease              NO  V3-V4
    ## 4 periodontal disease              NO  V3-V4
    ## 5 periodontal disease              NO  V3-V4
    ## 6 periodontal disease              NO  V3-V4
    ##                                                                                             title Disease_study PCR_primers
    ## 1 The Clinical Potential of Oral Microbiota as a Screening Tool for Oral Squamous Cell Carcinomas           YES   341F-806R
    ## 2 The Clinical Potential of Oral Microbiota as a Screening Tool for Oral Squamous Cell Carcinomas           YES   341F-806R
    ## 3 The Clinical Potential of Oral Microbiota as a Screening Tool for Oral Squamous Cell Carcinomas           YES   341F-806R
    ## 4 The Clinical Potential of Oral Microbiota as a Screening Tool for Oral Squamous Cell Carcinomas           YES   341F-806R
    ## 5 The Clinical Potential of Oral Microbiota as a Screening Tool for Oral Squamous Cell Carcinomas           YES   341F-806R
    ## 6 The Clinical Potential of Oral Microbiota as a Screening Tool for Oral Squamous Cell Carcinomas           YES   341F-806R

# Step 4: Summarize and Validate the Metadata

    print(paste("Number of unique runs:", length(unique(metadata$Run))))

    ## [1] "Number of unique runs: 30"

    print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))

    ## [1] "Number of unique sample names: 30"

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
