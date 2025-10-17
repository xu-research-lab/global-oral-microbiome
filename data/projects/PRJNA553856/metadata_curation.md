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

    setwd("C:/Users/zyc/Desktop/oral_project/data/projects/PRJNA553856")
    initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
    head(initial_meta)

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample              BioSampleModel    Bytes         Center.Name Collection_Date
    ## 1 SRR9723153    RNA-Seq        239 24514949 PRJNA553856 SAMN12252306 Metagenome or environmental 11592741 NANCHANG UNIVERSITY      2018-10-15
    ## 2 SRR9723154    RNA-Seq        239 25160457 PRJNA553856 SAMN12252306 Metagenome or environmental 12021804 NANCHANG UNIVERSITY      2018-10-15
    ## 3 SRR9723155    RNA-Seq        239 12452111 PRJNA553856 SAMN12252306 Metagenome or environmental  5992875 NANCHANG UNIVERSITY      2018-10-15
    ## 4 SRR9723156    RNA-Seq        240 12895693 PRJNA553856 SAMN12252306 Metagenome or environmental  6313538 NANCHANG UNIVERSITY      2018-10-15
    ## 5 SRR9723157    RNA-Seq        239  9453908 PRJNA553856 SAMN12252306 Metagenome or environmental  4550777 NANCHANG UNIVERSITY      2018-10-15
    ## 6 SRR9723158    RNA-Seq        240 11019127 PRJNA553856 SAMN12252306 Metagenome or environmental  5405864 NANCHANG UNIVERSITY      2018-10-15
    ##   Consent DATASTORE.filetype DATASTORE.provider                     DATASTORE.region Experiment geo_loc_name_country
    ## 1  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6480725                China
    ## 2  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6480724                China
    ## 3  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6480723                China
    ## 4  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6480722                China
    ## 5  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6480721                China
    ## 6  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6480720                China
    ##   geo_loc_name_country_continent geo_loc_name         HOST          Instrument isolation_source lat_lon Library.Name LibraryLayout
    ## 1                           Asia        China Homo sapiens Illumina HiSeq 2000   oral secretion missing         CBM4        SINGLE
    ## 2                           Asia        China Homo sapiens Illumina HiSeq 2000   oral secretion missing         CBM5        SINGLE
    ## 3                           Asia        China Homo sapiens Illumina HiSeq 2000   oral secretion missing         ABM4        SINGLE
    ## 4                           Asia        China Homo sapiens Illumina HiSeq 2000   oral secretion missing         ABM5        SINGLE
    ## 5                           Asia        China Homo sapiens Illumina HiSeq 2000   oral secretion missing         CBW3        SINGLE
    ## 6                           Asia        China Homo sapiens Illumina HiSeq 2000   oral secretion missing         CBW4        SINGLE
    ##   LibrarySelection LibrarySource              Organism Platform          ReleaseDate          create_date version
    ## 1              PCR       GENOMIC human oral metagenome ILLUMINA 2019-07-22T00:00:00Z 2019-07-22T09:10:00Z       1
    ## 2              PCR       GENOMIC human oral metagenome ILLUMINA 2019-07-22T00:00:00Z 2019-07-22T09:10:00Z       1
    ## 3              PCR       GENOMIC human oral metagenome ILLUMINA 2019-07-22T00:00:00Z 2019-07-22T09:10:00Z       1
    ## 4              PCR       GENOMIC human oral metagenome ILLUMINA 2019-07-22T00:00:00Z 2019-07-22T09:10:00Z       1
    ## 5              PCR       GENOMIC human oral metagenome ILLUMINA 2019-07-22T00:00:00Z 2019-07-22T09:11:00Z       1
    ## 6              PCR       GENOMIC human oral metagenome ILLUMINA 2019-07-22T00:00:00Z 2019-07-22T09:11:00Z       1
    ##                           Sample.Name SRA.Study
    ## 1 toothbrush microbe and oral microbe SRP215850
    ## 2 toothbrush microbe and oral microbe SRP215850
    ## 3 toothbrush microbe and oral microbe SRP215850
    ## 4 toothbrush microbe and oral microbe SRP215850
    ## 5 toothbrush microbe and oral microbe SRP215850
    ## 6 toothbrush microbe and oral microbe SRP215850

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
      sample_name = initial_meta$Library.Name,
      Age = NA,
      BMI = NA,
      Sex = NA,
      Platform = initial_meta$Instrument,
      Country = initial_meta$geo_loc_name_country,
      # Site = initial_meta$isolation_source,
      public_time = initial_meta$ReleaseDate,
      # Status = initial_meta$Sample.Name,
      stringsAsFactors = FALSE
    )

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

    metadata <- metadata %>% filter(str_detect(sample_name, "AO") | str_detect(sample_name, "CO"))


     # metadata <- metadata %>%
     #  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
    # metadata <-metadata%>%
    #   mutate(Status = case_when(Site == "mucous" ~ "OSCC",
    #                             Site == "dental plaque" ~ "OSCC",
    #                             Site == "tumor" ~ "OSCC",
    #                             TRUE ~ NA))

    metadata$Site <- "oral cavity"
    metadata$Status <- "HC"
    metadata$matched_samples <- "NO"
    metadata$Region <- "V4"
    metadata$title <- "Interaction of Oral and Toothbrush Microbiota Affects Oral Cavity Health"
    metadata$Disease_study <- "NO"
    metadata$PCR_primers <- "515F-806R"
    head(metadata)

    ##    BioProject        Run sample_name Age BMI Sex            Platform Country          public_time        Site Status matched_samples Region
    ## 1 PRJNA553856 SRR9723171        AOW4  NA  NA  NA Illumina HiSeq 2000   China 2019-07-22T00:00:00Z oral cavity     HC              NO     V4
    ## 2 PRJNA553856 SRR9723172        AOW3  NA  NA  NA Illumina HiSeq 2000   China 2019-07-22T00:00:00Z oral cavity     HC              NO     V4
    ## 3 PRJNA553856 SRR9723173        AOW2  NA  NA  NA Illumina HiSeq 2000   China 2019-07-22T00:00:00Z oral cavity     HC              NO     V4
    ## 4 PRJNA553856 SRR9723174        AOW1  NA  NA  NA Illumina HiSeq 2000   China 2019-07-22T00:00:00Z oral cavity     HC              NO     V4
    ## 5 PRJNA553856 SRR9723175        AOM3  NA  NA  NA Illumina HiSeq 2000   China 2019-07-22T00:00:00Z oral cavity     HC              NO     V4
    ## 6 PRJNA553856 SRR9723176        AOM2  NA  NA  NA Illumina HiSeq 2000   China 2019-07-22T00:00:00Z oral cavity     HC              NO     V4
    ##                                                                      title Disease_study PCR_primers
    ## 1 Interaction of Oral and Toothbrush Microbiota Affects Oral Cavity Health            NO   515F-806R
    ## 2 Interaction of Oral and Toothbrush Microbiota Affects Oral Cavity Health            NO   515F-806R
    ## 3 Interaction of Oral and Toothbrush Microbiota Affects Oral Cavity Health            NO   515F-806R
    ## 4 Interaction of Oral and Toothbrush Microbiota Affects Oral Cavity Health            NO   515F-806R
    ## 5 Interaction of Oral and Toothbrush Microbiota Affects Oral Cavity Health            NO   515F-806R
    ## 6 Interaction of Oral and Toothbrush Microbiota Affects Oral Cavity Health            NO   515F-806R

# Step 4: Summarize and Validate the Metadata

    print(paste("Number of unique runs:", length(unique(metadata$Run))))

    ## [1] "Number of unique runs: 20"

    print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))

    ## [1] "Number of unique sample names: 20"

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
