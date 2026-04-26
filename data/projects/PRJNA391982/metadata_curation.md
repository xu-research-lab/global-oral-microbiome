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

    setwd("C:/Users/zyc/Desktop/oral_project/data/projects/PRJNA391982")
    initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
    head(initial_meta)

    ##          Run Assay.Type AvgSpotLen   Bases  BioProject    BioSample              BioSampleModel   Bytes               Center.Name
    ## 1 SRR5758886   AMPLICON        472 4369596 PRJNA391982 SAMN07279873 Metagenome or environmental 2624008 THE OHIO STATE UNIVERSITY
    ## 2 SRR5758887   AMPLICON        464 3379245 PRJNA391982 SAMN07279876 Metagenome or environmental 2061569 THE OHIO STATE UNIVERSITY
    ## 3 SRR5758888   AMPLICON        445 2043533 PRJNA391982 SAMN07279877 Metagenome or environmental 1339117 THE OHIO STATE UNIVERSITY
    ## 4 SRR5758889   AMPLICON        473 1629283 PRJNA391982 SAMN07279916 Metagenome or environmental 1031222 THE OHIO STATE UNIVERSITY
    ## 5 SRR5758890   AMPLICON        430 1463810 PRJNA391982 SAMN07279917 Metagenome or environmental  903723 THE OHIO STATE UNIVERSITY
    ## 6 SRR5758891   AMPLICON        467 4764347 PRJNA391982 SAMN07279908 Metagenome or environmental 2849701 THE OHIO STATE UNIVERSITY
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider                     DATASTORE.region Experiment geo_loc_name_country
    ## 1   not available  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX2958844                  USA
    ## 2   not available  public   sra,fastq,run.zq         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX2958843                  USA
    ## 3   not available  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX2958842                  USA
    ## 4   not available  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX2958841                  USA
    ## 5   not available  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX2958840                  USA
    ## 6   not available  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX2958839                  USA
    ##   geo_loc_name_country_continent        geo_loc_name         HOST          Instrument isolation_source             lat_lon Library.Name
    ## 1                  North America USA: Columbus\\, OH Homo sapiens 454 GS FLX Titanium mother-permanent 39.9612 N 82.9988 W  239m saliva
    ## 2                  North America USA: Columbus\\, OH Homo sapiens 454 GS FLX Titanium  child_permanent 39.9612 N 82.9988 W  241c saliva
    ## 3                  North America USA: Columbus\\, OH Homo sapiens 454 GS FLX Titanium mother-permanent 39.9612 N 82.9988 W  241m saliva
    ## 4                  North America USA: Columbus\\, OH Homo sapiens 454 GS FLX Titanium  child_permanent 39.9612 N 82.9988 W  261c saliva
    ## 5                  North America USA: Columbus\\, OH Homo sapiens 454 GS FLX Titanium mother-permanent 39.9612 N 82.9988 W  261m saliva
    ## 6                  North America USA: Columbus\\, OH Homo sapiens 454 GS FLX Titanium child_predentate 39.9612 N 82.9988 W  257c mucosa
    ##   LibraryLayout LibrarySelection LibrarySource              Organism Platform          ReleaseDate          create_date version sample.
    ## 1        SINGLE           RANDOM   METAGENOMIC human oral metagenome    LS454 2018-07-31T00:00:00Z 2017-06-27T03:57:00Z       1      76
    ## 2        SINGLE           RANDOM   METAGENOMIC human oral metagenome    LS454 2018-07-31T00:00:00Z 2017-06-27T03:57:00Z       1      79
    ## 3        SINGLE           RANDOM   METAGENOMIC human oral metagenome    LS454 2018-07-31T00:00:00Z 2017-06-27T03:57:00Z       1      80
    ## 4        SINGLE           RANDOM   METAGENOMIC human oral metagenome    LS454 2018-07-31T00:00:00Z 2017-06-27T03:57:00Z       1     119
    ## 5        SINGLE           RANDOM   METAGENOMIC human oral metagenome    LS454 2018-07-31T00:00:00Z 2017-06-27T03:57:00Z       1     120
    ## 6        SINGLE           RANDOM   METAGENOMIC human oral metagenome    LS454 2018-07-31T00:00:00Z 2017-06-27T03:57:00Z       1     111
    ##   Sample.Name SRA.Study
    ## 1        239m SRP110543
    ## 2        241c SRP110543
    ## 3        241m SRP110543
    ## 4        261c SRP110543
    ## 5        261m SRP110543
    ## 6        257c SRP110543

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
      Site = initial_meta$Library.Name,
      public_time = initial_meta$ReleaseDate,
      # Status = initial_meta$Sample.Name,
      stringsAsFactors = FALSE
    )

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

    metadata <- metadata %>% mutate(Site = case_when(
      grepl("saliva",Site) ~ "saliva",
      grepl("mucosa",Site) ~ "mucosa"
    )) 

    # metadata <- metadata %>%
     #  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
    # metadata <-metadata%>%
    #   mutate(Status = case_when(Site == "mucous" ~ "OSCC",
    #                             Site == "dental plaque" ~ "OSCC",
    #                             Site == "tumor" ~ "OSCC",
    #                             TRUE ~ NA))

    metadata$Status <- "HC"
    metadata$matched_samples <- "NO"
    metadata$Region <- "V1-V3/V7-V9"
    metadata$title <- "Characterizing oral microbial communities across dentition states and colonization niches"
    metadata$Disease_study <- "NO"
    metadata$PCR_primers <- NA
    head(metadata)

    ##    BioProject        Run sample_name Age BMI Sex            Platform Country   Site          public_time Status matched_samples      Region
    ## 1 PRJNA391982 SRR5758886        239m  NA  NA  NA 454 GS FLX Titanium     USA saliva 2018-07-31T00:00:00Z     HC              NO V1-V3/V7-V9
    ## 2 PRJNA391982 SRR5758887        241c  NA  NA  NA 454 GS FLX Titanium     USA saliva 2018-07-31T00:00:00Z     HC              NO V1-V3/V7-V9
    ## 3 PRJNA391982 SRR5758888        241m  NA  NA  NA 454 GS FLX Titanium     USA saliva 2018-07-31T00:00:00Z     HC              NO V1-V3/V7-V9
    ## 4 PRJNA391982 SRR5758889        261c  NA  NA  NA 454 GS FLX Titanium     USA saliva 2018-07-31T00:00:00Z     HC              NO V1-V3/V7-V9
    ## 5 PRJNA391982 SRR5758890        261m  NA  NA  NA 454 GS FLX Titanium     USA saliva 2018-07-31T00:00:00Z     HC              NO V1-V3/V7-V9
    ## 6 PRJNA391982 SRR5758891        257c  NA  NA  NA 454 GS FLX Titanium     USA mucosa 2018-07-31T00:00:00Z     HC              NO V1-V3/V7-V9
    ##                                                                                       title Disease_study PCR_primers
    ## 1 Characterizing oral microbial communities across dentition states and colonization niches            NO          NA
    ## 2 Characterizing oral microbial communities across dentition states and colonization niches            NO          NA
    ## 3 Characterizing oral microbial communities across dentition states and colonization niches            NO          NA
    ## 4 Characterizing oral microbial communities across dentition states and colonization niches            NO          NA
    ## 5 Characterizing oral microbial communities across dentition states and colonization niches            NO          NA
    ## 6 Characterizing oral microbial communities across dentition states and colonization niches            NO          NA

# Step 4: Summarize and Validate the Metadata

    print(paste("Number of unique runs:", length(unique(metadata$Run))))

    ## [1] "Number of unique runs: 240"

    print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))

    ## [1] "Number of unique sample names: 120"

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
