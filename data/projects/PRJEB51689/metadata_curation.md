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

    setwd("C:/Users/zyc/Desktop/oral_project/data/projects/PRJEB51689")
    initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
    head(initial_meta)

    ##          Run Assay.Type AvgSpotLen    Bases BioProject     BioSample    Bytes             Center.Name Consent DATASTORE.filetype
    ## 1 ERR9289529   AMPLICON        601 79537025 PRJEB51689 SAMEA13421241 49839242 MEDICAL UNIVERSITY GRAZ  public   fastq,run.zq,sra
    ## 2 ERR9289530   AMPLICON        601 85097844 PRJEB51689 SAMEA13421242 53858318 MEDICAL UNIVERSITY GRAZ  public   fastq,run.zq,sra
    ## 3 ERR9289531   AMPLICON        601 99036522 PRJEB51689 SAMEA13421243 62022572 MEDICAL UNIVERSITY GRAZ  public   fastq,run.zq,sra
    ## 4 ERR9289532   AMPLICON        601 71033784 PRJEB51689 SAMEA13421244 44377889 MEDICAL UNIVERSITY GRAZ  public   fastq,run.zq,sra
    ## 5 ERR9289533   AMPLICON        601 64036401 PRJEB51689 SAMEA13421245 39909137 MEDICAL UNIVERSITY GRAZ  public   fastq,run.zq,sra
    ## 6 ERR9289534   AMPLICON        601 68113400 PRJEB51689 SAMEA13421246 43238365 MEDICAL UNIVERSITY GRAZ  public   fastq,run.zq,sra
    ##   DATASTORE.provider                     DATASTORE.region ENA.FIRST.PUBLIC..run. ENA.FIRST.PUBLIC ENA.LAST.UPDATE..run. ENA.LAST.UPDATE
    ## 1         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1             2023-01-03       2023-01-03            2023-01-03      2023-01-03
    ## 2         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1             2023-01-03       2023-01-03            2023-01-03      2023-01-03
    ## 3         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1             2023-01-03       2023-01-03            2023-01-03      2023-01-03
    ## 4         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1             2023-01-03       2023-01-03            2023-01-03      2023-01-03
    ## 5         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1             2023-01-03       2023-01-03            2023-01-03      2023-01-03
    ## 6         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1             2023-01-03       2023-01-03            2023-01-03      2023-01-03
    ##   Experiment   external_id       INSDC_center_name   INSDC_first_public    INSDC_last_update INSDC_status     Instrument Library.Name
    ## 1 ERX8831179 SAMEA13421241 MEDICAL UNIVERSITY GRAZ 2023-01-03T00:31:42Z 2023-01-03T00:31:42Z       public Illumina MiSeq          16S
    ## 2 ERX8831180 SAMEA13421242 MEDICAL UNIVERSITY GRAZ 2023-01-03T00:31:42Z 2023-01-03T00:31:42Z       public Illumina MiSeq          16S
    ## 3 ERX8831181 SAMEA13421243 MEDICAL UNIVERSITY GRAZ 2023-01-03T00:31:42Z 2023-01-03T00:31:42Z       public Illumina MiSeq          16S
    ## 4 ERX8831182 SAMEA13421244 MEDICAL UNIVERSITY GRAZ 2023-01-03T00:31:42Z 2023-01-03T00:31:42Z       public Illumina MiSeq          16S
    ## 5 ERX8831183 SAMEA13421245 MEDICAL UNIVERSITY GRAZ 2023-01-03T00:31:42Z 2023-01-03T00:31:42Z       public Illumina MiSeq          16S
    ## 6 ERX8831184 SAMEA13421246 MEDICAL UNIVERSITY GRAZ 2023-01-03T00:31:42Z 2023-01-03T00:31:42Z       public Illumina MiSeq          16S
    ##   LibraryLayout LibrarySelection LibrarySource             Organism Platform          ReleaseDate          create_date version   Sample.Name
    ## 1        PAIRED              PCR   METAGENOMIC uncultured bacterium ILLUMINA 2023-01-03T00:00:00Z 2023-01-10T02:30:00Z       1 SAMEA13421241
    ## 2        PAIRED              PCR   METAGENOMIC uncultured bacterium ILLUMINA 2023-01-03T00:00:00Z 2023-01-10T02:30:00Z       1 SAMEA13421242
    ## 3        PAIRED              PCR   METAGENOMIC uncultured bacterium ILLUMINA 2023-01-03T00:00:00Z 2023-01-10T02:30:00Z       1 SAMEA13421243
    ## 4        PAIRED              PCR   METAGENOMIC uncultured bacterium ILLUMINA 2023-01-03T00:00:00Z 2023-01-10T02:30:00Z       1 SAMEA13421244
    ## 5        PAIRED              PCR   METAGENOMIC uncultured bacterium ILLUMINA 2023-01-03T00:00:00Z 2023-01-10T02:30:00Z       1 SAMEA13421245
    ## 6        PAIRED              PCR   METAGENOMIC uncultured bacterium ILLUMINA 2023-01-03T00:00:00Z 2023-01-10T02:30:00Z       1 SAMEA13421246
    ##   Sample_name      scientific_name SRA.Study Submitter_Id      INSDC_center_alias
    ## 1    001Co124 uncultured bacterium ERP136342     001Co124 MEDICAL UNIVERSITY GRAZ
    ## 2    001Co185 uncultured bacterium ERP136342     001Co185 MEDICAL UNIVERSITY GRAZ
    ## 3    001Co367 uncultured bacterium ERP136342     001Co367 MEDICAL UNIVERSITY GRAZ
    ## 4     001GYC5 uncultured bacterium ERP136342      001GYC5 MEDICAL UNIVERSITY GRAZ
    ## 5     001GYD1 uncultured bacterium ERP136342      001GYD1 MEDICAL UNIVERSITY GRAZ
    ## 6   001GYpreR uncultured bacterium ERP136342    001GYpreR MEDICAL UNIVERSITY GRAZ

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
      sample_name = initial_meta$Sample_name,
      Age = NA,
      BMI = NA,
      Sex = NA,
      Platform = initial_meta$Instrument,
      # Country = initial_meta$geo_loc_name_country,
      # Site = initial_meta$isolation_source,
      public_time = initial_meta$ReleaseDate,
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

    metadata$Status <- "non-metastatic breast cancer"
    metadata$Site <- "saliva"
    metadata$Country <- "USA"
    metadata$matched_samples <- "NO"
    metadata$Region <- "V4-V5"
    metadata$title <- "Chemotherapy-associated oral microbiome changes in breast cancer patients"
    metadata$Disease_study <- "YES"
    metadata$PCR_primers <- "515F-926R"
    head(metadata)

    ##   BioProject        Run sample_name Age BMI Sex       Platform          public_time                       Status   Site Country
    ## 1 PRJEB51689 ERR9289529    001Co124  NA  NA  NA Illumina MiSeq 2023-01-03T00:00:00Z non-metastatic breast cancer saliva     USA
    ## 2 PRJEB51689 ERR9289530    001Co185  NA  NA  NA Illumina MiSeq 2023-01-03T00:00:00Z non-metastatic breast cancer saliva     USA
    ## 3 PRJEB51689 ERR9289531    001Co367  NA  NA  NA Illumina MiSeq 2023-01-03T00:00:00Z non-metastatic breast cancer saliva     USA
    ## 4 PRJEB51689 ERR9289532     001GYC5  NA  NA  NA Illumina MiSeq 2023-01-03T00:00:00Z non-metastatic breast cancer saliva     USA
    ## 5 PRJEB51689 ERR9289533     001GYD1  NA  NA  NA Illumina MiSeq 2023-01-03T00:00:00Z non-metastatic breast cancer saliva     USA
    ## 6 PRJEB51689 ERR9289534   001GYpreR  NA  NA  NA Illumina MiSeq 2023-01-03T00:00:00Z non-metastatic breast cancer saliva     USA
    ##   matched_samples Region                                                                     title Disease_study PCR_primers
    ## 1              NO  V4-V5 Chemotherapy-associated oral microbiome changes in breast cancer patients           YES   515F-926R
    ## 2              NO  V4-V5 Chemotherapy-associated oral microbiome changes in breast cancer patients           YES   515F-926R
    ## 3              NO  V4-V5 Chemotherapy-associated oral microbiome changes in breast cancer patients           YES   515F-926R
    ## 4              NO  V4-V5 Chemotherapy-associated oral microbiome changes in breast cancer patients           YES   515F-926R
    ## 5              NO  V4-V5 Chemotherapy-associated oral microbiome changes in breast cancer patients           YES   515F-926R
    ## 6              NO  V4-V5 Chemotherapy-associated oral microbiome changes in breast cancer patients           YES   515F-926R

# Step 4: Summarize and Validate the Metadata

    print(paste("Number of unique runs:", length(unique(metadata$Run))))

    ## [1] "Number of unique runs: 108"

    print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))

    ## [1] "Number of unique sample names: 108"

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
