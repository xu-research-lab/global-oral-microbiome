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

    setwd("C:/Users/zyc/Desktop/oral_project/data/projects/PRJNA555622")
    initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
    head(initial_meta)

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample              BioSampleModel    Bytes               Center.Name
    ## 1 SRR9727013   AMPLICON        600 25494719 PRJNA555622 SAMN12324795 Metagenome or environmental 18467971 UNIVERSITY OF CONNECTICUT
    ## 2 SRR9727014   AMPLICON        600 67829409 PRJNA555622 SAMN12324796 Metagenome or environmental 48716396 UNIVERSITY OF CONNECTICUT
    ## 3 SRR9727015   AMPLICON        600 63297026 PRJNA555622 SAMN12324793 Metagenome or environmental 45383418 UNIVERSITY OF CONNECTICUT
    ## 4 SRR9727016   AMPLICON        600 62271299 PRJNA555622 SAMN12324794 Metagenome or environmental 44601290 UNIVERSITY OF CONNECTICUT
    ## 5 SRR9727017   AMPLICON        600 55261546 PRJNA555622 SAMN12324791 Metagenome or environmental 39674354 UNIVERSITY OF CONNECTICUT
    ## 6 SRR9727018   AMPLICON        600 45523636 PRJNA555622 SAMN12324792 Metagenome or environmental 32674536 UNIVERSITY OF CONNECTICUT
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider                     DATASTORE.region Experiment geo_loc_name_country
    ## 1              NA  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6484551         uncalculated
    ## 2              NA  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6484550         uncalculated
    ## 3              NA  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6484549         uncalculated
    ## 4              NA  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6484548         uncalculated
    ## 5              NA  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6484547         uncalculated
    ## 6              NA  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6484546         uncalculated
    ##   geo_loc_name_country_continent geo_loc_name         HOST     Instrument lat_lon Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1                   uncalculated           NA Homo sapiens Illumina MiSeq      NA           C4        PAIRED              PCR   METAGENOMIC
    ## 2                   uncalculated           NA Homo sapiens Illumina MiSeq      NA           C5        PAIRED              PCR   METAGENOMIC
    ## 3                   uncalculated           NA Homo sapiens Illumina MiSeq      NA           C2        PAIRED              PCR   METAGENOMIC
    ## 4                   uncalculated           NA Homo sapiens Illumina MiSeq      NA           C3        PAIRED              PCR   METAGENOMIC
    ## 5                   uncalculated           NA Homo sapiens Illumina MiSeq      NA          B20        PAIRED              PCR   METAGENOMIC
    ## 6                   uncalculated           NA Homo sapiens Illumina MiSeq      NA           C1        PAIRED              PCR   METAGENOMIC
    ##                Organism Platform          ReleaseDate          create_date version Sample.Name source_material_id SRA.Study
    ## 1 human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-22T12:50:00Z       1          C4                 49 SRP215893
    ## 2 human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-22T12:50:00Z       1          C5                 50 SRP215893
    ## 3 human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-22T12:50:00Z       1          C2                 47 SRP215893
    ## 4 human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-22T12:50:00Z       1          C3                 48 SRP215893
    ## 5 human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-22T12:50:00Z       1         B20                 45 SRP215893
    ## 6 human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-22T12:50:00Z       1          C1                 46 SRP215893

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
      # Country = initial_meta$geo_loc_name_country,
      # Status = initial_meta$Host_disease,
      # Site = initial_meta$isolation_source,
      public_time = initial_meta$ReleaseDate,
      stringsAsFactors = FALSE
    )

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

     # metadata <- metadata %>%
     #  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
     # metadata <- metadata %>%mutate(Age = Age/12)

    metadata$Country = "America"
    metadata$Site = "Supragingival plaques"
    metadata$Status = NA
    metadata$matched_samples = "NO"
    metadata$Region = "V3-V4"
    metadata$title <- "Ethnic variation of oral microbiota in children"
    metadata$Disease_study <- "NO"
    metadata$PCR_primers <- "F:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG;R:GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC"
    head(metadata)

    ##    BioProject        Run sample_name Age BMI Sex       Platform          public_time Country                  Site Status matched_samples
    ## 1 PRJNA555622 SRR9727013          C4  NA  NA  NA Illumina MiSeq 2019-12-31T00:00:00Z America Supragingival plaques     NA              NO
    ## 2 PRJNA555622 SRR9727014          C5  NA  NA  NA Illumina MiSeq 2019-12-31T00:00:00Z America Supragingival plaques     NA              NO
    ## 3 PRJNA555622 SRR9727015          C2  NA  NA  NA Illumina MiSeq 2019-12-31T00:00:00Z America Supragingival plaques     NA              NO
    ## 4 PRJNA555622 SRR9727016          C3  NA  NA  NA Illumina MiSeq 2019-12-31T00:00:00Z America Supragingival plaques     NA              NO
    ## 5 PRJNA555622 SRR9727017         B20  NA  NA  NA Illumina MiSeq 2019-12-31T00:00:00Z America Supragingival plaques     NA              NO
    ## 6 PRJNA555622 SRR9727018          C1  NA  NA  NA Illumina MiSeq 2019-12-31T00:00:00Z America Supragingival plaques     NA              NO
    ##   Region                                           title Disease_study
    ## 1  V3-V4 Ethnic variation of oral microbiota in children            NO
    ## 2  V3-V4 Ethnic variation of oral microbiota in children            NO
    ## 3  V3-V4 Ethnic variation of oral microbiota in children            NO
    ## 4  V3-V4 Ethnic variation of oral microbiota in children            NO
    ## 5  V3-V4 Ethnic variation of oral microbiota in children            NO
    ## 6  V3-V4 Ethnic variation of oral microbiota in children            NO
    ##                                                                                                      PCR_primers
    ## 1 F:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG;R:GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC
    ## 2 F:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG;R:GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC
    ## 3 F:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG;R:GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC
    ## 4 F:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG;R:GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC
    ## 5 F:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG;R:GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC
    ## 6 F:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG;R:GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC

# Step 4: Summarize and Validate the Metadata

    print(paste("Number of unique runs:", length(unique(metadata$Run))))

    ## [1] "Number of unique runs: 96"

    print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))

    ## [1] "Number of unique sample names: 96"

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
              # plots$PCR_primers,
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
