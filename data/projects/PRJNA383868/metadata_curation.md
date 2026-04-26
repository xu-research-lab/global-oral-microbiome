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

    setwd("C:/Users/zyc/Desktop/oral_project/data/projects/PRJNA383868")
    initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
    head(initial_meta)

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample                       BioSampleModel    Bytes
    ## 1 SRR5467622   AMPLICON        502 12991258 PRJNA383868 SAMN06818328 MIMS.me,MIGS/MIMS/MIMARKS.human-oral  8392295
    ## 2 SRR5467623   AMPLICON        502 21643228 PRJNA383868 SAMN06818327 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 13831072
    ## 3 SRR5467624   AMPLICON        502 21135706 PRJNA383868 SAMN06818326 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 13491196
    ## 4 SRR5467625   AMPLICON        502  7511928 PRJNA383868 SAMN06818325 MIMS.me,MIGS/MIMS/MIMARKS.human-oral  4718731
    ## 5 SRR5467626   AMPLICON        502  6295582 PRJNA383868 SAMN06818324 MIMS.me,MIGS/MIMS/MIMARKS.human-oral  3970507
    ## 6 SRR5467627   AMPLICON        502 15582582 PRJNA383868 SAMN06818323 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 10037418
    ##                 Center.Name chem_administration Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1 J. CRAIG VENTER INSTITUTE      not applicable      2014-10-27  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2 J. CRAIG VENTER INSTITUTE      not applicable      2014-11-24  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3 J. CRAIG VENTER INSTITUTE      not applicable      2014-11-24  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4 J. CRAIG VENTER INSTITUTE      not applicable      2014-06-25  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5 J. CRAIG VENTER INSTITUTE      not applicable      2014-06-25  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6 J. CRAIG VENTER INSTITUTE      not applicable      2014-12-17  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region     env_biome   env_feature  env_material Experiment geo_loc_name_country
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 ENVO:01000418 ENVO:00002034 ENVO:01000156 SRX2752978            Australia
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 ENVO:01000418 ENVO:00002034 ENVO:01000156 SRX2752979            Australia
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 ENVO:01000418 ENVO:00002034 ENVO:01000156 SRX2752980            Australia
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 ENVO:01000418 ENVO:00002034 ENVO:01000156 SRX2752981            Australia
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 ENVO:01000418 ENVO:00002034 ENVO:01000156 SRX2752982            Australia
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 ENVO:01000418 ENVO:00002034 ENVO:01000156 SRX2752983            Australia
    ##   geo_loc_name_country_continent                geo_loc_name  Host_Age       Host_disease host_family_relationship
    ## 1                        Oceania Australia: Melbourne\\, VIC 78 months      dental caries           sibling 2028.2
    ## 2                        Oceania Australia: Melbourne\\, VIC 79 months      dental caries           sibling 2026.1
    ## 3                        Oceania Australia: Melbourne\\, VIC 79 months      dental caries           sibling 2026.2
    ## 4                        Oceania Australia: Melbourne\\, VIC 77 months dental caries free           sibling 2023.1
    ## 5                        Oceania Australia: Melbourne\\, VIC 77 months dental caries free           sibling 2023.2
    ## 6                        Oceania Australia: Melbourne\\, VIC 80 months dental caries free           sibling 2020.1
    ##           HOST host_sex host_subject_id host_tissue_sampled     Instrument isolation_source                 lat_lon
    ## 1 Homo sapiens   female          2028.1       dental plaque Illumina MiSeq            teeth 37.793742 S 144.94921 E
    ## 2 Homo sapiens   female          2026.2       dental plaque Illumina MiSeq            teeth 37.793742 S 144.94921 E
    ## 3 Homo sapiens   female          2026.1       dental plaque Illumina MiSeq            teeth 37.793742 S 144.94921 E
    ## 4 Homo sapiens     male          2023.2       dental plaque Illumina MiSeq            teeth 37.793742 S 144.94921 E
    ## 5 Homo sapiens     male          2023.1       dental plaque Illumina MiSeq            teeth 37.793742 S 144.94921 E
    ## 6 Homo sapiens     male          2020.2       dental plaque Illumina MiSeq            teeth 37.793742 S 144.94921 E
    ##     Library.Name LibraryLayout LibrarySelection LibrarySource nose_mouth_teeth_throat_disord              Organism
    ## 1 CAR_2028.1_RD1        PAIRED              PCR   METAGENOMIC                 not applicable human oral metagenome
    ## 2 CAR_2026.2_RD1        PAIRED              PCR   METAGENOMIC                          cleft human oral metagenome
    ## 3 CAR_2026.1_RD1        PAIRED              PCR   METAGENOMIC                 not applicable human oral metagenome
    ## 4 CAR_2023.2_RD1        PAIRED              PCR   METAGENOMIC                 not applicable human oral metagenome
    ## 5 CAR_2023.1_RD1        PAIRED              PCR   METAGENOMIC                 not applicable human oral metagenome
    ## 6 CAR_2020.2_RD1        PAIRED              PCR   METAGENOMIC                 not applicable human oral metagenome
    ##   Platform          ReleaseDate samp_collect_device    samp_store_temp samp_vol_we_dna_ext Sample.Name SRA.Study
    ## 1 ILLUMINA 2017-04-22T00:00:00Z          oral swabs -80 degree Celsius              750 uL  2028.1_RD1 SRP104798
    ## 2 ILLUMINA 2017-04-22T00:00:00Z          oral swabs -80 degree Celsius              750 uL  2026.2_RD1 SRP104798
    ## 3 ILLUMINA 2017-04-22T00:00:00Z          oral swabs -80 degree Celsius              750 uL  2026.1_RD1 SRP104798
    ## 4 ILLUMINA 2017-04-22T00:00:00Z          oral swabs -80 degree Celsius              750 uL  2023.2_RD1 SRP104798
    ## 5 ILLUMINA 2017-04-22T00:00:00Z          oral swabs -80 degree Celsius              750 uL  2023.1_RD1 SRP104798
    ## 6 ILLUMINA 2017-04-22T00:00:00Z          oral swabs -80 degree Celsius              750 uL  2020.2_RD1 SRP104798
    ##            create_date version  host_body_site host_body_product misc_param source_material_id
    ## 1 2017-04-22T00:07:00Z       1 right dentition                                                
    ## 2 2017-04-22T00:07:00Z       1 right dentition                                                
    ## 3 2017-04-22T00:07:00Z       1 right dentition                                                
    ## 4 2017-04-22T00:06:00Z       1 right dentition                                                
    ## 5 2017-04-22T00:07:00Z       1 right dentition                                                
    ## 6 2017-04-22T00:07:00Z       1 right dentition

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
      assaytype = initial_meta$Assay.Type,
      BioProject = initial_meta$BioProject,
      Run = initial_meta$Run,
      sample_name = initial_meta$Library.Name,
      Age = initial_meta$Host_Age, 
      BMI = NA,  
      Sex = initial_meta$host_sex,
      Region = "V4",
      Platform = initial_meta$Instrument,
      Country = initial_meta$geo_loc_name_country,
      Site = initial_meta$host_tissue_sampled,
      public_time = initial_meta$ReleaseDate,
      matched_samples = "NO",
      Status = initial_meta$Host_disease,
      stringsAsFactors = FALSE
    )

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

    metadata <- metadata %>%filter(assaytype == "AMPLICON")
    metadata <- metadata %>%mutate(Status = case_when(
      Status == "dental caries free" ~ "HC",
      TRUE ~ "dental caries"
    ))
    metadata <- metadata %>%
      mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
    metadata <- metadata %>%mutate(Age = as.numeric(sapply(str_split(Age, " "), function(x) head(x, 1))))%>%mutate(Age = Age/12)

    metadata$title <- "Host Genetic Control of the Oral Microbiome in Health and Disease"
    metadata$Disease_study <- "YES"
    metadata$PCR_primers <- "515F-806R"

    head(metadata)

    ##   assaytype  BioProject        Run    sample_name      Age BMI Sex Region       Platform   Country          Site
    ## 1  AMPLICON PRJNA383868 SRR5467622 CAR_2028.1_RD1 6.500000  NA   F     V4 Illumina MiSeq Australia dental plaque
    ## 2  AMPLICON PRJNA383868 SRR5467623 CAR_2026.2_RD1 6.583333  NA   F     V4 Illumina MiSeq Australia dental plaque
    ## 3  AMPLICON PRJNA383868 SRR5467624 CAR_2026.1_RD1 6.583333  NA   F     V4 Illumina MiSeq Australia dental plaque
    ## 4  AMPLICON PRJNA383868 SRR5467625 CAR_2023.2_RD1 6.416667  NA   M     V4 Illumina MiSeq Australia dental plaque
    ## 5  AMPLICON PRJNA383868 SRR5467626 CAR_2023.1_RD1 6.416667  NA   M     V4 Illumina MiSeq Australia dental plaque
    ## 6  AMPLICON PRJNA383868 SRR5467627 CAR_2020.2_RD1 6.666667  NA   M     V4 Illumina MiSeq Australia dental plaque
    ##            public_time matched_samples        Status
    ## 1 2017-04-22T00:00:00Z              NO dental caries
    ## 2 2017-04-22T00:00:00Z              NO dental caries
    ## 3 2017-04-22T00:00:00Z              NO dental caries
    ## 4 2017-04-22T00:00:00Z              NO            HC
    ## 5 2017-04-22T00:00:00Z              NO            HC
    ## 6 2017-04-22T00:00:00Z              NO            HC
    ##                                                               title Disease_study PCR_primers
    ## 1 Host Genetic Control of the Oral Microbiome in Health and Disease           YES   515F-806R
    ## 2 Host Genetic Control of the Oral Microbiome in Health and Disease           YES   515F-806R
    ## 3 Host Genetic Control of the Oral Microbiome in Health and Disease           YES   515F-806R
    ## 4 Host Genetic Control of the Oral Microbiome in Health and Disease           YES   515F-806R
    ## 5 Host Genetic Control of the Oral Microbiome in Health and Disease           YES   515F-806R
    ## 6 Host Genetic Control of the Oral Microbiome in Health and Disease           YES   515F-806R

# Step 4: Summarize and Validate the Metadata

    print(paste("Number of unique runs:", length(unique(metadata$Run))))

    ## [1] "Number of unique runs: 540"

    print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))

    ## [1] "Number of unique sample names: 540"

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
            theme(axis.text.x = element_text(angle = 30, hjust = 1))
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
    ggsave("metadata_summary.png", plot = p, width = 8, height = 8, dpi = 150)

    knitr::include_graphics("metadata_summary.png")

![](metadata_summary.png)

# Step 6: Finalize and Export the Metadata

    metadata <- metadata[, c(
      "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
      "Region", "Platform", "Country", "Site", "public_time", 
      "PCR_primers", "title", "Disease_study", "Status","matched_samples"
    )]

    write.csv(metadata, file = "metadata.csv", quote = FALSE)
