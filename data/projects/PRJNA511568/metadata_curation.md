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
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.3     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggplot2)
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

``` r
setwd("/home/dongbiao/word_embedding_microbiome/all_data/oral/global-oral-microbiome/data/projects/PRJNA511568")
initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR8497113   AMPLICON        473 30485239 PRJNA511568 SAMN10834315
    ## 2 SRR8497114   AMPLICON        473 29216085 PRJNA511568 SAMN10834316
    ## 3 SRR8497115   AMPLICON        473 31260563 PRJNA511568 SAMN10834317
    ## 4 SRR8497116   AMPLICON        475 34187713 PRJNA511568 SAMN10834318
    ## 5 SRR8497117   AMPLICON        474 35750629 PRJNA511568 SAMN10834311
    ## 6 SRR8497118   AMPLICON        474 33073169 PRJNA511568 SAMN10834312
    ##                BioSampleModel    Bytes            Center.Name Collection_Date
    ## 1 Metagenome or environmental 19587267 SUN YAT-SEN UNIVERSITY      2016-12-24
    ## 2 Metagenome or environmental 18722208 SUN YAT-SEN UNIVERSITY      2016-12-24
    ## 3 Metagenome or environmental 19895407 SUN YAT-SEN UNIVERSITY      2016-12-24
    ## 4 Metagenome or environmental 21617625 SUN YAT-SEN UNIVERSITY      2016-12-24
    ## 5 Metagenome or environmental 22708387 SUN YAT-SEN UNIVERSITY      2016-12-24
    ## 6 Metagenome or environmental 20928425 SUN YAT-SEN UNIVERSITY      2016-12-24
    ##   Consent DATASTORE.filetype DATASTORE.provider
    ## 1  public         run.zq,sra         gs,ncbi,s3
    ## 2  public         run.zq,sra         gs,ncbi,s3
    ## 3  public         run.zq,sra         gs,ncbi,s3
    ## 4  public         run.zq,sra         gs,ncbi,s3
    ## 5  public         run.zq,sra         gs,ncbi,s3
    ## 6  public         run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region Experiment geo_loc_name_country
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX5301186         uncalculated
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX5301185         uncalculated
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX5301184         uncalculated
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX5301183         uncalculated
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX5301182         uncalculated
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX5301181         uncalculated
    ##   geo_loc_name_country_continent   geo_loc_name              HOST
    ## 1                   uncalculated not applicable Homo sapiens No.5
    ## 2                   uncalculated not applicable Homo sapiens No.6
    ## 3                   uncalculated not applicable Homo sapiens No.7
    ## 4                   uncalculated not applicable Homo sapiens No.8
    ## 5                   uncalculated not applicable Homo sapiens No.1
    ## 6                   uncalculated not applicable Homo sapiens No.2
    ##       Instrument isolation_source        lat_lon                   Library.Name
    ## 1 Illumina MiSeq    dental plaque not applicable dental plaque 16S rRNA 2016 C5
    ## 2 Illumina MiSeq    dental plaque not applicable dental plaque 16S rRNA 2016 C6
    ## 3 Illumina MiSeq    dental plaque not applicable dental plaque 16S rRNA 2016 C7
    ## 4 Illumina MiSeq    dental plaque not applicable dental plaque 16S rRNA 2016 C8
    ## 5 Illumina MiSeq    dental plaque not applicable dental plaque 16S rRNA 2016 C1
    ## 6 Illumina MiSeq    dental plaque not applicable dental plaque 16S rRNA 2016 C2
    ##   LibraryLayout LibrarySelection LibrarySource              Organism Platform
    ## 1        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 2        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 3        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 4        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 5        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 6        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ##            ReleaseDate          create_date version Sample.Name SRA.Study
    ## 1 2020-02-28T00:00:00Z 2019-01-28T13:58:00Z       1          C5 SRP182703
    ## 2 2020-02-28T00:00:00Z 2019-01-28T13:58:00Z       1          C6 SRP182703
    ## 3 2020-02-28T00:00:00Z 2019-01-28T13:58:00Z       1          C7 SRP182703
    ## 4 2020-02-28T00:00:00Z 2019-01-28T13:58:00Z       1          C8 SRP182703
    ## 5 2020-02-28T00:00:00Z 2019-01-28T13:58:00Z       1          C1 SRP182703
    ## 6 2020-02-28T00:00:00Z 2019-01-28T13:58:00Z       1          C2 SRP182703

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
  Platform = initial_meta$Instrument,
  # Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$isolation_source,
  public_time = initial_meta$ReleaseDate,
  # Status = initial_meta$Sample.Name,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
 # metadata <- metadata %>%
 #  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
# metadata <-metadata%>%
#   mutate(Status = case_when(Site == "mucous" ~ "OSCC",
#                             Site == "dental plaque" ~ "OSCC",
#                             Site == "tumor" ~ "OSCC",
#                             TRUE ~ NA))
```

``` r
metadata$Country <- "China"
metadata$Status <- "HC"
metadata$matched_samples <- "NO"
metadata$Region <- "V4-V5"
metadata$title <- "Effects of a probiotic drink containing Lactobacillus casei strain Shirota on dental plaque microbiota"
metadata$Disease_study <- "NO"
metadata$PCR_primers <- "(forward, 5'-GTGACAGCAGCCGCGGTAA-3' and reverse, 5'-CCGTCAATTTATTTAAGTTT-3')"
head(metadata)
```

    ##    BioProject        Run sample_name Age BMI Sex       Platform          Site
    ## 1 PRJNA511568 SRR8497113          C5  NA  NA  NA Illumina MiSeq dental plaque
    ## 2 PRJNA511568 SRR8497114          C6  NA  NA  NA Illumina MiSeq dental plaque
    ## 3 PRJNA511568 SRR8497115          C7  NA  NA  NA Illumina MiSeq dental plaque
    ## 4 PRJNA511568 SRR8497116          C8  NA  NA  NA Illumina MiSeq dental plaque
    ## 5 PRJNA511568 SRR8497117          C1  NA  NA  NA Illumina MiSeq dental plaque
    ## 6 PRJNA511568 SRR8497118          C2  NA  NA  NA Illumina MiSeq dental plaque
    ##            public_time Country Status matched_samples Region
    ## 1 2020-02-28T00:00:00Z   China     HC              NO  V4-V5
    ## 2 2020-02-28T00:00:00Z   China     HC              NO  V4-V5
    ## 3 2020-02-28T00:00:00Z   China     HC              NO  V4-V5
    ## 4 2020-02-28T00:00:00Z   China     HC              NO  V4-V5
    ## 5 2020-02-28T00:00:00Z   China     HC              NO  V4-V5
    ## 6 2020-02-28T00:00:00Z   China     HC              NO  V4-V5
    ##                                                                                                    title
    ## 1 Effects of a probiotic drink containing Lactobacillus casei strain Shirota on dental plaque microbiota
    ## 2 Effects of a probiotic drink containing Lactobacillus casei strain Shirota on dental plaque microbiota
    ## 3 Effects of a probiotic drink containing Lactobacillus casei strain Shirota on dental plaque microbiota
    ## 4 Effects of a probiotic drink containing Lactobacillus casei strain Shirota on dental plaque microbiota
    ## 5 Effects of a probiotic drink containing Lactobacillus casei strain Shirota on dental plaque microbiota
    ## 6 Effects of a probiotic drink containing Lactobacillus casei strain Shirota on dental plaque microbiota
    ##   Disease_study
    ## 1            NO
    ## 2            NO
    ## 3            NO
    ## 4            NO
    ## 5            NO
    ## 6            NO
    ##                                                                    PCR_primers
    ## 1 (forward, 5'-GTGACAGCAGCCGCGGTAA-3' and reverse, 5'-CCGTCAATTTATTTAAGTTT-3')
    ## 2 (forward, 5'-GTGACAGCAGCCGCGGTAA-3' and reverse, 5'-CCGTCAATTTATTTAAGTTT-3')
    ## 3 (forward, 5'-GTGACAGCAGCCGCGGTAA-3' and reverse, 5'-CCGTCAATTTATTTAAGTTT-3')
    ## 4 (forward, 5'-GTGACAGCAGCCGCGGTAA-3' and reverse, 5'-CCGTCAATTTATTTAAGTTT-3')
    ## 5 (forward, 5'-GTGACAGCAGCCGCGGTAA-3' and reverse, 5'-CCGTCAATTTATTTAAGTTT-3')
    ## 6 (forward, 5'-GTGACAGCAGCCGCGGTAA-3' and reverse, 5'-CCGTCAATTTATTTAAGTTT-3')

# Step 4: Summarize and Validate the Metadata

``` r
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 30"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 30"

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
ggsave("Metadata_Summary.png", plot = p, width = 8, height = 6, dpi = 150)
```

``` r
knitr::include_graphics("Metadata_Summary.png")
```

<img src="Metadata_Summary.png" width="1200" /> \# Step 6: Finalize and
Export the Metadata

``` r
metadata <- metadata[, c(
  "BioProject", "Run", "sample_name", "Age", "BMI", "Sex",
  "Region", "Platform", "Country", "Site", "public_time",
  "PCR_primers", "title", "Disease_study", "Status","matched_samples"
)]



write.csv(metadata, file = "metadata.csv", quote = TRUE)
```
