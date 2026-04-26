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
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(cowplot)

# Set working directory (update this path to your local directory)
setwd("E:/bioproject/PRJNA1149521")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = T) 

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases   BioProject    BioSample
    ## 1 SRR30305021   AMPLICON        602 52281894 PRJNA1149521 SAMN43250877
    ## 2 SRR30305022   AMPLICON        602 33502504 PRJNA1149521 SAMN43250876
    ## 3 SRR30305023   AMPLICON        593 27260117 PRJNA1149521 SAMN43250651
    ## 4 SRR30305024   AMPLICON        602 49908810 PRJNA1149521 SAMN43250875
    ## 5 SRR30305025   AMPLICON        602 60962132 PRJNA1149521 SAMN43250874
    ## 6 SRR30305026   AMPLICON        602 50201984 PRJNA1149521 SAMN43250873
    ##                BioSampleModel    Bytes
    ## 1 Metagenome or environmental 29051769
    ## 2 Metagenome or environmental 18779359
    ## 3 Metagenome or environmental 13673684
    ## 4 Metagenome or environmental 27201363
    ## 5 Metagenome or environmental 31724859
    ## 6 Metagenome or environmental 27001123
    ##                                            Center.Name Collection_Date Consent
    ## 1 WEST CHINA SCHOOL OF STOMATOLOGY, SICHUAN UNIVERSITY      2023-07-08  public
    ## 2 WEST CHINA SCHOOL OF STOMATOLOGY, SICHUAN UNIVERSITY      2023-07-08  public
    ## 3 WEST CHINA SCHOOL OF STOMATOLOGY, SICHUAN UNIVERSITY      2023-07-08  public
    ## 4 WEST CHINA SCHOOL OF STOMATOLOGY, SICHUAN UNIVERSITY      2023-07-08  public
    ## 5 WEST CHINA SCHOOL OF STOMATOLOGY, SICHUAN UNIVERSITY      2023-07-08  public
    ## 6 WEST CHINA SCHOOL OF STOMATOLOGY, SICHUAN UNIVERSITY      2023-07-08  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##    Experiment geo_loc_name_country geo_loc_name_country_continent
    ## 1 SRX25765749                China                           Asia
    ## 2 SRX25765748                China                           Asia
    ## 3 SRX25765747                China                           Asia
    ## 4 SRX25765746                China                           Asia
    ## 5 SRX25765745                China                           Asia
    ## 6 SRX25765744                China                           Asia
    ##     geo_loc_name         HOST     Instrument      isolate isolation_source
    ## 1 China: Sichuan Homo sapiens Illumina MiSeq Human siliva           saliva
    ## 2 China: Sichuan Homo sapiens Illumina MiSeq Human siliva           saliva
    ## 3 China: Sichuan Homo sapiens Illumina MiSeq Human siliva           saliva
    ## 4 China: Sichuan Homo sapiens Illumina MiSeq Human siliva           saliva
    ## 5 China: Sichuan Homo sapiens Illumina MiSeq Human siliva           saliva
    ## 6 China: Sichuan Homo sapiens Illumina MiSeq Human siliva           saliva
    ##   lat_lon Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1 missing         CS64        PAIRED              PCR   METAGENOMIC
    ## 2 missing         CS63        PAIRED              PCR   METAGENOMIC
    ## 3 missing        CS121        PAIRED              PCR   METAGENOMIC
    ## 4 missing         CS62        PAIRED              PCR   METAGENOMIC
    ## 5 missing         CS61        PAIRED              PCR   METAGENOMIC
    ## 6 missing         CS60        PAIRED              PCR   METAGENOMIC
    ##                Organism Platform          ReleaseDate          create_date
    ## 1 human oral metagenome ILLUMINA 2024-09-02T00:00:00Z 2024-08-20T02:50:00Z
    ## 2 human oral metagenome ILLUMINA 2024-09-02T00:00:00Z 2024-08-20T02:50:00Z
    ## 3 human oral metagenome ILLUMINA 2024-09-02T00:00:00Z 2024-08-20T02:50:00Z
    ## 4 human oral metagenome ILLUMINA 2024-09-02T00:00:00Z 2024-08-20T02:50:00Z
    ## 5 human oral metagenome ILLUMINA 2024-09-02T00:00:00Z 2024-08-20T02:50:00Z
    ## 6 human oral metagenome ILLUMINA 2024-09-02T00:00:00Z 2024-08-20T02:49:00Z
    ##   version Sample.Name SRA.Study AGE BIOMATERIAL_PROVIDER sex tissue sample_type
    ## 1       1        CS64 SRP527429  NA                                      Saliva
    ## 2       1        CS63 SRP527429  NA                                      Saliva
    ## 3       1       CS121 SRP527429  NA                                      Saliva
    ## 4       1        CS62 SRP527429  NA                                      Saliva
    ## 5       1        CS61 SRP527429  NA                                      Saliva
    ## 6       1        CS60 SRP527429  NA                                      Saliva
    ##           strain
    ## 1 not applicable
    ## 2 not applicable
    ## 3 not applicable
    ## 4 not applicable
    ## 5 not applicable
    ## 6 not applicable

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
# Create a uniform metadata dataframe
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Sample.Name,
  Age = "18-50", 
  BMI = NA,  
  Sex = NA,
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$isolation_source,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)

# Remove entries non-oral samples
metadata <- metadata %>% 
  filter(Site != "missing" | is.na(Site))
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Region <- "V3-V4"
metadata$PCR_primers <- "338F-806R"
metadata$title <- "hTAS2R38 polymorphisms modulate oral microbiota and influence the prevalence and treatment outcome of halitosis"
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status <- "HC"

#extract subject id by remove prefix of sample_name
metadata$sample_name <- gsub("C|Co|CS", "", metadata$sample_name)

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##    1  1_1  1_2  1_3  1_4   10 10_1 10_2 10_3 10_4  100  101  102  103  104  105 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  106  107  108  109   11 11_1 11_2 11_3 11_4  110  111  112  113  114  115  116 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  117  118  119   12 12_1 12_2 12_3 12_4  120  121  122  123  124  125  126  127 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  128  129   13 13_1 13_2 13_3 13_4  130  131  132  133  134  135  136  137  138 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  139   14 14_1 14_2 14_3 14_4  140  141  142  143  144  145  146  147  148  149 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   15 15_1 15_2 15_3 15_4  150  151  152  153  154  155  156  157  158  159   16 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 16_1 16_2 16_3 16_4  160  161  162  163  164  165  166  167  168  169   17 17_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 17_2 17_3 17_4  170  171  172  173  174  175  176  177  178  179   18 18_1 18_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 18_3 18_4  180  181  182  183  184  185  186  187  188  189   19 19_1 19_2 19_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 19_4  190  191  192  193  194  195  196  197  198  199    2  2_1  2_2  2_3  2_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   20 20_1 20_2 20_3 20_4  200  201  202  203  204  205  206  207  208  209   21 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 21_1 21_2 21_3 21_4  210  211  212  213  214  215  216  217  218  219   22 22_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 22_2 22_3 22_4  220  221  222  223  224  225  226  227  228  229   23 23_1 23_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 23_3 23_4  230  231  232  233  234  235  236  237  238  239   24 24_1 24_2 24_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 24_4  240  241  242  243  244  245  246  247  248  249   25 25_1 25_2 25_3 25_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  250  251  252  253  254  255  256  257  258  259   26 26_1 26_2 26_3 26_4  260 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  261  262  263  264  265  266  267  268  269   27 27_1 27_2 27_3 27_4  270  271 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  272  273  274  275  276  277  278  279   28 28_1 28_2 28_3 28_4  280  281  282 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  283  284  285  286  287  288  289   29 29_1 29_2 29_3 29_4  290    3  3_1  3_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  3_3  3_4   30 30_1 30_2 30_3 30_4   31 31_1 31_2 31_3 31_4   32 32_1 32_2 32_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 32_4   33 33_1 33_2 33_3 33_4   34 34_1 34_2 34_3 34_4   35 35_1 35_2 35_3 35_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   36 36_1 36_2 36_3 36_4   37 37_1 37_2 37_3 37_4   38 38_1 38_2 38_3 38_4   39 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 39_1 39_2 39_3 39_4    4  4_1  4_2  4_3  4_4   40 40_1 40_2 40_3 40_4   41 41_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 41_2 41_3 41_4   42 42_1 42_2 42_3 42_4   43 43_1 43_2 43_3 43_4   44 44_1 44_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 44_3 44_4   45 45_1 45_2 45_3 45_4   46 46_1 46_2 46_3 46_4   47 47_1 47_2 47_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 47_4   48 48_1 48_2 48_3 48_4   49 49_1 49_2 49_3 49_4    5  5_1  5_2  5_3  5_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   50 50_1 50_2 50_3 50_4   51 51_1 51_2 51_3 51_4   52 52_1 52_2 52_3 52_4   53 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 53_1 53_2 53_3 53_4   54 54_1 54_2 54_3 54_4   55 55_1 55_2 55_3 55_4   56 56_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 56_2 56_3 56_4   57 57_1 57_2 57_3 57_4   58 58_1 58_2 58_3 58_4   59 59_1 59_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 59_3 59_4    6  6_1  6_2  6_3  6_4   60 60_1 60_2 60_3 60_4   61   62   63   64 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   65   66   67   68   69    7  7_1  7_2  7_3  7_4   70   71   72   73   74   75 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   76   77   78   79    8  8_1  8_2  8_3  8_4   80   81   82   83   84   85   86 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   87   88   89    9  9_1  9_2  9_3  9_4   90   91   92   93   94   95   96   97 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   98   99 
    ##    1    1

``` r
table(metadata$Site)
```

    ## 
    ## saliva 
    ##    530

``` r
# Calculate sampling frequency per subject
sample_counts <- table(metadata$sample_name)
cat("\nNumber of samples per subject:\n")
```

    ## 
    ## Number of samples per subject:

``` r
print(sample_counts)
```

    ## 
    ##    1  1_1  1_2  1_3  1_4   10 10_1 10_2 10_3 10_4  100  101  102  103  104  105 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  106  107  108  109   11 11_1 11_2 11_3 11_4  110  111  112  113  114  115  116 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  117  118  119   12 12_1 12_2 12_3 12_4  120  121  122  123  124  125  126  127 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  128  129   13 13_1 13_2 13_3 13_4  130  131  132  133  134  135  136  137  138 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  139   14 14_1 14_2 14_3 14_4  140  141  142  143  144  145  146  147  148  149 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   15 15_1 15_2 15_3 15_4  150  151  152  153  154  155  156  157  158  159   16 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 16_1 16_2 16_3 16_4  160  161  162  163  164  165  166  167  168  169   17 17_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 17_2 17_3 17_4  170  171  172  173  174  175  176  177  178  179   18 18_1 18_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 18_3 18_4  180  181  182  183  184  185  186  187  188  189   19 19_1 19_2 19_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 19_4  190  191  192  193  194  195  196  197  198  199    2  2_1  2_2  2_3  2_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   20 20_1 20_2 20_3 20_4  200  201  202  203  204  205  206  207  208  209   21 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 21_1 21_2 21_3 21_4  210  211  212  213  214  215  216  217  218  219   22 22_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 22_2 22_3 22_4  220  221  222  223  224  225  226  227  228  229   23 23_1 23_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 23_3 23_4  230  231  232  233  234  235  236  237  238  239   24 24_1 24_2 24_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 24_4  240  241  242  243  244  245  246  247  248  249   25 25_1 25_2 25_3 25_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  250  251  252  253  254  255  256  257  258  259   26 26_1 26_2 26_3 26_4  260 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  261  262  263  264  265  266  267  268  269   27 27_1 27_2 27_3 27_4  270  271 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  272  273  274  275  276  277  278  279   28 28_1 28_2 28_3 28_4  280  281  282 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  283  284  285  286  287  288  289   29 29_1 29_2 29_3 29_4  290    3  3_1  3_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  3_3  3_4   30 30_1 30_2 30_3 30_4   31 31_1 31_2 31_3 31_4   32 32_1 32_2 32_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 32_4   33 33_1 33_2 33_3 33_4   34 34_1 34_2 34_3 34_4   35 35_1 35_2 35_3 35_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   36 36_1 36_2 36_3 36_4   37 37_1 37_2 37_3 37_4   38 38_1 38_2 38_3 38_4   39 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 39_1 39_2 39_3 39_4    4  4_1  4_2  4_3  4_4   40 40_1 40_2 40_3 40_4   41 41_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 41_2 41_3 41_4   42 42_1 42_2 42_3 42_4   43 43_1 43_2 43_3 43_4   44 44_1 44_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 44_3 44_4   45 45_1 45_2 45_3 45_4   46 46_1 46_2 46_3 46_4   47 47_1 47_2 47_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 47_4   48 48_1 48_2 48_3 48_4   49 49_1 49_2 49_3 49_4    5  5_1  5_2  5_3  5_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   50 50_1 50_2 50_3 50_4   51 51_1 51_2 51_3 51_4   52 52_1 52_2 52_3 52_4   53 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 53_1 53_2 53_3 53_4   54 54_1 54_2 54_3 54_4   55 55_1 55_2 55_3 55_4   56 56_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 56_2 56_3 56_4   57 57_1 57_2 57_3 57_4   58 58_1 58_2 58_3 58_4   59 59_1 59_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 59_3 59_4    6  6_1  6_2  6_3  6_4   60 60_1 60_2 60_3 60_4   61   62   63   64 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   65   66   67   68   69    7  7_1  7_2  7_3  7_4   70   71   72   73   74   75 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   76   77   78   79    8  8_1  8_2  8_3  8_4   80   81   82   83   84   85   86 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   87   88   89    9  9_1  9_2  9_3  9_4   90   91   92   93   94   95   96   97 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   98   99 
    ##    1    1

``` r
# Create matched_samples column
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")

cat("\nMatched samples distribution:\n")
```

    ## 
    ## Matched samples distribution:

``` r
table(metadata$matched_samples)
```

    ## 
    ##  NO 
    ## 530

``` r
# Identify samples for removal: single samples excluding saliva and tongue coat
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & 
                     !(metadata$Site %in% c("saliva", "tongue coat"))

cat("\nNumber of samples to be removed:", sum(samples_to_remove), "\n")
```

    ## 
    ## Number of samples to be removed: 0

``` r
# Remove identified samples
metadata_clean <- metadata[!samples_to_remove, ]

cat("\nFinal sample distribution after cleaning:\n")
```

    ## 
    ## Final sample distribution after cleaning:

``` r
table(metadata_clean$sample_name)
```

    ## 
    ##    1  1_1  1_2  1_3  1_4   10 10_1 10_2 10_3 10_4  100  101  102  103  104  105 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  106  107  108  109   11 11_1 11_2 11_3 11_4  110  111  112  113  114  115  116 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  117  118  119   12 12_1 12_2 12_3 12_4  120  121  122  123  124  125  126  127 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  128  129   13 13_1 13_2 13_3 13_4  130  131  132  133  134  135  136  137  138 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  139   14 14_1 14_2 14_3 14_4  140  141  142  143  144  145  146  147  148  149 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   15 15_1 15_2 15_3 15_4  150  151  152  153  154  155  156  157  158  159   16 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 16_1 16_2 16_3 16_4  160  161  162  163  164  165  166  167  168  169   17 17_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 17_2 17_3 17_4  170  171  172  173  174  175  176  177  178  179   18 18_1 18_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 18_3 18_4  180  181  182  183  184  185  186  187  188  189   19 19_1 19_2 19_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 19_4  190  191  192  193  194  195  196  197  198  199    2  2_1  2_2  2_3  2_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   20 20_1 20_2 20_3 20_4  200  201  202  203  204  205  206  207  208  209   21 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 21_1 21_2 21_3 21_4  210  211  212  213  214  215  216  217  218  219   22 22_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 22_2 22_3 22_4  220  221  222  223  224  225  226  227  228  229   23 23_1 23_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 23_3 23_4  230  231  232  233  234  235  236  237  238  239   24 24_1 24_2 24_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 24_4  240  241  242  243  244  245  246  247  248  249   25 25_1 25_2 25_3 25_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  250  251  252  253  254  255  256  257  258  259   26 26_1 26_2 26_3 26_4  260 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  261  262  263  264  265  266  267  268  269   27 27_1 27_2 27_3 27_4  270  271 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  272  273  274  275  276  277  278  279   28 28_1 28_2 28_3 28_4  280  281  282 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  283  284  285  286  287  288  289   29 29_1 29_2 29_3 29_4  290    3  3_1  3_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##  3_3  3_4   30 30_1 30_2 30_3 30_4   31 31_1 31_2 31_3 31_4   32 32_1 32_2 32_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 32_4   33 33_1 33_2 33_3 33_4   34 34_1 34_2 34_3 34_4   35 35_1 35_2 35_3 35_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   36 36_1 36_2 36_3 36_4   37 37_1 37_2 37_3 37_4   38 38_1 38_2 38_3 38_4   39 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 39_1 39_2 39_3 39_4    4  4_1  4_2  4_3  4_4   40 40_1 40_2 40_3 40_4   41 41_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 41_2 41_3 41_4   42 42_1 42_2 42_3 42_4   43 43_1 43_2 43_3 43_4   44 44_1 44_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 44_3 44_4   45 45_1 45_2 45_3 45_4   46 46_1 46_2 46_3 46_4   47 47_1 47_2 47_3 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 47_4   48 48_1 48_2 48_3 48_4   49 49_1 49_2 49_3 49_4    5  5_1  5_2  5_3  5_4 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   50 50_1 50_2 50_3 50_4   51 51_1 51_2 51_3 51_4   52 52_1 52_2 52_3 52_4   53 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 53_1 53_2 53_3 53_4   54 54_1 54_2 54_3 54_4   55 55_1 55_2 55_3 55_4   56 56_1 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 56_2 56_3 56_4   57 57_1 57_2 57_3 57_4   58 58_1 58_2 58_3 58_4   59 59_1 59_2 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## 59_3 59_4    6  6_1  6_2  6_3  6_4   60 60_1 60_2 60_3 60_4   61   62   63   64 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   65   66   67   68   69    7  7_1  7_2  7_3  7_4   70   71   72   73   74   75 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   76   77   78   79    8  8_1  8_2  8_3  8_4   80   81   82   83   84   85   86 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   87   88   89    9  9_1  9_2  9_3  9_4   90   91   92   93   94   95   96   97 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ##   98   99 
    ##    1    1

``` r
table(metadata_clean$Site)
```

    ## 
    ## saliva 
    ##    530

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ##  NO 
    ## 530

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 530"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 530"

# Step 5: Visualize the Metadata

``` r
knitr::include_graphics("Metadata_Summary.png")
```

![](Metadata_Summary.png)<!-- -->

# Step 6: Finalize and Export the Metadata

``` r
# Reorder columns to match the desired format
metadata <- metadata[, c(
  "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
  "Region", "Platform", "Country", "Site", "public_time", 
  "PCR_primers", "title", "Disease_study", "Status", "matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##     BioProject         Run sample_name   Age BMI Sex Region       Platform
    ## 1 PRJNA1149521 SRR30305021          64 18-50  NA  NA  V3-V4 Illumina MiSeq
    ## 2 PRJNA1149521 SRR30305022          63 18-50  NA  NA  V3-V4 Illumina MiSeq
    ## 3 PRJNA1149521 SRR30305023         121 18-50  NA  NA  V3-V4 Illumina MiSeq
    ## 4 PRJNA1149521 SRR30305024          62 18-50  NA  NA  V3-V4 Illumina MiSeq
    ## 5 PRJNA1149521 SRR30305025          61 18-50  NA  NA  V3-V4 Illumina MiSeq
    ## 6 PRJNA1149521 SRR30305026          60 18-50  NA  NA  V3-V4 Illumina MiSeq
    ##   Country   Site          public_time PCR_primers
    ## 1   China saliva 2024-09-02T00:00:00Z   338F-806R
    ## 2   China saliva 2024-09-02T00:00:00Z   338F-806R
    ## 3   China saliva 2024-09-02T00:00:00Z   338F-806R
    ## 4   China saliva 2024-09-02T00:00:00Z   338F-806R
    ## 5   China saliva 2024-09-02T00:00:00Z   338F-806R
    ## 6   China saliva 2024-09-02T00:00:00Z   338F-806R
    ##                                                                                                             title
    ## 1 hTAS2R38 polymorphisms modulate oral microbiota and influence the prevalence and treatment outcome of halitosis
    ## 2 hTAS2R38 polymorphisms modulate oral microbiota and influence the prevalence and treatment outcome of halitosis
    ## 3 hTAS2R38 polymorphisms modulate oral microbiota and influence the prevalence and treatment outcome of halitosis
    ## 4 hTAS2R38 polymorphisms modulate oral microbiota and influence the prevalence and treatment outcome of halitosis
    ## 5 hTAS2R38 polymorphisms modulate oral microbiota and influence the prevalence and treatment outcome of halitosis
    ## 6 hTAS2R38 polymorphisms modulate oral microbiota and influence the prevalence and treatment outcome of halitosis
    ##   Disease_study Status matched_samples
    ## 1            NO     HC              NO
    ## 2            NO     HC              NO
    ## 3            NO     HC              NO
    ## 4            NO     HC              NO
    ## 5            NO     HC              NO
    ## 6            NO     HC              NO

## Notes for Students:

1.  **Metadata Accuracy**: Always cross-check the metadata with the
    original publication to ensure accuracy.
2.  **Manual Entries**: Columns like `PCR_primers`, and `Status` often
    require manual entry based on the publication.
3.  **Disease Status**: Be careful when assigning `Status`—some studies
    include both disease and healthy samples.
4.  **Visualization**: The plots help quickly identify the distribution
    of values in each column, which is useful for quality control.
5.  **Reproducibility**: Update the working directory and file paths to
    match your local environment.
