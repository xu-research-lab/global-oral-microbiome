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
setwd("E:/bioproject/PRJNA910326")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = T) 

# Remove entries non-oral samples
initial_meta <- initial_meta[initial_meta$Sample.Name != "DF", ]

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR22711069   AMPLICON        495 39050550 PRJNA910326 SAMN32122688
    ## 2 SRR22711070   AMPLICON        492 37252272 PRJNA910326 SAMN32122687
    ## 3 SRR22711071   AMPLICON        493 42630203 PRJNA910326 SAMN32122686
    ## 4 SRR22711072   AMPLICON        494 35275058 PRJNA910326 SAMN32122685
    ## 5 SRR22711073   AMPLICON        495 40858290 PRJNA910326 SAMN32122684
    ## 6 SRR22711074   AMPLICON        492 38685960 PRJNA910326 SAMN32122683
    ##                BioSampleModel    Bytes              Center.Name Collection_Date
    ## 1 Metagenome or environmental 12313205 SCHOOL OF BASIC MEDICINE      2021-08-01
    ## 2 Metagenome or environmental 11983424 SCHOOL OF BASIC MEDICINE      2021-08-01
    ## 3 Metagenome or environmental 13591788 SCHOOL OF BASIC MEDICINE      2021-08-01
    ## 4 Metagenome or environmental 11356048 SCHOOL OF BASIC MEDICINE      2021-08-01
    ## 5 Metagenome or environmental 13184345 SCHOOL OF BASIC MEDICINE      2021-08-01
    ## 6 Metagenome or environmental 12555840 SCHOOL OF BASIC MEDICINE      2021-08-01
    ##   Consent DATASTORE.filetype DATASTORE.provider
    ## 1  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region  Experiment geo_loc_name_country
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX18672881                China
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX18672880                China
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX18672879                China
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX18672878                China
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX18672877                China
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX18672876                China
    ##   geo_loc_name_country_continent   geo_loc_name         HOST
    ## 1                           Asia China:Shanghai Homo sapiens
    ## 2                           Asia China:Shanghai Homo sapiens
    ## 3                           Asia China:Shanghai Homo sapiens
    ## 4                           Asia China:Shanghai Homo sapiens
    ## 5                           Asia China:Shanghai Homo sapiens
    ## 6                           Asia China:Shanghai Homo sapiens
    ##              Instrument isolation_source          lat_lon Library.Name
    ## 1 Illumina NovaSeq 6000   Tongue coating 31.11 N 121.35 E        DT094
    ## 2 Illumina NovaSeq 6000   Tongue coating 31.11 N 121.35 E        DT093
    ## 3 Illumina NovaSeq 6000   Tongue coating 31.11 N 121.35 E        DT092
    ## 4 Illumina NovaSeq 6000   Tongue coating 31.11 N 121.35 E        DT091
    ## 5 Illumina NovaSeq 6000   Tongue coating 31.11 N 121.35 E        DT090
    ## 6 Illumina NovaSeq 6000   Tongue coating 31.11 N 121.35 E        DT089
    ##   LibraryLayout LibrarySelection LibrarySource              Organism Platform
    ## 1        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 2        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 3        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 4        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 5        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 6        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ##            ReleaseDate          create_date version Sample.Name SRA.Study
    ## 1 2023-02-02T00:00:00Z 2022-12-13T14:41:00Z       1       DT094 SRP412812
    ## 2 2023-02-02T00:00:00Z 2022-12-13T14:41:00Z       1       DT093 SRP412812
    ## 3 2023-02-02T00:00:00Z 2022-12-13T14:41:00Z       1       DT092 SRP412812
    ## 4 2023-02-02T00:00:00Z 2022-12-13T14:41:00Z       1       DT091 SRP412812
    ## 5 2023-02-02T00:00:00Z 2022-12-13T14:41:00Z       1       DT090 SRP412812
    ## 6 2023-02-02T00:00:00Z 2022-12-13T14:42:00Z       1       DT089 SRP412812

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
  Age = "18-75", 
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
  filter(Site != "gut" | is.na(Site))
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Region <- "V3-V4"
metadata$PCR_primers <- "338F-806R"
metadata$title <- "Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata <- metadata %>%
  mutate(
    Status = case_when(
      grepl("^DT", sample_name) ~ "T2DM",
      grepl("^HT", sample_name) ~ "HC",
      TRUE ~ NA)
    )

#extract subject id by remove prefix of sample_name
metadata$sample_name <- gsub("DT|HT", "", metadata$sample_name)

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##   001   002   003   004   005   006   007   008   009   010   011   012   013 
    ##     2     1     2     2     2     2     2     2     1     2     2     2     2 
    ##   014   015   016   017   018   019   020   021   022   023   024   025   026 
    ##     2     2     2     2     2     2     2     2     1     2     2     1     2 
    ##   027   028   029   030   031   032   033   034   035   036   037   038   039 
    ##     2     2     2     2     2     2     1     2     1     1     1     1     1 
    ##   040   041   042   043   045   046   047   048   050   051   052   053   054 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   055   057   058   061   062   063   064   065   066   067   068   069   070 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   071   073   074   075   076   077   079   082   083   085   086   087   088 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   089   090   091   092   093   094   095   096   097   098   099   100   101 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   103   104   105   106   107   109   110   111   112   113   114   115   116 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   118   119   120   122   123   124   125   126   127   128   129   131   132 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   133   134   135   136   137   138   139   140   143   144   145   147   148 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   150   152   154   155   156   157   159   160   161   162   163   165   166 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   167   168   169   170   171   172   173   174   175   176   177   178   179 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   181   182   183   185   186   187   188   191   192   193   194   195   197 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   198   199   200   201   202   203   205   206   207   208   209   210   212 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   213   215   216   217   218   219   225 DF001 DF003 DF004 DF005 DF006 DF007 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF008 DF010 DF011 DF012 DF013 DF014 DF015 DF016 DF017 DF018 DF019 DF020 DF021 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF023 DF024 DF026 DF027 DF029 DF030 DF031 DF032 DF035 DF036 DF037 DF038 DF039 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF041 DF042 DF043 DF045 DF046 DF048 DF050 DF051 DF052 DF053 DF058 DF061 DF062 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF063 DF064 DF066 DF067 DF068 DF069 DF070 DF071 DF073 DF074 DF075 DF079 DF083 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF085 DF086 DF090 DF091 DF092 DF094 DF095 DF096 DF098 DF100 DF106 DF107 DF109 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF111 DF112 DF113 DF119 DF122 DF129 DF131 DF132 DF133 DF134 DF135 DF136 DF137 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF138 DF140 DF145 DF147 DF148 DF150 DF152 DF153 DF156 DF157 DF160 DF161 DF163 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF165 DF166 DF168 DF172 DF173 DF174 DF175 DF176 DF179 DF181 DF182 DF183 DF186 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF187 DF191 DF193 DF194 DF197 DF199 DF200 DF201 DF202 DF203 DF206 DF207 DF208 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF210 DF212 DF213 DF217 DF218 DF219 
    ##     1     1     1     1     1     1

``` r
table(metadata$Site)
```

    ## 
    ##            Gut Tongue coating 
    ##            129            218

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
    ##   001   002   003   004   005   006   007   008   009   010   011   012   013 
    ##     2     1     2     2     2     2     2     2     1     2     2     2     2 
    ##   014   015   016   017   018   019   020   021   022   023   024   025   026 
    ##     2     2     2     2     2     2     2     2     1     2     2     1     2 
    ##   027   028   029   030   031   032   033   034   035   036   037   038   039 
    ##     2     2     2     2     2     2     1     2     1     1     1     1     1 
    ##   040   041   042   043   045   046   047   048   050   051   052   053   054 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   055   057   058   061   062   063   064   065   066   067   068   069   070 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   071   073   074   075   076   077   079   082   083   085   086   087   088 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   089   090   091   092   093   094   095   096   097   098   099   100   101 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   103   104   105   106   107   109   110   111   112   113   114   115   116 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   118   119   120   122   123   124   125   126   127   128   129   131   132 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   133   134   135   136   137   138   139   140   143   144   145   147   148 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   150   152   154   155   156   157   159   160   161   162   163   165   166 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   167   168   169   170   171   172   173   174   175   176   177   178   179 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   181   182   183   185   186   187   188   191   192   193   194   195   197 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   198   199   200   201   202   203   205   206   207   208   209   210   212 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   213   215   216   217   218   219   225 DF001 DF003 DF004 DF005 DF006 DF007 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF008 DF010 DF011 DF012 DF013 DF014 DF015 DF016 DF017 DF018 DF019 DF020 DF021 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF023 DF024 DF026 DF027 DF029 DF030 DF031 DF032 DF035 DF036 DF037 DF038 DF039 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF041 DF042 DF043 DF045 DF046 DF048 DF050 DF051 DF052 DF053 DF058 DF061 DF062 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF063 DF064 DF066 DF067 DF068 DF069 DF070 DF071 DF073 DF074 DF075 DF079 DF083 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF085 DF086 DF090 DF091 DF092 DF094 DF095 DF096 DF098 DF100 DF106 DF107 DF109 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF111 DF112 DF113 DF119 DF122 DF129 DF131 DF132 DF133 DF134 DF135 DF136 DF137 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF138 DF140 DF145 DF147 DF148 DF150 DF152 DF153 DF156 DF157 DF160 DF161 DF163 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF165 DF166 DF168 DF172 DF173 DF174 DF175 DF176 DF179 DF181 DF182 DF183 DF186 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF187 DF191 DF193 DF194 DF197 DF199 DF200 DF201 DF202 DF203 DF206 DF207 DF208 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## DF210 DF212 DF213 DF217 DF218 DF219 
    ##     1     1     1     1     1     1

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
    ##  NO YES 
    ## 289  58

``` r
# Identify samples for removal: single samples excluding saliva and tongue coat
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & 
                     !(metadata$Site %in% c("saliva", "tongue coat", "Tongue coating"))

cat("\nNumber of samples to be removed:", sum(samples_to_remove), "\n")
```

    ## 
    ## Number of samples to be removed: 129

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
    ## 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 
    ##   2   1   2   2   2   2   2   2   1   2   2   2   2   2   2   2   2   2   2   2 
    ## 021 022 023 024 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 
    ##   2   1   2   2   1   2   2   2   2   2   2   2   1   2   1   1   1   1   1   1 
    ## 041 042 043 045 046 047 048 050 051 052 053 054 055 057 058 061 062 063 064 065 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 066 067 068 069 070 071 073 074 075 076 077 079 082 083 085 086 087 088 089 090 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 091 092 093 094 095 096 097 098 099 100 101 103 104 105 106 107 109 110 111 112 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 113 114 115 116 118 119 120 122 123 124 125 126 127 128 129 131 132 133 134 135 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 136 137 138 139 140 143 144 145 147 148 150 152 154 155 156 157 159 160 161 162 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 163 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 181 182 183 185 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 186 187 188 191 192 193 194 195 197 198 199 200 201 202 203 205 206 207 208 209 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 210 212 213 215 216 217 218 219 225 
    ##   1   1   1   1   1   1   1   1   1

``` r
table(metadata_clean$Site)
```

    ## 
    ## Tongue coating 
    ##            218

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ##  NO YES 
    ## 160  58

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 218"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 189"

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
  "PCR_primers", "title", "Disease_study", "Status","matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##    BioProject         Run sample_name   Age BMI Sex Region
    ## 1 PRJNA910326 SRR22711069         094 18-75  NA  NA  V3-V4
    ## 2 PRJNA910326 SRR22711070         093 18-75  NA  NA  V3-V4
    ## 3 PRJNA910326 SRR22711071         092 18-75  NA  NA  V3-V4
    ## 4 PRJNA910326 SRR22711072         091 18-75  NA  NA  V3-V4
    ## 5 PRJNA910326 SRR22711073         090 18-75  NA  NA  V3-V4
    ## 6 PRJNA910326 SRR22711074         089 18-75  NA  NA  V3-V4
    ##                Platform Country           Site          public_time PCR_primers
    ## 1 Illumina NovaSeq 6000   China Tongue coating 2023-02-02T00:00:00Z   338F-806R
    ## 2 Illumina NovaSeq 6000   China Tongue coating 2023-02-02T00:00:00Z   338F-806R
    ## 3 Illumina NovaSeq 6000   China Tongue coating 2023-02-02T00:00:00Z   338F-806R
    ## 4 Illumina NovaSeq 6000   China Tongue coating 2023-02-02T00:00:00Z   338F-806R
    ## 5 Illumina NovaSeq 6000   China Tongue coating 2023-02-02T00:00:00Z   338F-806R
    ## 6 Illumina NovaSeq 6000   China Tongue coating 2023-02-02T00:00:00Z   338F-806R
    ##                                                                                                                                       title
    ## 1 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ## 2 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ## 3 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ## 4 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ## 5 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ## 6 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ##   Disease_study Status matched_samples
    ## 1           YES   T2DM              NO
    ## 2           YES   T2DM              NO
    ## 3           YES   T2DM              NO
    ## 4           YES   T2DM              NO
    ## 5           YES   T2DM              NO
    ## 6           YES   T2DM              NO

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
