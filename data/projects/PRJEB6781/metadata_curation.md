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
setwd("E:/bioproject_YL/PRJEB6781")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = T)

# Display the first few rows of the metadata
head(initial_meta)
```

    ##         Run AGE Assay.Type AvgSpotLen   Bases BioProject    BioSample  Bytes
    ## 1 ERR618776  35   AMPLICON        299  705597  PRJEB6781 SAMEA2672934 480172
    ## 2 ERR618777  35   AMPLICON        327 1146038  PRJEB6781 SAMEA2672935 741667
    ## 3 ERR618778  35   AMPLICON        302  695622  PRJEB6781 SAMEA2672936 475764
    ## 4 ERR618779  35   AMPLICON        307  933490  PRJEB6781 SAMEA2672937 629230
    ## 5 ERR618780  35   AMPLICON        316  913272  PRJEB6781 SAMEA2672938 617439
    ## 6 ERR618781  39   AMPLICON        321 1445040  PRJEB6781 SAMEA2672939 951180
    ##   Center.Name Collection_Date common_name Consent Couple_Id DATASTORE.filetype
    ## 1         TNO      2012-07-01       human  public  Couple 1   fastq,run.zq,sra
    ## 2         TNO      2012-07-01       human  public  Couple 1   fastq,run.zq,sra
    ## 3         TNO      2012-07-01       human  public  Couple 2   fastq,run.zq,sra
    ## 4         TNO      2012-07-01       human  public  Couple 3   fastq,run.zq,sra
    ## 5         TNO      2012-07-01       human  public  Couple 5   fastq,run.zq,sra
    ## 6         TNO      2012-07-01       human  public  Couple 6   fastq,run.zq,sra
    ##   DATASTORE.provider                         DATASTORE.region Donor.Acceptor
    ## 1     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1              D
    ## 2     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1              D
    ## 3     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1              D
    ## 4     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1              D
    ## 5     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1              D
    ## 6     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1              A
    ##   ENA.FIRST.PUBLIC..run. ENA.LAST.UPDATE..run. environment_.biome.
    ## 1             2014-09-17            2018-11-16    Human associated
    ## 2             2014-09-17            2018-11-16    Human associated
    ## 3             2014-09-17            2018-11-16    Human associated
    ## 4             2014-09-17            2018-11-16    Human associated
    ## 5             2014-09-17            2018-11-16    Human associated
    ## 6             2014-09-17            2018-11-16    Human associated
    ##   environment_.feature. environment_.material. Experiment  external_id
    ## 1                 Mouth                 Saliva  ERX575547 SAMEA2672934
    ## 2                 Mouth                 Saliva  ERX575548 SAMEA2672935
    ## 3                 Mouth                 Saliva  ERX575549 SAMEA2672936
    ## 4                 Mouth                 Saliva  ERX575550 SAMEA2672937
    ## 5                 Mouth                 Saliva  ERX575551 SAMEA2672938
    ## 6                 Mouth                 Saliva  ERX575552 SAMEA2672939
    ##   geo_loc_name_country geo_loc_name_country_continent
    ## 1          Netherlands                         Europe
    ## 2          Netherlands                         Europe
    ## 3          Netherlands                         Europe
    ## 4          Netherlands                         Europe
    ## 5          Netherlands                         Europe
    ## 6          Netherlands                         Europe
    ##   geographic_location_.country_and.or_sea. geographic_location_.latitude.
    ## 1                              Netherlands                       52.36609
    ## 2                              Netherlands                       52.36609
    ## 3                              Netherlands                       52.36609
    ## 4                              Netherlands                       52.36609
    ## 5                              Netherlands                       52.36609
    ## 6                              Netherlands                       52.36609
    ##   geographic_location_.longitude. host.associated_environmental_package
    ## 1                        4.916525                       host-associated
    ## 2                        4.916525                       host-associated
    ## 3                        4.916525                       host-associated
    ## 4                        4.916525                       host-associated
    ## 5                        4.916525                       host-associated
    ## 6                        4.916525                       host-associated
    ##   host_body_product host_sex INSDC_center_name   INSDC_first_public
    ## 1            Saliva   female               TNO 2014-07-30T11:25:28Z
    ## 2            Saliva   female               TNO 2014-09-17T17:04:10Z
    ## 3            Saliva   female               TNO 2014-09-17T17:04:10Z
    ## 4            Saliva   female               TNO 2014-09-17T17:04:10Z
    ## 5            Saliva   female               TNO 2014-09-17T17:04:10Z
    ## 6            Saliva   female               TNO 2014-09-17T17:04:10Z
    ##      INSDC_last_update INSDC_status          Instrument investigation_type
    ## 1 2016-10-21T09:39:01Z       public 454 GS FLX Titanium     mimarks-survey
    ## 2 2016-10-21T09:39:01Z       public 454 GS FLX Titanium     mimarks-survey
    ## 3 2016-10-21T09:39:01Z       public 454 GS FLX Titanium     mimarks-survey
    ## 4 2016-10-21T09:39:01Z       public 454 GS FLX Titanium     mimarks-survey
    ## 5 2016-10-21T09:39:01Z       public 454 GS FLX Titanium     mimarks-survey
    ## 6 2016-10-21T09:39:01Z       public 454 GS FLX Titanium     mimarks-survey
    ##   Library.Name LibraryLayout LibrarySelection LibrarySource     Organism
    ## 1  unspecified        SINGLE              PCR         OTHER Homo sapiens
    ## 2  unspecified        SINGLE              PCR         OTHER Homo sapiens
    ## 3  unspecified        SINGLE              PCR         OTHER Homo sapiens
    ## 4  unspecified        SINGLE              PCR         OTHER Homo sapiens
    ## 5  unspecified        SINGLE              PCR         OTHER Homo sapiens
    ## 6  unspecified        SINGLE              PCR         OTHER Homo sapiens
    ##   Platform                                    project_name          ReleaseDate
    ## 1    LS454 Shaping the oral microbiota by intimate kissing 2014-09-18T00:00:00Z
    ## 2    LS454 Shaping the oral microbiota by intimate kissing 2014-09-18T00:00:00Z
    ## 3    LS454 Shaping the oral microbiota by intimate kissing 2014-09-18T00:00:00Z
    ## 4    LS454 Shaping the oral microbiota by intimate kissing 2014-09-18T00:00:00Z
    ## 5    LS454 Shaping the oral microbiota by intimate kissing 2014-09-18T00:00:00Z
    ## 6    LS454 Shaping the oral microbiota by intimate kissing 2014-09-18T00:00:00Z
    ##   Reported_Kissing_with_partner_Frequency          create_date version
    ## 1                             2 or 3/week 2014-09-25T19:05:00Z       1
    ## 2                             2 or 3/week 2014-09-25T19:10:00Z       1
    ## 3                             2 or 3/week 2014-09-25T19:10:00Z       1
    ## 4                             2 or 3/week 2014-09-25T19:05:00Z       1
    ## 5                             2 or 3/week 2014-09-25T19:10:00Z       1
    ## 6                                  7/week 2014-09-25T19:05:00Z       1
    ##    Sample.Name      Sample_name sampling_time_point Sequencing_method SRA.Study
    ## 1 SAMEA2672934 Oral-Microbiota1                   1    MiSeq Illumina ERP006418
    ## 2 SAMEA2672935 Oral-Microbiota2                   1    MiSeq Illumina ERP006418
    ## 3 SAMEA2672936 Oral-Microbiota3                   2    MiSeq Illumina ERP006418
    ## 4 SAMEA2672937 Oral-Microbiota4                   3    MiSeq Illumina ERP006418
    ## 5 SAMEA2672938 Oral-Microbiota5                   3    MiSeq Illumina ERP006418
    ## 6 SAMEA2672939 Oral-Microbiota6                   1    MiSeq Illumina ERP006418
    ##       Submitter_Id Time_since_last_intimate_kissing_with_partner
    ## 1 Oral-Microbiota1                                   2 to 3 days
    ## 2 Oral-Microbiota2                                   3 to 3 days
    ## 3 Oral-Microbiota3                                   4 to 3 days
    ## 4 Oral-Microbiota4                                   5 to 3 days
    ## 5 Oral-Microbiota5                                   6 to 3 days
    ## 6 Oral-Microbiota6                                         2days
    ##   ENA_first_public ENA_last_update INSDC_center_alias ENA.FIRST.PUBLIC
    ## 1       2014-07-30      2016-10-21                TNO                 
    ## 2       2014-09-17      2016-10-21                TNO                 
    ## 3       2014-09-17      2016-10-21                TNO                 
    ## 4       2014-09-17      2016-10-21                TNO                 
    ## 5       2014-09-17      2016-10-21                TNO                 
    ## 6       2014-09-17      2016-10-21                TNO                 
    ##   ENA.LAST.UPDATE Scientific_Name
    ## 1                                
    ## 2                                
    ## 3                                
    ## 4                                
    ## 5                                
    ## 6

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
  sample_name = initial_meta$Submitter_Id,
  Age = initial_meta$AGE,
  BMI = NA,  
  Sex = initial_meta$host_sex,
  Region = "V5-V7",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$environment_.material.,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)

# Use abbreviations to indicate gender
metadata <- metadata %>%
  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "NA"
metadata$title <- "Shaping the oral microbiota through intimate kissing"
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status <- "HC"

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##   Oral-Microbiota1  Oral-Microbiota10 Oral-Microbiota100 Oral-Microbiota101 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota102 Oral-Microbiota103 Oral-Microbiota104 Oral-Microbiota105 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota106 Oral-Microbiota107 Oral-Microbiota108 Oral-Microbiota109 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota11 Oral-Microbiota110 Oral-Microbiota111 Oral-Microbiota112 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota113 Oral-Microbiota114 Oral-Microbiota115 Oral-Microbiota116 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota117 Oral-Microbiota118 Oral-Microbiota119  Oral-Microbiota12 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota120 Oral-Microbiota121 Oral-Microbiota122 Oral-Microbiota123 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota124 Oral-Microbiota125 Oral-Microbiota126 Oral-Microbiota127 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota128 Oral-Microbiota129  Oral-Microbiota13 Oral-Microbiota130 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota131 Oral-Microbiota132 Oral-Microbiota133 Oral-Microbiota134 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota135 Oral-Microbiota136 Oral-Microbiota137 Oral-Microbiota138 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota139  Oral-Microbiota14 Oral-Microbiota140 Oral-Microbiota141 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota142 Oral-Microbiota143 Oral-Microbiota144 Oral-Microbiota145 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota146 Oral-Microbiota147 Oral-Microbiota148 Oral-Microbiota149 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota15 Oral-Microbiota150 Oral-Microbiota151 Oral-Microbiota152 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota153 Oral-Microbiota154 Oral-Microbiota155 Oral-Microbiota156 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota157 Oral-Microbiota158 Oral-Microbiota159  Oral-Microbiota16 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota160 Oral-Microbiota161 Oral-Microbiota162 Oral-Microbiota163 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota164 Oral-Microbiota165 Oral-Microbiota166 Oral-Microbiota167 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota168 Oral-Microbiota169  Oral-Microbiota17 Oral-Microbiota170 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota171 Oral-Microbiota172 Oral-Microbiota173 Oral-Microbiota174 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota175 Oral-Microbiota176 Oral-Microbiota177 Oral-Microbiota178 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota179  Oral-Microbiota18 Oral-Microbiota180 Oral-Microbiota181 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota182 Oral-Microbiota183 Oral-Microbiota184 Oral-Microbiota185 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota186 Oral-Microbiota187 Oral-Microbiota188 Oral-Microbiota189 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota19 Oral-Microbiota190 Oral-Microbiota191 Oral-Microbiota192 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota193 Oral-Microbiota194 Oral-Microbiota195 Oral-Microbiota196 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota197 Oral-Microbiota198 Oral-Microbiota199   Oral-Microbiota2 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota20 Oral-Microbiota200 Oral-Microbiota201 Oral-Microbiota202 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota203 Oral-Microbiota204 Oral-Microbiota205 Oral-Microbiota206 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota207 Oral-Microbiota208 Oral-Microbiota209  Oral-Microbiota21 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota210 Oral-Microbiota211 Oral-Microbiota212 Oral-Microbiota213 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota214 Oral-Microbiota215 Oral-Microbiota216 Oral-Microbiota217 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota218 Oral-Microbiota219  Oral-Microbiota22 Oral-Microbiota220 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota221 Oral-Microbiota222 Oral-Microbiota223 Oral-Microbiota224 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota225 Oral-Microbiota226 Oral-Microbiota227 Oral-Microbiota228 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota229  Oral-Microbiota23 Oral-Microbiota230 Oral-Microbiota231 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota232 Oral-Microbiota233 Oral-Microbiota234 Oral-Microbiota235 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota236 Oral-Microbiota237 Oral-Microbiota238 Oral-Microbiota239 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota24 Oral-Microbiota240 Oral-Microbiota241 Oral-Microbiota242 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota243 Oral-Microbiota244 Oral-Microbiota245 Oral-Microbiota246 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota247 Oral-Microbiota248 Oral-Microbiota249  Oral-Microbiota25 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota250 Oral-Microbiota251 Oral-Microbiota252 Oral-Microbiota253 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota254 Oral-Microbiota255 Oral-Microbiota256 Oral-Microbiota257 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota258 Oral-Microbiota259  Oral-Microbiota26 Oral-Microbiota260 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota261 Oral-Microbiota262 Oral-Microbiota263 Oral-Microbiota264 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota265 Oral-Microbiota266 Oral-Microbiota267 Oral-Microbiota268 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota269  Oral-Microbiota27 Oral-Microbiota270 Oral-Microbiota271 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota272 Oral-Microbiota273 Oral-Microbiota274 Oral-Microbiota275 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota276 Oral-Microbiota277 Oral-Microbiota278 Oral-Microbiota279 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota28 Oral-Microbiota280 Oral-Microbiota281 Oral-Microbiota282 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota283 Oral-Microbiota284  Oral-Microbiota29   Oral-Microbiota3 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota30  Oral-Microbiota31  Oral-Microbiota32  Oral-Microbiota33 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota34  Oral-Microbiota35  Oral-Microbiota36  Oral-Microbiota37 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota38  Oral-Microbiota39   Oral-Microbiota4  Oral-Microbiota40 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota41  Oral-Microbiota42  Oral-Microbiota43  Oral-Microbiota44 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota45  Oral-Microbiota46  Oral-Microbiota47  Oral-Microbiota48 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota49   Oral-Microbiota5  Oral-Microbiota50  Oral-Microbiota51 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota52  Oral-Microbiota53  Oral-Microbiota54  Oral-Microbiota55 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota56  Oral-Microbiota57  Oral-Microbiota58  Oral-Microbiota59 
    ##                  1                  1                  1                  1 
    ##   Oral-Microbiota6  Oral-Microbiota60  Oral-Microbiota61  Oral-Microbiota62 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota63  Oral-Microbiota64  Oral-Microbiota65  Oral-Microbiota66 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota67  Oral-Microbiota68  Oral-Microbiota69   Oral-Microbiota7 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota70  Oral-Microbiota71  Oral-Microbiota72  Oral-Microbiota73 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota74  Oral-Microbiota75  Oral-Microbiota76  Oral-Microbiota77 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota78  Oral-Microbiota79   Oral-Microbiota8  Oral-Microbiota80 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota81  Oral-Microbiota82  Oral-Microbiota83  Oral-Microbiota84 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota85  Oral-Microbiota86  Oral-Microbiota87  Oral-Microbiota88 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota89   Oral-Microbiota9  Oral-Microbiota90  Oral-Microbiota91 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota92  Oral-Microbiota93  Oral-Microbiota94  Oral-Microbiota95 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota96  Oral-Microbiota97  Oral-Microbiota98  Oral-Microbiota99 
    ##                  1                  1                  1                  1

``` r
table(metadata$Site)
```

    ## 
    ##          Saliva Tongue scraping 
    ##             141             143

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
    ##   Oral-Microbiota1  Oral-Microbiota10 Oral-Microbiota100 Oral-Microbiota101 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota102 Oral-Microbiota103 Oral-Microbiota104 Oral-Microbiota105 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota106 Oral-Microbiota107 Oral-Microbiota108 Oral-Microbiota109 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota11 Oral-Microbiota110 Oral-Microbiota111 Oral-Microbiota112 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota113 Oral-Microbiota114 Oral-Microbiota115 Oral-Microbiota116 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota117 Oral-Microbiota118 Oral-Microbiota119  Oral-Microbiota12 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota120 Oral-Microbiota121 Oral-Microbiota122 Oral-Microbiota123 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota124 Oral-Microbiota125 Oral-Microbiota126 Oral-Microbiota127 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota128 Oral-Microbiota129  Oral-Microbiota13 Oral-Microbiota130 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota131 Oral-Microbiota132 Oral-Microbiota133 Oral-Microbiota134 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota135 Oral-Microbiota136 Oral-Microbiota137 Oral-Microbiota138 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota139  Oral-Microbiota14 Oral-Microbiota140 Oral-Microbiota141 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota142 Oral-Microbiota143 Oral-Microbiota144 Oral-Microbiota145 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota146 Oral-Microbiota147 Oral-Microbiota148 Oral-Microbiota149 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota15 Oral-Microbiota150 Oral-Microbiota151 Oral-Microbiota152 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota153 Oral-Microbiota154 Oral-Microbiota155 Oral-Microbiota156 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota157 Oral-Microbiota158 Oral-Microbiota159  Oral-Microbiota16 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota160 Oral-Microbiota161 Oral-Microbiota162 Oral-Microbiota163 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota164 Oral-Microbiota165 Oral-Microbiota166 Oral-Microbiota167 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota168 Oral-Microbiota169  Oral-Microbiota17 Oral-Microbiota170 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota171 Oral-Microbiota172 Oral-Microbiota173 Oral-Microbiota174 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota175 Oral-Microbiota176 Oral-Microbiota177 Oral-Microbiota178 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota179  Oral-Microbiota18 Oral-Microbiota180 Oral-Microbiota181 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota182 Oral-Microbiota183 Oral-Microbiota184 Oral-Microbiota185 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota186 Oral-Microbiota187 Oral-Microbiota188 Oral-Microbiota189 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota19 Oral-Microbiota190 Oral-Microbiota191 Oral-Microbiota192 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota193 Oral-Microbiota194 Oral-Microbiota195 Oral-Microbiota196 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota197 Oral-Microbiota198 Oral-Microbiota199   Oral-Microbiota2 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota20 Oral-Microbiota200 Oral-Microbiota201 Oral-Microbiota202 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota203 Oral-Microbiota204 Oral-Microbiota205 Oral-Microbiota206 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota207 Oral-Microbiota208 Oral-Microbiota209  Oral-Microbiota21 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota210 Oral-Microbiota211 Oral-Microbiota212 Oral-Microbiota213 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota214 Oral-Microbiota215 Oral-Microbiota216 Oral-Microbiota217 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota218 Oral-Microbiota219  Oral-Microbiota22 Oral-Microbiota220 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota221 Oral-Microbiota222 Oral-Microbiota223 Oral-Microbiota224 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota225 Oral-Microbiota226 Oral-Microbiota227 Oral-Microbiota228 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota229  Oral-Microbiota23 Oral-Microbiota230 Oral-Microbiota231 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota232 Oral-Microbiota233 Oral-Microbiota234 Oral-Microbiota235 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota236 Oral-Microbiota237 Oral-Microbiota238 Oral-Microbiota239 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota24 Oral-Microbiota240 Oral-Microbiota241 Oral-Microbiota242 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota243 Oral-Microbiota244 Oral-Microbiota245 Oral-Microbiota246 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota247 Oral-Microbiota248 Oral-Microbiota249  Oral-Microbiota25 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota250 Oral-Microbiota251 Oral-Microbiota252 Oral-Microbiota253 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota254 Oral-Microbiota255 Oral-Microbiota256 Oral-Microbiota257 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota258 Oral-Microbiota259  Oral-Microbiota26 Oral-Microbiota260 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota261 Oral-Microbiota262 Oral-Microbiota263 Oral-Microbiota264 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota265 Oral-Microbiota266 Oral-Microbiota267 Oral-Microbiota268 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota269  Oral-Microbiota27 Oral-Microbiota270 Oral-Microbiota271 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota272 Oral-Microbiota273 Oral-Microbiota274 Oral-Microbiota275 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota276 Oral-Microbiota277 Oral-Microbiota278 Oral-Microbiota279 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota28 Oral-Microbiota280 Oral-Microbiota281 Oral-Microbiota282 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota283 Oral-Microbiota284  Oral-Microbiota29   Oral-Microbiota3 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota30  Oral-Microbiota31  Oral-Microbiota32  Oral-Microbiota33 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota34  Oral-Microbiota35  Oral-Microbiota36  Oral-Microbiota37 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota38  Oral-Microbiota39   Oral-Microbiota4  Oral-Microbiota40 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota41  Oral-Microbiota42  Oral-Microbiota43  Oral-Microbiota44 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota45  Oral-Microbiota46  Oral-Microbiota47  Oral-Microbiota48 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota49   Oral-Microbiota5  Oral-Microbiota50  Oral-Microbiota51 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota52  Oral-Microbiota53  Oral-Microbiota54  Oral-Microbiota55 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota56  Oral-Microbiota57  Oral-Microbiota58  Oral-Microbiota59 
    ##                  1                  1                  1                  1 
    ##   Oral-Microbiota6  Oral-Microbiota60  Oral-Microbiota61  Oral-Microbiota62 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota63  Oral-Microbiota64  Oral-Microbiota65  Oral-Microbiota66 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota67  Oral-Microbiota68  Oral-Microbiota69   Oral-Microbiota7 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota70  Oral-Microbiota71  Oral-Microbiota72  Oral-Microbiota73 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota74  Oral-Microbiota75  Oral-Microbiota76  Oral-Microbiota77 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota78  Oral-Microbiota79   Oral-Microbiota8  Oral-Microbiota80 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota81  Oral-Microbiota82  Oral-Microbiota83  Oral-Microbiota84 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota85  Oral-Microbiota86  Oral-Microbiota87  Oral-Microbiota88 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota89   Oral-Microbiota9  Oral-Microbiota90  Oral-Microbiota91 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota92  Oral-Microbiota93  Oral-Microbiota94  Oral-Microbiota95 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota96  Oral-Microbiota97  Oral-Microbiota98  Oral-Microbiota99 
    ##                  1                  1                  1                  1

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
    ## 284

``` r
# Identify samples for removal: single samples excluding Saliva
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & 
                     !(metadata$Site %in% c("Saliva"))

cat("\nNumber of samples to be removed:", sum(samples_to_remove), "\n")
```

    ## 
    ## Number of samples to be removed: 143

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
    ##   Oral-Microbiota1  Oral-Microbiota10 Oral-Microbiota100 Oral-Microbiota101 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota102 Oral-Microbiota103 Oral-Microbiota104 Oral-Microbiota105 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota106 Oral-Microbiota107 Oral-Microbiota108 Oral-Microbiota109 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota11 Oral-Microbiota110 Oral-Microbiota111 Oral-Microbiota112 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota113 Oral-Microbiota114 Oral-Microbiota115 Oral-Microbiota116 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota117 Oral-Microbiota118 Oral-Microbiota119  Oral-Microbiota12 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota120 Oral-Microbiota121 Oral-Microbiota122 Oral-Microbiota123 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota124 Oral-Microbiota125 Oral-Microbiota126 Oral-Microbiota127 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota128 Oral-Microbiota129  Oral-Microbiota13 Oral-Microbiota130 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota131 Oral-Microbiota132 Oral-Microbiota133 Oral-Microbiota134 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota135 Oral-Microbiota136 Oral-Microbiota137 Oral-Microbiota138 
    ##                  1                  1                  1                  1 
    ## Oral-Microbiota139  Oral-Microbiota14 Oral-Microbiota140 Oral-Microbiota141 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota15  Oral-Microbiota16  Oral-Microbiota17  Oral-Microbiota18 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota19   Oral-Microbiota2  Oral-Microbiota20  Oral-Microbiota21 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota22  Oral-Microbiota23  Oral-Microbiota24  Oral-Microbiota25 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota26  Oral-Microbiota27  Oral-Microbiota28  Oral-Microbiota29 
    ##                  1                  1                  1                  1 
    ##   Oral-Microbiota3  Oral-Microbiota30  Oral-Microbiota31  Oral-Microbiota32 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota33  Oral-Microbiota34  Oral-Microbiota35  Oral-Microbiota36 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota37  Oral-Microbiota38  Oral-Microbiota39   Oral-Microbiota4 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota40  Oral-Microbiota41  Oral-Microbiota42  Oral-Microbiota43 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota44  Oral-Microbiota45  Oral-Microbiota46  Oral-Microbiota47 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota48  Oral-Microbiota49   Oral-Microbiota5  Oral-Microbiota50 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota51  Oral-Microbiota52  Oral-Microbiota53  Oral-Microbiota54 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota55  Oral-Microbiota56  Oral-Microbiota57  Oral-Microbiota58 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota59   Oral-Microbiota6  Oral-Microbiota60  Oral-Microbiota61 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota62  Oral-Microbiota63  Oral-Microbiota64  Oral-Microbiota65 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota66  Oral-Microbiota67  Oral-Microbiota68  Oral-Microbiota69 
    ##                  1                  1                  1                  1 
    ##   Oral-Microbiota7  Oral-Microbiota70  Oral-Microbiota71  Oral-Microbiota72 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota73  Oral-Microbiota74  Oral-Microbiota75  Oral-Microbiota76 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota77  Oral-Microbiota78  Oral-Microbiota79   Oral-Microbiota8 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota80  Oral-Microbiota81  Oral-Microbiota82  Oral-Microbiota83 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota84  Oral-Microbiota85  Oral-Microbiota86  Oral-Microbiota87 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota88  Oral-Microbiota89   Oral-Microbiota9  Oral-Microbiota90 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota91  Oral-Microbiota92  Oral-Microbiota93  Oral-Microbiota94 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota95  Oral-Microbiota96  Oral-Microbiota97  Oral-Microbiota98 
    ##                  1                  1                  1                  1 
    ##  Oral-Microbiota99 
    ##                  1

``` r
table(metadata_clean$Site)
```

    ## 
    ## Saliva 
    ##    141

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ##  NO 
    ## 141

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 141"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 141"

# Step 5: Visualize the Metadata

``` r
knitr::include_graphics("Metadata_Summary.png")
```

<img src="Metadata_Summary.png" width="1200" />

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

    ##   BioProject       Run      sample_name Age BMI Sex Region            Platform
    ## 1  PRJEB6781 ERR618776 Oral-Microbiota1  35  NA   F  V5-V7 454 GS FLX Titanium
    ## 2  PRJEB6781 ERR618777 Oral-Microbiota2  35  NA   F  V5-V7 454 GS FLX Titanium
    ## 3  PRJEB6781 ERR618778 Oral-Microbiota3  35  NA   F  V5-V7 454 GS FLX Titanium
    ## 4  PRJEB6781 ERR618779 Oral-Microbiota4  35  NA   F  V5-V7 454 GS FLX Titanium
    ## 5  PRJEB6781 ERR618780 Oral-Microbiota5  35  NA   F  V5-V7 454 GS FLX Titanium
    ## 6  PRJEB6781 ERR618781 Oral-Microbiota6  39  NA   F  V5-V7 454 GS FLX Titanium
    ##       Country   Site          public_time PCR_primers
    ## 1 Netherlands Saliva 2014-09-18T00:00:00Z          NA
    ## 2 Netherlands Saliva 2014-09-18T00:00:00Z          NA
    ## 3 Netherlands Saliva 2014-09-18T00:00:00Z          NA
    ## 4 Netherlands Saliva 2014-09-18T00:00:00Z          NA
    ## 5 Netherlands Saliva 2014-09-18T00:00:00Z          NA
    ## 6 Netherlands Saliva 2014-09-18T00:00:00Z          NA
    ##                                                  title Disease_study Status
    ## 1 Shaping the oral microbiota through intimate kissing            NO     HC
    ## 2 Shaping the oral microbiota through intimate kissing            NO     HC
    ## 3 Shaping the oral microbiota through intimate kissing            NO     HC
    ## 4 Shaping the oral microbiota through intimate kissing            NO     HC
    ## 5 Shaping the oral microbiota through intimate kissing            NO     HC
    ## 6 Shaping the oral microbiota through intimate kissing            NO     HC
    ##   matched_samples
    ## 1              NO
    ## 2              NO
    ## 3              NO
    ## 4              NO
    ## 5              NO
    ## 6              NO

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
