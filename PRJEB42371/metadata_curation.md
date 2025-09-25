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
setwd("E:/bioproject/PRJEB42371")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = T) 

# Remove entries non-oral samples
initial_meta <- initial_meta %>%
  filter(target_gene..exp. != "ITS",
         env_material != "sterile water")

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen barcode..exp.    Bases BioProject
    ## 1 ERR5050934   AMPLICON        250         AAAAA  8649606 PRJEB42371
    ## 2 ERR5050935   AMPLICON        250         AAAAA 12412209 PRJEB42371
    ## 3 ERR5050936   AMPLICON        250         AAAAA  5886944 PRJEB42371
    ## 4 ERR5050937   AMPLICON        250         AAAAA  7570303 PRJEB42371
    ## 5 ERR5050938   AMPLICON        250         AAAAA  5406358 PRJEB42371
    ## 6 ERR5050939   AMPLICON        250         AAAAA  7028778 PRJEB42371
    ##      BioSample   Bytes center_name..exp.
    ## 1 SAMEA7767094 6297508     U Puerto Rico
    ## 2 SAMEA7767095 9072165     U Puerto Rico
    ## 3 SAMEA7767096 3212096     U Puerto Rico
    ## 4 SAMEA7767097 4162448     U Puerto Rico
    ## 5 SAMEA7767098 3142792     U Puerto Rico
    ## 6 SAMEA7767099 3949993     U Puerto Rico
    ##                                          Center.Name
    ## 1 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ## 2 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ## 3 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ## 4 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ## 5 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ## 6 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ##     center_project_name..exp. collection_timestamp Consent     COUNTRY
    ## 1 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ## 2 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ## 3 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ## 4 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ## 5 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ## 6 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   dna_extracted elevation          empo_1            empo_2           empo_3
    ## 1          TRUE         8 Host-associated Animal-associated Animal secretion
    ## 2          TRUE         8 Host-associated Animal-associated Animal secretion
    ## 3          TRUE         8 Host-associated Animal-associated Animal secretion
    ## 4          TRUE         8 Host-associated Animal-associated Animal secretion
    ## 5          TRUE         8 Host-associated Animal-associated Animal secretion
    ## 6          TRUE         8 Host-associated Animal-associated Animal secretion
    ##   ENA.FIRST.PUBLIC..run. ENA.LAST.UPDATE..run.   env_biome
    ## 1              2021/1/13              2021/1/6 urban biome
    ## 2              2021/1/13              2021/1/6 urban biome
    ## 3              2021/1/13              2021/1/6 urban biome
    ## 4              2021/1/13              2021/1/6 urban biome
    ## 5              2021/1/13              2021/1/6 urban biome
    ## 6              2021/1/13              2021/1/6 urban biome
    ##                env_feature env_material env_package Experiment  external_id
    ## 1 human-associated habitat       saliva  human-oral ERX4857266 SAMEA7767094
    ## 2 human-associated habitat       saliva  human-oral ERX4857267 SAMEA7767095
    ## 3 human-associated habitat       saliva  human-oral ERX4857268 SAMEA7767096
    ## 4 human-associated habitat       saliva  human-oral ERX4857269 SAMEA7767097
    ## 5 human-associated habitat       saliva  human-oral ERX4857270 SAMEA7767098
    ## 6 human-associated habitat       saliva  human-oral ERX4857271 SAMEA7767099
    ##   geo_loc_name_country geo_loc_name_country_continent         geo_loc_name
    ## 1          Puerto Rico                  North America Puerto Rico:San Juan
    ## 2          Puerto Rico                  North America Puerto Rico:San Juan
    ## 3          Puerto Rico                  North America Puerto Rico:San Juan
    ## 4          Puerto Rico                  North America Puerto Rico:San Juan
    ## 5          Puerto Rico                  North America Puerto Rico:San Juan
    ## 6          Puerto Rico                  North America Puerto Rico:San Juan
    ##   host_age_units host_height_units host_subject_id host_weight_units
    ## 1          years                cm              O1                kg
    ## 2          years                cm             O10                kg
    ## 3          years                cm            O100                kg
    ## 4          years                cm            O101                kg
    ## 5          years                cm            O102                kg
    ## 6          years                cm            O103                kg
    ##                                          INSDC_center_name   INSDC_first_public
    ## 1 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ## 2 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ## 3 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ## 4 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ## 5 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ## 6 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ##      INSDC_last_update INSDC_status     Instrument latitude Library.Name
    ## 1 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837     13193.O1
    ## 2 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837    13193.O10
    ## 3 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837   13193.O100
    ## 4 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837   13193.O101
    ## 5 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837   13193.O102
    ## 6 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837   13193.O103
    ##   LibraryLayout LibrarySelection LibrarySource longitude
    ## 1        SINGLE              PCR   METAGENOMIC -66.07453
    ## 2        SINGLE              PCR   METAGENOMIC -66.07453
    ## 3        SINGLE              PCR   METAGENOMIC -66.07453
    ## 4        SINGLE              PCR   METAGENOMIC -66.07453
    ## 5        SINGLE              PCR   METAGENOMIC -66.07453
    ## 6        SINGLE              PCR   METAGENOMIC -66.07453
    ##                  Organism                         physical_specimen_location
    ## 1 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ## 2 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ## 3 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ## 4 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ## 5 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ## 6 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ##   physical_specimen_remaining Platform          ReleaseDate
    ## 1                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ## 2                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ## 3                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ## 4                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ## 5                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ## 6                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ##             run_center..exp.          create_date version
    ## 1 Louisiana State University 2021-01-28T10:32:00Z       1
    ## 2 Louisiana State University 2021-01-28T10:32:00Z       1
    ## 3 Louisiana State University 2021-01-28T10:31:00Z       1
    ## 4 Louisiana State University 2021-01-28T10:31:00Z       1
    ## 5 Louisiana State University 2021-01-28T10:31:00Z       1
    ## 6 Louisiana State University 2021-01-28T10:32:00Z       1
    ##           run_prefix..exp.  Sample.Name                Sample_name sample_type
    ## 1      GodO01_S175_L001_R1 SAMEA7767094   qiita_sid_13193:13193.O1      saliva
    ## 2      GodO10_S184_L001_R1 SAMEA7767095  qiita_sid_13193:13193.O10      saliva
    ## 3 VitO100_S221_L001_R1_001 SAMEA7767096 qiita_sid_13193:13193.O100      saliva
    ## 4 VitO101_S233_L001_R1_001 SAMEA7767097 qiita_sid_13193:13193.O101      saliva
    ## 5 VitO102_S245_L001_R1_001 SAMEA7767098 qiita_sid_13193:13193.O102      saliva
    ## 6 VitO103_S257_L001_R1_001 SAMEA7767099 qiita_sid_13193:13193.O103      saliva
    ##     sequencing_meth..exp. SRA.Study               Submitter_Id
    ## 1 Sequencing by synthesis ERP126217   qiita_sid_13193:13193.O1
    ## 2 Sequencing by synthesis ERP126217  qiita_sid_13193:13193.O10
    ## 3 Sequencing by synthesis ERP126217 qiita_sid_13193:13193.O100
    ## 4 Sequencing by synthesis ERP126217 qiita_sid_13193:13193.O101
    ## 5 Sequencing by synthesis ERP126217 qiita_sid_13193:13193.O102
    ## 6 Sequencing by synthesis ERP126217 qiita_sid_13193:13193.O103
    ##   target_gene..exp. target_subfragment..exp. age_category
    ## 1          16S rRNA                       V4        31_40
    ## 2          16S rRNA                       V4        41_49
    ## 3          16S rRNA                       V4        21_30
    ## 4          16S rRNA                       V4        31_40
    ## 5          16S rRNA                       V4        31_40
    ## 6          16S rRNA                       V4        21_30
    ##   alcohol_consumption_last_month     Host_age host_body_habitat
    ## 1                            yes not provided     UBERON:saliva
    ## 2                             no not provided     UBERON:saliva
    ## 3                            yes not provided     UBERON:saliva
    ## 4                            yes not provided     UBERON:saliva
    ## 5                            yes not provided     UBERON:saliva
    ## 6                            yes not provided     UBERON:saliva
    ##   host_body_mass_index host_body_product host_body_site host_common_name
    ## 1         not provided     UBERON:saliva  UBERON:saliva            human
    ## 2         not provided     UBERON:saliva  UBERON:saliva            human
    ## 3         not provided     UBERON:saliva  UBERON:saliva            human
    ## 4         not provided     UBERON:saliva  UBERON:saliva            human
    ## 5         not provided     UBERON:saliva  UBERON:saliva            human
    ## 6         not provided     UBERON:saliva  UBERON:saliva            human
    ##    host_height host_life_stage host_scientific_name host_taxid  host_weight
    ## 1 not provided    not provided         Homo sapiens       9606 not provided
    ## 2 not provided    not provided         Homo sapiens       9606 not provided
    ## 3 not provided    not provided         Homo sapiens       9606 not provided
    ## 4 not provided    not provided         Homo sapiens       9606 not provided
    ## 5 not provided    not provided         Homo sapiens       9606 not provided
    ## 6 not provided    not provided         Homo sapiens       9606 not provided
    ##   perio_3_categories    sex Smoking perio_4_categories ENA_first_public
    ## 1                 no female      no                 no        2021/1/13
    ## 2                 no female      no                 no        2021/1/13
    ## 3                 no   male      no                 no        2021/1/13
    ## 4                 no   male      no                 no        2021/1/13
    ## 5                 no   male      no                 no        2021/1/13
    ## 6                 no   male     yes                 no        2021/1/13
    ##   ena_last_update INSDC_center_alias oral_sex_practices        primer..exp.
    ## 1        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ## 2        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ## 3        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ## 4        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ## 5        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ## 6        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ##   linker..exp.                                  PCR_Primer..exp.
    ## 1           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ## 2           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ## 3           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ## 4           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ## 5           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ## 6           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ##   barcodesequence..exp. ena_first_public ena_last_update.1 Scientific_Name
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
  sample_name = initial_meta$host_subject_id,
  Age = initial_meta$age_category, 
  BMI = NA,  
  Sex = initial_meta$sex,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$Instrument,
  Country = initial_meta$COUNTRY,
  Site = initial_meta$env_material,
  public_time = initial_meta$ReleaseDate,
  Status = initial_meta$perio_4_categories,
  stringsAsFactors = FALSE
)

# Use abbreviations to indicate gender
# Remove the rows where the target_gene (exp) column equals "16S rRNA"
# At the same time, keep the original filter env_material != "sterile water"
metadata <- metadata %>%
  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "F:5'-GTGYCAGCMGCCGCGGTAA-3'; R:5'-GGACTACNVGGGTWTCTAAT-3'"
metadata$title <- "Ecological competition in the oral mycobiome of Hispanic adults living in Puerto Rico associates with periodontitis"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata <- metadata %>%
  mutate(
    Status = case_when(
      Status %in% c("moderate", "mild", "severe") ~ "PD",
      Status == "no" ~ "HC",
      Status == "" ~ NA_character_,
      TRUE ~ Status
    )
  )

#extract subject id by remove prefix of sample_name
metadata$sample_name <- gsub("O", "", metadata$sample_name)

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##   1  10 100 101 102 103 104 105 106 107 108 109  11 110 111 112 113 114 115 116 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 117 118 119  12 120 121 122 123 124 125 126 127 128 129  13 130 131 132 133 134 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 135  14  15  16  17  18  19   2  20  21  22  23  24  25  26  27  28  29   3  30 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  31  32  33  34  35  36  37  38  39   4  40  41  42  43  44  45  46  47  48  49 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##   5  50  51  52  53  54  55  56  57  58  59   6  60  61  62  63  64  65  66  67 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  68  69   7  70  71  72  73  74  75  76  77  78  79   8  80  81  82  83  84  85 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  86  87  88  89   9  90  91  92  93  94  95  96  97  98  99 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1

``` r
table(metadata$Site)
```

    ## 
    ## saliva 
    ##    135

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
    ##   1  10 100 101 102 103 104 105 106 107 108 109  11 110 111 112 113 114 115 116 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 117 118 119  12 120 121 122 123 124 125 126 127 128 129  13 130 131 132 133 134 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 135  14  15  16  17  18  19   2  20  21  22  23  24  25  26  27  28  29   3  30 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  31  32  33  34  35  36  37  38  39   4  40  41  42  43  44  45  46  47  48  49 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##   5  50  51  52  53  54  55  56  57  58  59   6  60  61  62  63  64  65  66  67 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  68  69   7  70  71  72  73  74  75  76  77  78  79   8  80  81  82  83  84  85 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  86  87  88  89   9  90  91  92  93  94  95  96  97  98  99 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1

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
    ## 135

``` r
# Identify samples for removal: single samples excluding saliva and tongue coat
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & 
                     !(metadata$Site %in% c("saliva", "tongue coat", "Tongue coating"))

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
    ##   1  10 100 101 102 103 104 105 106 107 108 109  11 110 111 112 113 114 115 116 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 117 118 119  12 120 121 122 123 124 125 126 127 128 129  13 130 131 132 133 134 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 135  14  15  16  17  18  19   2  20  21  22  23  24  25  26  27  28  29   3  30 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  31  32  33  34  35  36  37  38  39   4  40  41  42  43  44  45  46  47  48  49 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##   5  50  51  52  53  54  55  56  57  58  59   6  60  61  62  63  64  65  66  67 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  68  69   7  70  71  72  73  74  75  76  77  78  79   8  80  81  82  83  84  85 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  86  87  88  89   9  90  91  92  93  94  95  96  97  98  99 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1

``` r
table(metadata_clean$Site)
```

    ## 
    ## saliva 
    ##    135

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ##  NO 
    ## 135

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 135"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 135"

# Step 5: Visualize the Metadata

``` r
# Select columns for visualization
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
        theme(axis.text.x = element_text(angle = 15, hjust = 1))
    }
    plots[[col]] <- p
  }
  
  return(plots)
}

# Generate and arrange plots
plots <- create_simple_plots(plot_data)
# grid.arrange(grobs = plots, ncol = 3, top = "Metadata Summary", align = "hv")
plot_grid(plots$Age, plots$BMI, plots$Sex, plots$Region, plots$Platform, 
          plots$Country, plots$Site, plots$public_time, 
          ncol = 3, align = "hv", axis = "lb")
```

<img src="metadata_curation_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />
\# Step 6: Finalize and Export the Metadata

``` r
# Reorder columns to match the desired format
metadata <- metadata[, c(
  "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
  "Region", "Platform", "Country", "Site", "public_time", 
  "PCR_primers", "title", "Disease_study", "Status", "matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
```

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
