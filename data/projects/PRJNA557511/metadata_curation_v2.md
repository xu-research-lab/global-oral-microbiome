R Notebook: Metadata Processing Workflow
================

# Step 1: Obtain Initial Metadata from Each Project

This step sets the project directory, loads the raw SRA metadata, and
optionally merges supplement tables.

``` r
# Load libraries
library(dplyr)
library(ggplot2)
library(cowplot)

# (1) Setting working directory
# EDIT HERE: change project_id for each study
project_id <- "PRJNA557511" 

base_dir <- "E:/oral/26/metadata-2604/Oral-Metadata-YL"
setwd(file.path(base_dir, project_id))

# (2) Loading raw SRA metadata file
sra_file <- "SraRunTable.csv"
initial_meta <- read.csv(sra_file, header = TRUE)

# (3) OPTIONAL: If extra metadata (e.g., Table.x or supplemental materials) are available from the article, you may need to merge the file with initial_meta manually. Below is example code for merging. Depending on the specific case of each study, you can customize the code as needed.
supp_table <- read.csv("supplemental_table_v2.csv", header = TRUE)

mapping <- initial_meta %>%
  filter(grepl("^ZZ_ORAL_HC_|^HC__", Sample.Name)) %>% # Because there are 384 samples in initial_meta, which is inconsistent with the 132 samples provided by supplemental_table_v2.csv, the samples that do not exist in supplemental_table_v2.csv should be filtered out first. 
  mutate(number = sub("^ZZ_ORAL_HC_|^HC__", "HC_", Sample.Name)) %>% # In addition, The values in the 'Sample.Name' column of 'SraRunTable.csv' are in the HC__2 format (double underscore), while the values in the 'number' column of 'supplemental_table_v2.csv' are in the HC_2 format (single underscore). The single underscore needs to be replaced with a double underscore in order to merge 'supplemental_table_v2.csv' into 'SraRunTable.csv'.
  group_by(number) %>%
  arrange(desc(Sample.Name)) %>% 
  slice(1) %>% 
  ungroup() %>%
  select(number, Sample.Name) # The original sentence of the article is: After a rigorous diagnosis and exclusion procedures, 235 samples completed the Miseq sequencing, 44 patients with CKD and 88 healthy controls samples from Zhengzhou, China, and 59 patients with CKD and 44 HCs samples from Hangzhou, China, were included (Figure 1)

initial_meta <- supp_table %>%
  inner_join(mapping, by = "number") %>% 
  left_join(initial_meta, by = "Sample.Name") 
```

# Step 2: Data Overview

First, you need to get an overview of the data, which includes two
steps: reading the article and reading the initial metadata. This step
is very IMPORTANT!

``` r
# (1) Read the article in the following order: abstract, the beginning of the results section, and the methods. Understand the sample size, experimental design, and disease group overview. Clarify the sampling sites; for oral and gut samples, specify the exact sites like supragingival plaque, tongue, etc.

# (2) IMPORTANT: MUST go through each column of initial_meta
head(initial_meta, 10)
```

    ##    number age gender.male.0.female.1.   BMI hypertension..no.0.yes.1. WBC  RBC
    ## 1    HC_2  48                       1 22.77                         0 5.7 4.57
    ## 2    HC_3  38                       0 23.83                         0 7.8 5.28
    ## 3    HC_5  46                       0 23.74                         0 6.0 4.56
    ## 4    HC_6  51                       0 24.02                         0 5.9 4.64
    ## 5    HC_7  35                       0 24.26                         0 5.3 5.35
    ## 6    HC_8  54                       0 24.02                         0 6.1 4.50
    ## 7   HC_10  42                       0 21.60                         0 5.4 4.94
    ## 8   HC_12  41                       0 25.12                         0 9.9 5.51
    ## 9   HC_13  56                       1 23.56                         0 7.6 4.59
    ## 10  HC_14  54                       1 26.77                         0 5.8 4.44
    ##    X24hUTP  ALB  BUN Scr  UA    eGFR   TG calcium. phosphate.   Sample.Name
    ## 1        — 53.2 3.57  51 175 109.536 1.45     2.35       1.23         HC__2
    ## 2        — 52.3 5.76  83 325 102.581 0.91     2.21       1.51  ZZ_ORAL_HC_3
    ## 3        — 50.8 6.68  70 278 107.581 1.06     2.51       0.92  ZZ_ORAL_HC_5
    ## 4        — 46.8 3.85  68 303 105.113 1.46     2.01       1.10  ZZ_ORAL_HC_6
    ## 5        — 54.8 5.43  87 370  98.970 2.29     2.32       0.97         HC__7
    ## 6        — 45.4 6.19  76 302  98.322 0.78     2.07       1.23  ZZ_ORAL_HC_8
    ## 7        — 48.5 5.42  76 343 106.969 1.85     2.38       1.58        HC__10
    ## 8        — 47.8 7.65  85 349  97.591 1.72     2.41       1.43        HC__12
    ## 9        — 41.6 4.97  53 197 102.248 1.26     2.23       1.37 ZZ_ORAL_HC_13
    ## 10       — 46.9 3.88  59 227 100.099 2.11     2.06       1.20 ZZ_ORAL_HC_14
    ##            Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1  SRR10188230   AMPLICON        602 30177056 PRJNA557511 SAMN12826774
    ## 2  SRR11186065   AMPLICON        602 23638132 PRJNA557511 SAMN14209571
    ## 3  SRR11186026   AMPLICON        602 34639080 PRJNA557511 SAMN14209573
    ## 4  SRR11186025   AMPLICON        602 27324780 PRJNA557511 SAMN14209574
    ## 5  SRR10188225   AMPLICON        602 34017816 PRJNA557511 SAMN12826779
    ## 6  SRR11186024   AMPLICON        602 32211214 PRJNA557511 SAMN14209575
    ## 7  SRR10188222   AMPLICON        602 28163366 PRJNA557511 SAMN12826782
    ## 8  SRR10188219   AMPLICON        602 40264168 PRJNA557511 SAMN12826784
    ## 9  SRR11186023   AMPLICON        602 25553696 PRJNA557511 SAMN14209576
    ## 10 SRR11186022   AMPLICON        602 25400186 PRJNA557511 SAMN14209577
    ##                 BioSampleModel    Bytes
    ## 1  Metagenome or environmental 16491086
    ## 2  Metagenome or environmental 12978750
    ## 3  Metagenome or environmental 18722127
    ## 4  Metagenome or environmental 14995870
    ## 5  Metagenome or environmental 18481438
    ## 6  Metagenome or environmental 17541103
    ## 7  Metagenome or environmental 15316481
    ## 8  Metagenome or environmental 21562235
    ## 9  Metagenome or environmental 13866602
    ## 10 Metagenome or environmental 14101350
    ##                                              Center.Name Collection_Date
    ## 1  THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03
    ## 2  THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03
    ## 3  THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03
    ## 4  THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03
    ## 5  THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03
    ## 6  THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03
    ## 7  THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03
    ## 8  THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03
    ## 9  THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03
    ## 10 THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03
    ##    Consent DATASTORE.filetype DATASTORE.provider
    ## 1   public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2   public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3   public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4   public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5   public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6   public   fastq,run.zq,sra         gs,ncbi,s3
    ## 7   public   fastq,run.zq,sra         gs,ncbi,s3
    ## 8   public   fastq,run.zq,sra         gs,ncbi,s3
    ## 9   public   fastq,run.zq,sra         gs,ncbi,s3
    ## 10  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                        DATASTORE.region Experiment geo_loc_name_country
    ## 1  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908761                China
    ## 2  gs.us-east1,ncbi.public,s3.us-east-1 SRX7806332                China
    ## 3  gs.us-east1,ncbi.public,s3.us-east-1 SRX7806371                China
    ## 4  gs.us-east1,ncbi.public,s3.us-east-1 SRX7806372                China
    ## 5  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908766                China
    ## 6  gs.us-east1,ncbi.public,s3.us-east-1 SRX7806373                China
    ## 7  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908769                China
    ## 8  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908772                China
    ## 9  gs.us-east1,ncbi.public,s3.us-east-1 SRX7806374                China
    ## 10 gs.us-east1,ncbi.public,s3.us-east-1 SRX7806375                China
    ##    geo_loc_name_country_continent    geo_loc_name      HOST     Instrument
    ## 1                            Asia China:Zhengzhou Sample071 Illumina MiSeq
    ## 2                            Asia China:Zhengzhou      HC_3 Illumina MiSeq
    ## 3                            Asia China:Zhengzhou      HC_5 Illumina MiSeq
    ## 4                            Asia China:Zhengzhou      HC_6 Illumina MiSeq
    ## 5                            Asia China:Zhengzhou Sample076 Illumina MiSeq
    ## 6                            Asia China:Zhengzhou      HC_8 Illumina MiSeq
    ## 7                            Asia China:Zhengzhou Sample079 Illumina MiSeq
    ## 8                            Asia China:Zhengzhou Sample081 Illumina MiSeq
    ## 9                            Asia China:Zhengzhou     HC_13 Illumina MiSeq
    ## 10                           Asia China:Zhengzhou     HC_14 Illumina MiSeq
    ##    isolation_source          lat_lon  Library.Name LibraryLayout
    ## 1            saliva 34.44 N 113.42 E         HC__2        PAIRED
    ## 2            saliva 34.44 N 113.42 E  ZZ_ORAL_HC_3        PAIRED
    ## 3            saliva 34.44 N 113.42 E  ZZ_ORAL_HC_5        PAIRED
    ## 4            saliva 34.44 N 113.42 E  ZZ_ORAL_HC_6        PAIRED
    ## 5            saliva 34.44 N 113.42 E         HC__7        PAIRED
    ## 6            saliva 34.44 N 113.42 E  ZZ_ORAL_HC_8        PAIRED
    ## 7            saliva 34.44 N 113.42 E        HC__10        PAIRED
    ## 8            saliva 34.44 N 113.42 E        HC__12        PAIRED
    ## 9            saliva 34.44 N 113.42 E ZZ_ORAL_HC_13        PAIRED
    ## 10           saliva 34.44 N 113.42 E ZZ_ORAL_HC_14        PAIRED
    ##    LibrarySelection LibrarySource              Organism Platform
    ## 1               PCR         OTHER human oral metagenome ILLUMINA
    ## 2               PCR         OTHER human oral metagenome ILLUMINA
    ## 3               PCR         OTHER human oral metagenome ILLUMINA
    ## 4               PCR         OTHER human oral metagenome ILLUMINA
    ## 5               PCR         OTHER human oral metagenome ILLUMINA
    ## 6               PCR         OTHER human oral metagenome ILLUMINA
    ## 7               PCR         OTHER human oral metagenome ILLUMINA
    ## 8               PCR         OTHER human oral metagenome ILLUMINA
    ## 9               PCR         OTHER human oral metagenome ILLUMINA
    ## 10              PCR         OTHER human oral metagenome ILLUMINA
    ##             ReleaseDate          create_date version SRA.Study
    ## 1  2020-10-23T00:00:00Z 2019-09-26T00:29:00Z       1 SRP223308
    ## 2  2022-03-31T00:00:00Z 2020-02-26T13:17:00Z       1 SRP223308
    ## 3  2022-03-31T00:00:00Z 2020-02-26T13:15:00Z       1 SRP223308
    ## 4  2022-03-31T00:00:00Z 2020-02-26T13:15:00Z       1 SRP223308
    ## 5  2020-10-23T00:00:00Z 2019-09-26T00:29:00Z       1 SRP223308
    ## 6  2022-03-31T00:00:00Z 2020-02-26T13:15:00Z       1 SRP223308
    ## 7  2020-10-23T00:00:00Z 2019-09-26T00:29:00Z       1 SRP223308
    ## 8  2020-10-23T00:00:00Z 2019-09-26T00:29:00Z       1 SRP223308
    ## 9  2022-03-31T00:00:00Z 2020-02-26T13:15:00Z       1 SRP223308
    ## 10 2022-03-31T00:00:00Z 2020-02-26T13:14:00Z       1 SRP223308

``` r
# IMPORTANT: Check the unique values in the column of interest or query column.
lapply(initial_meta[c('isolation_source')], table) # site
```

    ## $isolation_source
    ## 
    ## saliva 
    ##    118

``` r
lapply(initial_meta[c('HOST', 'Organism')], table) # organism
```

    ## $HOST
    ## 
    ##    HC_103    HC_104    HC_105    HC_106    HC_108    HC_109    HC_110    HC_112 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_114    HC_115    HC_116    HC_117    HC_118    HC_119    HC_120    HC_122 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_126    HC_128    HC_129     HC_13    HC_130    HC_134    HC_135    HC_137 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_138     HC_14    HC_140    HC_142    HC_143    HC_145    HC_147     HC_18 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_20     HC_21     HC_22     HC_23     HC_24     HC_25     HC_26     HC_27 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_28     HC_29      HC_3     HC_30     HC_32     HC_33     HC_34     HC_35 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_38     HC_39     HC_40     HC_41     HC_42     HC_45     HC_47     HC_48 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_49      HC_5     HC_50     HC_53     HC_54     HC_55     HC_56     HC_57 
    ##         1         1         1         1         1         1         1         1 
    ##      HC_6     HC_60     HC_62     HC_65     HC_66     HC_68     HC_69     HC_70 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_71     HC_73     HC_78     HC_79      HC_8     HC_82     HC_83     HC_84 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_85     HC_87     HC_88     HC_90     HC_92     HC_93     HC_94     HC_95 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_96     HC_97     HC_98 Sample071 Sample076 Sample078 Sample079 Sample080 
    ##         1         1         1         1         1         1         1         1 
    ## Sample081 Sample084 Sample085 Sample100 Sample105 Sample106 Sample115 Sample120 
    ##         1         1         1         1         1         1         1         1 
    ## Sample121 Sample128 Sample130 Sample132 Sample133 Sample136 Sample143 Sample145 
    ##         1         1         1         1         1         1         1         1 
    ## Sample146 Sample150 Sample155 Sample158 Sample160 Sample168 
    ##         1         1         1         1         1         1 
    ## 
    ## $Organism
    ## 
    ## human oral metagenome 
    ##                   118

``` r
table(initial_meta$Sample.Name) # subject id
```

    ## 
    ##         HC__10         HC__11         HC__12         HC__15         HC__16 
    ##              1              1              1              1              1 
    ##          HC__2         HC__31         HC__36         HC__37         HC__46 
    ##              1              1              1              1              1 
    ##         HC__51         HC__52         HC__59         HC__61         HC__63 
    ##              1              1              1              1              1 
    ##         HC__64         HC__67          HC__7         HC__74         HC__76 
    ##              1              1              1              1              1 
    ##         HC__77         HC__81         HC__86         HC__89          HC__9 
    ##              1              1              1              1              1 
    ##         HC__91         HC__99 ZZ_ORAL_HC_103 ZZ_ORAL_HC_104 ZZ_ORAL_HC_105 
    ##              1              1              1              1              1 
    ## ZZ_ORAL_HC_106 ZZ_ORAL_HC_108 ZZ_ORAL_HC_109 ZZ_ORAL_HC_110 ZZ_ORAL_HC_112 
    ##              1              1              1              1              1 
    ## ZZ_ORAL_HC_114 ZZ_ORAL_HC_115 ZZ_ORAL_HC_116 ZZ_ORAL_HC_117 ZZ_ORAL_HC_118 
    ##              1              1              1              1              1 
    ## ZZ_ORAL_HC_119 ZZ_ORAL_HC_120 ZZ_ORAL_HC_122 ZZ_ORAL_HC_126 ZZ_ORAL_HC_128 
    ##              1              1              1              1              1 
    ## ZZ_ORAL_HC_129  ZZ_ORAL_HC_13 ZZ_ORAL_HC_130 ZZ_ORAL_HC_134 ZZ_ORAL_HC_135 
    ##              1              1              1              1              1 
    ## ZZ_ORAL_HC_137 ZZ_ORAL_HC_138  ZZ_ORAL_HC_14 ZZ_ORAL_HC_140 ZZ_ORAL_HC_142 
    ##              1              1              1              1              1 
    ## ZZ_ORAL_HC_143 ZZ_ORAL_HC_145 ZZ_ORAL_HC_147  ZZ_ORAL_HC_18  ZZ_ORAL_HC_20 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_21  ZZ_ORAL_HC_22  ZZ_ORAL_HC_23  ZZ_ORAL_HC_24  ZZ_ORAL_HC_25 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_26  ZZ_ORAL_HC_27  ZZ_ORAL_HC_28  ZZ_ORAL_HC_29   ZZ_ORAL_HC_3 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_30  ZZ_ORAL_HC_32  ZZ_ORAL_HC_33  ZZ_ORAL_HC_34  ZZ_ORAL_HC_35 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_38  ZZ_ORAL_HC_39  ZZ_ORAL_HC_40  ZZ_ORAL_HC_41  ZZ_ORAL_HC_42 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_45  ZZ_ORAL_HC_47  ZZ_ORAL_HC_48  ZZ_ORAL_HC_49   ZZ_ORAL_HC_5 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_50  ZZ_ORAL_HC_53  ZZ_ORAL_HC_54  ZZ_ORAL_HC_55  ZZ_ORAL_HC_56 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_57   ZZ_ORAL_HC_6  ZZ_ORAL_HC_60  ZZ_ORAL_HC_62  ZZ_ORAL_HC_65 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_66  ZZ_ORAL_HC_68  ZZ_ORAL_HC_69  ZZ_ORAL_HC_70  ZZ_ORAL_HC_71 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_73  ZZ_ORAL_HC_78  ZZ_ORAL_HC_79   ZZ_ORAL_HC_8  ZZ_ORAL_HC_82 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_83  ZZ_ORAL_HC_84  ZZ_ORAL_HC_85  ZZ_ORAL_HC_87  ZZ_ORAL_HC_88 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_90  ZZ_ORAL_HC_92  ZZ_ORAL_HC_93  ZZ_ORAL_HC_94  ZZ_ORAL_HC_95 
    ##              1              1              1              1              1 
    ##  ZZ_ORAL_HC_96  ZZ_ORAL_HC_97  ZZ_ORAL_HC_98 
    ##              1              1              1

``` r
# table(initial_meta$sequencing_meth..exp.) # output is pyrosequencing
table(initial_meta$Instrument) # output is Illumina. according to the methods, Illumina is correct.
```

    ## 
    ## Illumina MiSeq 
    ##            118

``` r
# table(initial_meta$target_subfragment..exp.)
```

# Step 3: Removing Control and Non-human Samples

First, we need to remove control samples(including both negative
controls and positive controls) and non-human samples.

``` r
# (1) EDIT HERE: Parameter setup
control_column <- "" # column name in initial_meta
# unique(initial_meta[[control_column]])
control_keywords <- c() # Enter exact keywords.  

species_column <- "" # column name in initial_meta
unique(initial_meta[[species_column]])
```

    ## NULL

``` r
non_human_keywords <- c() # Enter exact keywords

# (2) Remove controls and non-human samples. DO NOT need to modify code in this part, just run it.
n_samples <- nrow(initial_meta)
# Control evaluation
if (length(control_keywords) > 0) {
  if (is.null(control_column) || control_column == "") {
    stop("Error: control_keywords provided but control_column is missing.")
  }
  
  if (!(control_column %in% colnames(initial_meta))) {
    stop(sprintf("Error: Column '%s' not found in initial_meta.", control_column))
  }
  
  control_pattern <- paste(control_keywords, collapse = "|")
  control_text <- initial_meta[[control_column]]
  is_control <- grepl(control_pattern, control_text)
  
} else {
  is_control <- rep(FALSE, n_samples)
}

# Non-human evaluation
if (length(non_human_keywords) > 0) {
  if (is.null(species_column) || species_column == "") {
    stop("Error: non_human_keywords provided but species_column is missing.")
  }
  
  if (!(species_column %in% colnames(initial_meta))) {
    stop(sprintf("Error: Column '%s' not found in initial_meta.", species_column))
  }
  
  non_human_pattern <- paste(non_human_keywords, collapse = "|")
  species_text <- initial_meta[[species_column]]
  is_non_human <- grepl(non_human_pattern, species_text)
} else {
  is_non_human <- rep(FALSE, n_samples)
}

# Execute removal
samples_to_remove <- is_control | is_non_human
initial_meta <- initial_meta[!samples_to_remove, ]

cat("Initial samples counts:",n_samples,"\n")
```

    ## Initial samples counts: 118

``` r
cat("Removed control samples:", sum(is_control), "\n")
```

    ## Removed control samples: 0

``` r
cat("Removed non-human samples:", sum(is_non_human), "\n")
```

    ## Removed non-human samples: 0

``` r
cat("Total samples removed ", sum(samples_to_remove), "\n")
```

    ## Total samples removed  0

``` r
cat("Remaining samples after Step 2:", nrow(initial_meta), "\n")
```

    ## Remaining samples after Step 2: 118

# Step 4: Creating a Uniform Metadata Template

This step converts the raw SRA table into the standard metadata format
used in our meta-analysis. The uniform metadata should include the
following columns: 1. **bioproject**: Project ID from SRA or other
database e.g. `PRJNA1234`, `CNP1234`, ‘OEP1234’, ‘qiita_1234’ 2.
**sample_id**: SRA sample/run ID. Unique id for each sample!
e.g. `DRR1234567`, `ERR1234567` and `SRR1234567`. If no such column
exists, can use column like sample_name. 3. **subject_id**: subject id
for each person if available. e.g. `patient_01`, `Host_03_before`,
`subject12` 4. **age**: Age of the subject if available. e.g. `35`, `42`
\| not `adult`, `chilren` 5. **age_info**: If the age is not explicitly
stated but approximate age range information is available, mark it in
this column. e.g., `18-55`, `children`, `infant`, `elder` 6. **bmi**:
Body Mass Index if available. e.g. `22.5`, `18.9` \| not `normal`,
`overweight` 7. **sex**: Sex of the subject. Mandatory to write as M/F
\| not `female`, `woman` 8. **sequencing_region**: Targeted 16S rRNA
region. e.g. `V3-V4`, `V4` 9. **platform**: Sequencing platform.
e.g. `Illumina MiSeq`, `454 GS FLX Titanium`, `Ion Torrent` 10.
**country**: Country of sample origin. e.g. `China`, `Japan` \| not city
name like `Beijing` 11. **site**: Sampling site. For oral and gut sites,
specify the exact sites as detailed as possible. For example, if it is
known to be supragingival plaque, do not record it as ‘dental plaque’,
unless the original text provides only a vague description.
e.g. `saliva`, `tongue coating` \| not broad term like `oral` 12.
**site_group**: only accepts one of following three values: ‘oral’,
‘gut’, or ‘other’. 13. **public_time**: Data release date.
e.g. `2021-06-22` \| not `June 2021` 14. **pcr_primers**: PCR primers
used. e.g. `515F-806R`, `27F-338R` \| not primer sequences 15.
**title**: Title of the associated publication. 16. **disease_study**:
Whether the study involves disease samples. If it’s not clear, can
record as NA. Only three values accepted: ‘yes’, ‘no’ or NA 17.
**status**: Disease status for each sample if available, else record as
NA. Only accept following abbreviations, Write full names for all other
diseases. HC: healthy control CD: Crohn’s disease CRC: colorectal cancer
COPD: chronic obstructive pulmonary disease HIV: human immunodeficiency
virus ICU: intensive care unit OSCC: oral squamous cell carcinoma RA:
rheumatoid arthritis UC: ulcerative colitis 18. **matched_samples**:
Whether the subject has matched samples. Only accept ‘yes’ or ‘no’. 19.
**is_time_series**: Whether the study is a time series study. Only
accept ‘yes’ or ‘no’.

``` r
# Extract columns from initial_meta as follows: (1) No modification needed → extract directly. (2) Modification needed → create a new column with the prefix raw_ and mark the original column as “pending”. (3) Need to find or determine based on the article.
metadata <- data.frame(
  bioproject = initial_meta$BioProject,
  sample_id = initial_meta$Run,
  subject_id = initial_meta$Sample.Name,
  age = initial_meta$age,
  age_info = NA,
  bmi = initial_meta$BMI,
  raw_sex = initial_meta$gender.male.0.female.1., # if need modify, then create a column start with 'raw_'
  sex = "pending", # and mark as "pending"
  sequencing_region = "V3-V4",
  platform = "Illumina MiSeq", # extract from methods section, NOT initial_meta
  country = initial_meta$geo_loc_name_country,
  site = initial_meta$isolation_source,
  site_group = "pending", # modify later
  public_time = initial_meta$ReleaseDate,
  pcr_primers = "341F-805R", # Format: xxF-xxR
  title = "Characteristics of human oral microbiome and its non-invasive diagnostic value in chronic kidney disease", # extract from article
  disease_study = "yes", # determine according to the article
  status = "HC", # the original sentence of the article is: After a rigorous diagnosis and exclusion procedures, 235 samples completed the Miseq sequencing, 44 patients with CKD and 88 healthy controls samples from Zhengzhou, China, and 59 patients with CKD and 44 HCs samples from Hangzhou, China, were included (Figure 1).
  matched_samples = "pending", # modify later
  is_time_series = "no", # determine according to the article
  stringsAsFactors = FALSE
)
metadata$public_time <- as.Date(metadata$public_time)
```

# Step 5: Manually solving pending issues

This step involves two parts: (1) handle the columns with the `raw_`
prefix, and next (2) handle the `site_group` and `matched_samples`
columns.

``` r
# (1) Handle the columns with the `raw_` prefix.
# In this part, you need to customize the code based on the situation of the raw_ columns to obtain the content of the target columns. IMPORTANT: Add comments if the modification is not straightforward!!
metadata <- metadata %>%
  mutate(sex = case_when(
    raw_sex == "0" ~ "M", # 0 value represents male
    raw_sex == "1" ~ "F", # 1 value represents female
    TRUE ~ NA_character_
  )) # According to the 'gender.male.0.female.1.' column of initial_meta
```

``` r
# (2) Filling the 'site_group' and 'matched_samples' columns
# editing site_group
unique(metadata$site)
```

    ## [1] "saliva"

``` r
# EDIT HERE: add oral sites
Oral_site = unique(metadata$site)
# EDIT HERE: add gut sites
Gut_site = c()
# just run
metadata <- metadata %>%
  mutate(site_group = case_when(
    site %in% Oral_site ~ "oral",
    site %in% Gut_site ~ "gut",
    TRUE ~ "other" 
  ))

# editing matched_samples, just run.
metadata$matched_samples <- ifelse(
  ave(metadata$subject_id, metadata$subject_id, FUN = length) > 1, 
  "yes", 
  "no"
)
```

``` r
# (3) check if there is pending value in metadata. If there is pending value, go on modify the metadata.
# Function to check if a column contains "pending"
has_pending <- function(col) {
  if (is.character(col) || is.factor(col)) {
    any(col == "pending", na.rm = TRUE)
  } else {
    FALSE
  }
}

# Check for any pending value
any_pending_value <- any(sapply(metadata, has_pending))

if (!any_pending_value) {
  message("Good! No cell contains the value 'pending'.")
} else {
  pending_cols <- names(metadata)[sapply(metadata, has_pending)]
  message("Columns containing 'pending': ", paste(pending_cols, collapse = ", "))
  message("Please continue to resolve these pending columns!")
}
```

    ## Good! No cell contains the value 'pending'.

# Step 6: Visualize the Metadata

This step makes a quick set of plots for checking the distribution of
key metadata fields.

``` r
# Select columns for visualization, just run
plot_data <- metadata %>%
  select(-bioproject, -sample_id, -subject_id, -title)

create_simple_plots <- function(data) {
  plots <- list()
  continuous_vars <- c("age", "bmi")

  for (col in colnames(data)) {
    if (all(is.na(data[[col]]))) {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "this column is NA", size = 6) +
        labs(title = col) +
        theme_void()
    } else if (col %in% continuous_vars && is.numeric(data[[col]])) {
      p <- ggplot(data, aes(x = .data[[col]])) +
        geom_histogram(fill = "blue", alpha = 0.6, bins = 20, na.rm = TRUE) +
        labs(title = col, x = col, y = "Count") +
        theme_minimal()
    } else {
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
p <- plot_grid(plots$age,
               plots$age_info,
               plots$bmi,
               plots$sex,
               plots$site,
               plots$site_group,
               plots$sequencing_region,
               plots$platform,
               plots$country,
               plots$pcr_primers,
               plots$status,
               plots$matched_samples,
               ncol = 3, align = "hv", axis = "lb")
ggsave("Metadata_Summary_v2.png", plot = p, width = 8, height = 6, dpi = 150)
```

``` r
# Insert Metadata_Summary.png into the current Rmd output
knitr::include_graphics("Metadata_Summary_v2.png")
```

<img src="Metadata_Summary_v2.png" alt="" width="1200" />

# Step 7: Sort and export the metadata

This step saves the cleaned metadata table as the final `metadata.csv`
file.

``` r
# Reorder the data frame
fixed_cols <- c("bioproject", "sample_id", "subject_id", "age", "age_info",
                "bmi", "sex", "platform", "sequencing_region", "pcr_primers", 
                "country", "site", "site_group", "disease_study", "status", 
                "title", "public_time", "matched_samples", "is_time_series")  
metadata <- metadata %>%
  select(all_of(fixed_cols), everything())

# save metadata 
# head(metadata)
write.csv(metadata, file = "metadata_v2.csv", quote = TRUE)
```
