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
supp_table <- read.csv("supplemental_table_v2.csv", header = TRUE) # Supplemental_table_v2.csv is not from this article Oral microbiome characteristics in patients with autoimmune hepatitis. It comes from another article Characteristics of human oral microbiome and its non-invasive diagnostic value in chronic kidney disease. Because Characteristics of human oral microbiome and its non-invasive diagnostic value in chronic kidney disease this article is also used by PRJNA557511 (Article: The raw Illumina read data for all samples were deposited in the European Bioinformatics Institute European Nucleotide Archive database under the accession number PRJNA557511(HC group) and PRJNA649074(CKD group).). According to supplemental_table_v2.csv, it can correspond to a part of the samples in SraRunTable.csv of PRJNA557511.
initial_meta <- initial_meta %>%
  mutate(match_key = sub("^ZZ_ORAL_HC_", "HC_", Sample.Name)) %>%
  left_join(supp_table, by = c("match_key" = "number")) %>%
  select(-match_key) # There are 99 samples with the prefix of HC__ and 69 samples with the prefix of AIH__ in the Sample.Name column of initial_metadata, which are not mentioned in the article. According to the prefix part, they are classified into HC group and autoimmune hepatitis group in the status column. 
# In the Sample.Name column, there are 68 values prefixed with ZZ_ORAL_AIH, 136 values prefixed with ZZ_ORAL_HC, and 12 values prefixed with HZ_ORAL_hzAIH_.These samples are clearly mentioned in the article. The original text is: we collected 204 saliva samples including 68 AIH patients and 136 healthy controls and performed microbial MiSeq sequencing after screening. Moreover, we collected samples of 12 AIH patients from Hangzhou for cross-regional validation.
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

    ##            Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1  SRR10188199   AMPLICON        602 37266810 PRJNA557511 SAMN12826803
    ## 2  SRR10188200   AMPLICON        602 28186844 PRJNA557511 SAMN12826802
    ## 3  SRR10188201   AMPLICON        602 34899746 PRJNA557511 SAMN12826801
    ## 4  SRR10188202   AMPLICON        602 55844530 PRJNA557511 SAMN12826800
    ## 5  SRR10188203   AMPLICON        602 31148082 PRJNA557511 SAMN12826799
    ## 6  SRR10188204   AMPLICON        602 35811776 PRJNA557511 SAMN12826798
    ## 7  SRR10188205   AMPLICON        602 32898096 PRJNA557511 SAMN12826797
    ## 8  SRR10188206   AMPLICON        602 53542482 PRJNA557511 SAMN12826796
    ## 9  SRR10188207   AMPLICON        602 53799536 PRJNA557511 SAMN12826795
    ## 10 SRR10188208   AMPLICON        602 28433062 PRJNA557511 SAMN12826794
    ##                 BioSampleModel    Bytes
    ## 1  Metagenome or environmental 20236127
    ## 2  Metagenome or environmental 15594036
    ## 3  Metagenome or environmental 19169390
    ## 4  Metagenome or environmental 30563896
    ## 5  Metagenome or environmental 17239379
    ## 6  Metagenome or environmental 19661010
    ## 7  Metagenome or environmental 17894644
    ## 8  Metagenome or environmental 28940868
    ## 9  Metagenome or environmental 29365700
    ## 10 Metagenome or environmental 15573194
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
    ## 1  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908792                China
    ## 2  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908791                China
    ## 3  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908790                China
    ## 4  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908789                China
    ## 5  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908788                China
    ## 6  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908787                China
    ## 7  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908786                China
    ## 8  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908785                China
    ## 9  gs.us-east1,ncbi.public,s3.us-east-1 SRX6908784                China
    ## 10 gs.us-east1,ncbi.public,s3.us-east-1 SRX6908783                China
    ##    geo_loc_name_country_continent    geo_loc_name      HOST     Instrument
    ## 1                            Asia China:Zhengzhou Sample100 Illumina MiSeq
    ## 2                            Asia China:Zhengzhou Sample099 Illumina MiSeq
    ## 3                            Asia China:Zhengzhou Sample098 Illumina MiSeq
    ## 4                            Asia China:Zhengzhou Sample097 Illumina MiSeq
    ## 5                            Asia China:Zhengzhou Sample096 Illumina MiSeq
    ## 6                            Asia China:Zhengzhou Sample095 Illumina MiSeq
    ## 7                            Asia China:Zhengzhou Sample094 Illumina MiSeq
    ## 8                            Asia China:Zhengzhou Sample093 Illumina MiSeq
    ## 9                            Asia China:Zhengzhou Sample092 Illumina MiSeq
    ## 10                           Asia China:Zhengzhou Sample091 Illumina MiSeq
    ##    isolation_source          lat_lon Library.Name LibraryLayout
    ## 1            saliva 34.44 N 113.42 E       HC__31        PAIRED
    ## 2            saliva 34.44 N 113.42 E       HC__30        PAIRED
    ## 3            saliva 34.44 N 113.42 E       HC__29        PAIRED
    ## 4            saliva 34.44 N 113.42 E       HC__28        PAIRED
    ## 5            saliva 34.44 N 113.42 E       HC__27        PAIRED
    ## 6            saliva 34.44 N 113.42 E       HC__26        PAIRED
    ## 7            saliva 34.44 N 113.42 E       HC__25        PAIRED
    ## 8            saliva 34.44 N 113.42 E       HC__24        PAIRED
    ## 9            saliva 34.44 N 113.42 E       HC__23        PAIRED
    ## 10           saliva 34.44 N 113.42 E       HC__22        PAIRED
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
    ##             ReleaseDate          create_date version Sample.Name SRA.Study age
    ## 1  2020-10-23T00:00:00Z 2019-09-26T00:28:00Z       1      HC__31 SRP223308  NA
    ## 2  2020-10-23T00:00:00Z 2019-09-26T00:28:00Z       1      HC__30 SRP223308  NA
    ## 3  2020-10-23T00:00:00Z 2019-09-26T00:28:00Z       1      HC__29 SRP223308  NA
    ## 4  2020-10-23T00:00:00Z 2019-09-26T00:28:00Z       1      HC__28 SRP223308  NA
    ## 5  2020-10-23T00:00:00Z 2019-09-26T00:28:00Z       1      HC__27 SRP223308  NA
    ## 6  2020-10-23T00:00:00Z 2019-09-26T00:28:00Z       1      HC__26 SRP223308  NA
    ## 7  2020-10-23T00:00:00Z 2019-09-26T00:28:00Z       1      HC__25 SRP223308  NA
    ## 8  2020-10-23T00:00:00Z 2019-09-26T00:29:00Z       1      HC__24 SRP223308  NA
    ## 9  2020-10-23T00:00:00Z 2019-09-26T00:29:00Z       1      HC__23 SRP223308  NA
    ## 10 2020-10-23T00:00:00Z 2019-09-26T00:28:00Z       1      HC__22 SRP223308  NA
    ##    gender.male.0.female.1. BMI hypertension..no.0.yes.1. WBC RBC X24hUTP ALB
    ## 1                       NA  NA                        NA  NA  NA    <NA>  NA
    ## 2                       NA  NA                        NA  NA  NA    <NA>  NA
    ## 3                       NA  NA                        NA  NA  NA    <NA>  NA
    ## 4                       NA  NA                        NA  NA  NA    <NA>  NA
    ## 5                       NA  NA                        NA  NA  NA    <NA>  NA
    ## 6                       NA  NA                        NA  NA  NA    <NA>  NA
    ## 7                       NA  NA                        NA  NA  NA    <NA>  NA
    ## 8                       NA  NA                        NA  NA  NA    <NA>  NA
    ## 9                       NA  NA                        NA  NA  NA    <NA>  NA
    ## 10                      NA  NA                        NA  NA  NA    <NA>  NA
    ##    BUN Scr UA eGFR TG calcium. phosphate.
    ## 1   NA  NA NA   NA NA       NA         NA
    ## 2   NA  NA NA   NA NA       NA         NA
    ## 3   NA  NA NA   NA NA       NA         NA
    ## 4   NA  NA NA   NA NA       NA         NA
    ## 5   NA  NA NA   NA NA       NA         NA
    ## 6   NA  NA NA   NA NA       NA         NA
    ## 7   NA  NA NA   NA NA       NA         NA
    ## 8   NA  NA NA   NA NA       NA         NA
    ## 9   NA  NA NA   NA NA       NA         NA
    ## 10  NA  NA NA   NA NA       NA         NA

``` r
# IMPORTANT: Check the unique values in the column of interest or query column.
lapply(initial_meta[c('isolation_source')], table) # site
```

    ## $isolation_source
    ## 
    ## saliva 
    ##    384

``` r
lapply(initial_meta[c('HOST', 'Organism')], table) # organism
```

    ## $HOST
    ## 
    ##     AIH_1    AIH_10    AIH_11    AIH_12    AIH_13    AIH_15    AIH_16    AIH_17 
    ##         1         1         1         1         1         1         1         1 
    ##    AIH_18    AIH_19     AIH_2    AIH_20    AIH_21    AIH_22    AIH_23    AIH_24 
    ##         1         1         1         1         1         1         1         1 
    ##    AIH_25    AIH_26    AIH_27    AIH_28    AIH_29     AIH_3    AIH_30    AIH_31 
    ##         1         1         1         1         1         1         1         1 
    ##    AIH_32    AIH_33    AIH_34    AIH_35    AIH_36    AIH_37    AIH_38    AIH_39 
    ##         1         1         1         1         1         1         1         1 
    ##     AIH_4    AIH_40    AIH_41    AIH_42    AIH_43    AIH_44    AIH_45    AIH_46 
    ##         1         1         1         1         1         1         1         1 
    ##    AIH_47    AIH_48    AIH_49     AIH_5    AIH_50    AIH_51    AIH_52    AIH_53 
    ##         1         1         1         1         1         1         1         1 
    ##    AIH_54    AIH_55    AIH_56    AIH_57    AIH_58    AIH_59     AIH_6    AIH_60 
    ##         1         1         1         1         1         1         1         1 
    ##    AIH_61    AIH_62    AIH_63    AIH_64    AIH_65    AIH_66    AIH_67    AIH_68 
    ##         1         1         1         1         1         1         1         1 
    ##    AIH_69     AIH_7     AIH_8     AIH_9      HC_1    HC_103    HC_104    HC_105 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_106    HC_108    HC_109    HC_110    HC_111    HC_112    HC_114    HC_115 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_116    HC_117    HC_118    HC_119    HC_120    HC_122    HC_126    HC_128 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_129     HC_13    HC_130    HC_131    HC_134    HC_135    HC_136    HC_137 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_138    HC_139     HC_14    HC_140    HC_142    HC_143    HC_144    HC_145 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_147    HC_148    HC_149    HC_150    HC_152    HC_153    HC_154    HC_156 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_157    HC_159    HC_160    HC_161    HC_162    HC_163    HC_164    HC_166 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_167    HC_168    HC_169     HC_17    HC_171    HC_173    HC_174    HC_176 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_177    HC_179     HC_18    HC_180    HC_181    HC_182    HC_183    HC_184 
    ##         1         1         1         1         1         1         1         1 
    ##    HC_185    HC_186     HC_19     HC_20     HC_21     HC_22     HC_23     HC_24 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_25     HC_26     HC_27     HC_28     HC_29      HC_3     HC_30     HC_32 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_33     HC_34     HC_35     HC_38     HC_39      HC_4     HC_40     HC_41 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_42     HC_43     HC_45     HC_47     HC_48     HC_49      HC_5     HC_50 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_53     HC_54     HC_55     HC_56     HC_57     HC_58      HC_6     HC_60 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_62     HC_65     HC_66     HC_68     HC_69     HC_70     HC_71     HC_72 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_73     HC_75     HC_78     HC_79      HC_8     HC_80     HC_82     HC_83 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_84     HC_85     HC_87     HC_88     HC_90     HC_92     HC_93     HC_94 
    ##         1         1         1         1         1         1         1         1 
    ##     HC_95     HC_96     HC_97     HC_98   hzAIH_1  hzAIH_10  hzAIH_11  hzAIH_12 
    ##         1         1         1         1         1         1         1         1 
    ##   hzAIH_2   hzAIH_3   hzAIH_4   hzAIH_5   hzAIH_6   hzAIH_7   hzAIH_8   hzAIH_9 
    ##         1         1         1         1         1         1         1         1 
    ## Sample001 Sample002 Sample003 Sample004 Sample005 Sample006 Sample007 Sample008 
    ##         1         1         1         1         1         1         1         1 
    ## Sample009 Sample010 Sample011 Sample012 Sample013 Sample014 Sample015 Sample016 
    ##         1         1         1         1         1         1         1         1 
    ## Sample017 Sample018 Sample019 Sample020 Sample021 Sample022 Sample023 Sample024 
    ##         1         1         1         1         1         1         1         1 
    ## Sample025 Sample026 Sample027 Sample028 Sample029 Sample030 Sample031 Sample032 
    ##         1         1         1         1         1         1         1         1 
    ## Sample033 Sample034 Sample035 Sample036 Sample037 Sample038 Sample039 Sample040 
    ##         1         1         1         1         1         1         1         1 
    ## Sample041 Sample042 Sample043 Sample044 Sample045 Sample046 Sample047 Sample048 
    ##         1         1         1         1         1         1         1         1 
    ## Sample049 Sample050 Sample051 Sample052 Sample053 Sample054 Sample055 Sample056 
    ##         1         1         1         1         1         1         1         1 
    ## Sample057 Sample058 Sample059 Sample060 Sample061 Sample062 Sample063 Sample064 
    ##         1         1         1         1         1         1         1         1 
    ## Sample065 Sample066 Sample067 Sample068 Sample069 Sample070 Sample071 Sample072 
    ##         1         1         1         1         1         1         1         1 
    ## Sample073 Sample074 Sample075 Sample076 Sample077 Sample078 Sample079 Sample080 
    ##         1         1         1         1         1         1         1         1 
    ## Sample081 Sample082 Sample083 Sample084 Sample085 Sample086 Sample087 Sample088 
    ##         1         1         1         1         1         1         1         1 
    ## Sample089 Sample090 Sample091 Sample092 Sample093 Sample094 Sample095 Sample096 
    ##         1         1         1         1         1         1         1         1 
    ## Sample097 Sample098 Sample099 Sample100 Sample101 Sample102 Sample103 Sample104 
    ##         1         1         1         1         1         1         1         1 
    ## Sample105 Sample106 Sample107 Sample108 Sample109 Sample110 Sample111 Sample112 
    ##         1         1         1         1         1         1         1         1 
    ## Sample113 Sample114 Sample115 Sample116 Sample117 Sample118 Sample119 Sample120 
    ##         1         1         1         1         1         1         1         1 
    ## Sample121 Sample122 Sample123 Sample124 Sample125 Sample126 Sample127 Sample128 
    ##         1         1         1         1         1         1         1         1 
    ## Sample129 Sample130 Sample131 Sample132 Sample133 Sample134 Sample135 Sample136 
    ##         1         1         1         1         1         1         1         1 
    ## Sample137 Sample138 Sample139 Sample140 Sample141 Sample142 Sample143 Sample144 
    ##         1         1         1         1         1         1         1         1 
    ## Sample145 Sample146 Sample147 Sample148 Sample149 Sample150 Sample151 Sample152 
    ##         1         1         1         1         1         1         1         1 
    ## Sample153 Sample154 Sample155 Sample156 Sample157 Sample158 Sample159 Sample160 
    ##         1         1         1         1         1         1         1         1 
    ## Sample161 Sample162 Sample163 Sample164 Sample165 Sample166 Sample167 Sample168 
    ##         1         1         1         1         1         1         1         1 
    ## 
    ## $Organism
    ## 
    ## human oral metagenome 
    ##                   384

``` r
table(initial_meta$Sample.Name) # subject id
```

    ## 
    ##           AIH__1          AIH__10          AIH__11          AIH__12 
    ##                1                1                1                1 
    ##          AIH__13          AIH__14          AIH__15          AIH__16 
    ##                1                1                1                1 
    ##          AIH__17          AIH__18          AIH__19           AIH__2 
    ##                1                1                1                1 
    ##          AIH__20          AIH__21          AIH__22          AIH__23 
    ##                1                1                1                1 
    ##          AIH__24          AIH__25          AIH__26          AIH__27 
    ##                1                1                1                1 
    ##          AIH__28          AIH__29           AIH__3          AIH__30 
    ##                1                1                1                1 
    ##          AIH__31          AIH__32          AIH__33          AIH__34 
    ##                1                1                1                1 
    ##          AIH__35          AIH__36          AIH__37          AIH__38 
    ##                1                1                1                1 
    ##          AIH__39           AIH__4          AIH__40          AIH__41 
    ##                1                1                1                1 
    ##          AIH__42          AIH__43          AIH__44          AIH__45 
    ##                1                1                1                1 
    ##          AIH__46          AIH__47          AIH__48          AIH__49 
    ##                1                1                1                1 
    ##           AIH__5          AIH__50          AIH__51          AIH__52 
    ##                1                1                1                1 
    ##          AIH__53          AIH__54          AIH__55          AIH__56 
    ##                1                1                1                1 
    ##          AIH__57          AIH__58          AIH__59           AIH__6 
    ##                1                1                1                1 
    ##          AIH__60          AIH__61          AIH__62          AIH__63 
    ##                1                1                1                1 
    ##          AIH__64          AIH__65          AIH__66          AIH__67 
    ##                1                1                1                1 
    ##          AIH__68          AIH__69           AIH__7           AIH__8 
    ##                1                1                1                1 
    ##           AIH__9            HC__1           HC__10           HC__11 
    ##                1                1                1                1 
    ##           HC__12           HC__13           HC__14           HC__15 
    ##                1                1                1                1 
    ##           HC__16           HC__17           HC__18           HC__19 
    ##                1                1                1                1 
    ##            HC__2           HC__20           HC__21           HC__22 
    ##                1                1                1                1 
    ##           HC__23           HC__24           HC__25           HC__26 
    ##                1                1                1                1 
    ##           HC__27           HC__28           HC__29            HC__3 
    ##                1                1                1                1 
    ##           HC__30           HC__31           HC__32           HC__33 
    ##                1                1                1                1 
    ##           HC__34           HC__35           HC__36           HC__37 
    ##                1                1                1                1 
    ##           HC__38           HC__39            HC__4           HC__40 
    ##                1                1                1                1 
    ##           HC__41           HC__42           HC__43           HC__44 
    ##                1                1                1                1 
    ##           HC__45           HC__46           HC__47           HC__48 
    ##                1                1                1                1 
    ##           HC__49            HC__5           HC__50           HC__51 
    ##                1                1                1                1 
    ##           HC__52           HC__53           HC__54           HC__55 
    ##                1                1                1                1 
    ##           HC__56           HC__57           HC__58           HC__59 
    ##                1                1                1                1 
    ##            HC__6           HC__60           HC__61           HC__62 
    ##                1                1                1                1 
    ##           HC__63           HC__64           HC__65           HC__66 
    ##                1                1                1                1 
    ##           HC__67           HC__68           HC__69            HC__7 
    ##                1                1                1                1 
    ##           HC__70           HC__71           HC__72           HC__73 
    ##                1                1                1                1 
    ##           HC__74           HC__75           HC__76           HC__77 
    ##                1                1                1                1 
    ##           HC__78           HC__79            HC__8           HC__80 
    ##                1                1                1                1 
    ##           HC__81           HC__82           HC__83           HC__84 
    ##                1                1                1                1 
    ##           HC__85           HC__86           HC__87           HC__88 
    ##                1                1                1                1 
    ##           HC__89            HC__9           HC__90           HC__91 
    ##                1                1                1                1 
    ##           HC__92           HC__93           HC__94           HC__95 
    ##                1                1                1                1 
    ##           HC__96           HC__97           HC__98           HC__99 
    ##                1                1                1                1 
    ##  HZ_ORAL_hzAIH_1 HZ_ORAL_hzAIH_10 HZ_ORAL_hzAIH_11 HZ_ORAL_hzAIH_12 
    ##                1                1                1                1 
    ##  HZ_ORAL_hzAIH_2  HZ_ORAL_hzAIH_3  HZ_ORAL_hzAIH_4  HZ_ORAL_hzAIH_5 
    ##                1                1                1                1 
    ##  HZ_ORAL_hzAIH_6  HZ_ORAL_hzAIH_7  HZ_ORAL_hzAIH_8  HZ_ORAL_hzAIH_9 
    ##                1                1                1                1 
    ##    ZZ_ORAL_AIH_1   ZZ_ORAL_AIH_10   ZZ_ORAL_AIH_11   ZZ_ORAL_AIH_12 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_13   ZZ_ORAL_AIH_15   ZZ_ORAL_AIH_16   ZZ_ORAL_AIH_17 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_18   ZZ_ORAL_AIH_19    ZZ_ORAL_AIH_2   ZZ_ORAL_AIH_20 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_21   ZZ_ORAL_AIH_22   ZZ_ORAL_AIH_23   ZZ_ORAL_AIH_24 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_25   ZZ_ORAL_AIH_26   ZZ_ORAL_AIH_27   ZZ_ORAL_AIH_28 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_29    ZZ_ORAL_AIH_3   ZZ_ORAL_AIH_30   ZZ_ORAL_AIH_31 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_32   ZZ_ORAL_AIH_33   ZZ_ORAL_AIH_34   ZZ_ORAL_AIH_35 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_36   ZZ_ORAL_AIH_37   ZZ_ORAL_AIH_38   ZZ_ORAL_AIH_39 
    ##                1                1                1                1 
    ##    ZZ_ORAL_AIH_4   ZZ_ORAL_AIH_40   ZZ_ORAL_AIH_41   ZZ_ORAL_AIH_42 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_43   ZZ_ORAL_AIH_44   ZZ_ORAL_AIH_45   ZZ_ORAL_AIH_46 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_47   ZZ_ORAL_AIH_48   ZZ_ORAL_AIH_49    ZZ_ORAL_AIH_5 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_50   ZZ_ORAL_AIH_51   ZZ_ORAL_AIH_52   ZZ_ORAL_AIH_53 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_54   ZZ_ORAL_AIH_55   ZZ_ORAL_AIH_56   ZZ_ORAL_AIH_57 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_58   ZZ_ORAL_AIH_59    ZZ_ORAL_AIH_6   ZZ_ORAL_AIH_60 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_61   ZZ_ORAL_AIH_62   ZZ_ORAL_AIH_63   ZZ_ORAL_AIH_64 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_65   ZZ_ORAL_AIH_66   ZZ_ORAL_AIH_67   ZZ_ORAL_AIH_68 
    ##                1                1                1                1 
    ##   ZZ_ORAL_AIH_69    ZZ_ORAL_AIH_7    ZZ_ORAL_AIH_8    ZZ_ORAL_AIH_9 
    ##                1                1                1                1 
    ##     ZZ_ORAL_HC_1   ZZ_ORAL_HC_103   ZZ_ORAL_HC_104   ZZ_ORAL_HC_105 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_106   ZZ_ORAL_HC_108   ZZ_ORAL_HC_109   ZZ_ORAL_HC_110 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_111   ZZ_ORAL_HC_112   ZZ_ORAL_HC_114   ZZ_ORAL_HC_115 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_116   ZZ_ORAL_HC_117   ZZ_ORAL_HC_118   ZZ_ORAL_HC_119 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_120   ZZ_ORAL_HC_122   ZZ_ORAL_HC_126   ZZ_ORAL_HC_128 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_129    ZZ_ORAL_HC_13   ZZ_ORAL_HC_130   ZZ_ORAL_HC_131 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_134   ZZ_ORAL_HC_135   ZZ_ORAL_HC_136   ZZ_ORAL_HC_137 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_138   ZZ_ORAL_HC_139    ZZ_ORAL_HC_14   ZZ_ORAL_HC_140 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_142   ZZ_ORAL_HC_143   ZZ_ORAL_HC_144   ZZ_ORAL_HC_145 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_147   ZZ_ORAL_HC_148   ZZ_ORAL_HC_149   ZZ_ORAL_HC_150 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_152   ZZ_ORAL_HC_153   ZZ_ORAL_HC_154   ZZ_ORAL_HC_156 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_157   ZZ_ORAL_HC_159   ZZ_ORAL_HC_160   ZZ_ORAL_HC_161 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_162   ZZ_ORAL_HC_163   ZZ_ORAL_HC_164   ZZ_ORAL_HC_166 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_167   ZZ_ORAL_HC_168   ZZ_ORAL_HC_169    ZZ_ORAL_HC_17 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_171   ZZ_ORAL_HC_173   ZZ_ORAL_HC_174   ZZ_ORAL_HC_176 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_177   ZZ_ORAL_HC_179    ZZ_ORAL_HC_18   ZZ_ORAL_HC_180 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_181   ZZ_ORAL_HC_182   ZZ_ORAL_HC_183   ZZ_ORAL_HC_184 
    ##                1                1                1                1 
    ##   ZZ_ORAL_HC_185   ZZ_ORAL_HC_186    ZZ_ORAL_HC_19    ZZ_ORAL_HC_20 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_21    ZZ_ORAL_HC_22    ZZ_ORAL_HC_23    ZZ_ORAL_HC_24 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_25    ZZ_ORAL_HC_26    ZZ_ORAL_HC_27    ZZ_ORAL_HC_28 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_29     ZZ_ORAL_HC_3    ZZ_ORAL_HC_30    ZZ_ORAL_HC_32 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_33    ZZ_ORAL_HC_34    ZZ_ORAL_HC_35    ZZ_ORAL_HC_38 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_39     ZZ_ORAL_HC_4    ZZ_ORAL_HC_40    ZZ_ORAL_HC_41 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_42    ZZ_ORAL_HC_43    ZZ_ORAL_HC_45    ZZ_ORAL_HC_47 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_48    ZZ_ORAL_HC_49     ZZ_ORAL_HC_5    ZZ_ORAL_HC_50 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_53    ZZ_ORAL_HC_54    ZZ_ORAL_HC_55    ZZ_ORAL_HC_56 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_57    ZZ_ORAL_HC_58     ZZ_ORAL_HC_6    ZZ_ORAL_HC_60 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_62    ZZ_ORAL_HC_65    ZZ_ORAL_HC_66    ZZ_ORAL_HC_68 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_69    ZZ_ORAL_HC_70    ZZ_ORAL_HC_71    ZZ_ORAL_HC_72 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_73    ZZ_ORAL_HC_75    ZZ_ORAL_HC_78    ZZ_ORAL_HC_79 
    ##                1                1                1                1 
    ##     ZZ_ORAL_HC_8    ZZ_ORAL_HC_80    ZZ_ORAL_HC_82    ZZ_ORAL_HC_83 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_84    ZZ_ORAL_HC_85    ZZ_ORAL_HC_87    ZZ_ORAL_HC_88 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_90    ZZ_ORAL_HC_92    ZZ_ORAL_HC_93    ZZ_ORAL_HC_94 
    ##                1                1                1                1 
    ##    ZZ_ORAL_HC_95    ZZ_ORAL_HC_96    ZZ_ORAL_HC_97    ZZ_ORAL_HC_98 
    ##                1                1                1                1

``` r
# table(initial_meta$sequencing_meth..exp.) # output is pyrosequencing
table(initial_meta$Instrument) # output is Illumina. according to the methods, Illumina is correct.
```

    ## 
    ## Illumina MiSeq 
    ##            384

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

    ## Initial samples counts: 384

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

    ## Remaining samples after Step 2: 384

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
  sequencing_region = NA,
  platform = "Illumina MiSeq", # extract from methods section, NOT initial_meta
  country = initial_meta$geo_loc_name_country,
  site = initial_meta$isolation_source,
  site_group = "pending", # modify later
  public_time = initial_meta$ReleaseDate,
  pcr_primers = "F:5'-ACTCCTACGGGAGGCAGCA-3';R:5'-GGACTACHVGGGTWTCTAAT-3'", # Format: xxF-xxR
  title = "Oral microbiome characteristics in patients with autoimmune hepatitis", # extract from article
  disease_study = "yes", # determine according to the article
  raw_status = initial_meta$Sample.Name, 
  status = "pending", # modify later
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

metadata <- metadata %>%
  mutate(status = case_when(
    grepl("AIH", raw_status) ~ "autoimmune hepatitis", # 
    grepl("HC", raw_status) ~ "HC",
    TRUE ~ NA_character_
  )) # according to the 'Sample.Name' column of initial_meta
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
