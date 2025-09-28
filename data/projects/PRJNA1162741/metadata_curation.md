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
setwd("D:/studylearning/董彪师兄项目/任务/修改20250922/PRJNA1162741")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases   BioProject    BioSample
    ## 1 SRR30719636   AMPLICON        602 82330122 PRJNA1162741 SAMN43824388
    ## 2 SRR30719637   AMPLICON        602 72602404 PRJNA1162741 SAMN43824387
    ## 3 SRR30719638   AMPLICON        602 80581914 PRJNA1162741 SAMN43824386
    ## 4 SRR30719639   AMPLICON        602 90572706 PRJNA1162741 SAMN43824309
    ## 5 SRR30719640   AMPLICON        602 82547444 PRJNA1162741 SAMN43824308
    ## 6 SRR30719641   AMPLICON        602 93978822 PRJNA1162741 SAMN43824307
    ##                                BioSampleModel    Bytes
    ## 1 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 26574022
    ## 2 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 23521568
    ## 3 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 26105623
    ## 4 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 29458724
    ## 5 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 26546782
    ## 6 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 30597385
    ##                        Center.Name Collection_Date Consent DATASTORE.filetype
    ## 1 COMPLUTENSE UNIVERSITY OF MADRID      2022-05-19  public   fastq,run.zq,sra
    ## 2 COMPLUTENSE UNIVERSITY OF MADRID      2022-06-23  public   fastq,run.zq,sra
    ## 3 COMPLUTENSE UNIVERSITY OF MADRID      2021-05-31  public   fastq,run.zq,sra
    ## 4 COMPLUTENSE UNIVERSITY OF MADRID      2022-03-11  public   fastq,run.zq,sra
    ## 5 COMPLUTENSE UNIVERSITY OF MADRID      2022-05-26  public   fastq,run.zq,sra
    ## 6 COMPLUTENSE UNIVERSITY OF MADRID      2022-03-17  public   fastq,run.zq,sra
    ##   DATASTORE.provider                     DATASTORE.region env_broad_scale
    ## 1         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1         missing
    ## 2         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1         missing
    ## 3         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1         missing
    ## 4         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1         missing
    ## 5         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1         missing
    ## 6         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1         missing
    ##   env_local_scale env_medium  Experiment geo_loc_name_country
    ## 1           mouth     saliva SRX26123114                Spain
    ## 2           mouth     saliva SRX26123113                Spain
    ## 3           mouth     saliva SRX26123024                Spain
    ## 4           mouth     saliva SRX26123111                Spain
    ## 5           mouth     saliva SRX26123110                Spain
    ## 6           mouth     saliva SRX26123109                Spain
    ##   geo_loc_name_country_continent geo_loc_name         Host host_subject_id
    ## 1                         Europe        Spain Homo sapiens       sample289
    ## 2                         Europe        Spain Homo sapiens       sample288
    ## 3                         Europe        Spain Homo sapiens       sample287
    ## 4                         Europe        Spain Homo sapiens       sample210
    ## 5                         Europe        Spain Homo sapiens       sample209
    ## 6                         Europe        Spain Homo sapiens       sample208
    ##     Instrument        lat_lon Library.Name LibraryLayout LibrarySelection
    ## 1 NextSeq 2000 40.46 N 3.74 W  PBT-129_001        PAIRED              PCR
    ## 2 NextSeq 2000 40.46 N 3.74 W  PBT-128_001        PAIRED              PCR
    ## 3 NextSeq 2000 40.46 N 3.74 W  PBT-127_001        PAIRED              PCR
    ## 4 NextSeq 2000 40.46 N 3.74 W  PBT-028_001        PAIRED              PCR
    ## 5 NextSeq 2000 40.46 N 3.74 W  PBT-026_001        PAIRED              PCR
    ## 6 NextSeq 2000 40.46 N 3.74 W  PBT-025_001        PAIRED              PCR
    ##   LibrarySource              Organism Platform          ReleaseDate
    ## 1   METAGENOMIC human oral metagenome ILLUMINA 2024-10-08T00:00:00Z
    ## 2   METAGENOMIC human oral metagenome ILLUMINA 2024-10-08T00:00:00Z
    ## 3   METAGENOMIC human oral metagenome ILLUMINA 2024-10-08T00:00:00Z
    ## 4   METAGENOMIC human oral metagenome ILLUMINA 2024-10-08T00:00:00Z
    ## 5   METAGENOMIC human oral metagenome ILLUMINA 2024-10-08T00:00:00Z
    ## 6   METAGENOMIC human oral metagenome ILLUMINA 2024-10-08T00:00:00Z
    ##            create_date version
    ## 1 2024-09-19T09:15:00Z       1
    ## 2 2024-09-19T09:15:00Z       1
    ## 3 2024-09-19T09:15:00Z       1
    ## 4 2024-09-19T09:15:00Z       1
    ## 5 2024-09-19T09:15:00Z       1
    ## 6 2024-09-19T09:15:00Z       1
    ##                                                   samp_collect_device
    ## 1 DNA/RNA Shield Safe Collect Saliva Collection tubes (Zymo Research)
    ## 2 DNA/RNA Shield Safe Collect Saliva Collection tubes (Zymo Research)
    ## 3 DNA/RNA Shield Safe Collect Saliva Collection tubes (Zymo Research)
    ## 4 DNA/RNA Shield Safe Collect Saliva Collection tubes (Zymo Research)
    ## 5 DNA/RNA Shield Safe Collect Saliva Collection tubes (Zymo Research)
    ## 6 DNA/RNA Shield Safe Collect Saliva Collection tubes (Zymo Research)
    ##                                samp_mat_process Sample.Name SRA.Study Host_Age
    ## 1 ZymoBIOMICS DNA Miniprep Kit\\, Zymo Research       PT129 SRP533546       44
    ## 2 ZymoBIOMICS DNA Miniprep Kit\\, Zymo Research       PT128 SRP533546       51
    ## 3 ZymoBIOMICS DNA Miniprep Kit\\, Zymo Research       PT127 SRP533546       43
    ## 4 ZymoBIOMICS DNA Miniprep Kit\\, Zymo Research       PT028 SRP533546       23
    ## 5 ZymoBIOMICS DNA Miniprep Kit\\, Zymo Research       PT026 SRP533546       51
    ## 6 ZymoBIOMICS DNA Miniprep Kit\\, Zymo Research       PT025 SRP533546       30
    ##   host_body_mass_index host_sex
    ## 1                24.22   female
    ## 2                17.63   female
    ## 3                23.05   female
    ## 4                18.37   female
    ## 5                21.74     male
    ## 6                22.65   female

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

``` r
# Create a uniform metadata dataframe
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Sample.Name,
  Age = initial_meta$Host_Age,
  BMI = initial_meta$host_body_mass_index,
  Sex=  initial_meta$host_sex,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$env_medium,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- NA_character_  # Update based on the publication
metadata$title <- "Probing the oral-brain connection: oral microbiome patterns in a large community cohort with anxiety、depression and trauma symptoms and periodontal outcomes"
metadata$Disease_study <- "YES"
metadata <- metadata %>%
  mutate(Sex = case_when(
    grepl("f", Sex) ~ "F",
    grepl("m", Sex) ~ "M"
  ))
metadata <- metadata %>%
  mutate(Status = case_when(
    grepl("PT", sample_name) ~ "PTSD",
    TRUE ~ "HC"
  ))

#extract subject id by remove prefix of sample_name
#metadata$sample_name <- gsub("T|CR|CS|C|N|S", "", metadata$sample_name)

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##   P001   P002   P003   P004   P005   P006   P007   P008   P009   P011   P013 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P014   P015   P016   P017   P018   P019   P020   P022   P023   P024   P025 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P026   P027   P028   P029   P030   P031   P032   P033   P034   P035   P036 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P037   P038   P040   P041   P042   P043   P044   P045   P046   P047   P048 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P049   P050   P051   P052   P053   P054   P055   P056   P057   P059   P060 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P062   P063   P064   P065   P066   P068   P069   P070   P071   P072   P075 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P076   P077   P079   P080   P081   P082   P083   P086   P087   P088   P089 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P090   P091   P093   P094   P095   P096   P097   P098   P099   P100   P101 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P102   P103   P104   P105   P106   P107   P108   P109   P110   P111   P112 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P114   P115   P116   P117   P118   P119   P123   P125   P126   P127   P128 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P129   P130   P132   P133   P134   P135   P141   P142   P143   P145   P146 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P147   P150   P151   P152   P154   P155   P156   P158   P160   P162   P163 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P164   P165   P166   P167   P169   P170   P171   P172   P174   P175   P176 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P178   P179   P180   P181   P182   P183   P184   P185   P186   P187   P188 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P190   P191   P192   P195   P196   P197   P199   P200   P201   P202   P203 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P204   P205   P206   P207   P209   P210   P211   P214   P215   P216   P217 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P218   P219   P220   P221   P222   P223   P224   P225   P228   P229  PT001 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT002  PT003  PT004  PT005  PT007  PT008  PT010  PT011  PT012  PT013  PT014 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT015  PT016  PT017  PT019  PT020  PT021  PT022  PT023  PT024  PT025  PT026 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT028  PT029  PT031  PT032  PT033  PT034  PT035  PT038  PT039  PT040  PT042 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT043  PT044  PT045  PT046  PT047  PT048  PT049  PT050  PT051  PT053  PT054 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT055  PT056  PT058  PT059  PT060  PT061  PT062  PT063  PT064  PT065  PT066 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT067  PT068  PT069  PT070  PT071  PT072  PT073  PT075  PT076  PT077  PT078 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT079  PT080  PT081  PT083  PT084  PT086  PT088  PT089  PT090  PT092  PT095 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT096  PT097  PT098  PT099  PT100  PT102  PT103  PT105  PT107  PT108  PT109 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT111  PT112  PT114  PT116  PT117  PT119  PT120  PT123  PT124  PT125  PT126 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT127  PT128  PT129  PT130  PT131  PT135  PT136  PT138  PT140  PT141  PT143 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT145  PT146  PT147  PT148  PT149  PT150  PT152  PT153  PT154  PT156  PT158 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT160  PT163  PT164  PT165  PT166  PT168  PT169  PT170  PT171  PT172  PT174 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT175  PT177  PT178  PT180  PT182  PT185  PT186  PT187  PT190  PT192  PT194 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT195  PT196  PT197  PT198  PT199  PT200  PT201  PT202  PT203  PT205  PT206 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT207  PT209  PT211  PT212  PT213  PT215  PT217  PT218  PT220  PT224  PT226 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT227  PT228  PT234  PT235  PT236  PT238  PT239  PT240  PT241  PT242  PT243 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT244  PT246  PT247  PT248  PT249  PT250  PT251  PT252  PT254  PT255  PT256 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT257  PT258  PT259  PT261  PT262  PT263  PT264  PT265  PT266  PT269  PT270 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT273  PT275  PT277  PT278  PT279  PT280  PT282  PT283  PT284  PT285  PT286 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT287  PT288  PT289  PT290  PT293  PT294  PT295  PT297  PT298  PT299  PT301 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT303  PT304  PT305  PT306  PT307  PT309  PT311  PT312  PT313  PT315  PT316 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT320  PT321 PTB001 PTB002 PTB003 PTB004 PTB006 PTB007 PTB011 PTB012 PTB013 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTB014 PTB016 PTB018 PTB020 PTB021 PTB022 PTB023 PTB025 PTB027 PTB028 PTB029 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTO002 PTO003 PTO004 PTO005 PTO006 PTO009 PTO010 PTO011 PTO012 PTO013 PTO014 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTO016 PTO017 PTO018 PTO022 PTO023 PTO024 PTV001 PTV002 PTV003 PTV004 PTV005 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTV006 PTV007 PTV009 PTV014 PTV015 PTV017 PTV018 PTV020 PTV021 
    ##      1      1      1      1      1      1      1      1      1

``` r
table(metadata$Site)
```

    ## 
    ## saliva 
    ##    471

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
    ##   P001   P002   P003   P004   P005   P006   P007   P008   P009   P011   P013 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P014   P015   P016   P017   P018   P019   P020   P022   P023   P024   P025 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P026   P027   P028   P029   P030   P031   P032   P033   P034   P035   P036 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P037   P038   P040   P041   P042   P043   P044   P045   P046   P047   P048 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P049   P050   P051   P052   P053   P054   P055   P056   P057   P059   P060 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P062   P063   P064   P065   P066   P068   P069   P070   P071   P072   P075 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P076   P077   P079   P080   P081   P082   P083   P086   P087   P088   P089 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P090   P091   P093   P094   P095   P096   P097   P098   P099   P100   P101 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P102   P103   P104   P105   P106   P107   P108   P109   P110   P111   P112 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P114   P115   P116   P117   P118   P119   P123   P125   P126   P127   P128 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P129   P130   P132   P133   P134   P135   P141   P142   P143   P145   P146 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P147   P150   P151   P152   P154   P155   P156   P158   P160   P162   P163 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P164   P165   P166   P167   P169   P170   P171   P172   P174   P175   P176 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P178   P179   P180   P181   P182   P183   P184   P185   P186   P187   P188 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P190   P191   P192   P195   P196   P197   P199   P200   P201   P202   P203 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P204   P205   P206   P207   P209   P210   P211   P214   P215   P216   P217 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P218   P219   P220   P221   P222   P223   P224   P225   P228   P229  PT001 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT002  PT003  PT004  PT005  PT007  PT008  PT010  PT011  PT012  PT013  PT014 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT015  PT016  PT017  PT019  PT020  PT021  PT022  PT023  PT024  PT025  PT026 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT028  PT029  PT031  PT032  PT033  PT034  PT035  PT038  PT039  PT040  PT042 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT043  PT044  PT045  PT046  PT047  PT048  PT049  PT050  PT051  PT053  PT054 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT055  PT056  PT058  PT059  PT060  PT061  PT062  PT063  PT064  PT065  PT066 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT067  PT068  PT069  PT070  PT071  PT072  PT073  PT075  PT076  PT077  PT078 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT079  PT080  PT081  PT083  PT084  PT086  PT088  PT089  PT090  PT092  PT095 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT096  PT097  PT098  PT099  PT100  PT102  PT103  PT105  PT107  PT108  PT109 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT111  PT112  PT114  PT116  PT117  PT119  PT120  PT123  PT124  PT125  PT126 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT127  PT128  PT129  PT130  PT131  PT135  PT136  PT138  PT140  PT141  PT143 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT145  PT146  PT147  PT148  PT149  PT150  PT152  PT153  PT154  PT156  PT158 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT160  PT163  PT164  PT165  PT166  PT168  PT169  PT170  PT171  PT172  PT174 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT175  PT177  PT178  PT180  PT182  PT185  PT186  PT187  PT190  PT192  PT194 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT195  PT196  PT197  PT198  PT199  PT200  PT201  PT202  PT203  PT205  PT206 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT207  PT209  PT211  PT212  PT213  PT215  PT217  PT218  PT220  PT224  PT226 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT227  PT228  PT234  PT235  PT236  PT238  PT239  PT240  PT241  PT242  PT243 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT244  PT246  PT247  PT248  PT249  PT250  PT251  PT252  PT254  PT255  PT256 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT257  PT258  PT259  PT261  PT262  PT263  PT264  PT265  PT266  PT269  PT270 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT273  PT275  PT277  PT278  PT279  PT280  PT282  PT283  PT284  PT285  PT286 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT287  PT288  PT289  PT290  PT293  PT294  PT295  PT297  PT298  PT299  PT301 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT303  PT304  PT305  PT306  PT307  PT309  PT311  PT312  PT313  PT315  PT316 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT320  PT321 PTB001 PTB002 PTB003 PTB004 PTB006 PTB007 PTB011 PTB012 PTB013 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTB014 PTB016 PTB018 PTB020 PTB021 PTB022 PTB023 PTB025 PTB027 PTB028 PTB029 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTO002 PTO003 PTO004 PTO005 PTO006 PTO009 PTO010 PTO011 PTO012 PTO013 PTO014 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTO016 PTO017 PTO018 PTO022 PTO023 PTO024 PTV001 PTV002 PTV003 PTV004 PTV005 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTV006 PTV007 PTV009 PTV014 PTV015 PTV017 PTV018 PTV020 PTV021 
    ##      1      1      1      1      1      1      1      1      1

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
    ## 471

``` r
# Identify samples for removal: single samples excluding saliva
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & !(metadata$Site %in% c("saliva"))

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
    ##   P001   P002   P003   P004   P005   P006   P007   P008   P009   P011   P013 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P014   P015   P016   P017   P018   P019   P020   P022   P023   P024   P025 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P026   P027   P028   P029   P030   P031   P032   P033   P034   P035   P036 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P037   P038   P040   P041   P042   P043   P044   P045   P046   P047   P048 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P049   P050   P051   P052   P053   P054   P055   P056   P057   P059   P060 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P062   P063   P064   P065   P066   P068   P069   P070   P071   P072   P075 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P076   P077   P079   P080   P081   P082   P083   P086   P087   P088   P089 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P090   P091   P093   P094   P095   P096   P097   P098   P099   P100   P101 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P102   P103   P104   P105   P106   P107   P108   P109   P110   P111   P112 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P114   P115   P116   P117   P118   P119   P123   P125   P126   P127   P128 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P129   P130   P132   P133   P134   P135   P141   P142   P143   P145   P146 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P147   P150   P151   P152   P154   P155   P156   P158   P160   P162   P163 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P164   P165   P166   P167   P169   P170   P171   P172   P174   P175   P176 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P178   P179   P180   P181   P182   P183   P184   P185   P186   P187   P188 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P190   P191   P192   P195   P196   P197   P199   P200   P201   P202   P203 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P204   P205   P206   P207   P209   P210   P211   P214   P215   P216   P217 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##   P218   P219   P220   P221   P222   P223   P224   P225   P228   P229  PT001 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT002  PT003  PT004  PT005  PT007  PT008  PT010  PT011  PT012  PT013  PT014 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT015  PT016  PT017  PT019  PT020  PT021  PT022  PT023  PT024  PT025  PT026 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT028  PT029  PT031  PT032  PT033  PT034  PT035  PT038  PT039  PT040  PT042 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT043  PT044  PT045  PT046  PT047  PT048  PT049  PT050  PT051  PT053  PT054 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT055  PT056  PT058  PT059  PT060  PT061  PT062  PT063  PT064  PT065  PT066 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT067  PT068  PT069  PT070  PT071  PT072  PT073  PT075  PT076  PT077  PT078 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT079  PT080  PT081  PT083  PT084  PT086  PT088  PT089  PT090  PT092  PT095 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT096  PT097  PT098  PT099  PT100  PT102  PT103  PT105  PT107  PT108  PT109 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT111  PT112  PT114  PT116  PT117  PT119  PT120  PT123  PT124  PT125  PT126 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT127  PT128  PT129  PT130  PT131  PT135  PT136  PT138  PT140  PT141  PT143 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT145  PT146  PT147  PT148  PT149  PT150  PT152  PT153  PT154  PT156  PT158 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT160  PT163  PT164  PT165  PT166  PT168  PT169  PT170  PT171  PT172  PT174 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT175  PT177  PT178  PT180  PT182  PT185  PT186  PT187  PT190  PT192  PT194 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT195  PT196  PT197  PT198  PT199  PT200  PT201  PT202  PT203  PT205  PT206 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT207  PT209  PT211  PT212  PT213  PT215  PT217  PT218  PT220  PT224  PT226 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT227  PT228  PT234  PT235  PT236  PT238  PT239  PT240  PT241  PT242  PT243 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT244  PT246  PT247  PT248  PT249  PT250  PT251  PT252  PT254  PT255  PT256 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT257  PT258  PT259  PT261  PT262  PT263  PT264  PT265  PT266  PT269  PT270 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT273  PT275  PT277  PT278  PT279  PT280  PT282  PT283  PT284  PT285  PT286 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT287  PT288  PT289  PT290  PT293  PT294  PT295  PT297  PT298  PT299  PT301 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT303  PT304  PT305  PT306  PT307  PT309  PT311  PT312  PT313  PT315  PT316 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  PT320  PT321 PTB001 PTB002 PTB003 PTB004 PTB006 PTB007 PTB011 PTB012 PTB013 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTB014 PTB016 PTB018 PTB020 PTB021 PTB022 PTB023 PTB025 PTB027 PTB028 PTB029 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTO002 PTO003 PTO004 PTO005 PTO006 PTO009 PTO010 PTO011 PTO012 PTO013 PTO014 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTO016 PTO017 PTO018 PTO022 PTO023 PTO024 PTV001 PTV002 PTV003 PTV004 PTV005 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## PTV006 PTV007 PTV009 PTV014 PTV015 PTV017 PTV018 PTV020 PTV021 
    ##      1      1      1      1      1      1      1      1      1

``` r
table(metadata_clean$Site)
```

    ## 
    ## saliva 
    ##    471

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ##  NO 
    ## 471

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 471"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 471"

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
  "PCR_primers", "title", "Disease_study", "Status","matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##     BioProject         Run sample_name Age   BMI Sex Region     Platform
    ## 1 PRJNA1162741 SRR30719636       PT129  44 24.22   F  V3-V4 NextSeq 2000
    ## 2 PRJNA1162741 SRR30719637       PT128  51 17.63   F  V3-V4 NextSeq 2000
    ## 3 PRJNA1162741 SRR30719638       PT127  43 23.05   F  V3-V4 NextSeq 2000
    ## 4 PRJNA1162741 SRR30719639       PT028  23 18.37   F  V3-V4 NextSeq 2000
    ## 5 PRJNA1162741 SRR30719640       PT026  51 21.74   M  V3-V4 NextSeq 2000
    ## 6 PRJNA1162741 SRR30719641       PT025  30 22.65   F  V3-V4 NextSeq 2000
    ##   Country   Site          public_time PCR_primers
    ## 1   Spain saliva 2024-10-08T00:00:00Z        <NA>
    ## 2   Spain saliva 2024-10-08T00:00:00Z        <NA>
    ## 3   Spain saliva 2024-10-08T00:00:00Z        <NA>
    ## 4   Spain saliva 2024-10-08T00:00:00Z        <NA>
    ## 5   Spain saliva 2024-10-08T00:00:00Z        <NA>
    ## 6   Spain saliva 2024-10-08T00:00:00Z        <NA>
    ##                                                                                                                                                           title
    ## 1 Probing the oral-brain connection: oral microbiome patterns in a large community cohort with anxiety、depression and trauma symptoms and periodontal outcomes
    ## 2 Probing the oral-brain connection: oral microbiome patterns in a large community cohort with anxiety、depression and trauma symptoms and periodontal outcomes
    ## 3 Probing the oral-brain connection: oral microbiome patterns in a large community cohort with anxiety、depression and trauma symptoms and periodontal outcomes
    ## 4 Probing the oral-brain connection: oral microbiome patterns in a large community cohort with anxiety、depression and trauma symptoms and periodontal outcomes
    ## 5 Probing the oral-brain connection: oral microbiome patterns in a large community cohort with anxiety、depression and trauma symptoms and periodontal outcomes
    ## 6 Probing the oral-brain connection: oral microbiome patterns in a large community cohort with anxiety、depression and trauma symptoms and periodontal outcomes
    ##   Disease_study Status matched_samples
    ## 1           YES   PTSD              NO
    ## 2           YES   PTSD              NO
    ## 3           YES   PTSD              NO
    ## 4           YES   PTSD              NO
    ## 5           YES   PTSD              NO
    ## 6           YES   PTSD              NO

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
