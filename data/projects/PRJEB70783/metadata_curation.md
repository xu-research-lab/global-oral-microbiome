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
setwd("E:/bioproject/PRJEB70783")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = T) 

# Update the cancer_type and ibd_type columns
initial_meta <- initial_meta %>%
  mutate(across(c(cancer_type, ibd_type), ~ na_if(.x, "")))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases BioProject      BioSample
    ## 1 ERR12514043   AMPLICON        600 31759483 PRJEB70783 SAMEA115126056
    ## 2 ERR12514044   AMPLICON        601  6772967 PRJEB70783 SAMEA115126057
    ## 3 ERR12514045   AMPLICON        601   956368 PRJEB70783 SAMEA115126058
    ## 4 ERR12514047   AMPLICON        601 36057426 PRJEB70783 SAMEA115126060
    ## 5 ERR12514048   AMPLICON        601 21265578 PRJEB70783 SAMEA115126061
    ## 6 ERR12514050   AMPLICON        601 26443160 PRJEB70783 SAMEA115126063
    ##   broad.scale_environmental_context    Bytes          Center.Name
    ## 1                        microbiome 18817319 DALHOUSIE UNIVERSITY
    ## 2                        microbiome  4422348 DALHOUSIE UNIVERSITY
    ## 3                        microbiome   736218 DALHOUSIE UNIVERSITY
    ## 4                        microbiome 23198621 DALHOUSIE UNIVERSITY
    ## 5                        microbiome 14197068 DALHOUSIE UNIVERSITY
    ## 6                        microbiome 17848012 DALHOUSIE UNIVERSITY
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1       2009/2015  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2       2009/2015  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3       2009/2015  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4       2009/2015  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5       2009/2015  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6       2009/2015  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region declercq_etal_2021_medication
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1                             N
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1                             Y
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1                             Y
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1                             Y
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1                             Y
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1                             Y
    ##   declercq_etal_age_frailty declercq_wright_etal_ibd ENA.FIRST.PUBLIC..run.
    ## 1                         N                        N             2024-01-22
    ## 2                         Y                        N             2024-01-22
    ## 3                         Y                        N             2024-01-22
    ## 4                         Y                        N             2024-01-22
    ## 5                         Y                        Y             2024-01-22
    ## 6                         Y                        N             2024-01-22
    ##   ENA.LAST.UPDATE..run.     environmental_medium  Experiment
    ## 1            2024-01-22 microbiome in human oral ERX11889573
    ## 2            2024-01-22 microbiome in human oral ERX11889574
    ## 3            2024-01-22 microbiome in human oral ERX11889575
    ## 4            2024-01-22 microbiome in human oral ERX11889577
    ## 5            2024-01-22 microbiome in human oral ERX11889578
    ## 6            2024-01-22 microbiome in human oral ERX11889580
    ##   geo_loc_name_country geo_loc_name_country_continent
    ## 1               Canada                  North America
    ## 2               Canada                  North America
    ## 3               Canada                  North America
    ## 4               Canada                  North America
    ## 5               Canada                  North America
    ## 6               Canada                  North America
    ##   geographic_location_.country_and.or_sea. geographic_location_.latitude.
    ## 1                                   Canada                       44.61455
    ## 2                                   Canada                       44.61455
    ## 3                                   Canada                       44.61455
    ## 4                                   Canada                       44.61455
    ## 5                                   Canada                       44.61455
    ## 6                                   Canada                       44.61455
    ##   geographic_location_.longitude.     Instrument Library.Name LibraryLayout
    ## 1                       -63.55807 Illumina MiSeq    100000735        PAIRED
    ## 2                       -63.55807 Illumina MiSeq    100002335        PAIRED
    ## 3                       -63.55807 Illumina MiSeq    100002535        PAIRED
    ## 4                       -63.55807 Illumina MiSeq    100003235        PAIRED
    ## 5                       -63.55807 Illumina MiSeq    100003935        PAIRED
    ## 6                       -63.55807 Illumina MiSeq    100004535        PAIRED
    ##   LibrarySelection LibrarySource local_environmental_context
    ## 1              PCR   METAGENOMIC            human microbiome
    ## 2              PCR   METAGENOMIC            human microbiome
    ## 3              PCR   METAGENOMIC            human microbiome
    ## 4              PCR   METAGENOMIC            human microbiome
    ## 5              PCR   METAGENOMIC            human microbiome
    ## 6              PCR   METAGENOMIC            human microbiome
    ##   nearing_etal_2020_healthy_oral_microbiome nearing_etal_2023_cancer
    ## 1                                         N                        Y
    ## 2                                         Y                        N
    ## 3                                         Y                        N
    ## 4                                         Y                        N
    ## 5                                         Y                        N
    ## 6                                         Y                        N
    ##                  Organism Platform      primers project_name
    ## 1 human saliva metagenome ILLUMINA 515FBï¿½926R   PRJEB70783
    ## 2 human saliva metagenome ILLUMINA 515FBï¿½926R   PRJEB70783
    ## 3 human saliva metagenome ILLUMINA 515FBï¿½926R   PRJEB70783
    ## 4 human saliva metagenome ILLUMINA 515FBï¿½926R   PRJEB70783
    ## 5 human saliva metagenome ILLUMINA   515FB?926R   PRJEB70783
    ## 6 human saliva metagenome ILLUMINA 515FBï¿½926R   PRJEB70783
    ##            ReleaseDate          create_date version    Sample.Name
    ## 1 2024-01-23T00:00:00Z 2024-01-25T10:29:00Z       1      100000735
    ## 2 2024-01-23T00:00:00Z 2024-01-26T22:43:00Z       1      100002335
    ## 3 2024-01-23T00:00:00Z 2024-01-26T22:41:00Z       1      100002535
    ## 4 2024-01-23T00:00:00Z 2024-01-26T17:23:00Z       1      100003235
    ## 5 2024-01-23T00:00:00Z 2024-01-26T22:43:00Z       1 SAMEA115126061
    ## 6 2024-01-23T00:00:00Z 2024-01-26T17:21:00Z       1      100004535
    ##           Scientific_Name Sequencing_method SRA.Study                organism
    ## 1 human saliva metagenome             MiSeq ERP155686 human saliva metagenome
    ## 2 human saliva metagenome             MiSeq ERP155686 human saliva metagenome
    ## 3 human saliva metagenome             MiSeq ERP155686 human saliva metagenome
    ## 4 human saliva metagenome             MiSeq ERP155686 human saliva metagenome
    ## 5 human saliva metagenome             MiSeq ERP155686                        
    ## 6 human saliva metagenome             MiSeq ERP155686 human saliva metagenome
    ##   cancer_case_control cancer_retrospective_prospective cancer_type
    ## 1             Control                      prospective      Breast
    ## 2                                                             <NA>
    ## 3                                                             <NA>
    ## 4                                                             <NA>
    ## 5                                                             <NA>
    ## 6                                                             <NA>
    ##   IBD_case_control ibd_type     ena_first_public      ena_last_update
    ## 1                      <NA>                                          
    ## 2                      <NA>                                          
    ## 3                      <NA>                                          
    ## 4                      <NA>                                          
    ## 5             Case   Crohns 2024-01-22T20:33:36Z 2024-01-22T20:33:36Z
    ## 6                      <NA>                                          
    ##      external_id    INSDC_center_name   INSDC_first_public    INSDC_last_update
    ## 1                                                                              
    ## 2                                                                              
    ## 3                                                                              
    ## 4                                                                              
    ## 5 SAMEA115126061 Dalhousie University 2024-01-22T20:33:36Z 2024-01-22T20:33:36Z
    ## 6                                                                              
    ##   INSDC_status Sample_name Submitter_Id
    ## 1                                      
    ## 2                                      
    ## 3                                      
    ## 4                                      
    ## 5       public   100003935    100003935
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
  sample_name = initial_meta$Sample.Name,
  Age = "35-70", 
  BMI = NA, 
  Sex = NA, 
  Region = initial_meta$Organism,
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$Organism,
  public_time = initial_meta$ReleaseDate,
  cancer_type = initial_meta$cancer_type,
  ibd_type = initial_meta$ibd_type,
  matched_samples = NA,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "515F-926R"
metadata$title <- "Oral microbial signatures associated with age and frailty in Canadian adults"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata <- metadata %>%
  mutate(
    Status = case_when(
      cancer_type == "Breast" ~ "Breast cancer",
      cancer_type == "Colon" ~ "Colon cancer",
      cancer_type == "Prostate" ~ "Prostate cancer",
      ibd_type == "Crohns" ~ "Crohns disease",
      ibd_type == "UC"~ "Ulcerative colitis",
      ibd_type == "Crohns_UC" ~ "Crohns disease_Ul",
      TRUE ~ "H")
    )

#extract subject id by remove prefix of sample_name
metadata$sample_name <- gsub("10", "", metadata$sample_name)

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##        0000135        0000735        0002335        0002535        0002935 
    ##              1              1              1              1              1 
    ##        0003235        0004035        0004535        0005135        0005235 
    ##              1              1              1              1              1 
    ##        0005435        0005735        0006235        0006935        0007135 
    ##              1              1              1              1              1 
    ##        0007635        0008335        0009635        0012935        0013735 
    ##              1              1              1              1              1 
    ##        0014135        0014235        0014435        0014635        0014935 
    ##              1              1              1              1              1 
    ##        0015235        0015435        0015535        0015835        0016235 
    ##              1              1              1              1              1 
    ##        0016435        0016535        0017835        0017935        0018635 
    ##              1              1              1              1              1 
    ##        0019535      0020535-2        0020735        0020935        0022235 
    ##              1              1              1              1              1 
    ##        0022335        0023235        0023535        0023935      0025335-2 
    ##              1              1              1              1              1 
    ##        0025435      0025835-2        0027035        0027235        0027335 
    ##              1              1              1              1              1 
    ##      0027535-2        0027635      0028435-2        0029235        0029435 
    ##              1              1              1              1              1 
    ##        0029935        0030635        0032635        0032835        0032935 
    ##              1              1              1              1              1 
    ##        0034635        0035535        0038235        0039235        0045735 
    ##              1              1              1              1              1 
    ##        0046235        0046335        0046635        0047235        0047835 
    ##              1              1              1              1              1 
    ##        0048035        0048335        0048635        0049835        0050435 
    ##              1              1              1              1              1 
    ##        0050535        0050735        0051235        0051835        0053535 
    ##              1              1              1              1              1 
    ##        0053635        0054435        0055235        0056135        0056435 
    ##              1              1              1              1              1 
    ##        0057235        0057535        0058435        0058535        0058935 
    ##              1              1              1              1              1 
    ##        0059035        0059635        0060135        0060235        0062535 
    ##              1              1              1              1              1 
    ##        0063035      0063135-2        0063935        0064135        0064635 
    ##              1              1              1              1              1 
    ##        0065835        0066235        0067235        0068635        0070435 
    ##              1              1              1              1              1 
    ##        0072435        0072735      0073835-2        0073935        0074535 
    ##              1              1              1              1              1 
    ##        0074735        0074935        0076135        0076935        0078035 
    ##              1              1              1              1              1 
    ##        0078835        0079135        0079835        0080235        0080835 
    ##              1              1              1              1              1 
    ##        0082035        0084635        0085135        0085435        0089135 
    ##              1              1              1              1              1 
    ##        0089435        0090035        0090235        0092235        0092535 
    ##              1              1              1              1              1 
    ##        0095035        0097335        0098835        0099635        0113135 
    ##              1              1              1              1              1 
    ##        0113435        0113535        0113735        0115135        0117135 
    ##              1              1              1              1              1 
    ##        0117235        0118535        0119335        0120135        0120435 
    ##              1              1              1              1              1 
    ##        0121535        0122735        0123435        0123535        0124135 
    ##              1              1              1              1              1 
    ##        0126435        0126635        0127135        0127235        0128835 
    ##              1              1              1              1              1 
    ##        0136635        0136935        0137735        0138335        0139435 
    ##              1              1              1              1              1 
    ##        0139735        0140935        0141235        0141435        0141935 
    ##              1              1              1              1              1 
    ##        0142635        0144435        0146735        0146935        0147035 
    ##              1              1              1              1              1 
    ##        0147535        0147835        0151935        0153535        0153935 
    ##              1              1              1              1              1 
    ##        0154835        0156735        0157235        0160535      0160935-2 
    ##              1              1              1              1              1 
    ##        0161635        0161935        0162035        0163935        0164035 
    ##              1              1              1              1              1 
    ##        0164235        0164435        0164635        0165535        0165935 
    ##              1              1              1              1              1 
    ##        0166135      0166235-2        0166335        0166835        0167935 
    ##              1              1              1              1              1 
    ##        0168235        0168835        0169835        0170135        0170235 
    ##              1              1              1              1              1 
    ##        0171335      0178935-2        0182835        0182935      0186835-2 
    ##              1              1              1              1              1 
    ##        0192635        0194535        0195835        0196135        0200835 
    ##              1              1              1              1              1 
    ##        0201235        0202135        0203535        0203835        0206535 
    ##              1              1              1              1              1 
    ##        0220135        0220235        0220835        0222835        0223535 
    ##              1              1              1              1              1 
    ##        0225635        0228235      0230735-2        0231935        0233435 
    ##              1              1              1              1              1 
    ##        0234435        0236935        0238235        0240135        0241235 
    ##              1              1              1              1              1 
    ##        0241735        0241935        0242035        0242635          02435 
    ##              1              1              1              1              1 
    ##        0243935        0245135        0245635        0253435          02535 
    ##              1              1              1              1              1 
    ##        0253835        0254335        0254835        0255235        0260835 
    ##              1              1              1              1              1 
    ##        0261735        0262335        0263235        0264035        0266335 
    ##              1              1              1              1              1 
    ##        0267935        0268035        0268235        0268735        0269835 
    ##              1              1              1              1              1 
    ##        0272935        0274135        0277935        0279535        0291735 
    ##              1              1              1              1              1 
    ##        0299535        0299835        0300135        0303935        0306235 
    ##              1              1              1              1              1 
    ##        0306635        0306735        0306835        0308035        0309235 
    ##              1              1              1              1              1 
    ##        0315735        0316935        0318135        0318435        0323935 
    ##              1              1              1              1              1 
    ##        0325435        0329635        0332235        0332735        0334635 
    ##              1              1              1              1              1 
    ##        0335635        0336235        0336435        0336935        0339235 
    ##              1              1              1              1              1 
    ##        0341835        0342535        0342635        0343135        0343835 
    ##              1              1              1              1              1 
    ##        0344835        0347935        0348235        0350635        0352535 
    ##              1              1              1              1              1 
    ##        0352935          03535        0354235        0356135        0358035 
    ##              1              1              1              1              1 
    ##        0360335        0364935        0366535        0367035        0367135 
    ##              1              1              1              1              1 
    ##        0367935        0371235        0371535        0378535        0379335 
    ##              1              1              1              1              1 
    ##        0379635        0380935        0381735        0382135        0383335 
    ##              1              1              1              1              1 
    ##        0383835        0383935        0385235        0386035        0386135 
    ##              1              1              1              1              1 
    ##        0387935        0389235        0389335        0391335        0391435 
    ##              1              1              1              1              1 
    ##        0391735        0392335        0393435          03935        0394335 
    ##              1              1              1              1              1 
    ##        0395535        0396135        0397735        0398735        0399135 
    ##              1              1              1              1              1 
    ##        0399835        0402635        0403035        0403635        0406135 
    ##              1              1              1              1              1 
    ##        0406235        0406335        0408335        0408735        0409435 
    ##              1              1              1              1              1 
    ##          04135        0414335        0414835        0417735        0417935 
    ##              1              1              1              1              1 
    ##        0419335        0422835        0425235        0426535        0427735 
    ##              1              1              1              1              1 
    ##        0430335        0431135        0432235        0432735        0432935 
    ##              1              1              1              1              1 
    ##        0434735        0439235        0439635        0440135        0443735 
    ##              1              1              1              1              1 
    ##        0443935        0444935        0446035        0446135        0446235 
    ##              1              1              1              1              1 
    ##        0448135        0448735        0449335        0450035        0450935 
    ##              1              1              1              1              1 
    ##        0451235        0451335        0452035        0453335          04535 
    ##              1              1              1              1              1 
    ##        0453535        0453635        0454035        0456535        0460035 
    ##              1              1              1              1              1 
    ##        0460135        0461535        0463335        0463535        0464235 
    ##              1              1              1              1              1 
    ##        0465335        0466235        0467435        0467635        0467835 
    ##              1              1              1              1              1 
    ##        0468035        0468335        0468835        0469635        0470835 
    ##              1              1              1              1              1 
    ##        0470935        0471135        0473135        0473335        0473635 
    ##              1              1              1              1              1 
    ##        0475435        0475535        0476535        0476835      0477535-2 
    ##              1              1              1              1              1 
    ##      0477835-2        0479535        0479635        0480135        0480535 
    ##              1              1              1              1              1 
    ##      0480735-2        0481535        0481935        0483435        0483835 
    ##              1              1              1              1              1 
    ##        0485035        0485635        0485835        0486235      0489435-2 
    ##              1              1              1              1              1 
    ##      0490435-2        0490535        0491335        0492135        0494035 
    ##              1              1              1              1              1 
    ##        0494435        0494735        0494835        0496135        0496735 
    ##              1              1              1              1              1 
    ##        0496935        0500835        0501135        0501835        0502535 
    ##              1              1              1              1              1 
    ##        0503635        0505735        0506235        0507835        0509835 
    ##              1              1              1              1              1 
    ##        0511135        0512035      0513035-2        0513335        0515135 
    ##              1              1              1              1              1 
    ##        0517335        0518835        0519235        0524435        0525835 
    ##              1              1              1              1              1 
    ##        0527635        0528435        0528735        0531335          05435 
    ##              1              1              1              1              1 
    ##        0561135        0563435        0563735        0564135        0565235 
    ##              1              1              1              1              1 
    ##        0565435        0566335        0566435        0566535        0566835 
    ##              1              1              1              1              1 
    ##        0567435        0567835        0569235        0570035      0570635-2 
    ##              1              1              1              1              1 
    ##        0572435        0572835        0573735        0575835        0576735 
    ##              1              1              1              1              1 
    ##        0577235        0577335        0578435        0579635        0579935 
    ##              1              1              1              1              1 
    ##        0580135        0580735        0580835        0582235        0582335 
    ##              1              1              1              1              1 
    ##        0582435          05835        0583635        0583735        0584035 
    ##              1              1              1              1              1 
    ##        0584435        0584735        0585635        0586735        0587635 
    ##              1              1              1              1              1 
    ##        0587835      0588635-2        0588735        0589235        0589435 
    ##              1              1              1              1              1 
    ##        0589535        0589735        0591235        0591335        0592235 
    ##              1              1              1              1              1 
    ##        0593235        0594035        0594135        0595535        0595835 
    ##              1              1              1              1              1 
    ##        0595935        0597835        0597935        0598035        0598935 
    ##              1              1              1              1              1 
    ##        0599335        0599435        0599535        0599835        0601135 
    ##              1              1              1              1              1 
    ##        0602835        0603435        0603635        0604235        0604635 
    ##              1              1              1              1              1 
    ##        0604935        0605535        0607035        0607935        0608135 
    ##              1              1              1              1              1 
    ##        0608235        0609535        0631635        0636135        0637835 
    ##              1              1              1              1              1 
    ##        0638035        0638135        0638735        0638835        0638935 
    ##              1              1              1              1              1 
    ##        0639535        0639735        0640035        0641235        0641435 
    ##              1              1              1              1              1 
    ##          06735        0680835        0681435        0681935        0682435 
    ##              1              1              1              1              1 
    ##          06835        0683735        0685035        0685135        0685835 
    ##              2              1              1              1              1 
    ##        0687935        0691935      0692735-2     0693835-2W        0695335 
    ##              1              1              1              1              1 
    ##        0695735        0695935        0696535        0696835        0697835 
    ##              1              1              1              1              1 
    ##        0700535        0702835        0703535        0704235      0705635-2 
    ##              1              1              1              1              1 
    ##        0705935        0706235      0706535-2        0707535        0707735 
    ##              1              1              1              1              1 
    ##        0708935      0711435-2        0712235        0713435          07135 
    ##              1              1              1              1              1 
    ##        0713735        0713835      0715135-2        0715935        0718935 
    ##              1              1              1              1              1 
    ##          07235        0723535        0729435        0730535        0730935 
    ##              1              1              1              1              1 
    ##        0731935        0732135        0732335        0733235        0733935 
    ##              1              1              1              1              1 
    ##        0734235        0735035        0735435        0735535        0736535 
    ##              1              1              1              1              1 
    ##        0740635        0742635        0743335          07435        0743635 
    ##              1              1              1              1              1 
    ##        0743835        0744235        0744735        0745035        0746535 
    ##              1              1              1              1              1 
    ##        0747035        0749535        0749835        0752035        0752535 
    ##              1              1              1              1              1 
    ##        0752735          07535        0753935        0755435        0755635 
    ##              1              1              1              1              1 
    ##        0756235        0756735        0757835        0759735        0759835 
    ##              1              1              1              1              1 
    ##        0763635        0765935        0766835        0768635        0769235 
    ##              1              1              1              1              1 
    ##        0769535        0769735        0771135        0771835        0774835 
    ##              1              1              1              1              1 
    ##        0774935        0775035        0776335        0776435     0778235-2W 
    ##              1              1              1              1              1 
    ##        0778735        0787735        0790335        0790935        0791135 
    ##              1              1              1              1              1 
    ##        0792035        0792935        0793235        0793335     0793835-2W 
    ##              1              1              1              1              1 
    ##        0796535        0797735        0797835        0798935        0800735 
    ##              1              1              1              1              1 
    ##        0800835        0801335        0802435        0804035        0805135 
    ##              1              1              1              1              1 
    ##        0807235        0807935        0809735        0811935          08135 
    ##              1              1              1              1              1 
    ##        0814335        0815535        0816935        0818835        0819135 
    ##              1              1              1              1              1 
    ##        0825635        0829735        0830535        0831735        0834435 
    ##              1              1              1              1              1 
    ##      0835235-2        0835935        0837235        0837735        0837835 
    ##              1              1              1              1              1 
    ##        0838435        0840435        0841835        0841935        0842035 
    ##              1              1              1              1              1 
    ##        0842335        0842835          08435        0843635        0844335 
    ##              1              1              1              1              1 
    ##        0844435        0847435        0850935        0852035        0852535 
    ##              1              1              1              1              1 
    ##        0852835        0853235        0853335        0853735        0854535 
    ##              1              1              1              1              1 
    ##        0855035        0859135        0860735        0868235        0872435 
    ##              1              1              1              1              1 
    ##          08735        0873635        0873735        0875835        0876135 
    ##              1              1              1              1              1 
    ##        0876535        0879535        0881235        0883535        0883835 
    ##              1              1              1              1              1 
    ##        0885335      0886035-2        0886435        0887035        0887135 
    ##              1              1              1              1              1 
    ##        0887835        0887935        0888735        0889235        0889435 
    ##              1              1              1              1              1 
    ##      0891235-2      0893735-2        0893935        0945035        0946235 
    ##              1              1              1              1              1 
    ##        0948735        0960035        0960235        0960735        0961235 
    ##              1              1              1              1              1 
    ##        0963435          09635        0964035        0964235        0964735 
    ##              1              1              1              1              1 
    ##        0967735        0969835        0973335      0974035-2        0975235 
    ##              1              1              1              1              1 
    ##        0976935        0986635        0986735        0988435        0989235 
    ##              1              1              1              1              1 
    ##        0990935        0991335        0992835        0994135        0995935 
    ##              1              1              1              1              1 
    ##        0996135        0996935        0997135        0997235        0998135 
    ##              1              1              1              1              1 
    ##        0998535        0998935        0999035        10135-2          10235 
    ##              1              1              1              1              1 
    ##          11035        1111135        1112435        1112935        1113335 
    ##              1              1              1              1              1 
    ##          11135        1114335        1114735        1119135        1119235 
    ##              2              1              1              1              1 
    ##        1120235      1120335-2        1120935        1121335        1121435 
    ##              1              1              1              1              1 
    ##        1122935        1123535        1159335        1159835        1160835 
    ##              1              1              1              1              1 
    ##        1161435        1164235        1165335        1166435        1170635 
    ##              1              1              1              1              1 
    ##        1170835        1173135          11735        12135-2      1233635-2 
    ##              1              1              1              1              1 
    ##        1234335        1234735        1237135        1237535        1237935 
    ##              1              1              1              1              1 
    ##        1240935        1241835        1241935        1242635        1243035 
    ##              1              1              1              1              1 
    ##        1247035        1248435        1249635          12535        1256335 
    ##              1              1              1              1              1 
    ##        1260935        1262335        1262435        1264135        1265035 
    ##              1              1              1              1              1 
    ##        1265235        1266635        1268435        1268535        1274335 
    ##              1              1              1              1              1 
    ##        1275435        1278235        1280735        1281235        1282235 
    ##              1              1              1              1              1 
    ##        1282435        1285035        1286435        1288335        1290735 
    ##              1              1              1              1              1 
    ##        1291935      1294535-2        1295335      1295735-2      1296535-2 
    ##              1              1              1              1              1 
    ##      1296935-2        1302735        1303635        1304735        1321135 
    ##              1              1              1              1              1 
    ##        1324935        1334035        1336635        1341635        1342635 
    ##              1              1              1              1              1 
    ##        1344035        1349335        1350135        1350435        1352235 
    ##              1              1              1              1              1 
    ##        1352835        1354835        1355135        1357035        1357635 
    ##              1              1              1              1              1 
    ##        1359835        1361835        1369035        1370535        1371735 
    ##              1              1              1              1              1 
    ##        1372235          13735        1379435        1379535      1380335-2 
    ##              1              1              1              1              1 
    ##        1381135        1381335      1381835-2        1382735        1383035 
    ##              1              1              1              1              1 
    ##      1383435-2        1387135        1391135        1399135        1399335 
    ##              1              1              1              1              1 
    ##      1399835-2        1401635        1401735        1401935      1402135-2 
    ##              1              1              1              1              1 
    ##        1402335          14035          14135          14335        1493835 
    ##              1              1              1              1              1 
    ##        1494035        1494135        1497235        1504635        1504735 
    ##              1              1              1              1              1 
    ##        1505335        1505535        1505835        1506335        1506535 
    ##              1              1              1              1              1 
    ##        1506635        1506935        1507435        1508635        1509735 
    ##              1              1              1              1              1 
    ##        1511135        1511935        1512435        1513135        1514335 
    ##              1              1              1              1              1 
    ##        1514435        1515335        1515535        1516435        1517035 
    ##              1              1              1              1              1 
    ##        1517735        1517835        1518435        1520235        1520935 
    ##              1              1              1              1              1 
    ##        1521135        1522435          15235        1524335        1524535 
    ##              1              1              1              1              1 
    ##        1525235        1527635        1529335        1530735        1531635 
    ##              1              1              1              1              1 
    ##        1532835        1533335          15335        1535435        1535835 
    ##              1              1              1              1              1 
    ##        1536335        1536835        1537335        1538235        1538435 
    ##              1              1              1              1              1 
    ##        1546235        1547535        1547935        1549035        1551135 
    ##              1              1              1              1              1 
    ##        1555735        1558435        1558635        1559335        1560035 
    ##              1              1              1              1              1 
    ##        1560135      1561635-2        1564635        1564935        1567135 
    ##              1              1              1              1              1 
    ##        1567335        1570335        1570535        1573235          15735 
    ##              1              1              1              1              2 
    ##        1573735        1574835        1576435        1576535        1576835 
    ##              1              1              1              1              1 
    ##        1579335        1581135          15835        1584435        1588135 
    ##              1              1              1              1              1 
    ##        1591535        1593235          15935        1593635        1594035 
    ##              1              1              1              1              1 
    ##        1634535        1635835        1637135        1637635        1638735 
    ##              1              1              1              1              1 
    ##        1639235        1640335        1642135        1643035        1643835 
    ##              1              1              1              1              1 
    ##        1644535        1645235        1648035        1669135        1670335 
    ##              1              1              1              1              1 
    ##        1670535        1672435        1673735        1673835        1691435 
    ##              1              1              1              1              1 
    ##        1691535        1691935        1692735        1699535        1700635 
    ##              1              1              1              1              1 
    ##        1700835        1701735      1702135-2          17035        1711635 
    ##              1              1              1              1              1 
    ##        1712135        1712335        1712935          17135        1713535 
    ##              1              1              1              1              1 
    ##        1714735        1716135          17535        1779732        1779932 
    ##              1              1              1              1              1 
    ##        1780132        1780832        1781632        1781832        1781932 
    ##              1              1              1              1              1 
    ##        1782232        1783232        1784132        1784432          17935 
    ##              1              1              1              1              1 
    ##        1859435        1859635        1861635        1862335        1862535 
    ##              1              1              1              1              1 
    ##          18635        1865235        1865535        1866735        1867135 
    ##              1              1              1              1              1 
    ##        1867335        1867935        1868035          18735      1875235-2 
    ##              1              1              1              1              1 
    ##        1877435      1877935-2        1880035        1881135        1882235 
    ##              1              1              1              1              1 
    ##        1882435        1882635        1883835      1884435-2        1884835 
    ##              1              1              1              1              1 
    ##        1887635        1887735        1889035        1889335        1890235 
    ##              1              1              1              1              1 
    ##        1891735        1894235        1894635        1895135        1897335 
    ##              1              1              1              1              1 
    ##        1898335        1899035        1899235        1901335        1901735 
    ##              1              1              1              1              1 
    ##          19035        1904035        1904935        1905035        1906335 
    ##              1              1              1              1              1 
    ##        1906435        1907035        1909635        1912335        1912635 
    ##              1              1              1              1              1 
    ##        1912935        1913135        1913735        1914535        1914835 
    ##              1              1              1              1              1 
    ##        1915235        1916435        1917635        1918535        1919135 
    ##              1              1              1              1              1 
    ##        1920135        1920335      1921235-2        1921535      1921835-2 
    ##              1              1              1              1              1 
    ##     1922035-2W      1922135-2        1922935      1923635-2      1923935-2 
    ##              1              1              1              1              1 
    ##        1924335        1926135        1928235        1929835        1932035 
    ##              1              1              1              1              1 
    ##        1932635        1932835        1933035        1933635        1933835 
    ##              1              1              1              1              1 
    ##        1934035        1935535        1938235      1939335-2        1943235 
    ##              1              1              1              1              1 
    ##          19435        1944035        1944435        1946435      1947335-2 
    ##              2              1              1              1              1 
    ##        1947835        1948335        1948835        1949035        1952635 
    ##              1              1              1              1              1 
    ##      1953135-2          19535        1956235        1960235        1962635 
    ##              1              1              1              1              1 
    ##        1963435          19635        1970935        1972635        1972835 
    ##              1              1              1              1              1 
    ##        1973435        1974135        1974335        1974535        1977035 
    ##              1              1              1              1              1 
    ##        1977935        1980435        1981135        1981335        1981835 
    ##              1              1              1              1              1 
    ##        1982235        1983635        1984635        1984835        1985635 
    ##              1              1              1              1              1 
    ##        1985935        1986035        1986235        1987435        1987635 
    ##              1              1              1              1              1 
    ##        1988535        1989335        1991735        1993035          19935 
    ##              1              1              1              1              1 
    ##        1994635        1995035        1995835        1997335        1997435 
    ##              1              1              1              1              1 
    ##        1997735        1998235        1999835        2000935        2001935 
    ##              1              1              1              1              1 
    ##        2003035          20035        2003635        2003835        2003935 
    ##              1              1              1              1              1 
    ##      2004335-2        2005035        2006335        2006435        2008035 
    ##              1              1              1              1              1 
    ##          20135        2013935        2014135        2014935        2015135 
    ##              1              1              1              1              1 
    ##        2015335        2015435        2016735        2017135        2017735 
    ##              1              1              1              1              1 
    ##        2018035        2019835        2020235        2020335        2021435 
    ##              1              1              1              1              1 
    ##        2023035          20235        2025635        2026035        2026235 
    ##              1              1              1              1              1 
    ##        2026535        2027535        2029335        2029635        2030035 
    ##              1              1              1              1              1 
    ##        2030335        2031335        2031735        2032335        2032935 
    ##              1              1              1              1              1 
    ##        2033035        2033235        2033435        20335-2        2036335 
    ##              1              1              1              1              1 
    ##        2037035        2037535        2038735        2039935        2040535 
    ##              1              1              1              1              1 
    ##        2040635        2040735        2042935        2046735        2049535 
    ##              1              1              1              1              1 
    ##        2050035        2050735        2050835          20535        2057035 
    ##              1              1              1              2              1 
    ##        2058235        2063535        2064035        2064735        2066435 
    ##              1              1              1              1              1 
    ##        2066935        2068435        2068535        2068935        2069735 
    ##              1              1              1              1              1 
    ##        2069835        2070635        2076335      2076635-2      2077035-2 
    ##              1              1              1              1              1 
    ##      2077235-2        2078235        2078635        2081235        2082235 
    ##              1              1              1              1              1 
    ##        2084235        2084735        2087935        2088235        2088435 
    ##              1              1              1              1              1 
    ##        2088735        2088835        2089235        2089935        2090035 
    ##              1              1              1              1              1 
    ##        2090635        2090935        2091235      2091335-2        2091635 
    ##              1              1              1              1              1 
    ##        2096335        2097035        2097335        2099635        2099735 
    ##              1              1              1              1              1 
    ##        2111135      2112535-2        2112835        2113235        2115335 
    ##              1              1              1              1              1 
    ##        2116435        2117235        2118735        2119135        2119335 
    ##              1              1              1              1              1 
    ##          21235        2123935        2125035        2125135      2125535-2 
    ##              1              1              1              1              1 
    ##        2125635        2125935        2126035        2127535        2129035 
    ##              1              1              1              1              1 
    ##        2129335        2129735        2131135        2132335        2133335 
    ##              1              1              1              1              1 
    ##        2133535        2136635        2137335        2139035        2139935 
    ##              1              1              1              1              1 
    ##        2140135        2140635        2141635        2142135        2145335 
    ##              1              1              1              1              1 
    ##        2145935        2147835        2148035        2148235        2148835 
    ##              1              1              1              1              1 
    ##        2149035        2149435        2150335        2151335        2152135 
    ##              1              1              1              1              1 
    ##        2153135        21535-2        2154035      2155135-2      2155335-2 
    ##              1              1              1              1              1 
    ##        2156735      2157535-2        2157635        2158335      2158835-2 
    ##              1              1              1              1              1 
    ##        2159635        2160135        2160435        2160635        2161635 
    ##              1              1              1              1              1 
    ##        2161735        2163435        2166535        2166835        2170035 
    ##              1              1              1              1              1 
    ##        2171735        2172335          21735        2174135        2174735 
    ##              1              1              1              1              1 
    ##        2175735        2177635        2178135        2179335        2181435 
    ##              1              1              1              1              1 
    ##        2182035        2183735        2186235        2187235        2187335 
    ##              1              1              1              1              1 
    ##        2189035        2190035        2192035      2193835-2      2196135-2 
    ##              1              1              1              1              1 
    ##        2198435        2198535        2198735        2199935        2200635 
    ##              1              1              1              1              1 
    ##        2202135        2206735        2208035        2208735        2209535 
    ##              1              1              1              1              1 
    ##        2213135          22135        2216135        2217635        2218635 
    ##              1              1              1              1              1 
    ##        2222035        2223635        2224435        2224735        2225335 
    ##              1              1              1              1              1 
    ##        2226135        2226735        2226935        2227435        2230035 
    ##              1              1              1              1              1 
    ##        2230635        2233135        2234935        2236035        2236235 
    ##              1              1              1              1              1 
    ##        2237035        2237835        2238935        2243935        2244335 
    ##              1              1              1              1              1 
    ##        2244735        2246535        2246635        2247335        2248935 
    ##              1              1              1              1              1 
    ##        2249935        2250935        2252035          22535        2253535 
    ##              1              1              1              1              1 
    ##        2253935        2258435      2258535-2        2258935      2263135-2 
    ##              1              1              1              1              1 
    ##        2264635        2265335        2267135        2267335        2268735 
    ##              1              1              1              1              1 
    ##        2269035        2269635        2270135        2270335        2271235 
    ##              1              1              1              1              1 
    ##        2271535        2271835        2272435        2272535        2272935 
    ##              1              1              1              1              1 
    ##        2273235          22735        2273535        2275935        2276435 
    ##              1              1              1              1              1 
    ##        2277935        2278035        2279835        2281335        2282535 
    ##              1              1              1              1              1 
    ##        2283735        2284135        2285435        2286335        2286735 
    ##              1              1              1              1              1 
    ##        2287035        2287135        2288235        2288435        2288735 
    ##              1              1              1              1              1 
    ##        2291735        2292535        2294235        2295835        2296435 
    ##              1              1              1              1              1 
    ##        2299135        2299935        2300235        2300335        2300535 
    ##              1              1              1              1              1 
    ##        2302035        2302335        2302935        2303835        2303935 
    ##              1              1              1              1              1 
    ##        2304535        2307635        2311535        2315135        2319135 
    ##              1              1              1              1              1 
    ##        2319335        2321835        2322135          23235        2324535 
    ##              1              1              1              2              1 
    ##        2326335        2327535        2329135        2329935        2330635 
    ##              1              1              1              1              1 
    ##        2331635        2334335        2334735        2337635        2339235 
    ##              1              1              1              1              1 
    ##        2339435        2341335        2341635        2343635        2346535 
    ##              1              1              1              1              1 
    ##        2348535        2349535        2349635        2350835        2350935 
    ##              1              1              1              1              1 
    ##        2352035        2352935        2353035        2353335        2353435 
    ##              1              1              1              1              1 
    ##        2353535        2353635        2355535        2357735        2358335 
    ##              1              1              1              1              1 
    ##        2361235        2362435          23635        2363535        2364635 
    ##              1              1              1              1              1 
    ##        2365035        2366635        2366735        2368635        2369835 
    ##              1              1              1              1              1 
    ##        2372335        2373235        2374435        2375335        2376035 
    ##              1              1              1              1              1 
    ##        2377635        2380135      2380935-2        2381135        2382835 
    ##              1              1              1              1              1 
    ##          23835        2384135        2384235        2384535        2386135 
    ##              2              1              1              1              1 
    ##        2386935        2389035        2389135      2389535-2          23935 
    ##              1              1              1              1              2 
    ##        2394535        2396835        2400035      2400835-2        2401335 
    ##              1              1              1              1              1 
    ##      2401635-2     2401735-2W        2403335        24035-2        2403535 
    ##              1              1              1              1              1 
    ##        2404035        2404135        2404235        2405635        2412135 
    ##              1              1              1              1              1 
    ##        2412235      2413135-2          24335          24835        2489235 
    ##              1              1              2              2              1 
    ##        2490535        2493235        2493335        2495535        2496535 
    ##              1              1              1              1              1 
    ##        2496735        2497135        2497235        2497835        2498235 
    ##              1              1              1              1              1 
    ##        2499235        2500935        2502135        2502435        2503335 
    ##              1              1              1              1              1 
    ##          25035        2509935        2512335        2512735        2513435 
    ##              1              1              1              1              1 
    ##          25135        2514935        2517335        2518235        2524135 
    ##              1              1              1              1              1 
    ##        2525335        2525835        2525935        2527935        2529935 
    ##              1              1              1              1              1 
    ##        2535735        2536035        2545135        2546335        2546935 
    ##              1              1              1              1              1 
    ##        2549135        2550435        2550735        2551135        2551935 
    ##              1              1              1              1              1 
    ##        2552735          25535     2559435-2W      2559835-2        2560135 
    ##              1              2              1              1              1 
    ##        2560335        2560435        2560935      2562935-2        2564835 
    ##              1              1              1              1              1 
    ##        2565135        2565235        2565435        2566235        2566735 
    ##              1              1              1              1              1 
    ##        2566835        2566935        2567135      2568135-2        2569335 
    ##              1              1              1              1              1 
    ##        2569635        2571135        2571735          25735     2579835-2W 
    ##              1              1              1              1              1 
    ##        2580135        2581435        2583035        2584735        2585335 
    ##              1              1              1              1              1 
    ##        2588335        2588735        2589735        2589835        2590135 
    ##              1              1              1              1              1 
    ##        2592035        2592435        2592635        2604635        2607535 
    ##              1              1              1              1              1 
    ##        2607635        2608935      2615035-2        2615535        2615635 
    ##              1              1              1              1              1 
    ##      2615935-2        2618235        2629335        2630135        2630635 
    ##              1              1              1              1              1 
    ##        2631535          26335        2640435        2649735        2650235 
    ##              1              1              1              1              1 
    ##        2650735        2651335        2659535        2659735        2660335 
    ##              1              1              1              1              1 
    ##        2662435        2662635          27735          29435          30335 
    ##              1              1              1              1              1 
    ##          32235          33035          33635          35935          36035 
    ##              1              1              1              1              1 
    ##          36335          37435          38735        41235-2          43135 
    ##              1              1              1              1              1 
    ##          43635          96535          96735          97235          97535 
    ##              1              1              1              1              1 
    ##        99635-2          99835 SAMEA115126061 SAMEA115126171 SAMEA115126306 
    ##              1              1              1              1              1 
    ## SAMEA115126326 SAMEA115126491 SAMEA115126500 SAMEA115126520 SAMEA115126613 
    ##              1              1              1              1              1 
    ## SAMEA115126615 SAMEA115126635   SAMEA1151267 SAMEA115126717 SAMEA115126775 
    ##              1              1              2              1              1 
    ## SAMEA115126863 SAMEA115126870 SAMEA115126888 SAMEA115127046 SAMEA115127059 
    ##              1              1              1              1              1 
    ## SAMEA115127165 SAMEA115127331 SAMEA115127484 SAMEA115127502 SAMEA115127533 
    ##              1              1              1              1              1 
    ## SAMEA115127719 
    ##              1

``` r
table(metadata$Site)
```

    ## 
    ## human saliva metagenome 
    ##                    1698

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
    ##        0000135        0000735        0002335        0002535        0002935 
    ##              1              1              1              1              1 
    ##        0003235        0004035        0004535        0005135        0005235 
    ##              1              1              1              1              1 
    ##        0005435        0005735        0006235        0006935        0007135 
    ##              1              1              1              1              1 
    ##        0007635        0008335        0009635        0012935        0013735 
    ##              1              1              1              1              1 
    ##        0014135        0014235        0014435        0014635        0014935 
    ##              1              1              1              1              1 
    ##        0015235        0015435        0015535        0015835        0016235 
    ##              1              1              1              1              1 
    ##        0016435        0016535        0017835        0017935        0018635 
    ##              1              1              1              1              1 
    ##        0019535      0020535-2        0020735        0020935        0022235 
    ##              1              1              1              1              1 
    ##        0022335        0023235        0023535        0023935      0025335-2 
    ##              1              1              1              1              1 
    ##        0025435      0025835-2        0027035        0027235        0027335 
    ##              1              1              1              1              1 
    ##      0027535-2        0027635      0028435-2        0029235        0029435 
    ##              1              1              1              1              1 
    ##        0029935        0030635        0032635        0032835        0032935 
    ##              1              1              1              1              1 
    ##        0034635        0035535        0038235        0039235        0045735 
    ##              1              1              1              1              1 
    ##        0046235        0046335        0046635        0047235        0047835 
    ##              1              1              1              1              1 
    ##        0048035        0048335        0048635        0049835        0050435 
    ##              1              1              1              1              1 
    ##        0050535        0050735        0051235        0051835        0053535 
    ##              1              1              1              1              1 
    ##        0053635        0054435        0055235        0056135        0056435 
    ##              1              1              1              1              1 
    ##        0057235        0057535        0058435        0058535        0058935 
    ##              1              1              1              1              1 
    ##        0059035        0059635        0060135        0060235        0062535 
    ##              1              1              1              1              1 
    ##        0063035      0063135-2        0063935        0064135        0064635 
    ##              1              1              1              1              1 
    ##        0065835        0066235        0067235        0068635        0070435 
    ##              1              1              1              1              1 
    ##        0072435        0072735      0073835-2        0073935        0074535 
    ##              1              1              1              1              1 
    ##        0074735        0074935        0076135        0076935        0078035 
    ##              1              1              1              1              1 
    ##        0078835        0079135        0079835        0080235        0080835 
    ##              1              1              1              1              1 
    ##        0082035        0084635        0085135        0085435        0089135 
    ##              1              1              1              1              1 
    ##        0089435        0090035        0090235        0092235        0092535 
    ##              1              1              1              1              1 
    ##        0095035        0097335        0098835        0099635        0113135 
    ##              1              1              1              1              1 
    ##        0113435        0113535        0113735        0115135        0117135 
    ##              1              1              1              1              1 
    ##        0117235        0118535        0119335        0120135        0120435 
    ##              1              1              1              1              1 
    ##        0121535        0122735        0123435        0123535        0124135 
    ##              1              1              1              1              1 
    ##        0126435        0126635        0127135        0127235        0128835 
    ##              1              1              1              1              1 
    ##        0136635        0136935        0137735        0138335        0139435 
    ##              1              1              1              1              1 
    ##        0139735        0140935        0141235        0141435        0141935 
    ##              1              1              1              1              1 
    ##        0142635        0144435        0146735        0146935        0147035 
    ##              1              1              1              1              1 
    ##        0147535        0147835        0151935        0153535        0153935 
    ##              1              1              1              1              1 
    ##        0154835        0156735        0157235        0160535      0160935-2 
    ##              1              1              1              1              1 
    ##        0161635        0161935        0162035        0163935        0164035 
    ##              1              1              1              1              1 
    ##        0164235        0164435        0164635        0165535        0165935 
    ##              1              1              1              1              1 
    ##        0166135      0166235-2        0166335        0166835        0167935 
    ##              1              1              1              1              1 
    ##        0168235        0168835        0169835        0170135        0170235 
    ##              1              1              1              1              1 
    ##        0171335      0178935-2        0182835        0182935      0186835-2 
    ##              1              1              1              1              1 
    ##        0192635        0194535        0195835        0196135        0200835 
    ##              1              1              1              1              1 
    ##        0201235        0202135        0203535        0203835        0206535 
    ##              1              1              1              1              1 
    ##        0220135        0220235        0220835        0222835        0223535 
    ##              1              1              1              1              1 
    ##        0225635        0228235      0230735-2        0231935        0233435 
    ##              1              1              1              1              1 
    ##        0234435        0236935        0238235        0240135        0241235 
    ##              1              1              1              1              1 
    ##        0241735        0241935        0242035        0242635          02435 
    ##              1              1              1              1              1 
    ##        0243935        0245135        0245635        0253435          02535 
    ##              1              1              1              1              1 
    ##        0253835        0254335        0254835        0255235        0260835 
    ##              1              1              1              1              1 
    ##        0261735        0262335        0263235        0264035        0266335 
    ##              1              1              1              1              1 
    ##        0267935        0268035        0268235        0268735        0269835 
    ##              1              1              1              1              1 
    ##        0272935        0274135        0277935        0279535        0291735 
    ##              1              1              1              1              1 
    ##        0299535        0299835        0300135        0303935        0306235 
    ##              1              1              1              1              1 
    ##        0306635        0306735        0306835        0308035        0309235 
    ##              1              1              1              1              1 
    ##        0315735        0316935        0318135        0318435        0323935 
    ##              1              1              1              1              1 
    ##        0325435        0329635        0332235        0332735        0334635 
    ##              1              1              1              1              1 
    ##        0335635        0336235        0336435        0336935        0339235 
    ##              1              1              1              1              1 
    ##        0341835        0342535        0342635        0343135        0343835 
    ##              1              1              1              1              1 
    ##        0344835        0347935        0348235        0350635        0352535 
    ##              1              1              1              1              1 
    ##        0352935          03535        0354235        0356135        0358035 
    ##              1              1              1              1              1 
    ##        0360335        0364935        0366535        0367035        0367135 
    ##              1              1              1              1              1 
    ##        0367935        0371235        0371535        0378535        0379335 
    ##              1              1              1              1              1 
    ##        0379635        0380935        0381735        0382135        0383335 
    ##              1              1              1              1              1 
    ##        0383835        0383935        0385235        0386035        0386135 
    ##              1              1              1              1              1 
    ##        0387935        0389235        0389335        0391335        0391435 
    ##              1              1              1              1              1 
    ##        0391735        0392335        0393435          03935        0394335 
    ##              1              1              1              1              1 
    ##        0395535        0396135        0397735        0398735        0399135 
    ##              1              1              1              1              1 
    ##        0399835        0402635        0403035        0403635        0406135 
    ##              1              1              1              1              1 
    ##        0406235        0406335        0408335        0408735        0409435 
    ##              1              1              1              1              1 
    ##          04135        0414335        0414835        0417735        0417935 
    ##              1              1              1              1              1 
    ##        0419335        0422835        0425235        0426535        0427735 
    ##              1              1              1              1              1 
    ##        0430335        0431135        0432235        0432735        0432935 
    ##              1              1              1              1              1 
    ##        0434735        0439235        0439635        0440135        0443735 
    ##              1              1              1              1              1 
    ##        0443935        0444935        0446035        0446135        0446235 
    ##              1              1              1              1              1 
    ##        0448135        0448735        0449335        0450035        0450935 
    ##              1              1              1              1              1 
    ##        0451235        0451335        0452035        0453335          04535 
    ##              1              1              1              1              1 
    ##        0453535        0453635        0454035        0456535        0460035 
    ##              1              1              1              1              1 
    ##        0460135        0461535        0463335        0463535        0464235 
    ##              1              1              1              1              1 
    ##        0465335        0466235        0467435        0467635        0467835 
    ##              1              1              1              1              1 
    ##        0468035        0468335        0468835        0469635        0470835 
    ##              1              1              1              1              1 
    ##        0470935        0471135        0473135        0473335        0473635 
    ##              1              1              1              1              1 
    ##        0475435        0475535        0476535        0476835      0477535-2 
    ##              1              1              1              1              1 
    ##      0477835-2        0479535        0479635        0480135        0480535 
    ##              1              1              1              1              1 
    ##      0480735-2        0481535        0481935        0483435        0483835 
    ##              1              1              1              1              1 
    ##        0485035        0485635        0485835        0486235      0489435-2 
    ##              1              1              1              1              1 
    ##      0490435-2        0490535        0491335        0492135        0494035 
    ##              1              1              1              1              1 
    ##        0494435        0494735        0494835        0496135        0496735 
    ##              1              1              1              1              1 
    ##        0496935        0500835        0501135        0501835        0502535 
    ##              1              1              1              1              1 
    ##        0503635        0505735        0506235        0507835        0509835 
    ##              1              1              1              1              1 
    ##        0511135        0512035      0513035-2        0513335        0515135 
    ##              1              1              1              1              1 
    ##        0517335        0518835        0519235        0524435        0525835 
    ##              1              1              1              1              1 
    ##        0527635        0528435        0528735        0531335          05435 
    ##              1              1              1              1              1 
    ##        0561135        0563435        0563735        0564135        0565235 
    ##              1              1              1              1              1 
    ##        0565435        0566335        0566435        0566535        0566835 
    ##              1              1              1              1              1 
    ##        0567435        0567835        0569235        0570035      0570635-2 
    ##              1              1              1              1              1 
    ##        0572435        0572835        0573735        0575835        0576735 
    ##              1              1              1              1              1 
    ##        0577235        0577335        0578435        0579635        0579935 
    ##              1              1              1              1              1 
    ##        0580135        0580735        0580835        0582235        0582335 
    ##              1              1              1              1              1 
    ##        0582435          05835        0583635        0583735        0584035 
    ##              1              1              1              1              1 
    ##        0584435        0584735        0585635        0586735        0587635 
    ##              1              1              1              1              1 
    ##        0587835      0588635-2        0588735        0589235        0589435 
    ##              1              1              1              1              1 
    ##        0589535        0589735        0591235        0591335        0592235 
    ##              1              1              1              1              1 
    ##        0593235        0594035        0594135        0595535        0595835 
    ##              1              1              1              1              1 
    ##        0595935        0597835        0597935        0598035        0598935 
    ##              1              1              1              1              1 
    ##        0599335        0599435        0599535        0599835        0601135 
    ##              1              1              1              1              1 
    ##        0602835        0603435        0603635        0604235        0604635 
    ##              1              1              1              1              1 
    ##        0604935        0605535        0607035        0607935        0608135 
    ##              1              1              1              1              1 
    ##        0608235        0609535        0631635        0636135        0637835 
    ##              1              1              1              1              1 
    ##        0638035        0638135        0638735        0638835        0638935 
    ##              1              1              1              1              1 
    ##        0639535        0639735        0640035        0641235        0641435 
    ##              1              1              1              1              1 
    ##          06735        0680835        0681435        0681935        0682435 
    ##              1              1              1              1              1 
    ##          06835        0683735        0685035        0685135        0685835 
    ##              2              1              1              1              1 
    ##        0687935        0691935      0692735-2     0693835-2W        0695335 
    ##              1              1              1              1              1 
    ##        0695735        0695935        0696535        0696835        0697835 
    ##              1              1              1              1              1 
    ##        0700535        0702835        0703535        0704235      0705635-2 
    ##              1              1              1              1              1 
    ##        0705935        0706235      0706535-2        0707535        0707735 
    ##              1              1              1              1              1 
    ##        0708935      0711435-2        0712235        0713435          07135 
    ##              1              1              1              1              1 
    ##        0713735        0713835      0715135-2        0715935        0718935 
    ##              1              1              1              1              1 
    ##          07235        0723535        0729435        0730535        0730935 
    ##              1              1              1              1              1 
    ##        0731935        0732135        0732335        0733235        0733935 
    ##              1              1              1              1              1 
    ##        0734235        0735035        0735435        0735535        0736535 
    ##              1              1              1              1              1 
    ##        0740635        0742635        0743335          07435        0743635 
    ##              1              1              1              1              1 
    ##        0743835        0744235        0744735        0745035        0746535 
    ##              1              1              1              1              1 
    ##        0747035        0749535        0749835        0752035        0752535 
    ##              1              1              1              1              1 
    ##        0752735          07535        0753935        0755435        0755635 
    ##              1              1              1              1              1 
    ##        0756235        0756735        0757835        0759735        0759835 
    ##              1              1              1              1              1 
    ##        0763635        0765935        0766835        0768635        0769235 
    ##              1              1              1              1              1 
    ##        0769535        0769735        0771135        0771835        0774835 
    ##              1              1              1              1              1 
    ##        0774935        0775035        0776335        0776435     0778235-2W 
    ##              1              1              1              1              1 
    ##        0778735        0787735        0790335        0790935        0791135 
    ##              1              1              1              1              1 
    ##        0792035        0792935        0793235        0793335     0793835-2W 
    ##              1              1              1              1              1 
    ##        0796535        0797735        0797835        0798935        0800735 
    ##              1              1              1              1              1 
    ##        0800835        0801335        0802435        0804035        0805135 
    ##              1              1              1              1              1 
    ##        0807235        0807935        0809735        0811935          08135 
    ##              1              1              1              1              1 
    ##        0814335        0815535        0816935        0818835        0819135 
    ##              1              1              1              1              1 
    ##        0825635        0829735        0830535        0831735        0834435 
    ##              1              1              1              1              1 
    ##      0835235-2        0835935        0837235        0837735        0837835 
    ##              1              1              1              1              1 
    ##        0838435        0840435        0841835        0841935        0842035 
    ##              1              1              1              1              1 
    ##        0842335        0842835          08435        0843635        0844335 
    ##              1              1              1              1              1 
    ##        0844435        0847435        0850935        0852035        0852535 
    ##              1              1              1              1              1 
    ##        0852835        0853235        0853335        0853735        0854535 
    ##              1              1              1              1              1 
    ##        0855035        0859135        0860735        0868235        0872435 
    ##              1              1              1              1              1 
    ##          08735        0873635        0873735        0875835        0876135 
    ##              1              1              1              1              1 
    ##        0876535        0879535        0881235        0883535        0883835 
    ##              1              1              1              1              1 
    ##        0885335      0886035-2        0886435        0887035        0887135 
    ##              1              1              1              1              1 
    ##        0887835        0887935        0888735        0889235        0889435 
    ##              1              1              1              1              1 
    ##      0891235-2      0893735-2        0893935        0945035        0946235 
    ##              1              1              1              1              1 
    ##        0948735        0960035        0960235        0960735        0961235 
    ##              1              1              1              1              1 
    ##        0963435          09635        0964035        0964235        0964735 
    ##              1              1              1              1              1 
    ##        0967735        0969835        0973335      0974035-2        0975235 
    ##              1              1              1              1              1 
    ##        0976935        0986635        0986735        0988435        0989235 
    ##              1              1              1              1              1 
    ##        0990935        0991335        0992835        0994135        0995935 
    ##              1              1              1              1              1 
    ##        0996135        0996935        0997135        0997235        0998135 
    ##              1              1              1              1              1 
    ##        0998535        0998935        0999035        10135-2          10235 
    ##              1              1              1              1              1 
    ##          11035        1111135        1112435        1112935        1113335 
    ##              1              1              1              1              1 
    ##          11135        1114335        1114735        1119135        1119235 
    ##              2              1              1              1              1 
    ##        1120235      1120335-2        1120935        1121335        1121435 
    ##              1              1              1              1              1 
    ##        1122935        1123535        1159335        1159835        1160835 
    ##              1              1              1              1              1 
    ##        1161435        1164235        1165335        1166435        1170635 
    ##              1              1              1              1              1 
    ##        1170835        1173135          11735        12135-2      1233635-2 
    ##              1              1              1              1              1 
    ##        1234335        1234735        1237135        1237535        1237935 
    ##              1              1              1              1              1 
    ##        1240935        1241835        1241935        1242635        1243035 
    ##              1              1              1              1              1 
    ##        1247035        1248435        1249635          12535        1256335 
    ##              1              1              1              1              1 
    ##        1260935        1262335        1262435        1264135        1265035 
    ##              1              1              1              1              1 
    ##        1265235        1266635        1268435        1268535        1274335 
    ##              1              1              1              1              1 
    ##        1275435        1278235        1280735        1281235        1282235 
    ##              1              1              1              1              1 
    ##        1282435        1285035        1286435        1288335        1290735 
    ##              1              1              1              1              1 
    ##        1291935      1294535-2        1295335      1295735-2      1296535-2 
    ##              1              1              1              1              1 
    ##      1296935-2        1302735        1303635        1304735        1321135 
    ##              1              1              1              1              1 
    ##        1324935        1334035        1336635        1341635        1342635 
    ##              1              1              1              1              1 
    ##        1344035        1349335        1350135        1350435        1352235 
    ##              1              1              1              1              1 
    ##        1352835        1354835        1355135        1357035        1357635 
    ##              1              1              1              1              1 
    ##        1359835        1361835        1369035        1370535        1371735 
    ##              1              1              1              1              1 
    ##        1372235          13735        1379435        1379535      1380335-2 
    ##              1              1              1              1              1 
    ##        1381135        1381335      1381835-2        1382735        1383035 
    ##              1              1              1              1              1 
    ##      1383435-2        1387135        1391135        1399135        1399335 
    ##              1              1              1              1              1 
    ##      1399835-2        1401635        1401735        1401935      1402135-2 
    ##              1              1              1              1              1 
    ##        1402335          14035          14135          14335        1493835 
    ##              1              1              1              1              1 
    ##        1494035        1494135        1497235        1504635        1504735 
    ##              1              1              1              1              1 
    ##        1505335        1505535        1505835        1506335        1506535 
    ##              1              1              1              1              1 
    ##        1506635        1506935        1507435        1508635        1509735 
    ##              1              1              1              1              1 
    ##        1511135        1511935        1512435        1513135        1514335 
    ##              1              1              1              1              1 
    ##        1514435        1515335        1515535        1516435        1517035 
    ##              1              1              1              1              1 
    ##        1517735        1517835        1518435        1520235        1520935 
    ##              1              1              1              1              1 
    ##        1521135        1522435          15235        1524335        1524535 
    ##              1              1              1              1              1 
    ##        1525235        1527635        1529335        1530735        1531635 
    ##              1              1              1              1              1 
    ##        1532835        1533335          15335        1535435        1535835 
    ##              1              1              1              1              1 
    ##        1536335        1536835        1537335        1538235        1538435 
    ##              1              1              1              1              1 
    ##        1546235        1547535        1547935        1549035        1551135 
    ##              1              1              1              1              1 
    ##        1555735        1558435        1558635        1559335        1560035 
    ##              1              1              1              1              1 
    ##        1560135      1561635-2        1564635        1564935        1567135 
    ##              1              1              1              1              1 
    ##        1567335        1570335        1570535        1573235          15735 
    ##              1              1              1              1              2 
    ##        1573735        1574835        1576435        1576535        1576835 
    ##              1              1              1              1              1 
    ##        1579335        1581135          15835        1584435        1588135 
    ##              1              1              1              1              1 
    ##        1591535        1593235          15935        1593635        1594035 
    ##              1              1              1              1              1 
    ##        1634535        1635835        1637135        1637635        1638735 
    ##              1              1              1              1              1 
    ##        1639235        1640335        1642135        1643035        1643835 
    ##              1              1              1              1              1 
    ##        1644535        1645235        1648035        1669135        1670335 
    ##              1              1              1              1              1 
    ##        1670535        1672435        1673735        1673835        1691435 
    ##              1              1              1              1              1 
    ##        1691535        1691935        1692735        1699535        1700635 
    ##              1              1              1              1              1 
    ##        1700835        1701735      1702135-2          17035        1711635 
    ##              1              1              1              1              1 
    ##        1712135        1712335        1712935          17135        1713535 
    ##              1              1              1              1              1 
    ##        1714735        1716135          17535        1779732        1779932 
    ##              1              1              1              1              1 
    ##        1780132        1780832        1781632        1781832        1781932 
    ##              1              1              1              1              1 
    ##        1782232        1783232        1784132        1784432          17935 
    ##              1              1              1              1              1 
    ##        1859435        1859635        1861635        1862335        1862535 
    ##              1              1              1              1              1 
    ##          18635        1865235        1865535        1866735        1867135 
    ##              1              1              1              1              1 
    ##        1867335        1867935        1868035          18735      1875235-2 
    ##              1              1              1              1              1 
    ##        1877435      1877935-2        1880035        1881135        1882235 
    ##              1              1              1              1              1 
    ##        1882435        1882635        1883835      1884435-2        1884835 
    ##              1              1              1              1              1 
    ##        1887635        1887735        1889035        1889335        1890235 
    ##              1              1              1              1              1 
    ##        1891735        1894235        1894635        1895135        1897335 
    ##              1              1              1              1              1 
    ##        1898335        1899035        1899235        1901335        1901735 
    ##              1              1              1              1              1 
    ##          19035        1904035        1904935        1905035        1906335 
    ##              1              1              1              1              1 
    ##        1906435        1907035        1909635        1912335        1912635 
    ##              1              1              1              1              1 
    ##        1912935        1913135        1913735        1914535        1914835 
    ##              1              1              1              1              1 
    ##        1915235        1916435        1917635        1918535        1919135 
    ##              1              1              1              1              1 
    ##        1920135        1920335      1921235-2        1921535      1921835-2 
    ##              1              1              1              1              1 
    ##     1922035-2W      1922135-2        1922935      1923635-2      1923935-2 
    ##              1              1              1              1              1 
    ##        1924335        1926135        1928235        1929835        1932035 
    ##              1              1              1              1              1 
    ##        1932635        1932835        1933035        1933635        1933835 
    ##              1              1              1              1              1 
    ##        1934035        1935535        1938235      1939335-2        1943235 
    ##              1              1              1              1              1 
    ##          19435        1944035        1944435        1946435      1947335-2 
    ##              2              1              1              1              1 
    ##        1947835        1948335        1948835        1949035        1952635 
    ##              1              1              1              1              1 
    ##      1953135-2          19535        1956235        1960235        1962635 
    ##              1              1              1              1              1 
    ##        1963435          19635        1970935        1972635        1972835 
    ##              1              1              1              1              1 
    ##        1973435        1974135        1974335        1974535        1977035 
    ##              1              1              1              1              1 
    ##        1977935        1980435        1981135        1981335        1981835 
    ##              1              1              1              1              1 
    ##        1982235        1983635        1984635        1984835        1985635 
    ##              1              1              1              1              1 
    ##        1985935        1986035        1986235        1987435        1987635 
    ##              1              1              1              1              1 
    ##        1988535        1989335        1991735        1993035          19935 
    ##              1              1              1              1              1 
    ##        1994635        1995035        1995835        1997335        1997435 
    ##              1              1              1              1              1 
    ##        1997735        1998235        1999835        2000935        2001935 
    ##              1              1              1              1              1 
    ##        2003035          20035        2003635        2003835        2003935 
    ##              1              1              1              1              1 
    ##      2004335-2        2005035        2006335        2006435        2008035 
    ##              1              1              1              1              1 
    ##          20135        2013935        2014135        2014935        2015135 
    ##              1              1              1              1              1 
    ##        2015335        2015435        2016735        2017135        2017735 
    ##              1              1              1              1              1 
    ##        2018035        2019835        2020235        2020335        2021435 
    ##              1              1              1              1              1 
    ##        2023035          20235        2025635        2026035        2026235 
    ##              1              1              1              1              1 
    ##        2026535        2027535        2029335        2029635        2030035 
    ##              1              1              1              1              1 
    ##        2030335        2031335        2031735        2032335        2032935 
    ##              1              1              1              1              1 
    ##        2033035        2033235        2033435        20335-2        2036335 
    ##              1              1              1              1              1 
    ##        2037035        2037535        2038735        2039935        2040535 
    ##              1              1              1              1              1 
    ##        2040635        2040735        2042935        2046735        2049535 
    ##              1              1              1              1              1 
    ##        2050035        2050735        2050835          20535        2057035 
    ##              1              1              1              2              1 
    ##        2058235        2063535        2064035        2064735        2066435 
    ##              1              1              1              1              1 
    ##        2066935        2068435        2068535        2068935        2069735 
    ##              1              1              1              1              1 
    ##        2069835        2070635        2076335      2076635-2      2077035-2 
    ##              1              1              1              1              1 
    ##      2077235-2        2078235        2078635        2081235        2082235 
    ##              1              1              1              1              1 
    ##        2084235        2084735        2087935        2088235        2088435 
    ##              1              1              1              1              1 
    ##        2088735        2088835        2089235        2089935        2090035 
    ##              1              1              1              1              1 
    ##        2090635        2090935        2091235      2091335-2        2091635 
    ##              1              1              1              1              1 
    ##        2096335        2097035        2097335        2099635        2099735 
    ##              1              1              1              1              1 
    ##        2111135      2112535-2        2112835        2113235        2115335 
    ##              1              1              1              1              1 
    ##        2116435        2117235        2118735        2119135        2119335 
    ##              1              1              1              1              1 
    ##          21235        2123935        2125035        2125135      2125535-2 
    ##              1              1              1              1              1 
    ##        2125635        2125935        2126035        2127535        2129035 
    ##              1              1              1              1              1 
    ##        2129335        2129735        2131135        2132335        2133335 
    ##              1              1              1              1              1 
    ##        2133535        2136635        2137335        2139035        2139935 
    ##              1              1              1              1              1 
    ##        2140135        2140635        2141635        2142135        2145335 
    ##              1              1              1              1              1 
    ##        2145935        2147835        2148035        2148235        2148835 
    ##              1              1              1              1              1 
    ##        2149035        2149435        2150335        2151335        2152135 
    ##              1              1              1              1              1 
    ##        2153135        21535-2        2154035      2155135-2      2155335-2 
    ##              1              1              1              1              1 
    ##        2156735      2157535-2        2157635        2158335      2158835-2 
    ##              1              1              1              1              1 
    ##        2159635        2160135        2160435        2160635        2161635 
    ##              1              1              1              1              1 
    ##        2161735        2163435        2166535        2166835        2170035 
    ##              1              1              1              1              1 
    ##        2171735        2172335          21735        2174135        2174735 
    ##              1              1              1              1              1 
    ##        2175735        2177635        2178135        2179335        2181435 
    ##              1              1              1              1              1 
    ##        2182035        2183735        2186235        2187235        2187335 
    ##              1              1              1              1              1 
    ##        2189035        2190035        2192035      2193835-2      2196135-2 
    ##              1              1              1              1              1 
    ##        2198435        2198535        2198735        2199935        2200635 
    ##              1              1              1              1              1 
    ##        2202135        2206735        2208035        2208735        2209535 
    ##              1              1              1              1              1 
    ##        2213135          22135        2216135        2217635        2218635 
    ##              1              1              1              1              1 
    ##        2222035        2223635        2224435        2224735        2225335 
    ##              1              1              1              1              1 
    ##        2226135        2226735        2226935        2227435        2230035 
    ##              1              1              1              1              1 
    ##        2230635        2233135        2234935        2236035        2236235 
    ##              1              1              1              1              1 
    ##        2237035        2237835        2238935        2243935        2244335 
    ##              1              1              1              1              1 
    ##        2244735        2246535        2246635        2247335        2248935 
    ##              1              1              1              1              1 
    ##        2249935        2250935        2252035          22535        2253535 
    ##              1              1              1              1              1 
    ##        2253935        2258435      2258535-2        2258935      2263135-2 
    ##              1              1              1              1              1 
    ##        2264635        2265335        2267135        2267335        2268735 
    ##              1              1              1              1              1 
    ##        2269035        2269635        2270135        2270335        2271235 
    ##              1              1              1              1              1 
    ##        2271535        2271835        2272435        2272535        2272935 
    ##              1              1              1              1              1 
    ##        2273235          22735        2273535        2275935        2276435 
    ##              1              1              1              1              1 
    ##        2277935        2278035        2279835        2281335        2282535 
    ##              1              1              1              1              1 
    ##        2283735        2284135        2285435        2286335        2286735 
    ##              1              1              1              1              1 
    ##        2287035        2287135        2288235        2288435        2288735 
    ##              1              1              1              1              1 
    ##        2291735        2292535        2294235        2295835        2296435 
    ##              1              1              1              1              1 
    ##        2299135        2299935        2300235        2300335        2300535 
    ##              1              1              1              1              1 
    ##        2302035        2302335        2302935        2303835        2303935 
    ##              1              1              1              1              1 
    ##        2304535        2307635        2311535        2315135        2319135 
    ##              1              1              1              1              1 
    ##        2319335        2321835        2322135          23235        2324535 
    ##              1              1              1              2              1 
    ##        2326335        2327535        2329135        2329935        2330635 
    ##              1              1              1              1              1 
    ##        2331635        2334335        2334735        2337635        2339235 
    ##              1              1              1              1              1 
    ##        2339435        2341335        2341635        2343635        2346535 
    ##              1              1              1              1              1 
    ##        2348535        2349535        2349635        2350835        2350935 
    ##              1              1              1              1              1 
    ##        2352035        2352935        2353035        2353335        2353435 
    ##              1              1              1              1              1 
    ##        2353535        2353635        2355535        2357735        2358335 
    ##              1              1              1              1              1 
    ##        2361235        2362435          23635        2363535        2364635 
    ##              1              1              1              1              1 
    ##        2365035        2366635        2366735        2368635        2369835 
    ##              1              1              1              1              1 
    ##        2372335        2373235        2374435        2375335        2376035 
    ##              1              1              1              1              1 
    ##        2377635        2380135      2380935-2        2381135        2382835 
    ##              1              1              1              1              1 
    ##          23835        2384135        2384235        2384535        2386135 
    ##              2              1              1              1              1 
    ##        2386935        2389035        2389135      2389535-2          23935 
    ##              1              1              1              1              2 
    ##        2394535        2396835        2400035      2400835-2        2401335 
    ##              1              1              1              1              1 
    ##      2401635-2     2401735-2W        2403335        24035-2        2403535 
    ##              1              1              1              1              1 
    ##        2404035        2404135        2404235        2405635        2412135 
    ##              1              1              1              1              1 
    ##        2412235      2413135-2          24335          24835        2489235 
    ##              1              1              2              2              1 
    ##        2490535        2493235        2493335        2495535        2496535 
    ##              1              1              1              1              1 
    ##        2496735        2497135        2497235        2497835        2498235 
    ##              1              1              1              1              1 
    ##        2499235        2500935        2502135        2502435        2503335 
    ##              1              1              1              1              1 
    ##          25035        2509935        2512335        2512735        2513435 
    ##              1              1              1              1              1 
    ##          25135        2514935        2517335        2518235        2524135 
    ##              1              1              1              1              1 
    ##        2525335        2525835        2525935        2527935        2529935 
    ##              1              1              1              1              1 
    ##        2535735        2536035        2545135        2546335        2546935 
    ##              1              1              1              1              1 
    ##        2549135        2550435        2550735        2551135        2551935 
    ##              1              1              1              1              1 
    ##        2552735          25535     2559435-2W      2559835-2        2560135 
    ##              1              2              1              1              1 
    ##        2560335        2560435        2560935      2562935-2        2564835 
    ##              1              1              1              1              1 
    ##        2565135        2565235        2565435        2566235        2566735 
    ##              1              1              1              1              1 
    ##        2566835        2566935        2567135      2568135-2        2569335 
    ##              1              1              1              1              1 
    ##        2569635        2571135        2571735          25735     2579835-2W 
    ##              1              1              1              1              1 
    ##        2580135        2581435        2583035        2584735        2585335 
    ##              1              1              1              1              1 
    ##        2588335        2588735        2589735        2589835        2590135 
    ##              1              1              1              1              1 
    ##        2592035        2592435        2592635        2604635        2607535 
    ##              1              1              1              1              1 
    ##        2607635        2608935      2615035-2        2615535        2615635 
    ##              1              1              1              1              1 
    ##      2615935-2        2618235        2629335        2630135        2630635 
    ##              1              1              1              1              1 
    ##        2631535          26335        2640435        2649735        2650235 
    ##              1              1              1              1              1 
    ##        2650735        2651335        2659535        2659735        2660335 
    ##              1              1              1              1              1 
    ##        2662435        2662635          27735          29435          30335 
    ##              1              1              1              1              1 
    ##          32235          33035          33635          35935          36035 
    ##              1              1              1              1              1 
    ##          36335          37435          38735        41235-2          43135 
    ##              1              1              1              1              1 
    ##          43635          96535          96735          97235          97535 
    ##              1              1              1              1              1 
    ##        99635-2          99835 SAMEA115126061 SAMEA115126171 SAMEA115126306 
    ##              1              1              1              1              1 
    ## SAMEA115126326 SAMEA115126491 SAMEA115126500 SAMEA115126520 SAMEA115126613 
    ##              1              1              1              1              1 
    ## SAMEA115126615 SAMEA115126635   SAMEA1151267 SAMEA115126717 SAMEA115126775 
    ##              1              1              2              1              1 
    ## SAMEA115126863 SAMEA115126870 SAMEA115126888 SAMEA115127046 SAMEA115127059 
    ##              1              1              1              1              1 
    ## SAMEA115127165 SAMEA115127331 SAMEA115127484 SAMEA115127502 SAMEA115127533 
    ##              1              1              1              1              1 
    ## SAMEA115127719 
    ##              1

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
    ##   NO  YES 
    ## 1674   24

``` r
# Identify samples for removal: single samples excluding saliva and tongue coat
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & 
                     !(metadata$Site %in% c("saliva", "tongue coat", "Tongue coating"))

cat("\nNumber of samples to be removed:", sum(samples_to_remove), "\n")
```

    ## 
    ## Number of samples to be removed: 1674

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
    ##        06835        11135        15735        19435        20535        23235 
    ##            2            2            2            2            2            2 
    ##        23835        23935        24335        24835        25535 SAMEA1151267 
    ##            2            2            2            2            2            2

``` r
table(metadata_clean$Site)
```

    ## 
    ## human saliva metagenome 
    ##                      24

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ## YES 
    ##  24

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 24"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 12"

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

    ##     BioProject         Run  sample_name   Age BMI Sex                  Region
    ## 32  PRJEB70783 ERR12514094 SAMEA1151267 35-70  NA  NA human saliva metagenome
    ## 342 PRJEB70783 ERR12514570        06835 35-70  NA  NA human saliva metagenome
    ## 353 PRJEB70783 ERR12514585        06835 35-70  NA  NA human saliva metagenome
    ## 485 PRJEB70783 ERR12514798        20535 35-70  NA  NA human saliva metagenome
    ## 488 PRJEB70783 ERR12514802        23935 35-70  NA  NA human saliva metagenome
    ## 489 PRJEB70783 ERR12514804        24835 35-70  NA  NA human saliva metagenome
    ##           Platform Country                    Site          public_time
    ## 32  Illumina MiSeq  Canada human saliva metagenome 2024-01-23T00:00:00Z
    ## 342 Illumina MiSeq  Canada human saliva metagenome 2024-01-23T00:00:00Z
    ## 353 Illumina MiSeq  Canada human saliva metagenome 2024-01-23T00:00:00Z
    ## 485 Illumina MiSeq  Canada human saliva metagenome 2024-01-23T00:00:00Z
    ## 488 Illumina MiSeq  Canada human saliva metagenome 2024-01-23T00:00:00Z
    ## 489 Illumina MiSeq  Canada human saliva metagenome 2024-01-23T00:00:00Z
    ##     PCR_primers
    ## 32    515F-926R
    ## 342   515F-926R
    ## 353   515F-926R
    ## 485   515F-926R
    ## 488   515F-926R
    ## 489   515F-926R
    ##                                                                            title
    ## 32  Oral microbial signatures associated with age and frailty in Canadian adults
    ## 342 Oral microbial signatures associated with age and frailty in Canadian adults
    ## 353 Oral microbial signatures associated with age and frailty in Canadian adults
    ## 485 Oral microbial signatures associated with age and frailty in Canadian adults
    ## 488 Oral microbial signatures associated with age and frailty in Canadian adults
    ## 489 Oral microbial signatures associated with age and frailty in Canadian adults
    ##     Disease_study             Status matched_samples
    ## 32            YES Ulcerative colitis             YES
    ## 342           YES Ulcerative colitis             YES
    ## 353           YES      Breast cancer             YES
    ## 485           YES Ulcerative colitis             YES
    ## 488           YES     Crohns disease             YES
    ## 489           YES Ulcerative colitis             YES

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
