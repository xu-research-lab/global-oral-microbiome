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

# Extract "human" and "metagenome"
initial_meta $Organism <- sub("human (\\w+) metagenome", "\\1", initial_meta$Organism)

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
    ##   nearing_etal_2020_healthy_oral_microbiome nearing_etal_2023_cancer Organism
    ## 1                                         N                        Y   saliva
    ## 2                                         Y                        N   saliva
    ## 3                                         Y                        N   saliva
    ## 4                                         Y                        N   saliva
    ## 5                                         Y                        N   saliva
    ## 6                                         Y                        N   saliva
    ##   Platform      primers project_name          ReleaseDate          create_date
    ## 1 ILLUMINA 515FBï¿½926R   PRJEB70783 2024-01-23T00:00:00Z 2024-01-25T10:29:00Z
    ## 2 ILLUMINA 515FBï¿½926R   PRJEB70783 2024-01-23T00:00:00Z 2024-01-26T22:43:00Z
    ## 3 ILLUMINA 515FBï¿½926R   PRJEB70783 2024-01-23T00:00:00Z 2024-01-26T22:41:00Z
    ## 4 ILLUMINA 515FBï¿½926R   PRJEB70783 2024-01-23T00:00:00Z 2024-01-26T17:23:00Z
    ## 5 ILLUMINA   515FB?926R   PRJEB70783 2024-01-23T00:00:00Z 2024-01-26T22:43:00Z
    ## 6 ILLUMINA 515FBï¿½926R   PRJEB70783 2024-01-23T00:00:00Z 2024-01-26T17:21:00Z
    ##   version    Sample.Name         Scientific_Name Sequencing_method SRA.Study
    ## 1       1      100000735 human saliva metagenome             MiSeq ERP155686
    ## 2       1      100002335 human saliva metagenome             MiSeq ERP155686
    ## 3       1      100002535 human saliva metagenome             MiSeq ERP155686
    ## 4       1      100003235 human saliva metagenome             MiSeq ERP155686
    ## 5       1 SAMEA115126061 human saliva metagenome             MiSeq ERP155686
    ## 6       1      100004535 human saliva metagenome             MiSeq ERP155686
    ##                  organism cancer_case_control cancer_retrospective_prospective
    ## 1 human saliva metagenome             Control                      prospective
    ## 2 human saliva metagenome                                                     
    ## 3 human saliva metagenome                                                     
    ## 4 human saliva metagenome                                                     
    ## 5                                                                             
    ## 6 human saliva metagenome                                                     
    ##   cancer_type IBD_case_control ibd_type     ena_first_public
    ## 1      Breast                      <NA>                     
    ## 2        <NA>                      <NA>                     
    ## 3        <NA>                      <NA>                     
    ## 4        <NA>                      <NA>                     
    ## 5        <NA>             Case   Crohns 2024-01-22T20:33:36Z
    ## 6        <NA>                      <NA>                     
    ##        ena_last_update    external_id    INSDC_center_name   INSDC_first_public
    ## 1                                                                              
    ## 2                                                                              
    ## 3                                                                              
    ## 4                                                                              
    ## 5 2024-01-22T20:33:36Z SAMEA115126061 Dalhousie University 2024-01-22T20:33:36Z
    ## 6                                                                              
    ##      INSDC_last_update INSDC_status Sample_name Submitter_Id
    ## 1                                                           
    ## 2                                                           
    ## 3                                                           
    ## 4                                                           
    ## 5 2024-01-22T20:33:36Z       public   100003935    100003935
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
  Region = "V4–V5",
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
      TRUE ~ "HC")
    )

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##      100000135      100000735      100002335      100002535      100002935 
    ##              1              1              1              1              1 
    ##      100003235      100004035      100004535      100005135      100005235 
    ##              1              1              1              1              1 
    ##      100005435      100005735      100006235      100006935      100007135 
    ##              1              1              1              1              1 
    ##      100007635      100008335      100009635      100012935      100013735 
    ##              1              1              1              1              1 
    ##      100014135      100014235      100014435      100014635      100014935 
    ##              1              1              1              1              1 
    ##      100015235      100015435      100015535      100015835      100016235 
    ##              1              1              1              1              1 
    ##      100016435      100016535      100017835      100017935      100018635 
    ##              1              1              1              1              1 
    ##      100019535    100020535-2      100020735      100020935      100022235 
    ##              1              1              1              1              1 
    ##      100022335      100023235      100023535      100023935    100025335-2 
    ##              1              1              1              1              1 
    ##      100025435    100025835-2      100027035      100027235      100027335 
    ##              1              1              1              1              1 
    ##    100027535-2      100027635    100028435-2      100029235      100029435 
    ##              1              1              1              1              1 
    ##      100029935      100030635      100032635      100032835      100032935 
    ##              1              1              1              1              1 
    ##      100034635      100035535      100038235      100039235      100045735 
    ##              1              1              1              1              1 
    ##      100046235      100046335      100046635      100047235      100047835 
    ##              1              1              1              1              1 
    ##      100048035      100048335      100048635      100049835      100050435 
    ##              1              1              1              1              1 
    ##      100050535      100050735      100051235      100051835      100053535 
    ##              1              1              1              1              1 
    ##      100053635      100054435      100055235      100056135      100056435 
    ##              1              1              1              1              1 
    ##      100057235      100057535      100058435      100058535      100058935 
    ##              1              1              1              1              1 
    ##      100059035      100059635      100060135      100060235      100062535 
    ##              1              1              1              1              1 
    ##      100063035    100063135-2      100063935      100064135      100064635 
    ##              1              1              1              1              1 
    ##      100065835      100066235      100067235      100068635      100070435 
    ##              1              1              1              1              1 
    ##      100072435      100072735    100073835-2      100073935      100074535 
    ##              1              1              1              1              1 
    ##      100074735      100074935      100076135      100076935      100078035 
    ##              1              1              1              1              1 
    ##      100078835      100079135      100079835      100080235      100080835 
    ##              1              1              1              1              1 
    ##      100082035      100084635      100085135      100085435      100089135 
    ##              1              1              1              1              1 
    ##      100089435      100090035      100090235      100092235      100092535 
    ##              1              1              1              1              1 
    ##      100095035      100097335      100098835      100099635      100102435 
    ##              1              1              1              1              1 
    ##      100102535      100106735      100107235      100107435      100107535 
    ##              1              1              1              1              1 
    ##      100108435      100113135      100113435      100113535      100113735 
    ##              1              1              1              1              1 
    ##      100115135      100117135      100117235      100118535      100119335 
    ##              1              1              1              1              1 
    ##      100120135      100120435      100121535      100122735      100123435 
    ##              1              1              1              1              1 
    ##      100123535      100124135      100126435      100126635      100127135 
    ##              1              1              1              1              1 
    ##      100127235      100128835      100136635      100136935      100137735 
    ##              1              1              1              1              1 
    ##      100138335      100139435      100139735      100140935      100141235 
    ##              1              1              1              1              1 
    ##      100141435      100141935      100142635      100144435      100146735 
    ##              1              1              1              1              1 
    ##      100146935      100147035      100147535      100147835      100151935 
    ##              1              1              1              1              1 
    ##      100153535      100153935      100154835      100156735      100157235 
    ##              1              1              1              1              1 
    ##      100160535    100160935-2      100161635      100161935      100162035 
    ##              1              1              1              1              1 
    ##      100163935      100164035      100164235      100164435      100164635 
    ##              1              1              1              1              1 
    ##      100165535      100165935      100166135    100166235-2      100166335 
    ##              1              1              1              1              1 
    ##      100166835      100167935      100168235      100168835      100169835 
    ##              1              1              1              1              1 
    ##      100170135      100170235      100171335    100178935-2      100182835 
    ##              1              1              1              1              1 
    ##      100182935    100186835-2      100192635      100194535      100195835 
    ##              1              1              1              1              1 
    ##      100196135      100200835      100201235      100202135      100203535 
    ##              1              1              1              1              1 
    ##      100203835      100206535      100220135      100220235      100220835 
    ##              1              1              1              1              1 
    ##      100222835      100223535      100225635      100228235    100230735-2 
    ##              1              1              1              1              1 
    ##      100231935      100233435      100234435      100236935      100238235 
    ##              1              1              1              1              1 
    ##      100240135      100241235      100241735      100241935      100242035 
    ##              1              1              1              1              1 
    ##      100242635      100243935      100245135      100245635      100253435 
    ##              1              1              1              1              1 
    ##      100253835      100254335      100254835      100255235      100260835 
    ##              1              1              1              1              1 
    ##      100261735      100262335      100263235      100264035      100266335 
    ##              1              1              1              1              1 
    ##      100267935      100268035      100268235      100268735      100269835 
    ##              1              1              1              1              1 
    ##      100272935      100274135      100277935      100279535      100291735 
    ##              1              1              1              1              1 
    ##      100299535      100299835      100300135      100303935      100306235 
    ##              1              1              1              1              1 
    ##      100306635      100306735      100306835      100308035      100309235 
    ##              1              1              1              1              1 
    ##      100315735      100316935      100318135      100318435      100323935 
    ##              1              1              1              1              1 
    ##      100325435      100329635      100332235      100332735      100334635 
    ##              1              1              1              1              1 
    ##      100335635      100336235      100336435      100336935      100339235 
    ##              1              1              1              1              1 
    ##      100341835      100342535      100342635      100343135      100343835 
    ##              1              1              1              1              1 
    ##      100344835      100347935      100348235      100350635      100352535 
    ##              1              1              1              1              1 
    ##      100352935      100354235      100356135      100358035      100360335 
    ##              1              1              1              1              1 
    ##      100364935      100366535      100367035      100367135      100367935 
    ##              1              1              1              1              1 
    ##      100371235      100371535      100378535      100379335      100379635 
    ##              1              1              1              1              1 
    ##      100380935      100381735      100382135      100383335      100383835 
    ##              1              1              1              1              1 
    ##      100383935      100385235      100386035      100386135      100387935 
    ##              1              1              1              1              1 
    ##      100389235      100389335      100391035      100391335      100391435 
    ##              1              1              1              1              1 
    ##      100391735      100392335      100393435      100394335      100395535 
    ##              1              1              1              1              1 
    ##      100396135      100397735      100398735      100399135      100399835 
    ##              1              1              1              1              1 
    ##      100402635      100403035      100403635      100406135      100406235 
    ##              1              1              1              1              1 
    ##      100406335      100408335      100408735      100409435      100410135 
    ##              1              1              1              1              1 
    ##      100414335      100414835      100417735      100417935      100419335 
    ##              1              1              1              1              1 
    ##      100422835      100425235      100426535      100427735      100430335 
    ##              1              1              1              1              1 
    ##      100431135      100432235      100432735      100432935      100434735 
    ##              1              1              1              1              1 
    ##      100439235      100439635      100440135      100443735      100443935 
    ##              1              1              1              1              1 
    ##      100444935      100446035      100446135      100446235      100448135 
    ##              1              1              1              1              1 
    ##      100448735      100449335      100450035      100450935      100451035 
    ##              1              1              1              1              1 
    ##      100451235      100451335      100452035      100453335      100453535 
    ##              1              1              1              1              1 
    ##      100453635      100454035      100456535      100460035      100460135 
    ##              1              1              1              1              1 
    ##      100461535      100463335      100463535      100464235      100465335 
    ##              1              1              1              1              1 
    ##      100466235      100467435      100467635      100467835      100468035 
    ##              1              1              1              1              1 
    ##      100468335      100468835      100469635      100470835      100470935 
    ##              1              1              1              1              1 
    ##      100471135      100473135      100473335      100473635      100475435 
    ##              1              1              1              1              1 
    ##      100475535      100476535      100476835    100477535-2    100477835-2 
    ##              1              1              1              1              1 
    ##      100479535      100479635      100480135      100480535    100480735-2 
    ##              1              1              1              1              1 
    ##      100481535      100481935      100483435      100483835      100485035 
    ##              1              1              1              1              1 
    ##      100485635      100485835      100486235    100489435-2    100490435-2 
    ##              1              1              1              1              1 
    ##      100490535      100491335      100492135      100494035      100494435 
    ##              1              1              1              1              1 
    ##      100494735      100494835      100496135      100496735      100496935 
    ##              1              1              1              1              1 
    ##      100500835      100501135      100501835      100502535      100503635 
    ##              1              1              1              1              1 
    ##      100505735      100506235      100507835      100509835      100511135 
    ##              1              1              1              1              1 
    ##      100512035    100513035-2      100513335      100515135      100517335 
    ##              1              1              1              1              1 
    ##      100518835      100519235      100524435      100525835      100527635 
    ##              1              1              1              1              1 
    ##      100528435      100528735      100531335      100561135      100563435 
    ##              1              1              1              1              1 
    ##      100563735      100564135      100565235      100565435      100566335 
    ##              1              1              1              1              1 
    ##      100566435      100566535      100566835      100567435      100567835 
    ##              1              1              1              1              1 
    ##      100569235      100570035    100570635-2      100572435      100572835 
    ##              1              1              1              1              1 
    ##      100573735      100575835      100576735      100577235      100577335 
    ##              1              1              1              1              1 
    ##      100578435      100579635      100579935      100580135      100580735 
    ##              1              1              1              1              1 
    ##      100580835      100581035      100582235      100582335      100582435 
    ##              1              1              1              1              1 
    ##      100583635      100583735      100584035      100584435      100584735 
    ##              1              1              1              1              1 
    ##      100585635      100586735      100587635      100587835    100588635-2 
    ##              1              1              1              1              1 
    ##      100588735      100589235      100589435      100589535      100589735 
    ##              1              1              1              1              1 
    ##      100591235      100591335      100592235      100593235      100594035 
    ##              1              1              1              1              1 
    ##      100594135      100595535      100595835      100595935      100597835 
    ##              1              1              1              1              1 
    ##      100597935      100598035      100598935      100599335      100599435 
    ##              1              1              1              1              1 
    ##      100599535      100599835      100601135      100602835      100603435 
    ##              1              1              1              1              1 
    ##      100603635      100604235      100604635      100604935      100605535 
    ##              1              1              1              1              1 
    ##      100607035      100607935      100608135      100608235      100609535 
    ##              1              1              1              1              1 
    ##      100610835      100631635      100636135      100637835      100638035 
    ##              1              1              1              1              1 
    ##      100638135      100638735      100638835      100638935      100639535 
    ##              1              1              1              1              1 
    ##      100639735      100640035      100641235      100641435      100680835 
    ##              1              1              1              1              1 
    ##      100681035      100681435      100681935      100682435      100683735 
    ##              1              1              1              1              1 
    ##      100685035      100685135      100685835      100687935      100691935 
    ##              1              1              1              1              1 
    ##    100692735-2   100693835-2W      100695335      100695735      100695935 
    ##              1              1              1              1              1 
    ##      100696535      100696835      100697835      100700535      100702835 
    ##              1              1              1              1              1 
    ##      100703535      100704235    100705635-2      100705935      100706235 
    ##              1              1              1              1              1 
    ##    100706535-2      100707535      100707735      100708935    100711435-2 
    ##              1              1              1              1              1 
    ##      100712235      100713435      100713735      100713835    100715135-2 
    ##              1              1              1              1              1 
    ##      100715935      100718935      100723535      100729435      100730535 
    ##              1              1              1              1              1 
    ##      100730935      100731935      100732135      100732335      100733235 
    ##              1              1              1              1              1 
    ##      100733935      100734235      100735035      100735435      100735535 
    ##              1              1              1              1              1 
    ##      100736535      100740635      100742635      100743335      100743635 
    ##              1              1              1              1              1 
    ##      100743835      100744235      100744735      100745035      100746535 
    ##              1              1              1              1              1 
    ##      100747035      100749535      100749835      100752035      100752535 
    ##              1              1              1              1              1 
    ##      100752735      100753935      100755435      100755635      100756235 
    ##              1              1              1              1              1 
    ##      100756735      100757835      100759735      100759835      100763635 
    ##              1              1              1              1              1 
    ##      100765935      100766835      100768635      100769235      100769535 
    ##              1              1              1              1              1 
    ##      100769735      100771135      100771835      100774835      100774935 
    ##              1              1              1              1              1 
    ##      100775035      100776335      100776435   100778235-2W      100778735 
    ##              1              1              1              1              1 
    ##      100787735      100790335      100790935      100791135      100792035 
    ##              1              1              1              1              1 
    ##      100792935      100793235      100793335   100793835-2W      100796535 
    ##              1              1              1              1              1 
    ##      100797735      100797835      100798935      100800735      100800835 
    ##              1              1              1              1              1 
    ##      100801335      100802435      100804035      100805135      100807235 
    ##              1              1              1              1              1 
    ##      100807935      100809735      100811035      100811935      100814335 
    ##              1              1              1              1              1 
    ##      100815535      100816935      100818835      100819135      100825635 
    ##              1              1              1              1              1 
    ##      100829735      100830535      100831735      100834435    100835235-2 
    ##              1              1              1              1              1 
    ##      100835935      100837235      100837735      100837835      100838435 
    ##              1              1              1              1              1 
    ##      100840435      100841835      100841935      100842035      100842335 
    ##              1              1              1              1              1 
    ##      100842835      100843635      100844335      100844435      100847435 
    ##              1              1              1              1              1 
    ##      100850935      100852035      100852535      100852835      100853235 
    ##              1              1              1              1              1 
    ##      100853335      100853735      100854535      100855035      100859135 
    ##              1              1              1              1              1 
    ##      100860735      100868235      100872435      100873635      100873735 
    ##              1              1              1              1              1 
    ##      100875835      100876135      100876535      100879535      100881235 
    ##              1              1              1              1              1 
    ##      100883535      100883835      100885335    100886035-2      100886435 
    ##              1              1              1              1              1 
    ##      100887035      100887135      100887835      100887935      100888735 
    ##              1              1              1              1              1 
    ##      100889235      100889435    100891235-2    100893735-2      100893935 
    ##              1              1              1              1              1 
    ##      100945035      100946235      100948735      100960035      100960235 
    ##              1              1              1              1              1 
    ##      100960735      100961035      100961235      100963435      100964035 
    ##              1              1              1              1              1 
    ##      100964235      100964735      100967735      100969835      100973335 
    ##              1              1              1              1              1 
    ##    100974035-2      100975235      100976935      100986635      100986735 
    ##              1              1              1              1              1 
    ##      100988435      100989235      100990935      100991335      100992835 
    ##              1              1              1              1              1 
    ##      100994135      100995935      100996135      100996935      100997135 
    ##              1              1              1              1              1 
    ##      100997235      100998135      100998535      100998935      100999035 
    ##              1              1              1              1              1 
    ##      101003535      101005435      101007135      101008735      101011735 
    ##              1              1              1              1              1 
    ##    101012135-2      101012535      101014135      101015935      101017535 
    ##              1              1              1              1              1 
    ##      101019535      101020535      101021735      101022135      101023835 
    ##              1              1              1              1              1 
    ##      101023935      101024335      101024835      101025735      101027735 
    ##              1              1              1              1              1 
    ##      101030335      101032235      101033035      101033635      101035935 
    ##              1              1              1              1              1 
    ##      101036035      101036335      101037435      101038735    101041235-2 
    ##              1              1              1              1              1 
    ##      101043135      101043635      101096535      101096735      101097235 
    ##              1              1              1              1              1 
    ##      101097535    101099635-2      101099835    101100135-2      101100235 
    ##              1              1              1              1              1 
    ##      101103735      101104335      101105235      101105735      101105835 
    ##              1              1              1              1              1 
    ##      101108635      101108735      101109635      101109935      101110035 
    ##              1              1              1              1              1 
    ##      101110135      101111035      101111135      101112435      101112935 
    ##              1              1              1              1              1 
    ##      101113335      101114335      101114735      101119135      101119235 
    ##              1              1              1              1              1 
    ##      101120235    101120335-2      101120935      101121335      101121435 
    ##              1              1              1              1              1 
    ##      101122935      101123535      101159335      101159835      101160835 
    ##              1              1              1              1              1 
    ##      101161435      101164235      101165335      101166435      101170635 
    ##              1              1              1              1              1 
    ##      101170835      101173135    101233635-2      101234335      101234735 
    ##              1              1              1              1              1 
    ##      101237135      101237535      101237935      101240935      101241835 
    ##              1              1              1              1              1 
    ##      101241935      101242635      101243035      101247035      101248435 
    ##              1              1              1              1              1 
    ##      101249635      101256335      101260935      101262335      101262435 
    ##              1              1              1              1              1 
    ##      101264135      101265035      101265235      101266635      101268435 
    ##              1              1              1              1              1 
    ##      101268535      101274335      101275435      101278235      101280735 
    ##              1              1              1              1              1 
    ##      101281235      101282235      101282435      101285035      101286435 
    ##              1              1              1              1              1 
    ##      101288335      101290735      101291935    101294535-2      101295335 
    ##              1              1              1              1              1 
    ##    101295735-2    101296535-2    101296935-2      101302735      101303635 
    ##              1              1              1              1              1 
    ##      101304735      101321135      101324935      101334035      101336635 
    ##              1              1              1              1              1 
    ##      101341635      101342635      101344035      101349335      101350135 
    ##              1              1              1              1              1 
    ##      101350435      101352235      101352835      101354835      101355135 
    ##              1              1              1              1              1 
    ##      101357035      101357635      101359835      101361835      101369035 
    ##              1              1              1              1              1 
    ##      101370535      101371735      101372235      101379435      101379535 
    ##              1              1              1              1              1 
    ##    101380335-2      101381135      101381335    101381835-2      101382735 
    ##              1              1              1              1              1 
    ##      101383035    101383435-2      101387135      101391135      101399135 
    ##              1              1              1              1              1 
    ##      101399335    101399835-2      101401035      101401635      101401735 
    ##              1              1              1              1              1 
    ##      101401935    101402135-2      101402335      101493835      101494035 
    ##              1              1              1              1              1 
    ##      101494135      101497235      101504635      101504735      101505335 
    ##              1              1              1              1              1 
    ##      101505535      101505835      101506335      101506535      101506635 
    ##              1              1              1              1              1 
    ##      101506935      101507435      101508635      101509735      101510735 
    ##              1              1              1              1              1 
    ##      101511135      101511935      101512435      101513135      101514335 
    ##              1              1              1              1              1 
    ##      101514435      101515335      101515535      101516435      101517035 
    ##              1              1              1              1              1 
    ##      101517735      101517835      101518435      101520235      101520935 
    ##              1              1              1              1              1 
    ##      101521135      101522435      101524335      101524535      101525235 
    ##              1              1              1              1              1 
    ##      101527635      101529335      101530735      101531035      101531635 
    ##              1              1              1              1              1 
    ##      101532835      101533335      101535435      101535835      101536335 
    ##              1              1              1              1              1 
    ##      101536835      101537335      101538235      101538435      101546235 
    ##              1              1              1              1              1 
    ##      101547535      101547935      101549035      101551135      101555735 
    ##              1              1              1              1              1 
    ##      101558435      101558635      101559335      101560035      101560135 
    ##              1              1              1              1              1 
    ##    101561635-2      101564635      101564935      101567135      101567335 
    ##              1              1              1              1              1 
    ##      101570335      101570535      101573235      101573735      101574835 
    ##              1              1              1              1              1 
    ##      101576435      101576535      101576835      101579335      101581135 
    ##              1              1              1              1              1 
    ##      101584435      101588135      101591535      101593235      101593635 
    ##              1              1              1              1              1 
    ##      101594035      101634535      101635835      101637135      101637635 
    ##              1              1              1              1              1 
    ##      101638735      101639235      101640335      101642135      101643035 
    ##              1              1              1              1              1 
    ##      101643835      101644535      101645235      101648035      101669135 
    ##              1              1              1              1              1 
    ##      101670335      101670535      101672435      101673735      101673835 
    ##              1              1              1              1              1 
    ##      101691435      101691535      101691935      101692735      101699535 
    ##              1              1              1              1              1 
    ##      101700635      101700835      101701035      101701735    101702135-2 
    ##              1              1              1              1              1 
    ##      101710135      101710935      101711635      101712135      101712335 
    ##              1              1              1              1              1 
    ##      101712935      101713535      101714735      101716135      101779732 
    ##              1              1              1              1              1 
    ##      101779932      101780132      101780832      101781632      101781832 
    ##              1              1              1              1              1 
    ##      101781932      101782232      101783232      101784132      101784432 
    ##              1              1              1              1              1 
    ##      101859435      101859635      101861635      101862335      101862535 
    ##              1              1              1              1              1 
    ##      101865235      101865535      101866735      101867135      101867335 
    ##              1              1              1              1              1 
    ##      101867935      101868035    101875235-2      101877435    101877935-2 
    ##              1              1              1              1              1 
    ##      101880035      101881135      101882235      101882435      101882635 
    ##              1              1              1              1              1 
    ##      101883835    101884435-2      101884835      101887635      101887735 
    ##              1              1              1              1              1 
    ##      101889035      101889335      101890235      101891735      101894235 
    ##              1              1              1              1              1 
    ##      101894635      101895135      101897335      101898335      101899035 
    ##              1              1              1              1              1 
    ##      101899235      101901035      101901335      101901735      101904035 
    ##              1              1              1              1              1 
    ##      101904935      101905035      101906335      101906435      101907035 
    ##              1              1              1              1              1 
    ##      101909635      101910435      101912335      101912635      101912935 
    ##              1              1              1              1              1 
    ##      101913135      101913735      101914535      101914835      101915235 
    ##              1              1              1              1              1 
    ##      101916435      101917635      101918535      101919135      101920135 
    ##              1              1              1              1              1 
    ##      101920335    101921235-2      101921535    101921835-2   101922035-2W 
    ##              1              1              1              1              1 
    ##    101922135-2      101922935    101923635-2    101923935-2      101924335 
    ##              1              1              1              1              1 
    ##      101926135      101928235      101929835      101932035      101932635 
    ##              1              1              1              1              1 
    ##      101932835      101933035      101933635      101933835      101934035 
    ##              1              1              1              1              1 
    ##      101935535      101938235    101939335-2      101941035      101943235 
    ##              1              1              1              1              1 
    ##      101944035      101944435      101946435    101947335-2      101947835 
    ##              1              1              1              1              1 
    ##      101948335      101948835      101949035      101952635    101953135-2 
    ##              1              1              1              1              1 
    ##      101956235      101960235      101962635      101963435      101970935 
    ##              1              1              1              1              1 
    ##      101972635      101972835      101973435      101974135      101974335 
    ##              1              1              1              1              1 
    ##      101974535      101977035      101977935      101980435      101981135 
    ##              1              1              1              1              1 
    ##      101981335      101981835      101982235      101983635      101984635 
    ##              1              1              1              1              1 
    ##      101984835      101985635      101985935      101986035      101986235 
    ##              1              1              1              1              1 
    ##      101987435      101987635      101988535      101989335      101991735 
    ##              1              1              1              1              1 
    ##      101993035      101994635      101995035      101995835      101997335 
    ##              1              1              1              1              1 
    ##      101997435      101997735      101998235      101999835      102000935 
    ##              1              1              1              1              1 
    ##      102001935      102003035      102003635      102003835      102003935 
    ##              1              1              1              1              1 
    ##    102004335-2      102005035      102006335      102006435      102008035 
    ##              1              1              1              1              1 
    ##      102010035      102010235    102010335-2      102013935      102014135 
    ##              1              1              1              1              1 
    ##      102014935      102015135      102015335      102015435      102016735 
    ##              1              1              1              1              1 
    ##      102017135      102017735      102018035      102019835      102020235 
    ##              1              1              1              1              1 
    ##      102020335      102021435      102023035      102025635      102026035 
    ##              1              1              1              1              1 
    ##      102026235      102026535      102027535      102029335      102029635 
    ##              1              1              1              1              1 
    ##      102030035      102030335      102031335      102031735      102032335 
    ##              1              1              1              1              1 
    ##      102032935      102033035      102033235      102033435      102036335 
    ##              1              1              1              1              1 
    ##      102037035      102037535      102038735      102039935      102040535 
    ##              1              1              1              1              1 
    ##      102040635      102040735      102042935      102046735      102049535 
    ##              1              1              1              1              1 
    ##      102050035      102050735      102050835      102057035      102058235 
    ##              1              1              1              1              1 
    ##      102063535      102064035      102064735      102066435      102066935 
    ##              1              1              1              1              1 
    ##      102068435      102068535      102068935      102069735      102069835 
    ##              1              1              1              1              1 
    ##      102070635      102076335    102076635-2    102077035-2    102077235-2 
    ##              1              1              1              1              1 
    ##      102078235      102078635      102081235      102082235      102084235 
    ##              1              1              1              1              1 
    ##      102084735      102087935      102088235      102088435      102088735 
    ##              1              1              1              1              1 
    ##      102088835      102089235      102089935      102090035      102090635 
    ##              1              1              1              1              1 
    ##      102090935      102091235    102091335-2      102091635      102096335 
    ##              1              1              1              1              1 
    ##      102097035      102097335      102099635      102099735      102100135 
    ##              1              1              1              1              1 
    ##      102100535      102101235      102103235      102103835      102103935 
    ##              1              1              1              1              1 
    ##      102104335      102104835      102106335      102109435      102111135 
    ##              1              1              1              1              1 
    ##    102112535-2      102112835      102113235      102115335      102116435 
    ##              1              1              1              1              1 
    ##      102117235      102118735      102119135      102119335      102123935 
    ##              1              1              1              1              1 
    ##      102125035      102125135    102125535-2      102125635      102125935 
    ##              1              1              1              1              1 
    ##      102126035      102127535      102129035      102129335      102129735 
    ##              1              1              1              1              1 
    ##      102131135      102132335      102133335      102133535      102136635 
    ##              1              1              1              1              1 
    ##      102137335      102139035      102139935      102140135      102140635 
    ##              1              1              1              1              1 
    ##      102141635      102142135      102145335      102145935      102147835 
    ##              1              1              1              1              1 
    ##      102148035      102148235      102148835      102149035      102149435 
    ##              1              1              1              1              1 
    ##      102150335    102151035-2      102151335      102152135      102153135 
    ##              1              1              1              1              1 
    ##      102154035    102155135-2    102155335-2      102156735    102157535-2 
    ##              1              1              1              1              1 
    ##      102157635      102158335    102158835-2      102159635      102160135 
    ##              1              1              1              1              1 
    ##      102160435      102160635      102161635      102161735      102163435 
    ##              1              1              1              1              1 
    ##      102166535      102166835      102170035      102171735      102172335 
    ##              1              1              1              1              1 
    ##      102174135      102174735      102175735      102177635      102178135 
    ##              1              1              1              1              1 
    ##      102179335      102181435      102182035      102183735      102186235 
    ##              1              1              1              1              1 
    ##      102187235      102187335      102189035      102190035      102192035 
    ##              1              1              1              1              1 
    ##    102193835-2    102196135-2      102198435      102198535      102198735 
    ##              1              1              1              1              1 
    ##      102199935      102200635      102202135      102206735      102208035 
    ##              1              1              1              1              1 
    ##      102208735      102209535      102210535      102213135      102216135 
    ##              1              1              1              1              1 
    ##      102217635      102218635      102222035      102223635      102224435 
    ##              1              1              1              1              1 
    ##      102224735      102225335      102226135      102226735      102226935 
    ##              1              1              1              1              1 
    ##      102227435      102230035      102230635      102233135      102234935 
    ##              1              1              1              1              1 
    ##      102236035      102236235      102237035      102237835      102238935 
    ##              1              1              1              1              1 
    ##      102243935      102244335      102244735      102246535      102246635 
    ##              1              1              1              1              1 
    ##      102247335      102248935      102249935      102250935      102252035 
    ##              1              1              1              1              1 
    ##      102253535      102253935      102258435    102258535-2      102258935 
    ##              1              1              1              1              1 
    ##    102263135-2      102264635      102265335      102267135      102267335 
    ##              1              1              1              1              1 
    ##      102268735      102269035      102269635      102270135      102270335 
    ##              1              1              1              1              1 
    ##      102271035      102271235      102271535      102271835      102272435 
    ##              1              1              1              1              1 
    ##      102272535      102272935      102273235      102273535      102275935 
    ##              1              1              1              1              1 
    ##      102276435      102277935      102278035      102279835      102281335 
    ##              1              1              1              1              1 
    ##      102282535      102283735      102284135      102285435      102286335 
    ##              1              1              1              1              1 
    ##      102286735      102287035      102287135      102288235      102288435 
    ##              1              1              1              1              1 
    ##      102288735      102291735      102292535      102294235      102295835 
    ##              1              1              1              1              1 
    ##      102296435      102299135      102299935      102300235      102300335 
    ##              1              1              1              1              1 
    ##      102300535      102302035      102302335      102302935      102303835 
    ##              1              1              1              1              1 
    ##      102303935      102304535      102307635      102311535      102315135 
    ##              1              1              1              1              1 
    ##      102319135      102319335      102321035      102321835      102322135 
    ##              1              1              1              1              1 
    ##      102324535      102326335      102327535      102329135      102329935 
    ##              1              1              1              1              1 
    ##      102330635      102331635      102334335      102334735      102337635 
    ##              1              1              1              1              1 
    ##      102339235      102339435      102341335      102341635      102343635 
    ##              1              1              1              1              1 
    ##      102346535      102348535      102349535      102349635      102350835 
    ##              1              1              1              1              1 
    ##      102350935      102352035      102352935      102353035      102353335 
    ##              1              1              1              1              1 
    ##      102353435      102353535      102353635      102355535      102357735 
    ##              1              1              1              1              1 
    ##      102358335      102361035      102361235      102362435      102363535 
    ##              1              1              1              1              1 
    ##      102364635      102365035      102366635      102366735      102368635 
    ##              1              1              1              1              1 
    ##      102369835      102372335      102373235      102374435      102375335 
    ##              1              1              1              1              1 
    ##      102376035      102377635      102380135    102380935-2      102381135 
    ##              1              1              1              1              1 
    ##      102382835      102384135      102384235      102384535      102386135 
    ##              1              1              1              1              1 
    ##      102386935      102389035      102389135    102389535-2      102394535 
    ##              1              1              1              1              1 
    ##      102396835      102400035    102400835-2    102401035-2      102401335 
    ##              1              1              1              1              1 
    ##    102401635-2   102401735-2W      102403335      102403535      102404035 
    ##              1              1              1              1              1 
    ##      102404135      102404235      102405635      102412135      102412235 
    ##              1              1              1              1              1 
    ##    102413135-2      102489235      102490535      102493235      102493335 
    ##              1              1              1              1              1 
    ##      102495535      102496535      102496735      102497135      102497235 
    ##              1              1              1              1              1 
    ##      102497835      102498235      102499235      102500935      102502135 
    ##              1              1              1              1              1 
    ##      102502435      102503335      102509935      102510035      102510135 
    ##              1              1              1              1              1 
    ##      102510535      102512335      102512735      102513435      102514935 
    ##              1              1              1              1              1 
    ##      102517335      102518235      102524135      102525335      102525835 
    ##              1              1              1              1              1 
    ##      102525935      102527935      102529935      102535735      102536035 
    ##              1              1              1              1              1 
    ##      102545135      102546335      102546935      102549135      102550435 
    ##              1              1              1              1              1 
    ##      102550735      102551035      102551135      102551935      102552735 
    ##              1              1              1              1              1 
    ##   102559435-2W    102559835-2      102560135      102560335      102560435 
    ##              1              1              1              1              1 
    ##      102560935    102562935-2      102564835      102565135      102565235 
    ##              1              1              1              1              1 
    ##      102565435      102566235      102566735      102566835      102566935 
    ##              1              1              1              1              1 
    ##      102567135    102568135-2      102569335      102569635      102571135 
    ##              1              1              1              1              1 
    ##      102571735   102579835-2W      102580135      102581435      102583035 
    ##              1              1              1              1              1 
    ##      102584735      102585335      102588335      102588735      102589735 
    ##              1              1              1              1              1 
    ##      102589835      102590135      102592035      102592435      102592635 
    ##              1              1              1              1              1 
    ##      102604635      102607535      102607635      102608935    102615035-2 
    ##              1              1              1              1              1 
    ##      102615535      102615635    102615935-2      102618235      102629335 
    ##              1              1              1              1              1 
    ##      102630135      102630635      102631535      102640435      102649735 
    ##              1              1              1              1              1 
    ##      102650235      102650735      102651335      102659535      102659735 
    ##              1              1              1              1              1 
    ##      102660335      102662435      102662635 SAMEA115126061 SAMEA115126107 
    ##              1              1              1              1              1 
    ## SAMEA115126171 SAMEA115126306 SAMEA115126326 SAMEA115126491 SAMEA115126500 
    ##              1              1              1              1              1 
    ## SAMEA115126520 SAMEA115126613 SAMEA115126615 SAMEA115126635 SAMEA115126710 
    ##              1              1              1              1              1 
    ## SAMEA115126717 SAMEA115126775 SAMEA115126863 SAMEA115126870 SAMEA115126888 
    ##              1              1              1              1              1 
    ## SAMEA115127046 SAMEA115127059 SAMEA115127165 SAMEA115127331 SAMEA115127484 
    ##              1              1              1              1              1 
    ## SAMEA115127502 SAMEA115127533 SAMEA115127719 
    ##              1              1              1

``` r
table(metadata$Site)
```

    ## 
    ## saliva 
    ##   1698

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
    ##      100000135      100000735      100002335      100002535      100002935 
    ##              1              1              1              1              1 
    ##      100003235      100004035      100004535      100005135      100005235 
    ##              1              1              1              1              1 
    ##      100005435      100005735      100006235      100006935      100007135 
    ##              1              1              1              1              1 
    ##      100007635      100008335      100009635      100012935      100013735 
    ##              1              1              1              1              1 
    ##      100014135      100014235      100014435      100014635      100014935 
    ##              1              1              1              1              1 
    ##      100015235      100015435      100015535      100015835      100016235 
    ##              1              1              1              1              1 
    ##      100016435      100016535      100017835      100017935      100018635 
    ##              1              1              1              1              1 
    ##      100019535    100020535-2      100020735      100020935      100022235 
    ##              1              1              1              1              1 
    ##      100022335      100023235      100023535      100023935    100025335-2 
    ##              1              1              1              1              1 
    ##      100025435    100025835-2      100027035      100027235      100027335 
    ##              1              1              1              1              1 
    ##    100027535-2      100027635    100028435-2      100029235      100029435 
    ##              1              1              1              1              1 
    ##      100029935      100030635      100032635      100032835      100032935 
    ##              1              1              1              1              1 
    ##      100034635      100035535      100038235      100039235      100045735 
    ##              1              1              1              1              1 
    ##      100046235      100046335      100046635      100047235      100047835 
    ##              1              1              1              1              1 
    ##      100048035      100048335      100048635      100049835      100050435 
    ##              1              1              1              1              1 
    ##      100050535      100050735      100051235      100051835      100053535 
    ##              1              1              1              1              1 
    ##      100053635      100054435      100055235      100056135      100056435 
    ##              1              1              1              1              1 
    ##      100057235      100057535      100058435      100058535      100058935 
    ##              1              1              1              1              1 
    ##      100059035      100059635      100060135      100060235      100062535 
    ##              1              1              1              1              1 
    ##      100063035    100063135-2      100063935      100064135      100064635 
    ##              1              1              1              1              1 
    ##      100065835      100066235      100067235      100068635      100070435 
    ##              1              1              1              1              1 
    ##      100072435      100072735    100073835-2      100073935      100074535 
    ##              1              1              1              1              1 
    ##      100074735      100074935      100076135      100076935      100078035 
    ##              1              1              1              1              1 
    ##      100078835      100079135      100079835      100080235      100080835 
    ##              1              1              1              1              1 
    ##      100082035      100084635      100085135      100085435      100089135 
    ##              1              1              1              1              1 
    ##      100089435      100090035      100090235      100092235      100092535 
    ##              1              1              1              1              1 
    ##      100095035      100097335      100098835      100099635      100102435 
    ##              1              1              1              1              1 
    ##      100102535      100106735      100107235      100107435      100107535 
    ##              1              1              1              1              1 
    ##      100108435      100113135      100113435      100113535      100113735 
    ##              1              1              1              1              1 
    ##      100115135      100117135      100117235      100118535      100119335 
    ##              1              1              1              1              1 
    ##      100120135      100120435      100121535      100122735      100123435 
    ##              1              1              1              1              1 
    ##      100123535      100124135      100126435      100126635      100127135 
    ##              1              1              1              1              1 
    ##      100127235      100128835      100136635      100136935      100137735 
    ##              1              1              1              1              1 
    ##      100138335      100139435      100139735      100140935      100141235 
    ##              1              1              1              1              1 
    ##      100141435      100141935      100142635      100144435      100146735 
    ##              1              1              1              1              1 
    ##      100146935      100147035      100147535      100147835      100151935 
    ##              1              1              1              1              1 
    ##      100153535      100153935      100154835      100156735      100157235 
    ##              1              1              1              1              1 
    ##      100160535    100160935-2      100161635      100161935      100162035 
    ##              1              1              1              1              1 
    ##      100163935      100164035      100164235      100164435      100164635 
    ##              1              1              1              1              1 
    ##      100165535      100165935      100166135    100166235-2      100166335 
    ##              1              1              1              1              1 
    ##      100166835      100167935      100168235      100168835      100169835 
    ##              1              1              1              1              1 
    ##      100170135      100170235      100171335    100178935-2      100182835 
    ##              1              1              1              1              1 
    ##      100182935    100186835-2      100192635      100194535      100195835 
    ##              1              1              1              1              1 
    ##      100196135      100200835      100201235      100202135      100203535 
    ##              1              1              1              1              1 
    ##      100203835      100206535      100220135      100220235      100220835 
    ##              1              1              1              1              1 
    ##      100222835      100223535      100225635      100228235    100230735-2 
    ##              1              1              1              1              1 
    ##      100231935      100233435      100234435      100236935      100238235 
    ##              1              1              1              1              1 
    ##      100240135      100241235      100241735      100241935      100242035 
    ##              1              1              1              1              1 
    ##      100242635      100243935      100245135      100245635      100253435 
    ##              1              1              1              1              1 
    ##      100253835      100254335      100254835      100255235      100260835 
    ##              1              1              1              1              1 
    ##      100261735      100262335      100263235      100264035      100266335 
    ##              1              1              1              1              1 
    ##      100267935      100268035      100268235      100268735      100269835 
    ##              1              1              1              1              1 
    ##      100272935      100274135      100277935      100279535      100291735 
    ##              1              1              1              1              1 
    ##      100299535      100299835      100300135      100303935      100306235 
    ##              1              1              1              1              1 
    ##      100306635      100306735      100306835      100308035      100309235 
    ##              1              1              1              1              1 
    ##      100315735      100316935      100318135      100318435      100323935 
    ##              1              1              1              1              1 
    ##      100325435      100329635      100332235      100332735      100334635 
    ##              1              1              1              1              1 
    ##      100335635      100336235      100336435      100336935      100339235 
    ##              1              1              1              1              1 
    ##      100341835      100342535      100342635      100343135      100343835 
    ##              1              1              1              1              1 
    ##      100344835      100347935      100348235      100350635      100352535 
    ##              1              1              1              1              1 
    ##      100352935      100354235      100356135      100358035      100360335 
    ##              1              1              1              1              1 
    ##      100364935      100366535      100367035      100367135      100367935 
    ##              1              1              1              1              1 
    ##      100371235      100371535      100378535      100379335      100379635 
    ##              1              1              1              1              1 
    ##      100380935      100381735      100382135      100383335      100383835 
    ##              1              1              1              1              1 
    ##      100383935      100385235      100386035      100386135      100387935 
    ##              1              1              1              1              1 
    ##      100389235      100389335      100391035      100391335      100391435 
    ##              1              1              1              1              1 
    ##      100391735      100392335      100393435      100394335      100395535 
    ##              1              1              1              1              1 
    ##      100396135      100397735      100398735      100399135      100399835 
    ##              1              1              1              1              1 
    ##      100402635      100403035      100403635      100406135      100406235 
    ##              1              1              1              1              1 
    ##      100406335      100408335      100408735      100409435      100410135 
    ##              1              1              1              1              1 
    ##      100414335      100414835      100417735      100417935      100419335 
    ##              1              1              1              1              1 
    ##      100422835      100425235      100426535      100427735      100430335 
    ##              1              1              1              1              1 
    ##      100431135      100432235      100432735      100432935      100434735 
    ##              1              1              1              1              1 
    ##      100439235      100439635      100440135      100443735      100443935 
    ##              1              1              1              1              1 
    ##      100444935      100446035      100446135      100446235      100448135 
    ##              1              1              1              1              1 
    ##      100448735      100449335      100450035      100450935      100451035 
    ##              1              1              1              1              1 
    ##      100451235      100451335      100452035      100453335      100453535 
    ##              1              1              1              1              1 
    ##      100453635      100454035      100456535      100460035      100460135 
    ##              1              1              1              1              1 
    ##      100461535      100463335      100463535      100464235      100465335 
    ##              1              1              1              1              1 
    ##      100466235      100467435      100467635      100467835      100468035 
    ##              1              1              1              1              1 
    ##      100468335      100468835      100469635      100470835      100470935 
    ##              1              1              1              1              1 
    ##      100471135      100473135      100473335      100473635      100475435 
    ##              1              1              1              1              1 
    ##      100475535      100476535      100476835    100477535-2    100477835-2 
    ##              1              1              1              1              1 
    ##      100479535      100479635      100480135      100480535    100480735-2 
    ##              1              1              1              1              1 
    ##      100481535      100481935      100483435      100483835      100485035 
    ##              1              1              1              1              1 
    ##      100485635      100485835      100486235    100489435-2    100490435-2 
    ##              1              1              1              1              1 
    ##      100490535      100491335      100492135      100494035      100494435 
    ##              1              1              1              1              1 
    ##      100494735      100494835      100496135      100496735      100496935 
    ##              1              1              1              1              1 
    ##      100500835      100501135      100501835      100502535      100503635 
    ##              1              1              1              1              1 
    ##      100505735      100506235      100507835      100509835      100511135 
    ##              1              1              1              1              1 
    ##      100512035    100513035-2      100513335      100515135      100517335 
    ##              1              1              1              1              1 
    ##      100518835      100519235      100524435      100525835      100527635 
    ##              1              1              1              1              1 
    ##      100528435      100528735      100531335      100561135      100563435 
    ##              1              1              1              1              1 
    ##      100563735      100564135      100565235      100565435      100566335 
    ##              1              1              1              1              1 
    ##      100566435      100566535      100566835      100567435      100567835 
    ##              1              1              1              1              1 
    ##      100569235      100570035    100570635-2      100572435      100572835 
    ##              1              1              1              1              1 
    ##      100573735      100575835      100576735      100577235      100577335 
    ##              1              1              1              1              1 
    ##      100578435      100579635      100579935      100580135      100580735 
    ##              1              1              1              1              1 
    ##      100580835      100581035      100582235      100582335      100582435 
    ##              1              1              1              1              1 
    ##      100583635      100583735      100584035      100584435      100584735 
    ##              1              1              1              1              1 
    ##      100585635      100586735      100587635      100587835    100588635-2 
    ##              1              1              1              1              1 
    ##      100588735      100589235      100589435      100589535      100589735 
    ##              1              1              1              1              1 
    ##      100591235      100591335      100592235      100593235      100594035 
    ##              1              1              1              1              1 
    ##      100594135      100595535      100595835      100595935      100597835 
    ##              1              1              1              1              1 
    ##      100597935      100598035      100598935      100599335      100599435 
    ##              1              1              1              1              1 
    ##      100599535      100599835      100601135      100602835      100603435 
    ##              1              1              1              1              1 
    ##      100603635      100604235      100604635      100604935      100605535 
    ##              1              1              1              1              1 
    ##      100607035      100607935      100608135      100608235      100609535 
    ##              1              1              1              1              1 
    ##      100610835      100631635      100636135      100637835      100638035 
    ##              1              1              1              1              1 
    ##      100638135      100638735      100638835      100638935      100639535 
    ##              1              1              1              1              1 
    ##      100639735      100640035      100641235      100641435      100680835 
    ##              1              1              1              1              1 
    ##      100681035      100681435      100681935      100682435      100683735 
    ##              1              1              1              1              1 
    ##      100685035      100685135      100685835      100687935      100691935 
    ##              1              1              1              1              1 
    ##    100692735-2   100693835-2W      100695335      100695735      100695935 
    ##              1              1              1              1              1 
    ##      100696535      100696835      100697835      100700535      100702835 
    ##              1              1              1              1              1 
    ##      100703535      100704235    100705635-2      100705935      100706235 
    ##              1              1              1              1              1 
    ##    100706535-2      100707535      100707735      100708935    100711435-2 
    ##              1              1              1              1              1 
    ##      100712235      100713435      100713735      100713835    100715135-2 
    ##              1              1              1              1              1 
    ##      100715935      100718935      100723535      100729435      100730535 
    ##              1              1              1              1              1 
    ##      100730935      100731935      100732135      100732335      100733235 
    ##              1              1              1              1              1 
    ##      100733935      100734235      100735035      100735435      100735535 
    ##              1              1              1              1              1 
    ##      100736535      100740635      100742635      100743335      100743635 
    ##              1              1              1              1              1 
    ##      100743835      100744235      100744735      100745035      100746535 
    ##              1              1              1              1              1 
    ##      100747035      100749535      100749835      100752035      100752535 
    ##              1              1              1              1              1 
    ##      100752735      100753935      100755435      100755635      100756235 
    ##              1              1              1              1              1 
    ##      100756735      100757835      100759735      100759835      100763635 
    ##              1              1              1              1              1 
    ##      100765935      100766835      100768635      100769235      100769535 
    ##              1              1              1              1              1 
    ##      100769735      100771135      100771835      100774835      100774935 
    ##              1              1              1              1              1 
    ##      100775035      100776335      100776435   100778235-2W      100778735 
    ##              1              1              1              1              1 
    ##      100787735      100790335      100790935      100791135      100792035 
    ##              1              1              1              1              1 
    ##      100792935      100793235      100793335   100793835-2W      100796535 
    ##              1              1              1              1              1 
    ##      100797735      100797835      100798935      100800735      100800835 
    ##              1              1              1              1              1 
    ##      100801335      100802435      100804035      100805135      100807235 
    ##              1              1              1              1              1 
    ##      100807935      100809735      100811035      100811935      100814335 
    ##              1              1              1              1              1 
    ##      100815535      100816935      100818835      100819135      100825635 
    ##              1              1              1              1              1 
    ##      100829735      100830535      100831735      100834435    100835235-2 
    ##              1              1              1              1              1 
    ##      100835935      100837235      100837735      100837835      100838435 
    ##              1              1              1              1              1 
    ##      100840435      100841835      100841935      100842035      100842335 
    ##              1              1              1              1              1 
    ##      100842835      100843635      100844335      100844435      100847435 
    ##              1              1              1              1              1 
    ##      100850935      100852035      100852535      100852835      100853235 
    ##              1              1              1              1              1 
    ##      100853335      100853735      100854535      100855035      100859135 
    ##              1              1              1              1              1 
    ##      100860735      100868235      100872435      100873635      100873735 
    ##              1              1              1              1              1 
    ##      100875835      100876135      100876535      100879535      100881235 
    ##              1              1              1              1              1 
    ##      100883535      100883835      100885335    100886035-2      100886435 
    ##              1              1              1              1              1 
    ##      100887035      100887135      100887835      100887935      100888735 
    ##              1              1              1              1              1 
    ##      100889235      100889435    100891235-2    100893735-2      100893935 
    ##              1              1              1              1              1 
    ##      100945035      100946235      100948735      100960035      100960235 
    ##              1              1              1              1              1 
    ##      100960735      100961035      100961235      100963435      100964035 
    ##              1              1              1              1              1 
    ##      100964235      100964735      100967735      100969835      100973335 
    ##              1              1              1              1              1 
    ##    100974035-2      100975235      100976935      100986635      100986735 
    ##              1              1              1              1              1 
    ##      100988435      100989235      100990935      100991335      100992835 
    ##              1              1              1              1              1 
    ##      100994135      100995935      100996135      100996935      100997135 
    ##              1              1              1              1              1 
    ##      100997235      100998135      100998535      100998935      100999035 
    ##              1              1              1              1              1 
    ##      101003535      101005435      101007135      101008735      101011735 
    ##              1              1              1              1              1 
    ##    101012135-2      101012535      101014135      101015935      101017535 
    ##              1              1              1              1              1 
    ##      101019535      101020535      101021735      101022135      101023835 
    ##              1              1              1              1              1 
    ##      101023935      101024335      101024835      101025735      101027735 
    ##              1              1              1              1              1 
    ##      101030335      101032235      101033035      101033635      101035935 
    ##              1              1              1              1              1 
    ##      101036035      101036335      101037435      101038735    101041235-2 
    ##              1              1              1              1              1 
    ##      101043135      101043635      101096535      101096735      101097235 
    ##              1              1              1              1              1 
    ##      101097535    101099635-2      101099835    101100135-2      101100235 
    ##              1              1              1              1              1 
    ##      101103735      101104335      101105235      101105735      101105835 
    ##              1              1              1              1              1 
    ##      101108635      101108735      101109635      101109935      101110035 
    ##              1              1              1              1              1 
    ##      101110135      101111035      101111135      101112435      101112935 
    ##              1              1              1              1              1 
    ##      101113335      101114335      101114735      101119135      101119235 
    ##              1              1              1              1              1 
    ##      101120235    101120335-2      101120935      101121335      101121435 
    ##              1              1              1              1              1 
    ##      101122935      101123535      101159335      101159835      101160835 
    ##              1              1              1              1              1 
    ##      101161435      101164235      101165335      101166435      101170635 
    ##              1              1              1              1              1 
    ##      101170835      101173135    101233635-2      101234335      101234735 
    ##              1              1              1              1              1 
    ##      101237135      101237535      101237935      101240935      101241835 
    ##              1              1              1              1              1 
    ##      101241935      101242635      101243035      101247035      101248435 
    ##              1              1              1              1              1 
    ##      101249635      101256335      101260935      101262335      101262435 
    ##              1              1              1              1              1 
    ##      101264135      101265035      101265235      101266635      101268435 
    ##              1              1              1              1              1 
    ##      101268535      101274335      101275435      101278235      101280735 
    ##              1              1              1              1              1 
    ##      101281235      101282235      101282435      101285035      101286435 
    ##              1              1              1              1              1 
    ##      101288335      101290735      101291935    101294535-2      101295335 
    ##              1              1              1              1              1 
    ##    101295735-2    101296535-2    101296935-2      101302735      101303635 
    ##              1              1              1              1              1 
    ##      101304735      101321135      101324935      101334035      101336635 
    ##              1              1              1              1              1 
    ##      101341635      101342635      101344035      101349335      101350135 
    ##              1              1              1              1              1 
    ##      101350435      101352235      101352835      101354835      101355135 
    ##              1              1              1              1              1 
    ##      101357035      101357635      101359835      101361835      101369035 
    ##              1              1              1              1              1 
    ##      101370535      101371735      101372235      101379435      101379535 
    ##              1              1              1              1              1 
    ##    101380335-2      101381135      101381335    101381835-2      101382735 
    ##              1              1              1              1              1 
    ##      101383035    101383435-2      101387135      101391135      101399135 
    ##              1              1              1              1              1 
    ##      101399335    101399835-2      101401035      101401635      101401735 
    ##              1              1              1              1              1 
    ##      101401935    101402135-2      101402335      101493835      101494035 
    ##              1              1              1              1              1 
    ##      101494135      101497235      101504635      101504735      101505335 
    ##              1              1              1              1              1 
    ##      101505535      101505835      101506335      101506535      101506635 
    ##              1              1              1              1              1 
    ##      101506935      101507435      101508635      101509735      101510735 
    ##              1              1              1              1              1 
    ##      101511135      101511935      101512435      101513135      101514335 
    ##              1              1              1              1              1 
    ##      101514435      101515335      101515535      101516435      101517035 
    ##              1              1              1              1              1 
    ##      101517735      101517835      101518435      101520235      101520935 
    ##              1              1              1              1              1 
    ##      101521135      101522435      101524335      101524535      101525235 
    ##              1              1              1              1              1 
    ##      101527635      101529335      101530735      101531035      101531635 
    ##              1              1              1              1              1 
    ##      101532835      101533335      101535435      101535835      101536335 
    ##              1              1              1              1              1 
    ##      101536835      101537335      101538235      101538435      101546235 
    ##              1              1              1              1              1 
    ##      101547535      101547935      101549035      101551135      101555735 
    ##              1              1              1              1              1 
    ##      101558435      101558635      101559335      101560035      101560135 
    ##              1              1              1              1              1 
    ##    101561635-2      101564635      101564935      101567135      101567335 
    ##              1              1              1              1              1 
    ##      101570335      101570535      101573235      101573735      101574835 
    ##              1              1              1              1              1 
    ##      101576435      101576535      101576835      101579335      101581135 
    ##              1              1              1              1              1 
    ##      101584435      101588135      101591535      101593235      101593635 
    ##              1              1              1              1              1 
    ##      101594035      101634535      101635835      101637135      101637635 
    ##              1              1              1              1              1 
    ##      101638735      101639235      101640335      101642135      101643035 
    ##              1              1              1              1              1 
    ##      101643835      101644535      101645235      101648035      101669135 
    ##              1              1              1              1              1 
    ##      101670335      101670535      101672435      101673735      101673835 
    ##              1              1              1              1              1 
    ##      101691435      101691535      101691935      101692735      101699535 
    ##              1              1              1              1              1 
    ##      101700635      101700835      101701035      101701735    101702135-2 
    ##              1              1              1              1              1 
    ##      101710135      101710935      101711635      101712135      101712335 
    ##              1              1              1              1              1 
    ##      101712935      101713535      101714735      101716135      101779732 
    ##              1              1              1              1              1 
    ##      101779932      101780132      101780832      101781632      101781832 
    ##              1              1              1              1              1 
    ##      101781932      101782232      101783232      101784132      101784432 
    ##              1              1              1              1              1 
    ##      101859435      101859635      101861635      101862335      101862535 
    ##              1              1              1              1              1 
    ##      101865235      101865535      101866735      101867135      101867335 
    ##              1              1              1              1              1 
    ##      101867935      101868035    101875235-2      101877435    101877935-2 
    ##              1              1              1              1              1 
    ##      101880035      101881135      101882235      101882435      101882635 
    ##              1              1              1              1              1 
    ##      101883835    101884435-2      101884835      101887635      101887735 
    ##              1              1              1              1              1 
    ##      101889035      101889335      101890235      101891735      101894235 
    ##              1              1              1              1              1 
    ##      101894635      101895135      101897335      101898335      101899035 
    ##              1              1              1              1              1 
    ##      101899235      101901035      101901335      101901735      101904035 
    ##              1              1              1              1              1 
    ##      101904935      101905035      101906335      101906435      101907035 
    ##              1              1              1              1              1 
    ##      101909635      101910435      101912335      101912635      101912935 
    ##              1              1              1              1              1 
    ##      101913135      101913735      101914535      101914835      101915235 
    ##              1              1              1              1              1 
    ##      101916435      101917635      101918535      101919135      101920135 
    ##              1              1              1              1              1 
    ##      101920335    101921235-2      101921535    101921835-2   101922035-2W 
    ##              1              1              1              1              1 
    ##    101922135-2      101922935    101923635-2    101923935-2      101924335 
    ##              1              1              1              1              1 
    ##      101926135      101928235      101929835      101932035      101932635 
    ##              1              1              1              1              1 
    ##      101932835      101933035      101933635      101933835      101934035 
    ##              1              1              1              1              1 
    ##      101935535      101938235    101939335-2      101941035      101943235 
    ##              1              1              1              1              1 
    ##      101944035      101944435      101946435    101947335-2      101947835 
    ##              1              1              1              1              1 
    ##      101948335      101948835      101949035      101952635    101953135-2 
    ##              1              1              1              1              1 
    ##      101956235      101960235      101962635      101963435      101970935 
    ##              1              1              1              1              1 
    ##      101972635      101972835      101973435      101974135      101974335 
    ##              1              1              1              1              1 
    ##      101974535      101977035      101977935      101980435      101981135 
    ##              1              1              1              1              1 
    ##      101981335      101981835      101982235      101983635      101984635 
    ##              1              1              1              1              1 
    ##      101984835      101985635      101985935      101986035      101986235 
    ##              1              1              1              1              1 
    ##      101987435      101987635      101988535      101989335      101991735 
    ##              1              1              1              1              1 
    ##      101993035      101994635      101995035      101995835      101997335 
    ##              1              1              1              1              1 
    ##      101997435      101997735      101998235      101999835      102000935 
    ##              1              1              1              1              1 
    ##      102001935      102003035      102003635      102003835      102003935 
    ##              1              1              1              1              1 
    ##    102004335-2      102005035      102006335      102006435      102008035 
    ##              1              1              1              1              1 
    ##      102010035      102010235    102010335-2      102013935      102014135 
    ##              1              1              1              1              1 
    ##      102014935      102015135      102015335      102015435      102016735 
    ##              1              1              1              1              1 
    ##      102017135      102017735      102018035      102019835      102020235 
    ##              1              1              1              1              1 
    ##      102020335      102021435      102023035      102025635      102026035 
    ##              1              1              1              1              1 
    ##      102026235      102026535      102027535      102029335      102029635 
    ##              1              1              1              1              1 
    ##      102030035      102030335      102031335      102031735      102032335 
    ##              1              1              1              1              1 
    ##      102032935      102033035      102033235      102033435      102036335 
    ##              1              1              1              1              1 
    ##      102037035      102037535      102038735      102039935      102040535 
    ##              1              1              1              1              1 
    ##      102040635      102040735      102042935      102046735      102049535 
    ##              1              1              1              1              1 
    ##      102050035      102050735      102050835      102057035      102058235 
    ##              1              1              1              1              1 
    ##      102063535      102064035      102064735      102066435      102066935 
    ##              1              1              1              1              1 
    ##      102068435      102068535      102068935      102069735      102069835 
    ##              1              1              1              1              1 
    ##      102070635      102076335    102076635-2    102077035-2    102077235-2 
    ##              1              1              1              1              1 
    ##      102078235      102078635      102081235      102082235      102084235 
    ##              1              1              1              1              1 
    ##      102084735      102087935      102088235      102088435      102088735 
    ##              1              1              1              1              1 
    ##      102088835      102089235      102089935      102090035      102090635 
    ##              1              1              1              1              1 
    ##      102090935      102091235    102091335-2      102091635      102096335 
    ##              1              1              1              1              1 
    ##      102097035      102097335      102099635      102099735      102100135 
    ##              1              1              1              1              1 
    ##      102100535      102101235      102103235      102103835      102103935 
    ##              1              1              1              1              1 
    ##      102104335      102104835      102106335      102109435      102111135 
    ##              1              1              1              1              1 
    ##    102112535-2      102112835      102113235      102115335      102116435 
    ##              1              1              1              1              1 
    ##      102117235      102118735      102119135      102119335      102123935 
    ##              1              1              1              1              1 
    ##      102125035      102125135    102125535-2      102125635      102125935 
    ##              1              1              1              1              1 
    ##      102126035      102127535      102129035      102129335      102129735 
    ##              1              1              1              1              1 
    ##      102131135      102132335      102133335      102133535      102136635 
    ##              1              1              1              1              1 
    ##      102137335      102139035      102139935      102140135      102140635 
    ##              1              1              1              1              1 
    ##      102141635      102142135      102145335      102145935      102147835 
    ##              1              1              1              1              1 
    ##      102148035      102148235      102148835      102149035      102149435 
    ##              1              1              1              1              1 
    ##      102150335    102151035-2      102151335      102152135      102153135 
    ##              1              1              1              1              1 
    ##      102154035    102155135-2    102155335-2      102156735    102157535-2 
    ##              1              1              1              1              1 
    ##      102157635      102158335    102158835-2      102159635      102160135 
    ##              1              1              1              1              1 
    ##      102160435      102160635      102161635      102161735      102163435 
    ##              1              1              1              1              1 
    ##      102166535      102166835      102170035      102171735      102172335 
    ##              1              1              1              1              1 
    ##      102174135      102174735      102175735      102177635      102178135 
    ##              1              1              1              1              1 
    ##      102179335      102181435      102182035      102183735      102186235 
    ##              1              1              1              1              1 
    ##      102187235      102187335      102189035      102190035      102192035 
    ##              1              1              1              1              1 
    ##    102193835-2    102196135-2      102198435      102198535      102198735 
    ##              1              1              1              1              1 
    ##      102199935      102200635      102202135      102206735      102208035 
    ##              1              1              1              1              1 
    ##      102208735      102209535      102210535      102213135      102216135 
    ##              1              1              1              1              1 
    ##      102217635      102218635      102222035      102223635      102224435 
    ##              1              1              1              1              1 
    ##      102224735      102225335      102226135      102226735      102226935 
    ##              1              1              1              1              1 
    ##      102227435      102230035      102230635      102233135      102234935 
    ##              1              1              1              1              1 
    ##      102236035      102236235      102237035      102237835      102238935 
    ##              1              1              1              1              1 
    ##      102243935      102244335      102244735      102246535      102246635 
    ##              1              1              1              1              1 
    ##      102247335      102248935      102249935      102250935      102252035 
    ##              1              1              1              1              1 
    ##      102253535      102253935      102258435    102258535-2      102258935 
    ##              1              1              1              1              1 
    ##    102263135-2      102264635      102265335      102267135      102267335 
    ##              1              1              1              1              1 
    ##      102268735      102269035      102269635      102270135      102270335 
    ##              1              1              1              1              1 
    ##      102271035      102271235      102271535      102271835      102272435 
    ##              1              1              1              1              1 
    ##      102272535      102272935      102273235      102273535      102275935 
    ##              1              1              1              1              1 
    ##      102276435      102277935      102278035      102279835      102281335 
    ##              1              1              1              1              1 
    ##      102282535      102283735      102284135      102285435      102286335 
    ##              1              1              1              1              1 
    ##      102286735      102287035      102287135      102288235      102288435 
    ##              1              1              1              1              1 
    ##      102288735      102291735      102292535      102294235      102295835 
    ##              1              1              1              1              1 
    ##      102296435      102299135      102299935      102300235      102300335 
    ##              1              1              1              1              1 
    ##      102300535      102302035      102302335      102302935      102303835 
    ##              1              1              1              1              1 
    ##      102303935      102304535      102307635      102311535      102315135 
    ##              1              1              1              1              1 
    ##      102319135      102319335      102321035      102321835      102322135 
    ##              1              1              1              1              1 
    ##      102324535      102326335      102327535      102329135      102329935 
    ##              1              1              1              1              1 
    ##      102330635      102331635      102334335      102334735      102337635 
    ##              1              1              1              1              1 
    ##      102339235      102339435      102341335      102341635      102343635 
    ##              1              1              1              1              1 
    ##      102346535      102348535      102349535      102349635      102350835 
    ##              1              1              1              1              1 
    ##      102350935      102352035      102352935      102353035      102353335 
    ##              1              1              1              1              1 
    ##      102353435      102353535      102353635      102355535      102357735 
    ##              1              1              1              1              1 
    ##      102358335      102361035      102361235      102362435      102363535 
    ##              1              1              1              1              1 
    ##      102364635      102365035      102366635      102366735      102368635 
    ##              1              1              1              1              1 
    ##      102369835      102372335      102373235      102374435      102375335 
    ##              1              1              1              1              1 
    ##      102376035      102377635      102380135    102380935-2      102381135 
    ##              1              1              1              1              1 
    ##      102382835      102384135      102384235      102384535      102386135 
    ##              1              1              1              1              1 
    ##      102386935      102389035      102389135    102389535-2      102394535 
    ##              1              1              1              1              1 
    ##      102396835      102400035    102400835-2    102401035-2      102401335 
    ##              1              1              1              1              1 
    ##    102401635-2   102401735-2W      102403335      102403535      102404035 
    ##              1              1              1              1              1 
    ##      102404135      102404235      102405635      102412135      102412235 
    ##              1              1              1              1              1 
    ##    102413135-2      102489235      102490535      102493235      102493335 
    ##              1              1              1              1              1 
    ##      102495535      102496535      102496735      102497135      102497235 
    ##              1              1              1              1              1 
    ##      102497835      102498235      102499235      102500935      102502135 
    ##              1              1              1              1              1 
    ##      102502435      102503335      102509935      102510035      102510135 
    ##              1              1              1              1              1 
    ##      102510535      102512335      102512735      102513435      102514935 
    ##              1              1              1              1              1 
    ##      102517335      102518235      102524135      102525335      102525835 
    ##              1              1              1              1              1 
    ##      102525935      102527935      102529935      102535735      102536035 
    ##              1              1              1              1              1 
    ##      102545135      102546335      102546935      102549135      102550435 
    ##              1              1              1              1              1 
    ##      102550735      102551035      102551135      102551935      102552735 
    ##              1              1              1              1              1 
    ##   102559435-2W    102559835-2      102560135      102560335      102560435 
    ##              1              1              1              1              1 
    ##      102560935    102562935-2      102564835      102565135      102565235 
    ##              1              1              1              1              1 
    ##      102565435      102566235      102566735      102566835      102566935 
    ##              1              1              1              1              1 
    ##      102567135    102568135-2      102569335      102569635      102571135 
    ##              1              1              1              1              1 
    ##      102571735   102579835-2W      102580135      102581435      102583035 
    ##              1              1              1              1              1 
    ##      102584735      102585335      102588335      102588735      102589735 
    ##              1              1              1              1              1 
    ##      102589835      102590135      102592035      102592435      102592635 
    ##              1              1              1              1              1 
    ##      102604635      102607535      102607635      102608935    102615035-2 
    ##              1              1              1              1              1 
    ##      102615535      102615635    102615935-2      102618235      102629335 
    ##              1              1              1              1              1 
    ##      102630135      102630635      102631535      102640435      102649735 
    ##              1              1              1              1              1 
    ##      102650235      102650735      102651335      102659535      102659735 
    ##              1              1              1              1              1 
    ##      102660335      102662435      102662635 SAMEA115126061 SAMEA115126107 
    ##              1              1              1              1              1 
    ## SAMEA115126171 SAMEA115126306 SAMEA115126326 SAMEA115126491 SAMEA115126500 
    ##              1              1              1              1              1 
    ## SAMEA115126520 SAMEA115126613 SAMEA115126615 SAMEA115126635 SAMEA115126710 
    ##              1              1              1              1              1 
    ## SAMEA115126717 SAMEA115126775 SAMEA115126863 SAMEA115126870 SAMEA115126888 
    ##              1              1              1              1              1 
    ## SAMEA115127046 SAMEA115127059 SAMEA115127165 SAMEA115127331 SAMEA115127484 
    ##              1              1              1              1              1 
    ## SAMEA115127502 SAMEA115127533 SAMEA115127719 
    ##              1              1              1

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
    ##   NO 
    ## 1698

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
    ##      100000135      100000735      100002335      100002535      100002935 
    ##              1              1              1              1              1 
    ##      100003235      100004035      100004535      100005135      100005235 
    ##              1              1              1              1              1 
    ##      100005435      100005735      100006235      100006935      100007135 
    ##              1              1              1              1              1 
    ##      100007635      100008335      100009635      100012935      100013735 
    ##              1              1              1              1              1 
    ##      100014135      100014235      100014435      100014635      100014935 
    ##              1              1              1              1              1 
    ##      100015235      100015435      100015535      100015835      100016235 
    ##              1              1              1              1              1 
    ##      100016435      100016535      100017835      100017935      100018635 
    ##              1              1              1              1              1 
    ##      100019535    100020535-2      100020735      100020935      100022235 
    ##              1              1              1              1              1 
    ##      100022335      100023235      100023535      100023935    100025335-2 
    ##              1              1              1              1              1 
    ##      100025435    100025835-2      100027035      100027235      100027335 
    ##              1              1              1              1              1 
    ##    100027535-2      100027635    100028435-2      100029235      100029435 
    ##              1              1              1              1              1 
    ##      100029935      100030635      100032635      100032835      100032935 
    ##              1              1              1              1              1 
    ##      100034635      100035535      100038235      100039235      100045735 
    ##              1              1              1              1              1 
    ##      100046235      100046335      100046635      100047235      100047835 
    ##              1              1              1              1              1 
    ##      100048035      100048335      100048635      100049835      100050435 
    ##              1              1              1              1              1 
    ##      100050535      100050735      100051235      100051835      100053535 
    ##              1              1              1              1              1 
    ##      100053635      100054435      100055235      100056135      100056435 
    ##              1              1              1              1              1 
    ##      100057235      100057535      100058435      100058535      100058935 
    ##              1              1              1              1              1 
    ##      100059035      100059635      100060135      100060235      100062535 
    ##              1              1              1              1              1 
    ##      100063035    100063135-2      100063935      100064135      100064635 
    ##              1              1              1              1              1 
    ##      100065835      100066235      100067235      100068635      100070435 
    ##              1              1              1              1              1 
    ##      100072435      100072735    100073835-2      100073935      100074535 
    ##              1              1              1              1              1 
    ##      100074735      100074935      100076135      100076935      100078035 
    ##              1              1              1              1              1 
    ##      100078835      100079135      100079835      100080235      100080835 
    ##              1              1              1              1              1 
    ##      100082035      100084635      100085135      100085435      100089135 
    ##              1              1              1              1              1 
    ##      100089435      100090035      100090235      100092235      100092535 
    ##              1              1              1              1              1 
    ##      100095035      100097335      100098835      100099635      100102435 
    ##              1              1              1              1              1 
    ##      100102535      100106735      100107235      100107435      100107535 
    ##              1              1              1              1              1 
    ##      100108435      100113135      100113435      100113535      100113735 
    ##              1              1              1              1              1 
    ##      100115135      100117135      100117235      100118535      100119335 
    ##              1              1              1              1              1 
    ##      100120135      100120435      100121535      100122735      100123435 
    ##              1              1              1              1              1 
    ##      100123535      100124135      100126435      100126635      100127135 
    ##              1              1              1              1              1 
    ##      100127235      100128835      100136635      100136935      100137735 
    ##              1              1              1              1              1 
    ##      100138335      100139435      100139735      100140935      100141235 
    ##              1              1              1              1              1 
    ##      100141435      100141935      100142635      100144435      100146735 
    ##              1              1              1              1              1 
    ##      100146935      100147035      100147535      100147835      100151935 
    ##              1              1              1              1              1 
    ##      100153535      100153935      100154835      100156735      100157235 
    ##              1              1              1              1              1 
    ##      100160535    100160935-2      100161635      100161935      100162035 
    ##              1              1              1              1              1 
    ##      100163935      100164035      100164235      100164435      100164635 
    ##              1              1              1              1              1 
    ##      100165535      100165935      100166135    100166235-2      100166335 
    ##              1              1              1              1              1 
    ##      100166835      100167935      100168235      100168835      100169835 
    ##              1              1              1              1              1 
    ##      100170135      100170235      100171335    100178935-2      100182835 
    ##              1              1              1              1              1 
    ##      100182935    100186835-2      100192635      100194535      100195835 
    ##              1              1              1              1              1 
    ##      100196135      100200835      100201235      100202135      100203535 
    ##              1              1              1              1              1 
    ##      100203835      100206535      100220135      100220235      100220835 
    ##              1              1              1              1              1 
    ##      100222835      100223535      100225635      100228235    100230735-2 
    ##              1              1              1              1              1 
    ##      100231935      100233435      100234435      100236935      100238235 
    ##              1              1              1              1              1 
    ##      100240135      100241235      100241735      100241935      100242035 
    ##              1              1              1              1              1 
    ##      100242635      100243935      100245135      100245635      100253435 
    ##              1              1              1              1              1 
    ##      100253835      100254335      100254835      100255235      100260835 
    ##              1              1              1              1              1 
    ##      100261735      100262335      100263235      100264035      100266335 
    ##              1              1              1              1              1 
    ##      100267935      100268035      100268235      100268735      100269835 
    ##              1              1              1              1              1 
    ##      100272935      100274135      100277935      100279535      100291735 
    ##              1              1              1              1              1 
    ##      100299535      100299835      100300135      100303935      100306235 
    ##              1              1              1              1              1 
    ##      100306635      100306735      100306835      100308035      100309235 
    ##              1              1              1              1              1 
    ##      100315735      100316935      100318135      100318435      100323935 
    ##              1              1              1              1              1 
    ##      100325435      100329635      100332235      100332735      100334635 
    ##              1              1              1              1              1 
    ##      100335635      100336235      100336435      100336935      100339235 
    ##              1              1              1              1              1 
    ##      100341835      100342535      100342635      100343135      100343835 
    ##              1              1              1              1              1 
    ##      100344835      100347935      100348235      100350635      100352535 
    ##              1              1              1              1              1 
    ##      100352935      100354235      100356135      100358035      100360335 
    ##              1              1              1              1              1 
    ##      100364935      100366535      100367035      100367135      100367935 
    ##              1              1              1              1              1 
    ##      100371235      100371535      100378535      100379335      100379635 
    ##              1              1              1              1              1 
    ##      100380935      100381735      100382135      100383335      100383835 
    ##              1              1              1              1              1 
    ##      100383935      100385235      100386035      100386135      100387935 
    ##              1              1              1              1              1 
    ##      100389235      100389335      100391035      100391335      100391435 
    ##              1              1              1              1              1 
    ##      100391735      100392335      100393435      100394335      100395535 
    ##              1              1              1              1              1 
    ##      100396135      100397735      100398735      100399135      100399835 
    ##              1              1              1              1              1 
    ##      100402635      100403035      100403635      100406135      100406235 
    ##              1              1              1              1              1 
    ##      100406335      100408335      100408735      100409435      100410135 
    ##              1              1              1              1              1 
    ##      100414335      100414835      100417735      100417935      100419335 
    ##              1              1              1              1              1 
    ##      100422835      100425235      100426535      100427735      100430335 
    ##              1              1              1              1              1 
    ##      100431135      100432235      100432735      100432935      100434735 
    ##              1              1              1              1              1 
    ##      100439235      100439635      100440135      100443735      100443935 
    ##              1              1              1              1              1 
    ##      100444935      100446035      100446135      100446235      100448135 
    ##              1              1              1              1              1 
    ##      100448735      100449335      100450035      100450935      100451035 
    ##              1              1              1              1              1 
    ##      100451235      100451335      100452035      100453335      100453535 
    ##              1              1              1              1              1 
    ##      100453635      100454035      100456535      100460035      100460135 
    ##              1              1              1              1              1 
    ##      100461535      100463335      100463535      100464235      100465335 
    ##              1              1              1              1              1 
    ##      100466235      100467435      100467635      100467835      100468035 
    ##              1              1              1              1              1 
    ##      100468335      100468835      100469635      100470835      100470935 
    ##              1              1              1              1              1 
    ##      100471135      100473135      100473335      100473635      100475435 
    ##              1              1              1              1              1 
    ##      100475535      100476535      100476835    100477535-2    100477835-2 
    ##              1              1              1              1              1 
    ##      100479535      100479635      100480135      100480535    100480735-2 
    ##              1              1              1              1              1 
    ##      100481535      100481935      100483435      100483835      100485035 
    ##              1              1              1              1              1 
    ##      100485635      100485835      100486235    100489435-2    100490435-2 
    ##              1              1              1              1              1 
    ##      100490535      100491335      100492135      100494035      100494435 
    ##              1              1              1              1              1 
    ##      100494735      100494835      100496135      100496735      100496935 
    ##              1              1              1              1              1 
    ##      100500835      100501135      100501835      100502535      100503635 
    ##              1              1              1              1              1 
    ##      100505735      100506235      100507835      100509835      100511135 
    ##              1              1              1              1              1 
    ##      100512035    100513035-2      100513335      100515135      100517335 
    ##              1              1              1              1              1 
    ##      100518835      100519235      100524435      100525835      100527635 
    ##              1              1              1              1              1 
    ##      100528435      100528735      100531335      100561135      100563435 
    ##              1              1              1              1              1 
    ##      100563735      100564135      100565235      100565435      100566335 
    ##              1              1              1              1              1 
    ##      100566435      100566535      100566835      100567435      100567835 
    ##              1              1              1              1              1 
    ##      100569235      100570035    100570635-2      100572435      100572835 
    ##              1              1              1              1              1 
    ##      100573735      100575835      100576735      100577235      100577335 
    ##              1              1              1              1              1 
    ##      100578435      100579635      100579935      100580135      100580735 
    ##              1              1              1              1              1 
    ##      100580835      100581035      100582235      100582335      100582435 
    ##              1              1              1              1              1 
    ##      100583635      100583735      100584035      100584435      100584735 
    ##              1              1              1              1              1 
    ##      100585635      100586735      100587635      100587835    100588635-2 
    ##              1              1              1              1              1 
    ##      100588735      100589235      100589435      100589535      100589735 
    ##              1              1              1              1              1 
    ##      100591235      100591335      100592235      100593235      100594035 
    ##              1              1              1              1              1 
    ##      100594135      100595535      100595835      100595935      100597835 
    ##              1              1              1              1              1 
    ##      100597935      100598035      100598935      100599335      100599435 
    ##              1              1              1              1              1 
    ##      100599535      100599835      100601135      100602835      100603435 
    ##              1              1              1              1              1 
    ##      100603635      100604235      100604635      100604935      100605535 
    ##              1              1              1              1              1 
    ##      100607035      100607935      100608135      100608235      100609535 
    ##              1              1              1              1              1 
    ##      100610835      100631635      100636135      100637835      100638035 
    ##              1              1              1              1              1 
    ##      100638135      100638735      100638835      100638935      100639535 
    ##              1              1              1              1              1 
    ##      100639735      100640035      100641235      100641435      100680835 
    ##              1              1              1              1              1 
    ##      100681035      100681435      100681935      100682435      100683735 
    ##              1              1              1              1              1 
    ##      100685035      100685135      100685835      100687935      100691935 
    ##              1              1              1              1              1 
    ##    100692735-2   100693835-2W      100695335      100695735      100695935 
    ##              1              1              1              1              1 
    ##      100696535      100696835      100697835      100700535      100702835 
    ##              1              1              1              1              1 
    ##      100703535      100704235    100705635-2      100705935      100706235 
    ##              1              1              1              1              1 
    ##    100706535-2      100707535      100707735      100708935    100711435-2 
    ##              1              1              1              1              1 
    ##      100712235      100713435      100713735      100713835    100715135-2 
    ##              1              1              1              1              1 
    ##      100715935      100718935      100723535      100729435      100730535 
    ##              1              1              1              1              1 
    ##      100730935      100731935      100732135      100732335      100733235 
    ##              1              1              1              1              1 
    ##      100733935      100734235      100735035      100735435      100735535 
    ##              1              1              1              1              1 
    ##      100736535      100740635      100742635      100743335      100743635 
    ##              1              1              1              1              1 
    ##      100743835      100744235      100744735      100745035      100746535 
    ##              1              1              1              1              1 
    ##      100747035      100749535      100749835      100752035      100752535 
    ##              1              1              1              1              1 
    ##      100752735      100753935      100755435      100755635      100756235 
    ##              1              1              1              1              1 
    ##      100756735      100757835      100759735      100759835      100763635 
    ##              1              1              1              1              1 
    ##      100765935      100766835      100768635      100769235      100769535 
    ##              1              1              1              1              1 
    ##      100769735      100771135      100771835      100774835      100774935 
    ##              1              1              1              1              1 
    ##      100775035      100776335      100776435   100778235-2W      100778735 
    ##              1              1              1              1              1 
    ##      100787735      100790335      100790935      100791135      100792035 
    ##              1              1              1              1              1 
    ##      100792935      100793235      100793335   100793835-2W      100796535 
    ##              1              1              1              1              1 
    ##      100797735      100797835      100798935      100800735      100800835 
    ##              1              1              1              1              1 
    ##      100801335      100802435      100804035      100805135      100807235 
    ##              1              1              1              1              1 
    ##      100807935      100809735      100811035      100811935      100814335 
    ##              1              1              1              1              1 
    ##      100815535      100816935      100818835      100819135      100825635 
    ##              1              1              1              1              1 
    ##      100829735      100830535      100831735      100834435    100835235-2 
    ##              1              1              1              1              1 
    ##      100835935      100837235      100837735      100837835      100838435 
    ##              1              1              1              1              1 
    ##      100840435      100841835      100841935      100842035      100842335 
    ##              1              1              1              1              1 
    ##      100842835      100843635      100844335      100844435      100847435 
    ##              1              1              1              1              1 
    ##      100850935      100852035      100852535      100852835      100853235 
    ##              1              1              1              1              1 
    ##      100853335      100853735      100854535      100855035      100859135 
    ##              1              1              1              1              1 
    ##      100860735      100868235      100872435      100873635      100873735 
    ##              1              1              1              1              1 
    ##      100875835      100876135      100876535      100879535      100881235 
    ##              1              1              1              1              1 
    ##      100883535      100883835      100885335    100886035-2      100886435 
    ##              1              1              1              1              1 
    ##      100887035      100887135      100887835      100887935      100888735 
    ##              1              1              1              1              1 
    ##      100889235      100889435    100891235-2    100893735-2      100893935 
    ##              1              1              1              1              1 
    ##      100945035      100946235      100948735      100960035      100960235 
    ##              1              1              1              1              1 
    ##      100960735      100961035      100961235      100963435      100964035 
    ##              1              1              1              1              1 
    ##      100964235      100964735      100967735      100969835      100973335 
    ##              1              1              1              1              1 
    ##    100974035-2      100975235      100976935      100986635      100986735 
    ##              1              1              1              1              1 
    ##      100988435      100989235      100990935      100991335      100992835 
    ##              1              1              1              1              1 
    ##      100994135      100995935      100996135      100996935      100997135 
    ##              1              1              1              1              1 
    ##      100997235      100998135      100998535      100998935      100999035 
    ##              1              1              1              1              1 
    ##      101003535      101005435      101007135      101008735      101011735 
    ##              1              1              1              1              1 
    ##    101012135-2      101012535      101014135      101015935      101017535 
    ##              1              1              1              1              1 
    ##      101019535      101020535      101021735      101022135      101023835 
    ##              1              1              1              1              1 
    ##      101023935      101024335      101024835      101025735      101027735 
    ##              1              1              1              1              1 
    ##      101030335      101032235      101033035      101033635      101035935 
    ##              1              1              1              1              1 
    ##      101036035      101036335      101037435      101038735    101041235-2 
    ##              1              1              1              1              1 
    ##      101043135      101043635      101096535      101096735      101097235 
    ##              1              1              1              1              1 
    ##      101097535    101099635-2      101099835    101100135-2      101100235 
    ##              1              1              1              1              1 
    ##      101103735      101104335      101105235      101105735      101105835 
    ##              1              1              1              1              1 
    ##      101108635      101108735      101109635      101109935      101110035 
    ##              1              1              1              1              1 
    ##      101110135      101111035      101111135      101112435      101112935 
    ##              1              1              1              1              1 
    ##      101113335      101114335      101114735      101119135      101119235 
    ##              1              1              1              1              1 
    ##      101120235    101120335-2      101120935      101121335      101121435 
    ##              1              1              1              1              1 
    ##      101122935      101123535      101159335      101159835      101160835 
    ##              1              1              1              1              1 
    ##      101161435      101164235      101165335      101166435      101170635 
    ##              1              1              1              1              1 
    ##      101170835      101173135    101233635-2      101234335      101234735 
    ##              1              1              1              1              1 
    ##      101237135      101237535      101237935      101240935      101241835 
    ##              1              1              1              1              1 
    ##      101241935      101242635      101243035      101247035      101248435 
    ##              1              1              1              1              1 
    ##      101249635      101256335      101260935      101262335      101262435 
    ##              1              1              1              1              1 
    ##      101264135      101265035      101265235      101266635      101268435 
    ##              1              1              1              1              1 
    ##      101268535      101274335      101275435      101278235      101280735 
    ##              1              1              1              1              1 
    ##      101281235      101282235      101282435      101285035      101286435 
    ##              1              1              1              1              1 
    ##      101288335      101290735      101291935    101294535-2      101295335 
    ##              1              1              1              1              1 
    ##    101295735-2    101296535-2    101296935-2      101302735      101303635 
    ##              1              1              1              1              1 
    ##      101304735      101321135      101324935      101334035      101336635 
    ##              1              1              1              1              1 
    ##      101341635      101342635      101344035      101349335      101350135 
    ##              1              1              1              1              1 
    ##      101350435      101352235      101352835      101354835      101355135 
    ##              1              1              1              1              1 
    ##      101357035      101357635      101359835      101361835      101369035 
    ##              1              1              1              1              1 
    ##      101370535      101371735      101372235      101379435      101379535 
    ##              1              1              1              1              1 
    ##    101380335-2      101381135      101381335    101381835-2      101382735 
    ##              1              1              1              1              1 
    ##      101383035    101383435-2      101387135      101391135      101399135 
    ##              1              1              1              1              1 
    ##      101399335    101399835-2      101401035      101401635      101401735 
    ##              1              1              1              1              1 
    ##      101401935    101402135-2      101402335      101493835      101494035 
    ##              1              1              1              1              1 
    ##      101494135      101497235      101504635      101504735      101505335 
    ##              1              1              1              1              1 
    ##      101505535      101505835      101506335      101506535      101506635 
    ##              1              1              1              1              1 
    ##      101506935      101507435      101508635      101509735      101510735 
    ##              1              1              1              1              1 
    ##      101511135      101511935      101512435      101513135      101514335 
    ##              1              1              1              1              1 
    ##      101514435      101515335      101515535      101516435      101517035 
    ##              1              1              1              1              1 
    ##      101517735      101517835      101518435      101520235      101520935 
    ##              1              1              1              1              1 
    ##      101521135      101522435      101524335      101524535      101525235 
    ##              1              1              1              1              1 
    ##      101527635      101529335      101530735      101531035      101531635 
    ##              1              1              1              1              1 
    ##      101532835      101533335      101535435      101535835      101536335 
    ##              1              1              1              1              1 
    ##      101536835      101537335      101538235      101538435      101546235 
    ##              1              1              1              1              1 
    ##      101547535      101547935      101549035      101551135      101555735 
    ##              1              1              1              1              1 
    ##      101558435      101558635      101559335      101560035      101560135 
    ##              1              1              1              1              1 
    ##    101561635-2      101564635      101564935      101567135      101567335 
    ##              1              1              1              1              1 
    ##      101570335      101570535      101573235      101573735      101574835 
    ##              1              1              1              1              1 
    ##      101576435      101576535      101576835      101579335      101581135 
    ##              1              1              1              1              1 
    ##      101584435      101588135      101591535      101593235      101593635 
    ##              1              1              1              1              1 
    ##      101594035      101634535      101635835      101637135      101637635 
    ##              1              1              1              1              1 
    ##      101638735      101639235      101640335      101642135      101643035 
    ##              1              1              1              1              1 
    ##      101643835      101644535      101645235      101648035      101669135 
    ##              1              1              1              1              1 
    ##      101670335      101670535      101672435      101673735      101673835 
    ##              1              1              1              1              1 
    ##      101691435      101691535      101691935      101692735      101699535 
    ##              1              1              1              1              1 
    ##      101700635      101700835      101701035      101701735    101702135-2 
    ##              1              1              1              1              1 
    ##      101710135      101710935      101711635      101712135      101712335 
    ##              1              1              1              1              1 
    ##      101712935      101713535      101714735      101716135      101779732 
    ##              1              1              1              1              1 
    ##      101779932      101780132      101780832      101781632      101781832 
    ##              1              1              1              1              1 
    ##      101781932      101782232      101783232      101784132      101784432 
    ##              1              1              1              1              1 
    ##      101859435      101859635      101861635      101862335      101862535 
    ##              1              1              1              1              1 
    ##      101865235      101865535      101866735      101867135      101867335 
    ##              1              1              1              1              1 
    ##      101867935      101868035    101875235-2      101877435    101877935-2 
    ##              1              1              1              1              1 
    ##      101880035      101881135      101882235      101882435      101882635 
    ##              1              1              1              1              1 
    ##      101883835    101884435-2      101884835      101887635      101887735 
    ##              1              1              1              1              1 
    ##      101889035      101889335      101890235      101891735      101894235 
    ##              1              1              1              1              1 
    ##      101894635      101895135      101897335      101898335      101899035 
    ##              1              1              1              1              1 
    ##      101899235      101901035      101901335      101901735      101904035 
    ##              1              1              1              1              1 
    ##      101904935      101905035      101906335      101906435      101907035 
    ##              1              1              1              1              1 
    ##      101909635      101910435      101912335      101912635      101912935 
    ##              1              1              1              1              1 
    ##      101913135      101913735      101914535      101914835      101915235 
    ##              1              1              1              1              1 
    ##      101916435      101917635      101918535      101919135      101920135 
    ##              1              1              1              1              1 
    ##      101920335    101921235-2      101921535    101921835-2   101922035-2W 
    ##              1              1              1              1              1 
    ##    101922135-2      101922935    101923635-2    101923935-2      101924335 
    ##              1              1              1              1              1 
    ##      101926135      101928235      101929835      101932035      101932635 
    ##              1              1              1              1              1 
    ##      101932835      101933035      101933635      101933835      101934035 
    ##              1              1              1              1              1 
    ##      101935535      101938235    101939335-2      101941035      101943235 
    ##              1              1              1              1              1 
    ##      101944035      101944435      101946435    101947335-2      101947835 
    ##              1              1              1              1              1 
    ##      101948335      101948835      101949035      101952635    101953135-2 
    ##              1              1              1              1              1 
    ##      101956235      101960235      101962635      101963435      101970935 
    ##              1              1              1              1              1 
    ##      101972635      101972835      101973435      101974135      101974335 
    ##              1              1              1              1              1 
    ##      101974535      101977035      101977935      101980435      101981135 
    ##              1              1              1              1              1 
    ##      101981335      101981835      101982235      101983635      101984635 
    ##              1              1              1              1              1 
    ##      101984835      101985635      101985935      101986035      101986235 
    ##              1              1              1              1              1 
    ##      101987435      101987635      101988535      101989335      101991735 
    ##              1              1              1              1              1 
    ##      101993035      101994635      101995035      101995835      101997335 
    ##              1              1              1              1              1 
    ##      101997435      101997735      101998235      101999835      102000935 
    ##              1              1              1              1              1 
    ##      102001935      102003035      102003635      102003835      102003935 
    ##              1              1              1              1              1 
    ##    102004335-2      102005035      102006335      102006435      102008035 
    ##              1              1              1              1              1 
    ##      102010035      102010235    102010335-2      102013935      102014135 
    ##              1              1              1              1              1 
    ##      102014935      102015135      102015335      102015435      102016735 
    ##              1              1              1              1              1 
    ##      102017135      102017735      102018035      102019835      102020235 
    ##              1              1              1              1              1 
    ##      102020335      102021435      102023035      102025635      102026035 
    ##              1              1              1              1              1 
    ##      102026235      102026535      102027535      102029335      102029635 
    ##              1              1              1              1              1 
    ##      102030035      102030335      102031335      102031735      102032335 
    ##              1              1              1              1              1 
    ##      102032935      102033035      102033235      102033435      102036335 
    ##              1              1              1              1              1 
    ##      102037035      102037535      102038735      102039935      102040535 
    ##              1              1              1              1              1 
    ##      102040635      102040735      102042935      102046735      102049535 
    ##              1              1              1              1              1 
    ##      102050035      102050735      102050835      102057035      102058235 
    ##              1              1              1              1              1 
    ##      102063535      102064035      102064735      102066435      102066935 
    ##              1              1              1              1              1 
    ##      102068435      102068535      102068935      102069735      102069835 
    ##              1              1              1              1              1 
    ##      102070635      102076335    102076635-2    102077035-2    102077235-2 
    ##              1              1              1              1              1 
    ##      102078235      102078635      102081235      102082235      102084235 
    ##              1              1              1              1              1 
    ##      102084735      102087935      102088235      102088435      102088735 
    ##              1              1              1              1              1 
    ##      102088835      102089235      102089935      102090035      102090635 
    ##              1              1              1              1              1 
    ##      102090935      102091235    102091335-2      102091635      102096335 
    ##              1              1              1              1              1 
    ##      102097035      102097335      102099635      102099735      102100135 
    ##              1              1              1              1              1 
    ##      102100535      102101235      102103235      102103835      102103935 
    ##              1              1              1              1              1 
    ##      102104335      102104835      102106335      102109435      102111135 
    ##              1              1              1              1              1 
    ##    102112535-2      102112835      102113235      102115335      102116435 
    ##              1              1              1              1              1 
    ##      102117235      102118735      102119135      102119335      102123935 
    ##              1              1              1              1              1 
    ##      102125035      102125135    102125535-2      102125635      102125935 
    ##              1              1              1              1              1 
    ##      102126035      102127535      102129035      102129335      102129735 
    ##              1              1              1              1              1 
    ##      102131135      102132335      102133335      102133535      102136635 
    ##              1              1              1              1              1 
    ##      102137335      102139035      102139935      102140135      102140635 
    ##              1              1              1              1              1 
    ##      102141635      102142135      102145335      102145935      102147835 
    ##              1              1              1              1              1 
    ##      102148035      102148235      102148835      102149035      102149435 
    ##              1              1              1              1              1 
    ##      102150335    102151035-2      102151335      102152135      102153135 
    ##              1              1              1              1              1 
    ##      102154035    102155135-2    102155335-2      102156735    102157535-2 
    ##              1              1              1              1              1 
    ##      102157635      102158335    102158835-2      102159635      102160135 
    ##              1              1              1              1              1 
    ##      102160435      102160635      102161635      102161735      102163435 
    ##              1              1              1              1              1 
    ##      102166535      102166835      102170035      102171735      102172335 
    ##              1              1              1              1              1 
    ##      102174135      102174735      102175735      102177635      102178135 
    ##              1              1              1              1              1 
    ##      102179335      102181435      102182035      102183735      102186235 
    ##              1              1              1              1              1 
    ##      102187235      102187335      102189035      102190035      102192035 
    ##              1              1              1              1              1 
    ##    102193835-2    102196135-2      102198435      102198535      102198735 
    ##              1              1              1              1              1 
    ##      102199935      102200635      102202135      102206735      102208035 
    ##              1              1              1              1              1 
    ##      102208735      102209535      102210535      102213135      102216135 
    ##              1              1              1              1              1 
    ##      102217635      102218635      102222035      102223635      102224435 
    ##              1              1              1              1              1 
    ##      102224735      102225335      102226135      102226735      102226935 
    ##              1              1              1              1              1 
    ##      102227435      102230035      102230635      102233135      102234935 
    ##              1              1              1              1              1 
    ##      102236035      102236235      102237035      102237835      102238935 
    ##              1              1              1              1              1 
    ##      102243935      102244335      102244735      102246535      102246635 
    ##              1              1              1              1              1 
    ##      102247335      102248935      102249935      102250935      102252035 
    ##              1              1              1              1              1 
    ##      102253535      102253935      102258435    102258535-2      102258935 
    ##              1              1              1              1              1 
    ##    102263135-2      102264635      102265335      102267135      102267335 
    ##              1              1              1              1              1 
    ##      102268735      102269035      102269635      102270135      102270335 
    ##              1              1              1              1              1 
    ##      102271035      102271235      102271535      102271835      102272435 
    ##              1              1              1              1              1 
    ##      102272535      102272935      102273235      102273535      102275935 
    ##              1              1              1              1              1 
    ##      102276435      102277935      102278035      102279835      102281335 
    ##              1              1              1              1              1 
    ##      102282535      102283735      102284135      102285435      102286335 
    ##              1              1              1              1              1 
    ##      102286735      102287035      102287135      102288235      102288435 
    ##              1              1              1              1              1 
    ##      102288735      102291735      102292535      102294235      102295835 
    ##              1              1              1              1              1 
    ##      102296435      102299135      102299935      102300235      102300335 
    ##              1              1              1              1              1 
    ##      102300535      102302035      102302335      102302935      102303835 
    ##              1              1              1              1              1 
    ##      102303935      102304535      102307635      102311535      102315135 
    ##              1              1              1              1              1 
    ##      102319135      102319335      102321035      102321835      102322135 
    ##              1              1              1              1              1 
    ##      102324535      102326335      102327535      102329135      102329935 
    ##              1              1              1              1              1 
    ##      102330635      102331635      102334335      102334735      102337635 
    ##              1              1              1              1              1 
    ##      102339235      102339435      102341335      102341635      102343635 
    ##              1              1              1              1              1 
    ##      102346535      102348535      102349535      102349635      102350835 
    ##              1              1              1              1              1 
    ##      102350935      102352035      102352935      102353035      102353335 
    ##              1              1              1              1              1 
    ##      102353435      102353535      102353635      102355535      102357735 
    ##              1              1              1              1              1 
    ##      102358335      102361035      102361235      102362435      102363535 
    ##              1              1              1              1              1 
    ##      102364635      102365035      102366635      102366735      102368635 
    ##              1              1              1              1              1 
    ##      102369835      102372335      102373235      102374435      102375335 
    ##              1              1              1              1              1 
    ##      102376035      102377635      102380135    102380935-2      102381135 
    ##              1              1              1              1              1 
    ##      102382835      102384135      102384235      102384535      102386135 
    ##              1              1              1              1              1 
    ##      102386935      102389035      102389135    102389535-2      102394535 
    ##              1              1              1              1              1 
    ##      102396835      102400035    102400835-2    102401035-2      102401335 
    ##              1              1              1              1              1 
    ##    102401635-2   102401735-2W      102403335      102403535      102404035 
    ##              1              1              1              1              1 
    ##      102404135      102404235      102405635      102412135      102412235 
    ##              1              1              1              1              1 
    ##    102413135-2      102489235      102490535      102493235      102493335 
    ##              1              1              1              1              1 
    ##      102495535      102496535      102496735      102497135      102497235 
    ##              1              1              1              1              1 
    ##      102497835      102498235      102499235      102500935      102502135 
    ##              1              1              1              1              1 
    ##      102502435      102503335      102509935      102510035      102510135 
    ##              1              1              1              1              1 
    ##      102510535      102512335      102512735      102513435      102514935 
    ##              1              1              1              1              1 
    ##      102517335      102518235      102524135      102525335      102525835 
    ##              1              1              1              1              1 
    ##      102525935      102527935      102529935      102535735      102536035 
    ##              1              1              1              1              1 
    ##      102545135      102546335      102546935      102549135      102550435 
    ##              1              1              1              1              1 
    ##      102550735      102551035      102551135      102551935      102552735 
    ##              1              1              1              1              1 
    ##   102559435-2W    102559835-2      102560135      102560335      102560435 
    ##              1              1              1              1              1 
    ##      102560935    102562935-2      102564835      102565135      102565235 
    ##              1              1              1              1              1 
    ##      102565435      102566235      102566735      102566835      102566935 
    ##              1              1              1              1              1 
    ##      102567135    102568135-2      102569335      102569635      102571135 
    ##              1              1              1              1              1 
    ##      102571735   102579835-2W      102580135      102581435      102583035 
    ##              1              1              1              1              1 
    ##      102584735      102585335      102588335      102588735      102589735 
    ##              1              1              1              1              1 
    ##      102589835      102590135      102592035      102592435      102592635 
    ##              1              1              1              1              1 
    ##      102604635      102607535      102607635      102608935    102615035-2 
    ##              1              1              1              1              1 
    ##      102615535      102615635    102615935-2      102618235      102629335 
    ##              1              1              1              1              1 
    ##      102630135      102630635      102631535      102640435      102649735 
    ##              1              1              1              1              1 
    ##      102650235      102650735      102651335      102659535      102659735 
    ##              1              1              1              1              1 
    ##      102660335      102662435      102662635 SAMEA115126061 SAMEA115126107 
    ##              1              1              1              1              1 
    ## SAMEA115126171 SAMEA115126306 SAMEA115126326 SAMEA115126491 SAMEA115126500 
    ##              1              1              1              1              1 
    ## SAMEA115126520 SAMEA115126613 SAMEA115126615 SAMEA115126635 SAMEA115126710 
    ##              1              1              1              1              1 
    ## SAMEA115126717 SAMEA115126775 SAMEA115126863 SAMEA115126870 SAMEA115126888 
    ##              1              1              1              1              1 
    ## SAMEA115127046 SAMEA115127059 SAMEA115127165 SAMEA115127331 SAMEA115127484 
    ##              1              1              1              1              1 
    ## SAMEA115127502 SAMEA115127533 SAMEA115127719 
    ##              1              1              1

``` r
table(metadata_clean$Site)
```

    ## 
    ## saliva 
    ##   1698

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ##   NO 
    ## 1698

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 1698"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 1698"

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

    ##   BioProject         Run    sample_name   Age BMI Sex Region       Platform
    ## 1 PRJEB70783 ERR12514043      100000735 35-70  NA  NA  V4–V5 Illumina MiSeq
    ## 2 PRJEB70783 ERR12514044      100002335 35-70  NA  NA  V4–V5 Illumina MiSeq
    ## 3 PRJEB70783 ERR12514045      100002535 35-70  NA  NA  V4–V5 Illumina MiSeq
    ## 4 PRJEB70783 ERR12514047      100003235 35-70  NA  NA  V4–V5 Illumina MiSeq
    ## 5 PRJEB70783 ERR12514048 SAMEA115126061 35-70  NA  NA  V4–V5 Illumina MiSeq
    ## 6 PRJEB70783 ERR12514050      100004535 35-70  NA  NA  V4–V5 Illumina MiSeq
    ##   Country   Site          public_time PCR_primers
    ## 1  Canada saliva 2024-01-23T00:00:00Z   515F-926R
    ## 2  Canada saliva 2024-01-23T00:00:00Z   515F-926R
    ## 3  Canada saliva 2024-01-23T00:00:00Z   515F-926R
    ## 4  Canada saliva 2024-01-23T00:00:00Z   515F-926R
    ## 5  Canada saliva 2024-01-23T00:00:00Z   515F-926R
    ## 6  Canada saliva 2024-01-23T00:00:00Z   515F-926R
    ##                                                                          title
    ## 1 Oral microbial signatures associated with age and frailty in Canadian adults
    ## 2 Oral microbial signatures associated with age and frailty in Canadian adults
    ## 3 Oral microbial signatures associated with age and frailty in Canadian adults
    ## 4 Oral microbial signatures associated with age and frailty in Canadian adults
    ## 5 Oral microbial signatures associated with age and frailty in Canadian adults
    ## 6 Oral microbial signatures associated with age and frailty in Canadian adults
    ##   Disease_study         Status matched_samples
    ## 1           YES  Breast cancer              NO
    ## 2           YES             HC              NO
    ## 3           YES             HC              NO
    ## 4           YES             HC              NO
    ## 5           YES Crohns disease              NO
    ## 6           YES             HC              NO

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
