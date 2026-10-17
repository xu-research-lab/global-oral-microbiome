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
project_id <- "PRJNA767189" 

base_dir <- "E:/oral/26/metadata-2604/Oral-Metadata-YL-02"
setwd(file.path(base_dir, project_id))

# (2) Loading raw SRA metadata file
sra_file <- "SraRunTable.csv"
initial_meta <- read.csv(sra_file, header = TRUE)

# (3) OPTIONAL: If extra metadata (e.g., Table.x or supplemental materials) are available from the article, you may need to merge the file with initial_meta manually. Below is example code for merging. Depending on the specific case of each study, you can customize the code as needed.
# supp_table <- read.csv("supplemental_table.csv", header = TRUE)
# initial_meta <- left_join(initial_meta, supp_table, by = merge_col_name)
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

    ##            Run       Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1  SRR16107365 Targeted-Capture        500 25018368 PRJNA767189 SAMN21892326
    ## 2  SRR16107366 Targeted-Capture        500 19201268 PRJNA767189 SAMN21892941
    ## 3  SRR16107367 Targeted-Capture        501 25540887 PRJNA767189 SAMN21892325
    ## 4  SRR16107368 Targeted-Capture        501 24943649 PRJNA767189 SAMN21892324
    ## 5  SRR16107369 Targeted-Capture        500 24208732 PRJNA767189 SAMN21892189
    ## 6  SRR16107370 Targeted-Capture        500 20122220 PRJNA767189 SAMN21892323
    ## 7  SRR16107371 Targeted-Capture        500 22377940 PRJNA767189 SAMN21892322
    ## 8  SRR16107372 Targeted-Capture        500 26524244 PRJNA767189 SAMN21892321
    ## 9  SRR16107373 Targeted-Capture        500 26524500 PRJNA767189 SAMN21892320
    ## 10 SRR16107374 Targeted-Capture        501 22128081 PRJNA767189 SAMN21892319
    ##                 BioSampleModel    Bytes case_control Center.Name
    ## 1  Metagenome or environmental 14035832      Control         NCI
    ## 2  Metagenome or environmental 10276724         Case         NCI
    ## 3  Metagenome or environmental 14414969         Case         NCI
    ## 4  Metagenome or environmental 13894777         Case         NCI
    ## 5  Metagenome or environmental 13711242         Case         NCI
    ## 6  Metagenome or environmental 11198327         Case         NCI
    ## 7  Metagenome or environmental 12255838      Control         NCI
    ## 8  Metagenome or environmental 14726084         Case         NCI
    ## 9  Metagenome or environmental 14911026      Control         NCI
    ## 10 Metagenome or environmental 12971259         Case         NCI
    ##    Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 7               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 8               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 9               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 10              NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                        DATASTORE.region  Experiment geo_loc_name_country
    ## 1  gs.us-east1,ncbi.public,s3.us-east-1 SRX12393312         uncalculated
    ## 2  gs.us-east1,ncbi.public,s3.us-east-1 SRX12393886         uncalculated
    ## 3  gs.us-east1,ncbi.public,s3.us-east-1 SRX12393311         uncalculated
    ## 4  gs.us-east1,ncbi.public,s3.us-east-1 SRX12393310         uncalculated
    ## 5  gs.us-east1,ncbi.public,s3.us-east-1 SRX12393309         uncalculated
    ## 6  gs.us-east1,ncbi.public,s3.us-east-1 SRX12393308         uncalculated
    ## 7  gs.us-east1,ncbi.public,s3.us-east-1 SRX12393307         uncalculated
    ## 8  gs.us-east1,ncbi.public,s3.us-east-1 SRX12393306         uncalculated
    ## 9  gs.us-east1,ncbi.public,s3.us-east-1 SRX12393305         uncalculated
    ## 10 gs.us-east1,ncbi.public,s3.us-east-1 SRX12393304         uncalculated
    ##    geo_loc_name_country_continent geo_loc_name         HOST     Instrument
    ## 1                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 2                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 3                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 4                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 5                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 6                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 7                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 8                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 9                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 10                   uncalculated           NA Homo sapiens Illumina MiSeq
    ##    lat_lon lib_const_method Library.Name LibraryLayout LibrarySelection
    ## 1       NA       paired-end    GBHS15271        PAIRED              PCR
    ## 2       NA       paired-end    GBHS57749        PAIRED              PCR
    ## 3       NA       paired-end    GBHS79343        PAIRED              PCR
    ## 4       NA       paired-end    GBHS96121        PAIRED              PCR
    ## 5       NA       paired-end    GBHS64114        PAIRED              PCR
    ## 6       NA       paired-end    GBHS26790        PAIRED              PCR
    ## 7       NA       paired-end    GBHS86104        PAIRED              PCR
    ## 8       NA       paired-end    GBHS49272        PAIRED              PCR
    ## 9       NA       paired-end    GBHS65505        PAIRED              PCR
    ## 10      NA       paired-end    GBHS25078        PAIRED              PCR
    ##    LibrarySource material   Organism Platform          ReleaseDate
    ## 1        GENOMIC   Saliva metagenome ILLUMINA 2021-12-01T00:00:00Z
    ## 2        GENOMIC   Saliva metagenome ILLUMINA 2021-12-01T00:00:00Z
    ## 3        GENOMIC   Saliva metagenome ILLUMINA 2021-12-01T00:00:00Z
    ## 4        GENOMIC   Saliva metagenome ILLUMINA 2021-12-01T00:00:00Z
    ## 5        GENOMIC   Saliva metagenome ILLUMINA 2021-12-01T00:00:00Z
    ## 6        GENOMIC   Saliva metagenome ILLUMINA 2021-12-01T00:00:00Z
    ## 7        GENOMIC   Saliva metagenome ILLUMINA 2021-12-01T00:00:00Z
    ## 8        GENOMIC   Saliva metagenome ILLUMINA 2021-12-01T00:00:00Z
    ## 9        GENOMIC   Saliva metagenome ILLUMINA 2021-12-01T00:00:00Z
    ## 10       GENOMIC   Saliva metagenome ILLUMINA 2021-12-01T00:00:00Z
    ##             create_date version Sample.Name sample_type SampleID2
    ## 1  2021-09-29T01:48:00Z       1    SC070332     Control GBHS15271
    ## 2  2021-09-29T01:48:00Z       1    SC097882      Cancer GBHS57749
    ## 3  2021-09-29T01:48:00Z       1    SC073825      Cancer GBHS79343
    ## 4  2021-09-29T01:48:00Z       1    SC103740      Cancer GBHS96121
    ## 5  2021-09-29T01:48:00Z       1    SC072309      Cancer GBHS64114
    ## 6  2021-09-29T01:49:00Z       1    SC073422      Cancer GBHS26790
    ## 7  2021-09-29T01:48:00Z       1    SC070856     Control GBHS86104
    ## 8  2021-09-29T01:51:00Z       1    SC104113      Cancer GBHS49272
    ## 9  2021-09-29T01:48:00Z       1    SC097715     Control GBHS65505
    ## 10 2021-09-29T01:48:00Z       1    SC106121      Cancer GBHS25078
    ##       seq_methods SRA.Study target_gene target_subfragment CONSENSUS_DIAG
    ## 1  Illumina MiSeq SRP339204    16S rRNA                 V4               
    ## 2  Illumina MiSeq SRP339204    16S rRNA                 V4     1)Invasive
    ## 3  Illumina MiSeq SRP339204    16S rRNA                 V4     1)Invasive
    ## 4  Illumina MiSeq SRP339204    16S rRNA                 V4     1)Invasive
    ## 5  Illumina MiSeq SRP339204    16S rRNA                 V4     1)Invasive
    ## 6  Illumina MiSeq SRP339204    16S rRNA                 V4     1)Invasive
    ## 7  Illumina MiSeq SRP339204    16S rRNA                 V4               
    ## 8  Illumina MiSeq SRP339204    16S rRNA                 V4     1)Invasive
    ## 9  Illumina MiSeq SRP339204    16S rRNA                 V4               
    ## 10 Illumina MiSeq SRP339204    16S rRNA                 V4     1)Invasive

``` r
# IMPORTANT: Check the unique values in the column of interest or query column.
lapply(initial_meta[c('material')], table) # site
```

    ## $material
    ## 
    ## Saliva 
    ##    947

``` r
lapply(initial_meta[c('HOST', 'Organism', 'sample_type')], table) # organism
```

    ## $HOST
    ## 
    ## Homo sapiens 
    ##          947 
    ## 
    ## $Organism
    ## 
    ## metagenome 
    ##        947 
    ## 
    ## $sample_type
    ## 
    ##           Benign           Cancer          Control Repeated Control 
    ##              104              388              411               44

``` r
table(initial_meta$Sample.Name) # subject id
```

    ## 
    ##              SC000087              SC000091              SC000092 
    ##                     1                     1                     1 
    ##              SC000095              SC000096              SC000097 
    ##                     1                     1                     1 
    ##              SC000098              SC000287              SC000331 
    ##                     1                     1                     1 
    ##              SC000338              SC000363              SC000368 
    ##                     1                     1                     1 
    ##              SC000373              SC000395              SC000396 
    ##                     1                     1                     1 
    ##              SC000401              SC000413              SC000415 
    ##                     1                     1                     1 
    ##              SC000416              SC000417              SC000447 
    ##                     1                     1                     1 
    ##              SC000450              SC000453              SC000465 
    ##                     1                     1                     1 
    ##              SC000468              SC000469              SC000476 
    ##                     1                     1                     1 
    ##              SC000477              SC000478              SC000481 
    ##                     1                     1                     1 
    ##              SC000483              SC068964              SC068965 
    ##                     1                     1                     1 
    ##              SC068973              SC068976              SC068978 
    ##                     1                     1                     1 
    ##              SC069103              SC069104              SC069107 
    ##                     1                     1                     1 
    ##              SC069113              SC069121              SC069122 
    ##                     1                     1                     1 
    ##              SC069123              SC069125              SC069127 
    ##                     1                     1                     1 
    ##              SC069128              SC069129              SC069130 
    ##                     1                     1                     1 
    ##              SC069137              SC069139              SC069141 
    ##                     1                     1                     1 
    ##              SC069142              SC069145              SC069149 
    ##                     1                     1                     1 
    ##              SC069150              SC069154              SC069155 
    ##                     1                     1                     1 
    ##              SC069166              SC069174              SC069175 
    ##                     1                     1                     1 
    ##              SC069179              SC069181              SC069189 
    ##                     1                     1                     1 
    ##              SC069195              SC069196              SC069197 
    ##                     1                     1                     1 
    ##              SC069201              SC069203              SC069646 
    ##                     1                     1                     1 
    ##              SC069647              SC069648              SC069654 
    ##                     1                     1                     1 
    ##              SC069656              SC069660              SC069664 
    ##                     1                     1                     1 
    ##              SC069666              SC069669              SC069670 
    ##                     1                     1                     1 
    ##              SC069673              SC069674              SC069677 
    ##                     1                     1                     1 
    ##              SC069679              SC069681              SC069690 
    ##                     1                     1                     1 
    ##              SC069766              SC069769              SC069780 
    ##                     1                     1                     1 
    ##              SC069783              SC069786              SC069787 
    ##                     1                     1                     1 
    ##              SC069789              SC069790              SC069794 
    ##                     1                     1                     1 
    ##              SC069798              SC069925              SC069927 
    ##                     1                     1                     1 
    ##              SC069928              SC069932              SC069933 
    ##                     1                     1                     1 
    ##              SC069934              SC069936              SC069938 
    ##                     1                     1                     1 
    ##              SC069939              SC069942              SC069945 
    ##                     1                     1                     1 
    ##              SC069946              SC069947              SC069948 
    ##                     1                     1                     1 
    ## SC069950-PC30190-A-07 SC069950-PC30190-E-12              SC069951 
    ##                     1                     1                     1 
    ##              SC069952              SC069953              SC069955 
    ##                     1                     1                     1 
    ##              SC069956              SC069957              SC069958 
    ##                     1                     1                     1 
    ##              SC069959              SC069962              SC069966 
    ##                     1                     1                     1 
    ##              SC069973              SC069979              SC069980 
    ##                     1                     1                     1 
    ##              SC069983              SC069988              SC070061 
    ##                     1                     1                     1 
    ##              SC070065              SC070066              SC070067 
    ##                     1                     1                     1 
    ##              SC070072              SC070076              SC070079 
    ##                     1                     1                     1 
    ##              SC070080              SC070325              SC070326 
    ##                     1                     1                     1 
    ##              SC070327              SC070329              SC070332 
    ##                     1                     1                     1 
    ##              SC070338              SC070343              SC070856 
    ##                     1                     1                     1 
    ##              SC070859              SC070864              SC070890 
    ##                     1                     1                     1 
    ##              SC070979              SC070982              SC070985 
    ##                     1                     1                     1 
    ##              SC070989              SC070995              SC071056 
    ##                     1                     1                     1 
    ##              SC071057              SC071060              SC071061 
    ##                     1                     1                     1 
    ##              SC071063              SC071064              SC071065 
    ##                     1                     1                     1 
    ##              SC071066              SC071068              SC071073 
    ##                     1                     1                     1 
    ##              SC071074              SC071075              SC071081 
    ##                     1                     1                     1 
    ##              SC071085              SC071087              SC071090 
    ##                     1                     1                     1 
    ##              SC071093              SC071096              SC071097 
    ##                     1                     1                     1 
    ##              SC071100              SC071102              SC071103 
    ##                     1                     1                     1 
    ##              SC071104              SC071106              SC071107 
    ##                     1                     1                     1 
    ##              SC071109              SC071117              SC071121 
    ##                     1                     1                     1 
    ##              SC071122              SC071123              SC071126 
    ##                     1                     1                     1 
    ##              SC071128              SC071129              SC071136 
    ##                     1                     1                     1 
    ##              SC071138              SC071140              SC071142 
    ##                     1                     1                     1 
    ##              SC071145              SC071148              SC071150 
    ##                     1                     1                     1 
    ##              SC071151              SC071157              SC071160 
    ##                     1                     1                     1 
    ##              SC071164              SC071165              SC071167 
    ##                     1                     1                     1 
    ##              SC071330              SC071331              SC071334 
    ##                     1                     1                     1 
    ##              SC071338              SC071342              SC071345 
    ##                     1                     1                     1 
    ##              SC071347              SC071353              SC071358 
    ##                     1                     1                     1 
    ##              SC071359              SC071363              SC071364 
    ##                     1                     1                     1 
    ##              SC071612              SC071613              SC071614 
    ##                     1                     1                     1 
    ##              SC071621              SC071624              SC071625 
    ##                     1                     1                     1 
    ##              SC071627              SC071628              SC071629 
    ##                     1                     1                     1 
    ##              SC071630              SC071632              SC071633 
    ##                     1                     1                     1 
    ##              SC071635              SC071640              SC071641 
    ##                     1                     1                     1 
    ##              SC071644              SC071645              SC071647 
    ##                     1                     1                     1 
    ##              SC071651              SC071652              SC071654 
    ##                     1                     1                     1 
    ##              SC071656              SC071660              SC071662 
    ##                     1                     1                     1 
    ##              SC071665              SC071666              SC071667 
    ##                     1                     1                     1 
    ##              SC071668              SC071669              SC071671 
    ##                     1                     1                     1 
    ##              SC071672              SC071674              SC071676 
    ##                     1                     1                     1 
    ##              SC071689              SC071691              SC071697 
    ##                     1                     1                     1 
    ##              SC071698              SC071746              SC071748 
    ##                     1                     1                     1 
    ##              SC071845              SC071849              SC071852 
    ##                     1                     1                     1 
    ##              SC071853              SC071857              SC071858 
    ##                     1                     1                     1 
    ##              SC071865              SC071867              SC071872 
    ##                     1                     1                     1 
    ##              SC071897              SC071902              SC071904 
    ##                     1                     1                     1 
    ##              SC071908              SC071909              SC071912 
    ##                     1                     1                     1 
    ##              SC071925              SC071926              SC071928 
    ##                     1                     1                     1 
    ##              SC071931              SC071938              SC071940 
    ##                     1                     1                     1 
    ##              SC071942              SC071944              SC071946 
    ##                     1                     1                     1 
    ##              SC071951              SC071955              SC071957 
    ##                     1                     1                     1 
    ##              SC071969              SC071980              SC071991 
    ##                     1                     1                     1 
    ##              SC071994              SC071995              SC071999 
    ##                     1                     1                     1 
    ##              SC072001              SC072003              SC072004 
    ##                     1                     1                     1 
    ##              SC072010              SC072011              SC072016 
    ##                     1                     1                     1 
    ##              SC072017              SC072018              SC072020 
    ##                     1                     1                     1 
    ##              SC072304              SC072309              SC072333 
    ##                     1                     1                     1 
    ##              SC072337              SC072340              SC072343 
    ##                     1                     1                     1 
    ##              SC072345              SC072350              SC072359 
    ##                     1                     1                     1 
    ##              SC072373              SC072378              SC072380 
    ##                     1                     1                     1 
    ##              SC072381              SC072384              SC072385 
    ##                     1                     1                     1 
    ##              SC072391              SC073422              SC073441 
    ##                     1                     1                     1 
    ##              SC073528              SC073530              SC073562 
    ##                     1                     1                     1 
    ##              SC073565              SC073755              SC073759 
    ##                     1                     1                     1 
    ##              SC073761              SC073763              SC073766 
    ##                     1                     1                     1 
    ## SC073769-PC30188-A-12 SC073769-PC30188-B-05              SC073773 
    ##                     1                     1                     1 
    ##              SC073776              SC073809              SC073818 
    ##                     1                     1                     1 
    ##              SC073819              SC073825              SC073827 
    ##                     1                     1                     1 
    ##              SC073831              SC073832              SC073839 
    ##                     1                     1                     1 
    ##              SC073841              SC073844              SC073845 
    ##                     1                     1                     1 
    ##              SC073854              SC073856 SC073857-PC30212-A-11 
    ##                     1                     1                     1 
    ## SC073857-PC30212-B-08 SC073858-PC30214-E-11 SC073858-PC30214-G-04 
    ##                     1                     1                     1 
    ##              SC073859              SC073861 SC073862-PC30212-A-03 
    ##                     1                     1                     1 
    ## SC073862-PC30212-E-06              SC073863              SC073865 
    ##                     1                     1                     1 
    ##              SC073866              SC073868              SC073870 
    ##                     1                     1                     1 
    ##              SC073872              SC073873              SC073876 
    ##                     1                     1                     1 
    ##              SC073879              SC073880              SC073888 
    ##                     1                     1                     1 
    ##              SC073890              SC073892              SC073893 
    ##                     1                     1                     1 
    ##              SC073894              SC073902              SC073903 
    ##                     1                     1                     1 
    ##              SC073904              SC073905              SC073906 
    ##                     1                     1                     1 
    ##              SC073907              SC073911              SC073914 
    ##                     1                     1                     1 
    ##              SC073915              SC073916              SC073919 
    ##                     1                     1                     1 
    ##              SC073920              SC073921              SC073922 
    ##                     1                     1                     1 
    ##              SC073923              SC073927              SC073931 
    ##                     1                     1                     1 
    ##              SC073932              SC073979              SC073982 
    ##                     1                     1                     1 
    ##              SC073983 SC073986-PC31744-C-01 SC073986-PC31744-F-04 
    ##                     1                     1                     1 
    ##              SC073988              SC074042              SC074045 
    ##                     1                     1                     1 
    ##              SC074087              SC074089              SC074090 
    ##                     1                     1                     1 
    ##              SC074092              SC074093              SC074098 
    ##                     1                     1                     1 
    ##              SC074101              SC074102              SC074106 
    ##                     1                     1                     1 
    ##              SC074117              SC074118              SC074122 
    ##                     1                     1                     1 
    ##              SC097697              SC097702              SC097708 
    ##                     1                     1                     1 
    ##              SC097710              SC097711              SC097715 
    ##                     1                     1                     1 
    ##              SC097718              SC097722              SC097723 
    ##                     1                     1                     1 
    ##              SC097726              SC097727              SC097733 
    ##                     1                     1                     1 
    ##              SC097737              SC097740              SC097741 
    ##                     1                     1                     1 
    ##              SC097743 SC097745-PC30188-F-05 SC097745-PC30188-F-07 
    ##                     1                     1                     1 
    ##              SC097758              SC097762              SC097763 
    ##                     1                     1                     1 
    ##              SC097765              SC097766              SC097767 
    ##                     1                     1                     1 
    ##              SC097768              SC097770              SC097771 
    ##                     1                     1                     1 
    ##              SC097773              SC097774              SC097775 
    ##                     1                     1                     1 
    ##              SC097778              SC097783              SC097786 
    ##                     1                     1                     1 
    ##              SC097789              SC097791              SC097792 
    ##                     1                     1                     1 
    ##              SC097798              SC097799              SC097867 
    ##                     1                     1                     1 
    ##              SC097869              SC097870              SC097873 
    ##                     1                     1                     1 
    ##              SC097874              SC097877              SC097879 
    ##                     1                     1                     1 
    ##              SC097882              SC097883              SC097884 
    ##                     1                     1                     1 
    ##              SC097887              SC097888              SC097899 
    ##                     1                     1                     1 
    ##              SC097906              SC097944              SC097947 
    ##                     1                     1                     1 
    ##              SC097951              SC097953              SC097956 
    ##                     1                     1                     1 
    ##              SC097957              SC097958              SC097959 
    ##                     1                     1                     1 
    ##              SC097963              SC097964              SC098006 
    ##                     1                     1                     1 
    ##              SC098010              SC098012              SC098017 
    ##                     1                     1                     1 
    ##              SC098938              SC098940              SC098941 
    ##                     1                     1                     1 
    ##              SC098945              SC099233              SC099234 
    ##                     1                     1                     1 
    ##              SC099235              SC099240              SC100015 
    ##                     1                     1                     1 
    ##              SC100424              SC100426              SC100427 
    ##                     1                     1                     1 
    ##              SC100532              SC100541              SC100546 
    ##                     1                     1                     1 
    ##              SC100550              SC100552              SC100553 
    ##                     1                     1                     1 
    ##              SC100645              SC100646              SC100647 
    ##                     1                     1                     1 
    ##              SC100648              SC100650              SC100654 
    ##                     1                     1                     1 
    ##              SC100655              SC100656              SC100705 
    ##                     1                     1                     1 
    ##              SC100706              SC100709              SC100710 
    ##                     1                     1                     1 
    ##              SC100711              SC100716              SC100719 
    ##                     1                     1                     1 
    ##              SC100720              SC100721 SC100722-PC30214-B-01 
    ##                     1                     1                     1 
    ## SC100722-PC30214-B-08              SC100723              SC100732 
    ##                     1                     1                     1 
    ##              SC101622              SC101629              SC101632 
    ##                     1                     1                     1 
    ##              SC101633              SC101634              SC101639 
    ##                     1                     1                     1 
    ##              SC101642              SC101643              SC101648 
    ##                     1                     1                     1 
    ##              SC101650              SC101652              SC101653 
    ##                     1                     1                     1 
    ##              SC101654              SC101660              SC101665 
    ##                     1                     1                     1 
    ##              SC101666 SC101668-PC30190-B-09 SC101668-PC30190-F-06 
    ##                     1                     1                     1 
    ##              SC101685              SC101738              SC101811 
    ##                     1                     1                     1 
    ##              SC101815              SC101820              SC101821 
    ##                     1                     1                     1 
    ##              SC101824              SC101825              SC103538 
    ##                     1                     1                     1 
    ##              SC103539              SC103541              SC103544 
    ##                     1                     1                     1 
    ##              SC103579              SC103582              SC103583 
    ##                     1                     1                     1 
    ##              SC103584              SC103585              SC103602 
    ##                     1                     1                     1 
    ##              SC103603 SC103605-PC31763-B-04 SC103605-PC31763-F-10 
    ##                     1                     1                     1 
    ##              SC103608              SC103609              SC103610 
    ##                     1                     1                     1 
    ##              SC103611              SC103612              SC103613 
    ##                     1                     1                     1 
    ##              SC103627              SC103631              SC103634 
    ##                     1                     1                     1 
    ##              SC103635              SC103636              SC103658 
    ##                     1                     1                     1 
    ##              SC103661              SC103662              SC103663 
    ##                     1                     1                     1 
    ##              SC103664              SC103665              SC103666 
    ##                     1                     1                     1 
    ##              SC103667              SC103668              SC103669 
    ##                     1                     1                     1 
    ##              SC103670              SC103671              SC103672 
    ##                     1                     1                     1 
    ##              SC103675              SC103678              SC103679 
    ##                     1                     1                     1 
    ##              SC103683              SC103688              SC103690 
    ##                     1                     1                     1 
    ##              SC103740              SC103749              SC103751 
    ##                     1                     1                     1 
    ##              SC103753              SC103754              SC103756 
    ##                     1                     1                     1 
    ##              SC103757              SC103765              SC103766 
    ##                     1                     1                     1 
    ##              SC103768              SC103770              SC103772 
    ##                     1                     1                     1 
    ##              SC103775              SC103786              SC103787 
    ##                     1                     1                     1 
    ##              SC103790              SC103791              SC103792 
    ##                     1                     1                     1 
    ##              SC103793              SC103794              SC103796 
    ##                     1                     1                     1 
    ##              SC103797              SC103802              SC103804 
    ##                     1                     1                     1 
    ##              SC103805              SC103806              SC103807 
    ##                     1                     1                     1 
    ##              SC103808              SC103809              SC103810 
    ##                     1                     1                     1 
    ##              SC103813              SC103814              SC103817 
    ##                     1                     1                     1 
    ##              SC103819              SC103821              SC103822 
    ##                     1                     1                     1 
    ##              SC103830              SC103832              SC103918 
    ##                     1                     1                     1 
    ##              SC103920              SC103921              SC103924 
    ##                     1                     1                     1 
    ##              SC104108              SC104113              SC104116 
    ##                     1                     1                     1 
    ##              SC104121              SC104126              SC104127 
    ##                     1                     1                     1 
    ##              SC104129              SC104142              SC104143 
    ##                     1                     1                     1 
    ##              SC104324              SC104335              SC104336 
    ##                     1                     1                     1 
    ##              SC104339              SC104350              SC104351 
    ##                     1                     1                     1 
    ##              SC104363              SC104365              SC104367 
    ##                     1                     1                     1 
    ##              SC104369              SC104370              SC104371 
    ##                     1                     1                     1 
    ##              SC104373 SC104375-PC31744-E-08 SC104375-PC31744-F-05 
    ##                     1                     1                     1 
    ##              SC104376              SC104377              SC104379 
    ##                     1                     1                     1 
    ##              SC104380              SC104382              SC104384 
    ##                     1                     1                     1 
    ##              SC104386              SC104387              SC104390 
    ##                     1                     1                     1 
    ##              SC104397              SC104404              SC104405 
    ##                     1                     1                     1 
    ##              SC104408              SC104410              SC104412 
    ##                     1                     1                     1 
    ##              SC104423              SC104425              SC104430 
    ##                     1                     1                     1 
    ##              SC104433              SC104434              SC104435 
    ##                     1                     1                     1 
    ##              SC104437              SC104438              SC104439 
    ##                     1                     1                     1 
    ##              SC104445              SC104446              SC104450 
    ##                     1                     1                     1 
    ##              SC104451              SC104454              SC104455 
    ##                     1                     1                     1 
    ##              SC104458              SC104462              SC104483 
    ##                     1                     1                     1 
    ##              SC104484              SC104486              SC104666 
    ##                     1                     1                     1 
    ##              SC104669              SC104671              SC104796 
    ##                     1                     1                     1 
    ##              SC104798              SC104800              SC104801 
    ##                     1                     1                     1 
    ##              SC104803              SC104804              SC104806 
    ##                     1                     1                     1 
    ##              SC104807              SC104809              SC104810 
    ##                     1                     1                     1 
    ##              SC104812              SC104813              SC104815 
    ##                     1                     1                     1 
    ##              SC104816              SC106034              SC106035 
    ##                     1                     1                     1 
    ##              SC106048              SC106054              SC106055 
    ##                     1                     1                     1 
    ##              SC106064              SC106067              SC106071 
    ##                     1                     1                     1 
    ##              SC106073              SC106079              SC106080 
    ##                     1                     1                     1 
    ##              SC106081              SC106086              SC106087 
    ##                     1                     1                     1 
    ##              SC106088              SC106092              SC106093 
    ##                     1                     1                     1 
    ## SC106094-PC31762-A-08 SC106094-PC31762-D-04              SC106095 
    ##                     1                     1                     1 
    ##              SC106097              SC106098              SC106104 
    ##                     1                     1                     1 
    ##              SC106106              SC106113              SC106120 
    ##                     1                     1                     1 
    ##              SC106121              SC106125              SC106126 
    ##                     1                     1                     1 
    ##              SC106131              SC106133              SC106134 
    ##                     1                     1                     1 
    ##              SC106140              SC106141              SC106142 
    ##                     1                     1                     1 
    ##              SC106144              SC106145              SC106147 
    ##                     1                     1                     1 
    ##              SC106148              SC106149              SC106393 
    ##                     1                     1                     1 
    ##              SC106394              SC106396              SC106403 
    ##                     1                     1                     1 
    ##              SC106404              SC106407              SC106408 
    ##                     1                     1                     1 
    ##              SC106410 SC106414-PC31860-A-01 SC106414-PC31860-E-09 
    ##                     1                     1                     1 
    ##              SC106435              SC106437              SC106439 
    ##                     1                     1                     1 
    ##              SC106441              SC106442              SC106447 
    ##                     1                     1                     1 
    ##              SC106451              SC106452              SC106454 
    ##                     1                     1                     1 
    ##              SC106455              SC106456              SC106457 
    ##                     1                     1                     1 
    ##              SC106459              SC106461              SC106462 
    ##                     1                     1                     1 
    ##              SC106463              SC106464              SC106465 
    ##                     1                     1                     1 
    ##              SC106466              SC106467              SC106468 
    ##                     1                     1                     1 
    ##              SC106469              SC106470              SC106472 
    ##                     1                     1                     1 
    ##              SC106473              SC106476              SC106477 
    ##                     1                     1                     1 
    ## SC106478-PC31856-D-12 SC106478-PC31856-G-02              SC106481 
    ##                     1                     1                     1 
    ##              SC106482              SC106483              SC106486 
    ##                     1                     1                     1 
    ##              SC106488              SC106489              SC106492 
    ##                     1                     1                     1 
    ##              SC106493              SC106494              SC106495 
    ##                     1                     1                     1 
    ##              SC106496              SC106497              SC106498 
    ##                     1                     1                     1 
    ##              SC106504              SC106507 SC106631-PC31856-B-06 
    ##                     1                     1                     1 
    ## SC106631-PC31856-B-08              SC106632              SC106633 
    ##                     1                     1                     1 
    ##              SC106634              SC106635              SC106636 
    ##                     1                     1                     1 
    ##              SC106637              SC106638              SC106639 
    ##                     1                     1                     1 
    ##              SC106640              SC106642              SC106643 
    ##                     1                     1                     1 
    ##              SC106644              SC106645              SC106646 
    ##                     1                     1                     1 
    ##              SC106650              SC106651              SC106653 
    ##                     1                     1                     1 
    ##              SC106654              SC106655              SC106657 
    ##                     1                     1                     1 
    ##              SC106659              SC106660              SC106662 
    ##                     1                     1                     1 
    ##              SC106664              SC106665 SC106666-PC31745-A-10 
    ##                     1                     1                     1 
    ## SC106666-PC31745-D-09              SC106667              SC106668 
    ##                     1                     1                     1 
    ##              SC106669              SC106670              SC106671 
    ##                     1                     1                     1 
    ##              SC106676              SC106677              SC106687 
    ##                     1                     1                     1 
    ##              SC106689              SC106691              SC106692 
    ##                     1                     1                     1 
    ##              SC106694              SC106695              SC106697 
    ##                     1                     1                     1 
    ## SC106698-PC31763-B-05 SC106698-PC31763-C-06              SC106704 
    ##                     1                     1                     1 
    ##              SC106717              SC106718              SC106722 
    ##                     1                     1                     1 
    ##              SC106723              SC106726              SC106728 
    ##                     1                     1                     1 
    ##              SC106731              SC106734              SC106736 
    ##                     1                     1                     1 
    ##              SC106763              SC106765              SC106768 
    ##                     1                     1                     1 
    ##              SC106772              SC106774              SC106775 
    ##                     1                     1                     1 
    ##              SC106776              SC106778              SC106780 
    ##                     1                     1                     1 
    ##              SC106782              SC106787              SC106791 
    ##                     1                     1                     1 
    ##              SC106794              SC106838              SC106841 
    ##                     1                     1                     1 
    ##              SC106844              SC106846              SC106848 
    ##                     1                     1                     1 
    ##              SC107004              SC107006              SC107007 
    ##                     1                     1                     1 
    ##              SC107010              SC107011              SC107012 
    ##                     1                     1                     1 
    ## SC107014-PC31857-B-01 SC107014-PC31857-D-05              SC107015 
    ##                     1                     1                     1 
    ##              SC107016              SC107017              SC107019 
    ##                     1                     1                     1 
    ##              SC107022              SC107023              SC107024 
    ##                     1                     1                     1 
    ##              SC107025              SC107027              SC107028 
    ##                     1                     1                     1 
    ##              SC107029              SC107030              SC107031 
    ##                     1                     1                     1 
    ##              SC107033              SC107036              SC107039 
    ##                     1                     1                     1 
    ##              SC107040              SC107041              SC107042 
    ##                     1                     1                     1 
    ##              SC107043              SC107044              SC107049 
    ##                     1                     1                     1 
    ## SC107050-PC31860-A-05 SC107050-PC31860-F-10              SC107054 
    ##                     1                     1                     1 
    ## SC107056-PC31857-D-01 SC107056-PC31857-D-04              SC107058 
    ##                     1                     1                     1 
    ##              SC107061              SC107062              SC107065 
    ##                     1                     1                     1 
    ##              SC107066              SC107068              SC107069 
    ##                     1                     1                     1 
    ##              SC107071              SC107075              SC107077 
    ##                     1                     1                     1 
    ##              SC107084              SC107087              SC107089 
    ##                     1                     1                     1 
    ##              SC107090              SC107096              SC107097 
    ##                     1                     1                     1 
    ##              SC107098              SC107107              SC107115 
    ##                     1                     1                     1 
    ##              SC107116              SC107121              SC107123 
    ##                     1                     1                     1 
    ##              SC107131              SC107137              SC107138 
    ##                     1                     1                     1 
    ##              SC107146              SC107154              SC107249 
    ##                     1                     1                     1 
    ##              SC107257              SC107316              SC107320 
    ##                     1                     1                     1 
    ##              SC107322              SC107324              SC107330 
    ##                     1                     1                     1 
    ##              SC107331              SC107332 SC107334-PC31762-C-01 
    ##                     1                     1                     1 
    ## SC107334-PC31762-E-09              SC107335              SC107336 
    ##                     1                     1                     1 
    ##              SC107337              SC107339              SC107363 
    ##                     1                     1                     1 
    ##              SC107364 SC107365-PC31745-C-12 SC107365-PC31745-D-02 
    ##                     1                     1                     1 
    ##              SC107366              SC107368              SC107372 
    ##                     1                     1                     1 
    ##              SC107373              SC107381              SC107383 
    ##                     1                     1                     1 
    ##              SC107386              SC107415              SC107418 
    ##                     1                     1                     1 
    ##              SC107421              SC107427              SC107429 
    ##                     1                     1                     1 
    ##              SC107431              SC107433 
    ##                     1                     1

``` r
#table(initial_meta$sequencing_meth..exp.) # output is Sequencing by synthesis
table(initial_meta$Platform) # output is Illumina. according to the methods,  Illumina is correct.
```

    ## 
    ## ILLUMINA 
    ##      947

``` r
table(initial_meta$target_subfragment) # output is V4
```

    ## 
    ##  V4 
    ## 947

# Step 3: Removing Control and Non-human Samples

First, we need to remove control samples(including both negative
controls and positive controls) and non-human samples.

``` r
# (1) EDIT HERE: Parameter setup
control_column <- "sample_type" # column name in initial_meta
# unique(initial_meta[[control_column]])
control_keywords <- c("Repeated Control") # Enter exact keywords. Article: Within each PCR  batch, eight QC samples were also included: two oral artificial  communities,16 two water blanks, two robogut samples16 and two  duplicates of randomly selected control samples.

species_column <- "CONSENSUS_DIAG" # column name in initial_meta
unique(initial_meta[[species_column]])
```

    ## [1] ""                  "1)Invasive"        "3)Benign"         
    ## [4] "4)Other Case Type" "2)In-situ"

``` r
non_human_keywords <- c("2)In-situ", "4)Other Case Type") # Enter exact keywords. Article: Within each PCR  batch, eight QC samples were also included: two oral artificial  communities,16 two water blanks, two robogut samples16 and two  duplicates of randomly selected control samples.

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

    ## Initial samples counts: 947

``` r
cat("Removed control samples:", sum(is_control), "\n")
```

    ## Removed control samples: 44

``` r
cat("Removed non-human samples:", sum(is_non_human), "\n")
```

    ## Removed non-human samples: 5

``` r
cat("Total samples removed ", sum(samples_to_remove), "\n")
```

    ## Total samples removed  49

``` r
cat("Remaining samples after Step 2:", nrow(initial_meta), "\n")
```

    ## Remaining samples after Step 2: 898

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
  raw_age = NA,
  age = NA, 
  age_info = "18-74", # extract from methods section 
  bmi = NA,
  sex = "F", # extract from abstract section
  sequencing_region = "V4", # extract from methods section
  platform = "Illumina MiSeq", # extract from methods section, NOT initial_meta
  country = "Ghana",
  raw_site = initial_meta$material, 
  site = "pending", # modify later
  site_group = "pending", # modify later
  public_time = initial_meta$ReleaseDate,
  pcr_primers = NA, # Format: xxF-xxR
  title = "The oral microbiome and breast cancer and nonmalignant breast disease, and its relationship with the fecal microbiome in the Ghana breast health study", # extract from article
  disease_study = "yes", # determine according to the article
  raw_status = initial_meta$sample_type,
  status = "pending", # modify later 
  matched_samples = "pending", # modify later
  is_time_series = "no", 
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
  mutate(site = case_when(
    raw_site == "Saliva" ~ "saliva",
    TRUE ~ NA_character_
  )) # according to the 'material' column of initial_meta

metadata <- metadata %>%
  mutate(status = case_when(
    grepl("Cancer", raw_status) ~ "breast cancer",
    grepl("Benign", raw_status) ~ "nonmalignant breast disease",
    grepl("Control", raw_status) ~ "HC",
    TRUE ~ NA_character_
  )) # according to the 'sample_type' column of initial_meta
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
