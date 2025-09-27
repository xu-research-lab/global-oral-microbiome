```R
library(dplyr)
## visualization
library(ggplot2)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(cowplot)
```

    
    Attaching package: ‘dplyr’
    
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    
    
    Attaching package: ‘gridExtra’
    
    
    The following object is masked from ‘package:dplyr’:
    
        combine
    
    



```R
setwd("/home/liuxuan/metadata/PRJEB20213")
initial_meta <- read.csv("SraRunTable.csv", header = TRUE)
head(initial_meta)
colnames(initial_meta)
```


<table class="dataframe">
<caption>A data.frame: 6 × 54</caption>
<thead>
	<tr><th></th><th scope=col>Run</th><th scope=col>Assay.Type</th><th scope=col>AvgSpotLen</th><th scope=col>Bases</th><th scope=col>BioProject</th><th scope=col>BioSample</th><th scope=col>Bytes</th><th scope=col>Center.Name</th><th scope=col>Collection_Date</th><th scope=col>common_name</th><th scope=col>⋯</th><th scope=col>sample_name</th><th scope=col>Sequencing_method</th><th scope=col>SRA.Study</th><th scope=col>Submitter_Id</th><th scope=col>ENA_first_public</th><th scope=col>ENA_last_update</th><th scope=col>scientific_name</th><th scope=col>ENA_first_public.1</th><th scope=col>ENA_last_update.1</th><th scope=col>INSDC_center_alias</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>ERR1907465</td><td>AMPLICON</td><td>449</td><td> 7339626</td><td>PRJEB20213</td><td>SAMEA103955437</td><td> 3791943</td><td>AUSTRALIAN CENTRE FOR ECOGENOMICS</td><td>2014</td><td>human</td><td>⋯</td><td>F01.d01</td><td>Illumina MiSeq</td><td>ERP022351</td><td>F01.d01</td><td>2019-04-02T17:03:53Z</td><td>2017-04-03T00:34:28Z</td><td>Homo sapiens</td><td></td><td></td><td></td></tr>
	<tr><th scope=row>2</th><td>ERR1907466</td><td>AMPLICON</td><td>412</td><td> 4554387</td><td>PRJEB20213</td><td>SAMEA103955438</td><td> 2389775</td><td>AUSTRALIAN CENTRE FOR ECOGENOMICS</td><td>2014</td><td>human</td><td>⋯</td><td>F01.d02</td><td>Illumina MiSeq</td><td>ERP022351</td><td>F01.d02</td><td>2019-04-02T17:03:53Z</td><td>2017-04-03T00:34:28Z</td><td>Homo sapiens</td><td></td><td></td><td></td></tr>
	<tr><th scope=row>3</th><td>ERR1907467</td><td>AMPLICON</td><td>435</td><td>10169478</td><td>PRJEB20213</td><td>SAMEA103955439</td><td> 5250508</td><td>AUSTRALIAN CENTRE FOR ECOGENOMICS</td><td>2014</td><td>human</td><td>⋯</td><td>F01.d03</td><td>Illumina MiSeq</td><td>ERP022351</td><td>F01.d03</td><td>2019-04-02T17:03:53Z</td><td>2017-04-03T00:34:28Z</td><td>Homo sapiens</td><td></td><td></td><td></td></tr>
	<tr><th scope=row>4</th><td>ERR1907468</td><td>AMPLICON</td><td>457</td><td>24337195</td><td>PRJEB20213</td><td>SAMEA103955440</td><td>12063328</td><td>AUSTRALIAN CENTRE FOR ECOGENOMICS</td><td>2014</td><td>human</td><td>⋯</td><td>F01.d05</td><td>Illumina MiSeq</td><td>ERP022351</td><td>F01.d05</td><td>2019-04-02T17:03:53Z</td><td>2017-04-03T00:34:28Z</td><td>Homo sapiens</td><td></td><td></td><td></td></tr>
	<tr><th scope=row>5</th><td>ERR1907469</td><td>AMPLICON</td><td>457</td><td>12659926</td><td>PRJEB20213</td><td>SAMEA103955441</td><td> 6270490</td><td>AUSTRALIAN CENTRE FOR ECOGENOMICS</td><td>2014</td><td>human</td><td>⋯</td><td>F01.d07</td><td>Illumina MiSeq</td><td>ERP022351</td><td>F01.d07</td><td>2019-04-02T17:03:53Z</td><td>2017-04-03T00:34:28Z</td><td>Homo sapiens</td><td></td><td></td><td></td></tr>
	<tr><th scope=row>6</th><td>ERR1907470</td><td>AMPLICON</td><td>458</td><td> 6407084</td><td>PRJEB20213</td><td>SAMEA103955442</td><td> 3316637</td><td>AUSTRALIAN CENTRE FOR ECOGENOMICS</td><td>2014</td><td>human</td><td>⋯</td><td>F01.d08</td><td>Illumina MiSeq</td><td>ERP022351</td><td>F01.d08</td><td>2019-04-02T17:03:53Z</td><td>2017-04-03T00:34:28Z</td><td>Homo sapiens</td><td></td><td></td><td></td></tr>
</tbody>
</table>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Run'</li><li>'Assay.Type'</li><li>'AvgSpotLen'</li><li>'Bases'</li><li>'BioProject'</li><li>'BioSample'</li><li>'Bytes'</li><li>'Center.Name'</li><li>'Collection_Date'</li><li>'common_name'</li><li>'Consent'</li><li>'DATASTORE.filetype'</li><li>'DATASTORE.provider'</li><li>'DATASTORE.region'</li><li>'ENA.FIRST.PUBLIC..run.'</li><li>'ENA.LAST.UPDATE..run.'</li><li>'environment_.biome.'</li><li>'environment_.feature.'</li><li>'environment_.material.'</li><li>'Experiment'</li><li>'external_id'</li><li>'geo_loc_name_country'</li><li>'geo_loc_name_country_continent'</li><li>'geographic_location_.country_and.or_sea.'</li><li>'geographic_location_.latitude.'</li><li>'geographic_location_.longitude.'</li><li>'humam_oral_environmental_package'</li><li>'INSDC_center_name'</li><li>'INSDC_first_public'</li><li>'INSDC_last_update'</li><li>'INSDC_status'</li><li>'Instrument'</li><li>'Investigation_type'</li><li>'Library.Name'</li><li>'LibraryLayout'</li><li>'LibrarySelection'</li><li>'LibrarySource'</li><li>'Organism'</li><li>'Platform'</li><li>'project_name'</li><li>'ReleaseDate'</li><li>'create_date'</li><li>'version'</li><li>'Sample.Name'</li><li>'sample_name'</li><li>'Sequencing_method'</li><li>'SRA.Study'</li><li>'Submitter_Id'</li><li>'ENA_first_public'</li><li>'ENA_last_update'</li><li>'scientific_name'</li><li>'ENA_first_public.1'</li><li>'ENA_last_update.1'</li><li>'INSDC_center_alias'</li></ol>




```R

uniform_meta <- data.frame(
    BioProject = initial_meta$BioProject,
    Run = initial_meta$Run,
    sample_name = initial_meta$'sample_name',
    Age = NA,
    BMI = NA,
    Sex = NA,
    Region = "V3-V4",
    Platform = initial_meta$Platform,
    Country = initial_meta$'geographic_location_.country_and.or_sea.',
    Site = 'saliva',
    public_time = initial_meta$create_date,
    PCR_primers = 'F: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG; R: GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC',
    title = 'Oxidative stress in the oral cavity is driven by individual-specific bacterial communities',
    Disease_study = "NO",
    Status = NA,
    matched_samples = 'NO')
print(length(unique(uniform_meta$Run)))
print(length(unique(uniform_meta$sample_name)))
```

    [1] 551
    [1] 551



```R
plot_data <- uniform_meta %>%
  select(-BioProject, -Run, -sample_name, -title)

create_simple_plots <- function(data) {
  plots <- list()
  
  for (col in colnames(data)) {
    if (all(is.na(data[[col]]))) {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "this column is NA", size = 6) +
        labs(title = col) +
        theme_void()
    } else if (is.numeric(data[[col]])) {
      p <- ggplot(data, aes(x = .data[[col]])) +
        geom_histogram(fill = "blue", alpha = 0.6, bins = 20, na.rm = TRUE) +
        labs(title = col) +
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

plots <- create_simple_plots(plot_data)
p <- plot_grid(plots$Age, plots$BMI, plots$Sex, plots$Region, 
               plots$Platform, plots$Country, plots$Site,
               plots$public_time, plots$Disease_study,
               plots$Status, plots$matched_samples, 
               ncol = 3, align = "hv", axis = "lb")
p
ggsave("Metadata_Summary.png", plot = p, width = 8, height = 6, dpi = 150)
```


    
![png](metadata_curation_files/metadata_curation_3_0.png)
    



```R
knitr::include_graphics("Metadata_Summary.png")
```


    [1] "Metadata_Summary.png"
    attr(,"class")
    [1] "knit_image_paths" "knit_asis"       



```R
uniform_meta <- uniform_meta[, c("BioProject", "Run", "sample_name", "Age", "BMI", "Sex", "Region", "Platform", "Country", "Site", "public_time", "PCR_primers", "title", "Disease_study", "Status")]
write.csv(uniform_meta, file = "/home/liuxuan/metadata/PRJEB20213/metadata.csv", quote = F)
```
