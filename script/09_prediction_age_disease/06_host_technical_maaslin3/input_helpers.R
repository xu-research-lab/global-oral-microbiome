load_standard_metadata <- function(dedup = TRUE) {
  base <- "/beegfs/huangxiaochang/global_oral/new_analysis"
  path <- file.path(base, ifelse(dedup, "metadata_oral_qc_dedup.tsv", "metadata_oral_qc.tsv"))
  metadata <- read.table(path, header = TRUE, sep = "\t", fill = TRUE,
                         comment.char = "", stringsAsFactors = FALSE, quote = "")

  rename_if_present <- function(old, new) {
    if (old %in% colnames(metadata)) colnames(metadata)[colnames(metadata) == old] <<- new
  }
  rename_if_present("sample_id", "SampleID")
  rename_if_present("bioproject_id", "BioProject")
  rename_if_present("site_cleaned", "Site_Cleaned")
  rename_if_present("continent_group", "Geography")
  rename_if_present("status_cleaned", "Status_Cleaned")
  rename_if_present("sex_cleaned", "Sex")
  rename_if_present("age_cleaned_group", "Age_Group")
  rename_if_present("sequencing_region_cleaned", "Region_custom")

  for (col in c("SampleID", "BioProject", "Site_Cleaned", "Geography", "Status_Cleaned", "Sex", "Age_Group", "Region_custom")) {
    if (col %in% colnames(metadata)) {
      metadata[[col]] <- trimws(as.character(metadata[[col]]))
      metadata[[col]][metadata[[col]] %in% c("", "nan", "NA", "None")] <- NA
    }
  }

  site_map <- c("saliva" = "Saliva", "tongue" = "Tongue", "tooth" = "Tooth",
                "oral swab" = "Oral_swab", "oral mucosa" = "Oral Mucosa",
                "gingiva" = "Gingiva", "other" = "Others")
  metadata$Site_Cleaned <- ifelse(metadata$Site_Cleaned %in% names(site_map),
                                  site_map[metadata$Site_Cleaned], metadata$Site_Cleaned)
  metadata$Status_Cleaned[metadata$Status_Cleaned %in% c("HC", "Healthy")] <- "Healthy Control"
  metadata$Age_Group <- gsub("＞", ">", metadata$Age_Group)
  metadata
}

load_genus_relabund <- function() {
  path <- "/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3/genus_table_relabund.tsv"
  read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE,
             comment.char = "", stringsAsFactors = FALSE, quote = "")
}
