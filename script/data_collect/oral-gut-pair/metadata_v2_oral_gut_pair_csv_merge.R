library(dplyr)
library(readr)
library(purrr)
library(stringr)

# 1. Define input and output paths
project_list_input <- "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-v2/global-oral-microbiome/data/oral-gut-pair-v2-qc/project_oral_gut_pair_list_qc.txt"
base_project_dir <- "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-v2/global-oral-microbiome/data/projects/"

ref_cleaned_file <- "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-v2/results/metadata-v2-merge-cleaned.csv"

# Define main output path and output file names
output_dir <- "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-v2/oral-gut-pair/results/"
output_file <- "metadata_oral_gut_pair_v2_merge.csv"
study_output_file <- "metadata_study_oral_gut_pair_v2_merge.csv" 

# Define additional output path for the 4-part pair_id data
output_additional_file <- file.path(output_dir, "metadata_pair_v2_additional_oral_site_merge.csv")
# Define output path for the status_unknown NA values data
output_unknown_file <- file.path(output_dir, "metadata_pair_v2_status_unknown_merge.csv")

# Create directories automatically if they do not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
if (!dir.exists(dirname(output_additional_file))) {
  dir.create(dirname(output_additional_file), recursive = TRUE)
}
if (!dir.exists(dirname(output_unknown_file))) {
  dir.create(dirname(output_unknown_file), recursive = TRUE)
}

# 2. Read the project list and generate target CSV file paths
cat("Start reading the project list...\n")
project_ids <- readLines(project_list_input, warn = FALSE)
project_ids <- project_ids[trimws(project_ids) != ""] # Remove empty lines

# --- Manually append the missing projects to the list ---
project_ids <- unique(c(project_ids, "PRJNA847970", "PRJNA535518"))

# Construct full paths for all files, retaining only those that exist
all_files <- file.path(base_project_dir, project_ids, "metadata_oral_gut_pair_v2.csv")
valid_files <- all_files[file.exists(all_files)]

cat("Total number of files to be merged found: ", length(valid_files), " / ", length(project_ids), "\n")

# 3. Specify columns to keep, arranged in order
target_cols <- c(
  "pair_id", "pair_time_point", "pair_bioproject", "pair_type", "pair_sample_number", "pair_subject_number", 
  "bioproject", "sample_id", "subject_id", "age", "age_info", "bmi", "sex", "platform", 
  "sequencing_region", "pcr_primers", "country", "site", "site_group", "disease_study", "status", 
  "title", "public_time", "matched_samples", "is_time_series", "age_cleaned", "bmi_cleaned", 
  "age_cleaned_group", "age_info_cleaned", "bmi_cleaned_group", "sex_cleaned", "platform_cleaned", 
  "sequencing_region_cleaned", "country_cleaned", "continent_group", "site_cleaned", 
  "site_group_cleaned", "disease_study_cleaned", "status_cleaned", "is_time_series_cleaned"
)

# 4. Loop through and process each file, recording individual statuses
merged_data_list <- list()
additional_data_list <- list() 
unknown_status_data_list <- list() 

target_projects_for_na <- c("PRJNA706727", "PRJEB22984", "PRJNA326866", "PRJEB69530", "PRJNA847970", "PRJNA535518")

for (f in valid_files) {
  df <- read_csv(f, col_types = cols(.default = "c"), show_col_types = FALSE)
  proj_dir_name <- basename(dirname(f))
  
  if(nrow(df) > 0) {
    needs_recalc <- FALSE 
    
    # Filter 1: Extract NA status_cleaned
    if(proj_dir_name %in% target_projects_for_na) {
      is_na_status <- is.na(df$status_cleaned) | trimws(tolower(as.character(df$status_cleaned))) %in% c("na", "", "null")
      na_rows <- df[is_na_status, ]
      count_na <- nrow(na_rows)
      
      if(count_na > 0) {
         na_rows$pair_bioproject <- proj_dir_name
         
         available_cols_na <- target_cols[target_cols %in% names(na_rows)]
         df_na_filtered <- na_rows %>% select(all_of(available_cols_na))
         for(c in setdiff(target_cols, names(df_na_filtered))) {
            df_na_filtered[[c]] <- NA_character_
         }
         df_na_filtered <- df_na_filtered %>% select(all_of(target_cols))
         unknown_status_data_list[[f]] <- df_na_filtered
      }
      
      df <- df[!is_na_status, ]
      needs_recalc <- TRUE
    }
    
    # Filter 2: Filter out 4-part pair_id
    oral_sites_pattern <- "(saliva|tooth|oral[ _]mucosa|oral[ _]swab|tongue|gingiva|other)$"
    is_4part <- grepl(oral_sites_pattern, trimws(df$pair_id), ignore.case = TRUE)
    count_4part <- sum(is_4part, na.rm = TRUE)
    
    if (count_4part > 0) {
       df_4part <- df[is_4part, ]
       df <- df[!is_4part, ]
       needs_recalc <- TRUE
       
       available_cols_4part <- target_cols[target_cols %in% names(df_4part)]
       df_4part_filtered <- df_4part %>% select(all_of(available_cols_4part))
       for(c in setdiff(target_cols, names(df_4part_filtered))) {
          df_4part_filtered[[c]] <- NA_character_
       }
       df_4part_filtered <- df_4part_filtered %>% select(all_of(target_cols))
       additional_data_list[[f]] <- df_4part_filtered
    }
    
    # Recalculate metrics
    if(needs_recalc && nrow(df) > 0) {
      valid_pairs <- df %>% group_by(pair_id) %>% summarise(count = n(), .groups='drop') %>% filter(count == 2) %>% pull(pair_id)
      df <- df %>% filter(pair_id %in% valid_pairs)
      
      if(nrow(df) > 0) {
        new_pair_sample_number <- length(unique(df$pair_id))
        new_pair_subject_number <- length(unique(df$subject_id))
        
        subj_pairs_stats <- df %>% group_by(subject_id) %>% summarise(p_count = n_distinct(pair_id), .groups='drop')
        new_pair_type <- ifelse(max(subj_pairs_stats$p_count, na.rm = TRUE) > 1, "many pairs", "single pair")
        
        df <- df %>% mutate(
          pair_sample_number = as.character(new_pair_sample_number),
          pair_subject_number = as.character(new_pair_subject_number),
          pair_type = new_pair_type
        )
      }
    }
    
    # Select and merge main data columns
    if(nrow(df) > 0) {
      available_cols <- target_cols[target_cols %in% names(df)]
      df_filtered <- df %>% select(all_of(available_cols))
      
      for(c in setdiff(target_cols, names(df_filtered))) {
         df_filtered[[c]] <- NA_character_
      }
      df_filtered <- df_filtered %>% select(all_of(target_cols))
      merged_data_list[[f]] <- df_filtered
    }
  }
}

# Bind all filtered lists into complete DataFrames
merged_data <- bind_rows(merged_data_list)
additional_data <- bind_rows(additional_data_list)
unknown_status_data <- bind_rows(unknown_status_data_list)
cat("\n=> Initial filtering, separation, and merging of all files completed!\n\n")

# 5. Calculate time series and final formatting (before mapping)
final_data <- merged_data %>% 
  group_by(pair_bioproject) %>%
  mutate(
    unique_tp_count = n_distinct(pair_time_point, na.rm = TRUE),
    is_time_series_cleaned = ifelse(unique_tp_count >= 3, "yes", "no")
  ) %>%
  ungroup() %>%
  mutate(
    is_time_series_cleaned = ifelse(pair_bioproject == "qiita_11404", "no", is_time_series_cleaned)
  ) %>%
  select(-unique_tp_count) %>%
  select(all_of(target_cols)) 

final_additional_data <- additional_data
final_unknown_data <- unknown_status_data


cat("Loading reference cleaned data for specific column updates...\n")
ref_df <- read_csv(ref_cleaned_file, col_types = cols(.default = "c"), show_col_types = FALSE)

sync_cols <- c(
  "age_cleaned", "bmi_cleaned", "age_cleaned_group", "age_info_cleaned", 
  "bmi_cleaned_group", "sex_cleaned", "platform_cleaned", 
  "sequencing_region_cleaned", "country_cleaned", "continent_group", 
  "site_cleaned", "site_group_cleaned", "disease_study_cleaned", "status_cleaned"
)

ref_subset <- ref_df %>%
  select(sample_id, all_of(sync_cols)) %>%
  distinct(sample_id, .keep_all = TRUE)

update_with_ref <- function(target_df, ref_data, update_cols, final_col_order) {
  if(nrow(target_df) == 0) return(target_df)
  
  df_dropped <- target_df %>% select(-any_of(update_cols))
  df_merged <- df_dropped %>% left_join(ref_data, by = "sample_id")
  df_merged <- df_merged %>%
    mutate(
      age = as.numeric(age),
      bmi = as.numeric(bmi),
      age_cleaned = as.numeric(age_cleaned),
      bmi_cleaned = as.numeric(bmi_cleaned),
      pair_sample_number = as.numeric(pair_sample_number),
      pair_subject_number = as.numeric(pair_subject_number)
    )
    
  df_merged %>% select(all_of(final_col_order))
}

cat("Mapping columns for Final Data...\n")
final_data <- update_with_ref(final_data, ref_subset, sync_cols, target_cols)

if(nrow(final_additional_data) > 0) {
  cat("Mapping columns for Additional Data...\n")
  final_additional_data <- update_with_ref(final_additional_data, ref_subset, sync_cols, target_cols) %>%
    mutate(pair_sample_number = NA_real_, pair_subject_number = NA_real_) 
}

if(nrow(final_unknown_data) > 0) {
  cat("Mapping columns for Unknown Data...\n")
  final_unknown_data <- update_with_ref(final_unknown_data, ref_subset, sync_cols, target_cols)
}
# =====================================================================

# 7. Output files
write_csv(final_data, file.path(output_dir, output_file))
cat(sprintf("=> Successfully output main QC table: %s (Total %d rows)\n", output_file, nrow(final_data)))

if(nrow(final_additional_data) > 0) {
  write_csv(final_additional_data, output_additional_file)
  cat(sprintf("=> Successfully output additional QC table: %s (Total %d rows)\n", basename(output_additional_file), nrow(final_additional_data)))
}

if(nrow(final_unknown_data) > 0) {
  write_csv(final_unknown_data, output_unknown_file)
  cat(sprintf("=> Successfully output unknown status QC table: %s (Total %d rows)\n", basename(output_unknown_file), nrow(final_unknown_data)))
}

# ---------------------------------------------------------
# Extract specified 6 columns, deduplicate by pair_bioproject to generate Study-level CSV 
# ---------------------------------------------------------
study_cols <- c("pair_id", "pair_time_point", "pair_bioproject", "pair_type", 
                "pair_sample_number", "pair_subject_number")

study_data <- final_data %>%
  select(all_of(study_cols)) %>%
  distinct(pair_bioproject, .keep_all = TRUE)

write_csv(study_data, file.path(output_dir, study_output_file))

cat("\nMerge completed. All data and files have been processed, and results are exported to the target directory!\n")