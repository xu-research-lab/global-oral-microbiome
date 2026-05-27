library(dplyr)
library(readr)
library(purrr)
library(stringr)

# 1. Define input and output paths
project_list_input <- "D:/YL2025/2026/global-oral-github-YL-v2/global-oral-microbiome/data/project_oral_gut_pair_list.txt"
base_project_dir <- "D:/YL2025/2026/global-oral-github-YL-v2/global-oral-microbiome/data/projects"

# Define main output path and output file names
output_dir <- "E:/oral/26/metadata-2604/Oral-Metadata-merge-YL/pair/260518/pair-v2/results-0520/merge-pair-0527"
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
# First is pair_id, ending with status_cleaned and then the new is_time_series_cleaned
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
additional_data_list <- list() # Store 4-part pair_id data
unknown_status_data_list <- list() # Store NA status_cleaned data for specific projects

# Target projects for NA extraction and tracking vector
target_projects_for_na <- c("PRJEB11425", "PRJEB69530", "PRJNA847970", "PRJNA535518")

for (f in valid_files) {
  df <- read_csv(f, col_types = cols(.default = "c"), show_col_types = FALSE)
  # Output folder (project) name to simplify the log and guarantee strict matching
  proj_dir_name <- basename(dirname(f))
  
  if(nrow(df) > 0) {
    original_nrow <- nrow(df)
    needs_recalc <- FALSE # Flag to mark if recalculation of metrics is needed
    
    # ---------------------------------------------------------
    # Filter 1: Extract NA status_cleaned for target 4 projects based on FOLDER name
    # ---------------------------------------------------------
    if(proj_dir_name %in% target_projects_for_na) {
      # Identify NA values, explicit "NA" strings, or empty strings
      is_na_status <- is.na(df$status_cleaned) | trimws(tolower(as.character(df$status_cleaned))) %in% c("na", "", "null")
      na_rows <- df[is_na_status, ]
      count_na <- nrow(na_rows)
      
      if(count_na > 0) {
         # --- Force assign the folder name so they aren't blank in the CSV ---
         na_rows$pair_bioproject <- proj_dir_name
         
         # Standardize column structure for extracted NA rows
         available_cols_na <- target_cols[target_cols %in% names(na_rows)]
         df_na_filtered <- na_rows %>% select(all_of(available_cols_na))
         for(c in setdiff(target_cols, names(df_na_filtered))) {
            df_na_filtered[[c]] <- NA_character_
         }
         df_na_filtered <- df_na_filtered %>% select(all_of(target_cols))
         
         unknown_status_data_list[[f]] <- df_na_filtered
      }
      
      # Remove NA rows from main dataframe
      df <- df[!is_na_status, ]
      needs_recalc <- TRUE
    }
    
    # ---------------------------------------------------------
    # Filter 2: Filter out 4-part pair_id (ending with 7 specific keywords) to a separate table
    # ---------------------------------------------------------
    oral_sites_pattern <- "(saliva|tooth|oral[ _]mucosa|oral[ _]swab|tongue|gingiva|other)$"
    is_4part <- grepl(oral_sites_pattern, trimws(df$pair_id), ignore.case = TRUE)
    
    count_4part <- sum(is_4part, na.rm = TRUE)
    if (count_4part > 0) {
       # Extract 4-part rows
       df_4part <- df[is_4part, ]
       # Remove these rows from the main dataset
       df <- df[!is_4part, ]
       needs_recalc <- TRUE
       
       # Standardize column structure for 4-part data
       available_cols_4part <- target_cols[target_cols %in% names(df_4part)]
       df_4part_filtered <- df_4part %>% select(all_of(available_cols_4part))
       for(c in setdiff(target_cols, names(df_4part_filtered))) {
          df_4part_filtered[[c]] <- NA_character_
       }
       df_4part_filtered <- df_4part_filtered %>% select(all_of(target_cols))
       
       additional_data_list[[f]] <- df_4part_filtered
    }
    
    # ---------------------------------------------------------
    # Recalculate metrics
    # ---------------------------------------------------------
    if(needs_recalc && nrow(df) > 0) {
      # Keep pair_id that appear exactly twice
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
    
    # ---------------------------------------------------------
    # Select and merge main data columns
    # ---------------------------------------------------------
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

# 5. Restore numerical properties and arrange columns for the main table
final_data <- merged_data %>% 
  select(all_of(target_cols)) %>% 
  mutate(
    age                 = as.numeric(age),
    bmi                 = as.numeric(bmi),
    age_cleaned         = as.numeric(age_cleaned),
    bmi_cleaned         = as.numeric(bmi_cleaned),
    pair_sample_number  = as.numeric(pair_sample_number),
    pair_subject_number = as.numeric(pair_subject_number)
  )

# Calculate is_time_series_cleaned based on unique pair_time_point values
# AND explicitly enforce qiita_11404 to be "no"
final_data <- final_data %>%
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
  select(all_of(target_cols)) # Ensure precise final ordering

# Output the 4-part additional file
if(nrow(additional_data) > 0) {
  final_additional_data <- additional_data %>%
    mutate(
      age                 = as.numeric(age),
      bmi                 = as.numeric(bmi),
      age_cleaned         = as.numeric(age_cleaned),
      bmi_cleaned         = as.numeric(bmi_cleaned),
      # Set pair metrics to NA for additional files
      pair_sample_number  = NA_real_,
      pair_subject_number = NA_real_
    )
  write_csv(final_additional_data, output_additional_file)
  cat(sprintf("=> Successfully output the additional 4-part pair_id table: %s (Total %d rows, pair_sample_number and pair_subject_number set to NA)\n\n", output_additional_file, nrow(final_additional_data)))
}

# Output the unknown status NA file
if(nrow(unknown_status_data) > 0) {
  final_unknown_data <- unknown_status_data %>%
    mutate(
      age                 = as.numeric(age),
      bmi                 = as.numeric(bmi),
      age_cleaned         = as.numeric(age_cleaned),
      bmi_cleaned         = as.numeric(bmi_cleaned),
      pair_sample_number  = as.numeric(pair_sample_number),
      pair_subject_number = as.numeric(pair_subject_number)
    )
  write_csv(final_unknown_data, output_unknown_file)
  cat(sprintf("=> Successfully output the status_cleaned NA extracted table: %s (Total %d rows)\n\n", output_unknown_file, nrow(final_unknown_data)))
}

# Output the merged main CSV file
write_csv(final_data, file.path(output_dir, output_file))

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