# ==============================================================================
# 01_analysis_oral_site_maaslin3_array.R
# MaAsLin3 pairwise oral site comparisons - Array Job Version
#
# This script is designed for Slurm array jobs. Each task processes ONE
# oral site pair comparison.
#
# Usage: Rscript script.R <task_id>
#   task_id: 1-15 (one for each pairwise comparison, C(6,2) = 15)
# ==============================================================================

suppressPackageStartupMessages({
  library(qiime2R)
  library(maaslin3)
  library(dplyr)
  library(data.table)
  library(Matrix)
})

# ==========================================================================
# 0. Get task ID from command line
# ==========================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript script.R <task_id>")
}
task_id <- as.integer(args[1])

set.seed(20260615)

# ==========================================================================
# 1. Configuration
# ==========================================================================
base_dir      <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/03-site-differential-analysis"
input_dir     <- file.path(base_dir, "input")
output_dir    <- file.path(base_dir, "output")
results_dir   <- file.path(output_dir, "maaslin3_site_pairwise")

genus_qza     <- file.path(input_dir, "genus_table_all_qc.qza")
metadata_file <- file.path(input_dir, "metadata_oral_qc_update.tsv")

# --- ANALYSIS PARAMETERS ---
# Main effect: oral site (site_cleaned)
group_col     <- "site_cleaned"
# Covariate: continent (now a confounder, not main effect)
continent_col <- "continent_group"
# Filter: all oral samples
site_group_col <- "site_group_cleaned"

sample_id_col   <- "sample_id"
age_col         <- "age_cleaned"
sex_col         <- "sex_cleaned"
bmi_col         <- "bmi_cleaned"
project_col     <- "bioproject_id"
platform_col    <- "platform_cleaned"
seq_region_col  <- "sequencing_region_cleaned"
pcr_col         <- "pcr_primers"

# All anatomical oral sites (excluding "other" and empty)
site_order <- c("saliva", "tooth", "oral mucosa", "oral swab",
                "tongue", "gingiva")

min_samples_per_site <- 10

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
alpha <- 0.05

# Helper: sanitize site names for directory use (replace spaces)
sanitize <- function(x) gsub(" ", "_", x)

# ==========================================================================
# 2. Load and prepare data
# ==========================================================================
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cat("==========================================\n")
cat("MaAsLin3 Site Pairwise Array Task: ", task_id, "\n")
cat("Start: ", format(Sys.time()), "\n")
cat("==========================================\n\n")

cat("[1/4] Loading data...\n")

# Load genus abundance table
qza_obj <- read_qza(genus_qza)
genus_sparse <- qza_obj$data

# Load metadata
metadata_raw <- fread(metadata_file, sep = "\t", header = TRUE,
                      quote = "", fill = TRUE, check.names = FALSE,
                      data.table = FALSE)

# ==========================================================================
# 3. Filter and clean metadata
# ==========================================================================
cat("\n[2/4] Filtering metadata...\n")

required_cols <- c(sample_id_col, group_col, site_group_col, continent_col, project_col)
missing_cols <- setdiff(required_cols, names(metadata_raw))
if (length(missing_cols) > 0) {
  stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
}

metadata_clean <- metadata_raw %>%
  filter(!is.na(.data[[sample_id_col]]) & .data[[sample_id_col]] != "") %>%
  distinct(.data[[sample_id_col]], .keep_all = TRUE) %>%
  mutate(
    temp_age = suppressWarnings(as.numeric(as.character(.data[[age_col]]))),
    temp_bmi = suppressWarnings(as.numeric(as.character(.data[[bmi_col]]))),
    !!group_col := trimws(as.character(.data[[group_col]])),
    !!sex_col := trimws(as.character(.data[[sex_col]])),
    !!continent_col := trimws(as.character(.data[[continent_col]])),
    !!platform_col := trimws(as.character(.data[[platform_col]])),
    !!seq_region_col := trimws(as.character(.data[[seq_region_col]])),
    !!pcr_col := trimws(as.character(.data[[pcr_col]]))
  ) %>%
  mutate(
    !!bmi_col := case_when(
      is.na(temp_bmi) ~ "Unknown",
      temp_bmi < 18.5 ~ "Underweight",
      temp_bmi >= 18.5 & temp_bmi < 24.0 ~ "Normal",
      temp_bmi >= 24.0 & temp_bmi < 28.0 ~ "Overweight",
      temp_bmi >= 28.0 ~ "Obese",
      TRUE ~ "Unknown"
    ),
    !!age_col := case_when(
      is.na(temp_age) ~ "Unknown",
      temp_age < 18 ~ "Minor",
      temp_age >= 18 & temp_age < 45 ~ "Young_Adult",
      temp_age >= 45 & temp_age < 65 ~ "Middle_Aged",
      temp_age >= 65 ~ "Elderly",
      TRUE ~ "Unknown"
    )
  ) %>%
  mutate(
    !!bmi_col := factor(.data[[bmi_col]], levels = c("Normal", "Underweight", "Overweight", "Obese", "Unknown")),
    !!age_col := factor(.data[[age_col]], levels = c("Young_Adult", "Minor", "Middle_Aged", "Elderly", "Unknown"))
  ) %>%
  mutate(
    !!sex_col := ifelse(.data[[sex_col]] %in% c("", "NA"), "Unknown", .data[[sex_col]]),
    !!platform_col := ifelse(.data[[platform_col]] %in% c("", "NA"), "Unknown", .data[[platform_col]]),
    !!seq_region_col := ifelse(.data[[seq_region_col]] %in% c("", "NA"), "Unknown", .data[[seq_region_col]]),
    !!pcr_col := ifelse(.data[[pcr_col]] %in% c("", "NA"), "Unknown", .data[[pcr_col]])
  ) %>%
  filter(
    .data[[site_group_col]] == "Oral" &
    .data[[group_col]] != "" & !is.na(.data[[group_col]]) &
    .data[[group_col]] %in% site_order
  ) %>%
  select(-temp_age, -temp_bmi)

metadata_clean <- as.data.frame(metadata_clean)
rownames(metadata_clean) <- metadata_clean[[sample_id_col]]

# Find common samples
common_samples <- intersect(rownames(metadata_clean), colnames(genus_sparse))
if (length(common_samples) == 0) {
  stop("No common samples between metadata and genus table!")
}

metadata_clean <- metadata_clean[common_samples, ]
genus_sparse <- genus_sparse[, common_samples]

# Identify valid sites
site_counts <- table(metadata_clean[[group_col]])
valid_sites <- names(site_counts)[site_counts >= min_samples_per_site]
valid_sites <- intersect(site_order, valid_sites)

if (length(valid_sites) < 2) {
  stop("Need at least 2 sites with >= ", min_samples_per_site, " samples!")
}

cat("    Valid sites: ", paste(valid_sites, collapse = ", "), "\n", sep = "")
cat("    Sample counts per site:\n")
print(sort(site_counts[valid_sites], decreasing = TRUE))

# ==========================================================================
# 4. Generate all pairs and select the one for this task
# ==========================================================================
cat("\n[3/4] Identifying comparison for task ", task_id, "...\n", sep = "")

site_pairs <- combn(valid_sites, 2, simplify = FALSE)
n_pairs_total <- length(site_pairs)

if (task_id > n_pairs_total) {
  cat("Task ID", task_id, "exceeds total pairs (", n_pairs_total, ") - nothing to do. Exiting.\n", sep = " ")
  quit(save = "no", status = 0)
}

pair <- site_pairs[[task_id]]
site1 <- pair[1]
site2 <- pair[2]

comparison_name <- paste0(sanitize(site1), "_vs_", sanitize(site2))
output_subdir <- file.path(results_dir, comparison_name)

cat("    Task ", task_id, "/", n_pairs_total, ": ", comparison_name, "\n", sep = "")

# Check if already completed (for resumability)
done_file <- file.path(output_subdir, "analysis.done")
if (file.exists(done_file)) {
  cat("    Already completed. Skipping.\n")
  quit(save = "no", status = 0)
}

# ==========================================================================
# 5. Run MaAsLin3 for this pair
# ==========================================================================
cat("\n[4/4] Running MaAsLin3 analysis...\n")

# Subset samples for this pair
pair_samples <- metadata_clean[[group_col]] %in% c(site1, site2)
pair_meta <- metadata_clean[pair_samples, ]
pair_genus <- genus_sparse[, pair_meta[[sample_id_col]]]

# Free unused large objects before dense conversion to reduce memory pressure
rm(qza_obj, genus_sparse, metadata_raw, metadata_clean)
gc()

# Set site as factor with reference level
pair_meta[[group_col]] <- factor(pair_meta[[group_col]],
                                   levels = c(site1, site2))

# Prepare input for MaAsLin3
input_data <- as.data.frame(as.matrix(t(pair_genus)))
input_meta <- as.data.frame(pair_meta)
rownames(input_meta) <- input_meta[[sample_id_col]]
input_meta <- input_meta[rownames(input_data), ]

# --- Dynamic Covariate Filtering to Prevent Singular Fit ---
all_covars <- c(age_col, sex_col, platform_col, seq_region_col, bmi_col, continent_col)
valid_covars <- character(0)

for (cov in all_covars) {
  vals <- as.character(input_meta[[cov]])
  tab <- sort(table(vals), decreasing = TRUE)

  if (length(tab) < 2) {
    cat("    WARNING: ", cov, " has only 1 level, dropping from formula.\n", sep = "")
    next
  }

  cross_tab <- table(input_meta[[group_col]], input_meta[[cov]])
  if (any(colSums(cross_tab > 0) < 2)) {
      cat("    WARNING: ", cov, " has non-overlapping levels between sites, dropping to prevent Singular Fit.\n", sep = "")
      next
  }

  ref_level <- names(tab)[1]
  if (cov == bmi_col && "Normal" %in% names(tab)) ref_level <- "Normal"
  if (cov == age_col && "Young_Adult" %in% names(tab)) ref_level <- "Young_Adult"

  other_levels <- setdiff(names(tab), ref_level)
  input_meta[[cov]] <- factor(vals, levels = c(ref_level, other_levels))
  valid_covars <- c(valid_covars, cov)
  cat("    [OK] Kept ", cov, " (Ref: ", ref_level, ")\n", sep = "")
}

if (length(valid_covars) > 0) {
  base_formula <- paste0("~ ", group_col, " + ", paste(valid_covars, collapse = " + "))
} else {
  base_formula <- paste0("~ ", group_col)
}

project_tab <- table(input_meta[[group_col]], input_meta[[project_col]])
shared_projects <- sum(colSums(project_tab > 0) == 2)

if (shared_projects > 0) {
  formula_str <- paste0(base_formula, " + (1|", project_col, ")")
  cat("    Using Random Effect for BioProject (", shared_projects, " shared projects found).\n", sep="")
} else {
  formula_str <- base_formula
  cat("    WARNING: No shared BioProjects between sites. Falling back to fixed-effects only (LM).\n", sep="")
}

cat("    Formula: ", formula_str, "\n", sep = "")
cat("    Samples: ", site1, " (", sum(pair_meta[[group_col]] == site1), "), ",
    site2, " (", sum(pair_meta[[group_col]] == site2), ")\n", sep = "")

# Run MaAsLin3
if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)

tryCatch({
  fit_out <- maaslin3(
    input_data     = input_data,
    input_metadata = input_meta,
    output         = output_subdir,
    formula        = formula_str,
    reference      = paste0(group_col, ",", site1),
    normalization  = "TSS",
    transform      = "LOG",
    standardize    = TRUE,
    cores          = n_cores,
    min_prevalence = 0.05,
    plot_associations = TRUE,
    plot_summary      = TRUE
  )

  # Read and annotate results
  res_file <- file.path(output_subdir, "all_results.tsv")
  if (file.exists(res_file)) {
    res_df <- fread(res_file, sep = "\t", header = TRUE, data.table = FALSE)

    # Add comparison metadata
    res_df$Site1 <- site1
    res_df$Site2 <- site2
    res_df$Comparison <- comparison_name

    # Determine q-value column (MaAsLin3 outputs qval_individual / qval_joint)
    qcol <- if ("qval_individual" %in% names(res_df)) "qval_individual" else "qval"

    # Filter significant results: site_cleaned main effect only
    sig_df <- res_df[!is.na(res_df[[qcol]]) & res_df[[qcol]] < alpha, ]
    if ("metadata" %in% names(sig_df)) {
      sig_df <- sig_df[sig_df$metadata == group_col, ]
    }

    if (nrow(sig_df) > 0) {
      sig_file <- file.path(output_subdir, "significant_results.tsv")
      fwrite(sig_df, sig_file, sep = "\t", quote = FALSE)
    }

    n_sig <- nrow(sig_df)
    cat("    Significant taxa (", qcol, " < ", alpha, ", metadata=", group_col, "): ", n_sig, "\n", sep = "")

    # Save summary
    summary_file <- file.path(output_subdir, "summary.tsv")
    summary_df <- data.frame(
      Site1 = site1,
      Site2 = site2,
      Comparison = comparison_name,
      N_Significant = n_sig,
      stringsAsFactors = FALSE
    )
    fwrite(summary_df, summary_file, sep = "\t", quote = FALSE)
  }

  # Mark as completed
  writeLines("done", done_file)
  cat("\n    SUCCESS!\n")

}, error = function(e) {
  cat("    ERROR: ", e$message, "\n", sep = "")
  quit(save = "no", status = 1)
})

cat("\n==========================================\n")
cat("Task ", task_id, " Complete\n", sep = "")
cat("End: ", format(Sys.time()), "\n")
cat("==========================================\n")