# ==============================================================================
# 02_merge_maaslin3_site_results.R
# Merge results from all MaAsLin3 site pairwise array tasks
#
# This script collects and combines results from all pairwise site comparisons
# after the array job completes.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ==========================================================================
# Configuration
# ==========================================================================
base_dir     <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/03-site-differential-analysis"
output_dir   <- file.path(base_dir, "output")
results_dir  <- file.path(output_dir, "maaslin3_site_pairwise")
alpha <- 0.05
group_col <- "site_cleaned"

cat("==========================================\n")
cat("Merging MaAsLin3 Site Pairwise Results\n")
cat("Start: ", format(Sys.time()), "\n")
cat("==========================================\n\n")

# ==========================================================================
# 1. Find all completed comparisons
# ==========================================================================
cat("[1/3] Scanning for completed comparisons...\n")

comparison_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
comparison_names <- basename(comparison_dirs)

# Filter to directories with completed analysis
completed <- sapply(comparison_dirs, function(d) {
  file.exists(file.path(d, "analysis.done"))
})

comparison_dirs <- comparison_dirs[completed]
comparison_names <- comparison_names[completed]

n_completed <- length(comparison_dirs)
cat("    Found ", n_completed, " completed comparisons\n", sep = "")

if (n_completed == 0) {
  stop("No completed comparisons found!")
}

# ==========================================================================
# 2. Merge all results
# ==========================================================================
cat("\n[2/3] Merging results...\n")

all_results_list <- list()

for (i in seq_along(comparison_dirs)) {
  comp_dir <- comparison_dirs[i]
  comp_name <- comparison_names[i]

  cat("    [", i, "/", n_completed, "] ", comp_name, "\n", sep = "")

  res_file <- file.path(comp_dir, "all_results.tsv")
  if (file.exists(res_file)) {
    res_df <- fread(res_file, sep = "\t", header = TRUE, data.table = FALSE)

    # If the array R script forgot to attach Site1/2/Comparison, attach now
    # by parsing the directory name "<S1>_vs_<S2>".
    if (!"Comparison" %in% names(res_df)) res_df$Comparison <- comp_name
    if (!"Site1" %in% names(res_df) || !"Site2" %in% names(res_df)) {
      parts <- strsplit(comp_name, "_vs_", fixed = TRUE)[[1]]
      if (length(parts) == 2) {
        # Replace underscores back to spaces for site names
        res_df$Site1 <- gsub("_", " ", parts[1])
        res_df$Site2 <- gsub("_", " ", parts[2])
      }
    }
    all_results_list[[comp_name]] <- res_df
  }
}

# Combine all results
if (length(all_results_list) == 0) stop("No all_results.tsv files found!")

full_results <- bind_rows(all_results_list)
rownames(full_results) <- NULL

all_results_file <- file.path(results_dir, "all_pairwise_results.tsv")
fwrite(full_results, all_results_file, sep = "\t", quote = FALSE)
cat("\n    Saved: ", all_results_file, "\n", sep = "")

# ==========================================================================
# 2b. Recompute summary counts directly from full_results
# ==========================================================================
qcol <- if ("qval_individual" %in% names(full_results)) "qval_individual" else "qval"
cat("\n    Using q-value column: ", qcol, "\n", sep = "")

is_group_row <- if ("metadata" %in% names(full_results)) {
  full_results[["metadata"]] == group_col
} else {
  rep(TRUE, nrow(full_results))
}

sig_all <- full_results[is_group_row &
                          !is.na(full_results[[qcol]]) &
                          full_results[[qcol]] < alpha, ]

# Save the corrected significant-only table for downstream plots
sig_file <- file.path(results_dir, "significant_pairwise_results.tsv")
fwrite(sig_all, sig_file, sep = "\t", quote = FALSE)
cat("    Saved: ", sig_file, " (", nrow(sig_all), " significant rows)\n", sep = "")

summary_counts <- sig_all %>%
  group_by(Site1, Site2, Comparison) %>%
  summarise(N_Significant = n_distinct(feature), .groups = "drop") %>%
  as.data.frame()

# Make sure every comparison appears, even if it has 0 significant taxa
all_comparisons <- unique(data.frame(
  Comparison = sapply(all_results_list, function(d) d$Comparison[1]),
  Site1 = sapply(all_results_list, function(d) d$Site1[1]),
  Site2 = sapply(all_results_list, function(d) d$Site2[1]),
  stringsAsFactors = FALSE
))
summary_counts <- merge(all_comparisons, summary_counts,
                        by = c("Site1", "Site2", "Comparison"),
                        all.x = TRUE)
summary_counts$N_Significant[is.na(summary_counts$N_Significant)] <- 0

# Symmetric form (needed by bubble matrix)
summary_symmetric <- summary_counts
for (i in seq_len(nrow(summary_counts))) {
  summary_symmetric <- rbind(summary_symmetric, data.frame(
    Site1 = summary_counts$Site2[i],
    Site2 = summary_counts$Site1[i],
    Comparison = paste0(gsub(" ", "_", summary_counts$Site2[i]), "_vs_", gsub(" ", "_", summary_counts$Site1[i])),
    N_Significant = summary_counts$N_Significant[i],
    stringsAsFactors = FALSE
  ))
}

summary_file_out <- file.path(results_dir, "summary_diff_counts.tsv")
fwrite(summary_symmetric, summary_file_out, sep = "\t", quote = FALSE)
cat("    Saved: ", summary_file_out, "\n", sep = "")

# ==========================================================================
# 3. Print summary
# ==========================================================================
cat("\n[3/3] Summary statistics...\n")

cat("\n==========================================\n")
cat("Pairwise Significant Taxa Counts (", qcol, " < ", alpha, ", metadata=", group_col, "):\n", sep = "")
cat("==========================================\n")
print(summary_counts)

cat("\n    Total comparisons: ", nrow(summary_counts), "\n", sep = "")
cat("    Total associations tested: ", nrow(full_results), "\n", sep = "")
cat("    Mean diff taxa per comparison: ", round(mean(summary_counts$N_Significant), 1), "\n", sep = "")
cat("    Range: ", min(summary_counts$N_Significant), " - ", max(summary_counts$N_Significant), "\n", sep = "")

cat("\n==========================================\n")
cat("Merge Complete\n")
cat("End: ", format(Sys.time()), "\n")
cat("==========================================\n")