#!/usr/bin/env Rscript
# ==============================================================================
# 00_load_and_validate.R
# Load all HOMD uncapture analysis data and validate consistency
# ==============================================================================

library(jsonlite)
library(dplyr)
library(tidyr)

BASE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"

cat("=== Loading HOMD Uncapture Analysis Data ===\n\n")

# ------------------------------------------------------------------------------
# 1. collect_stats.json — overall ASV counts
# ------------------------------------------------------------------------------
cat("1. Loading collect_stats.json ...\n")
collect_stats <- fromJSON(file.path(BASE, "collect_stats.json"))
cat(sprintf("   Total ASVs: %s\n", format(collect_stats$total_asvs, big.mark = ",")))
cat(sprintf("   Derep sequences: %s\n", format(collect_stats$exact_derep_sequences, big.mark = ",")))
cat(sprintf("   BioProjects: %d\n", collect_stats$metadata_bioprojects))

# Build per-region summary
region_df <- bind_rows(lapply(names(collect_stats$per_region), function(r) {
  tibble(region = r,
         asvs = collect_stats$per_region[[r]]$asvs,
         derep = collect_stats$per_region[[r]]$derep_sequences)
}))
region_df <- region_df %>% arrange(desc(asvs))
cat(sprintf("   Top region: %s (%s ASVs)\n", region_df$region[1],
            format(region_df$asvs[1], big.mark = ",")))

# ------------------------------------------------------------------------------
# 2. homd_mapping_stats.json — HOMD capture vs uncaptured
# ------------------------------------------------------------------------------
cat("\n2. Loading homd_mapping_stats.json ...\n")
homd_stats <- fromJSON(file.path(BASE, "homd_mapping", "homd_mapping_stats.json"))
cat(sprintf("   Total derep: %s\n", format(homd_stats$total_derep_sequences, big.mark = ",")))
cat(sprintf("   HOMD captured: %s (%.1f%%)\n",
            format(homd_stats$homd_captured, big.mark = ","),
            100 * homd_stats$homd_captured / homd_stats$total_derep_sequences))
cat(sprintf("   HOMD uncaptured: %s (%.1f%%)\n",
            format(homd_stats$homd_uncaptured, big.mark = ","),
            100 * homd_stats$homd_uncaptured / homd_stats$total_derep_sequences))
cat("   Uncaptured breakdown:\n")
for (reason in names(homd_stats$by_reason)) {
  if (reason != "captured") {
    cat(sprintf("     - %s: %s\n", reason,
                format(homd_stats$by_reason[[reason]], big.mark = ",")))
  }
}

# ------------------------------------------------------------------------------
# 3. homd_mapping_per_region_stats.json — capture by region
# ------------------------------------------------------------------------------
cat("\n3. Loading homd_mapping_per_region_stats.json ...\n")
homd_region <- fromJSON(file.path(BASE, "homd_mapping", "homd_mapping_per_region_stats.json"))
homd_region_df <- bind_rows(lapply(names(homd_region), function(r) {
  tibble(region = r,
         total = homd_region[[r]]$total,
         captured = homd_region[[r]]$captured,
         uncaptured = homd_region[[r]]$uncaptured)
})) %>%
  mutate(capture_rate = captured / total * 100) %>%
  arrange(desc(total))

cat("   Top regions by ASV count:\n")
for (i in 1:min(8, nrow(homd_region_df))) {
  cat(sprintf("     %s: %s total, %.1f%% captured\n",
              homd_region_df$region[i],
              format(homd_region_df$total[i], big.mark = ","),
              homd_region_df$capture_rate[i]))
}

# ------------------------------------------------------------------------------
# 4. chimera_stats.json — chimera removal
# ------------------------------------------------------------------------------
cat("\n4. Loading chimera_stats.json ...\n")
chimera_stats <- fromJSON(file.path(BASE, "chimera_removal", "chimera_stats.json"))
cat(sprintf("   Total derep: %s\n", format(chimera_stats$total_derep_sequences, big.mark = ",")))
cat(sprintf("   Chimeras: %s (%.1f%%)\n",
            format(chimera_stats$chimera_flagged, big.mark = ","),
            chimera_stats$chimera_pct))
cat(sprintf("   HOMD uncaptured: %s\n", format(chimera_stats$homd_uncaptured_total, big.mark = ",")))
cat(sprintf("   Candidates post-chimera: %s (%.1f%%)\n",
            format(chimera_stats$candidates_post_chimera, big.mark = ","),
            chimera_stats$candidate_pct))

# ------------------------------------------------------------------------------
# 5. Candidate counts at each pipeline stage
# ------------------------------------------------------------------------------
cat("\n5. Pipeline stage counts ...\n")
n_derep <- collect_stats$exact_derep_sequences
n_homd_captured <- homd_stats$homd_captured
n_homd_uncaptured <- homd_stats$homd_uncaptured

# Read actual counts from files
candidates_all <- read.delim(file.path(BASE, "pipeline", "01_candidates", "candidates_all.tsv"),
                              header = TRUE, nrows = 1)
n_support_floor <- nrow(read.delim(file.path(BASE, "pipeline", "01_candidates", "candidates_all.tsv"),
                                    header = TRUE))

n_nonchimeric <- nrow(read.delim(file.path(BASE, "pipeline", "02_chimera", "candidates_all_nonchimeric.tsv"),
                                  header = TRUE))

cat(sprintf("   Dereplicated: %s\n", format(n_derep, big.mark = ",")))
cat(sprintf("   HOMD uncaptured: %s\n", format(n_homd_uncaptured, big.mark = ",")))
cat(sprintf("   After support floor: %s\n", format(n_support_floor, big.mark = ",")))
cat(sprintf("   After chimera removal: %s\n", format(n_nonchimeric, big.mark = ",")))

# ------------------------------------------------------------------------------
# 6. Candidate novel phylotypes (updated version, Jul 23)
# ------------------------------------------------------------------------------
cat("\n6. Loading candidate novel phylotypes (updated) ...\n")
phy_file <- file.path(BASE, "pipeline", "06_phylotypes_sepp", "phylotypes_candidate_novel_updated.tsv")
phy <- read.delim(phy_file, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("   Total candidate novel phylotypes: %d\n", nrow(phy)))
cat(sprintf("   Columns: %d\n", ncol(phy)))

# Check novelty class
if ("novelty_class" %in% names(phy)) {
  novelty_counts <- table(phy$novelty_class)
  cat("   Novelty classification:\n")
  for (cls in names(novelty_counts)) {
    cat(sprintf("     - %s: %d\n", cls, novelty_counts[cls]))
  }
} else {
  cat("   WARNING: novelty_class column not found in updated file\n")
  cat("   Falling back to phylotypes_novelty_classified.tsv ...\n")
  phy_novelty <- read.delim(file.path(BASE, "pipeline", "07_novelty_sepp",
                                       "phylotypes_novelty_classified.tsv"),
                             header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("   Loaded %d phylotypes from novelty file\n", nrow(phy_novelty)))
}

# Subject and Bioproject ranges
cat(sprintf("   Subject range: %d - %d\n", min(phy$max_subject_prevalence), max(phy$max_subject_prevalence)))
cat(sprintf("   Bioproject range: %d - %d\n", min(phy$max_n_bioprojects), max(phy$max_n_bioprojects)))
cat(sprintf("   Reads range: %s - %s\n",
            format(min(phy$sum_total_reads), big.mark = ","),
            format(max(phy$sum_total_reads), big.mark = ",")))

# Contamination info
if ("phylotype_contamination_concern" %in% names(phy)) {
  contam_counts <- table(phy$phylotype_contamination_concern)
  cat("   Contamination concern:\n")
  for (cls in names(contam_counts)) {
    cat(sprintf("     - %s: %d\n", cls, contam_counts[cls]))
  }
}

# ------------------------------------------------------------------------------
# 7. Integrated candidate annotations — contamination summary
# ------------------------------------------------------------------------------
cat("\n7. Loading integrated candidate annotations ...\n")
integ_file <- file.path(BASE, "pipeline", "05_integrated_candidates", "integrated_candidate_annotations.tsv")
# Read just the header and first row to check columns
integ_header <- read.delim(integ_file, header = TRUE, nrows = 0)
cat(sprintf("   Columns: %d\n", ncol(integ_header)))
cat(sprintf("   Column names: %s\n", paste(names(integ_header), collapse = ", ")))

# Read contamination concern level distribution
if ("contamination_concern_level" %in% names(integ_header)) {
  cat("   Reading contamination levels (this may take a moment) ...\n")
  col_types <- rep("NULL", ncol(integ_header))
  col_types[41] <- "character"  # contamination_concern_level is column 41
  integ_levels <- read.delim(integ_file, header = TRUE, stringsAsFactors = FALSE,
                             colClasses = col_types)
  level_counts <- table(integ_levels[[1]])
  cat("   Contamination concern level distribution:\n")
  for (cls in names(level_counts)) {
    cat(sprintf("     - %s: %s (%.1f%%)\n", cls,
                format(level_counts[cls], big.mark = ","),
                100 * level_counts[cls] / sum(level_counts)))
  }
}

# ------------------------------------------------------------------------------
# 8. Save validated data as RDS for downstream scripts
# ------------------------------------------------------------------------------
cat("\n8. Saving validated data as RDS ...\n")

validated_data <- list(
  # Counts
  n_derep = n_derep,
  n_homd_captured = n_homd_captured,
  n_homd_uncaptured = n_homd_uncaptured,
  n_support_floor = n_support_floor,
  n_nonchimeric = n_nonchimeric,
  n_projects = collect_stats$metadata_bioprojects,

  # Per-region HOMD mapping
  homd_region_df = homd_region_df,

  # Chimera stats
  chimera_pct = chimera_stats$chimera_pct,

  # Uncaptured reasons
  uncaptured_reasons = tibble(
    reason = names(homd_stats$by_reason),
    count = as.numeric(homd_stats$by_reason)
  ) %>% filter(reason != "captured"),

  # Phylotypes
  phylotypes = phy,

  # Overall stats
  collect_stats = collect_stats,
  homd_stats = homd_stats
)

saveRDS(validated_data, file.path(BASE, "pipeline", "analysis_report", "R", "validated_data.rds"))
cat("   Saved: pipeline/analysis_report/R/validated_data.rds\n")

cat("\n=== Data Validation Complete ===\n")
cat(sprintf("Final: %d candidate novel phylotypes from %d BioProjects\n",
            nrow(phy), collect_stats$metadata_bioprojects))
