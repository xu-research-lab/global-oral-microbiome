#!/bin/bash
#SBATCH --job-name=decontam
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/decontam_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/decontam_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=12:00:00

# ============================================================
# Phase 2 Step 4: decontam within each blank-bearing study
# Prevalence mode: identifies ASVs statistically enriched in controls
# ============================================================

set -e
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2
export R_PROFILE_USER="${HOME}/.Rprofile.qiime2"

BASE_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"
DECONTAM_DIR="${BASE_DIR}/decontam_results"
mkdir -p "${DECONTAM_DIR}"

echo "=== decontam within-study analysis ==="
echo "Started: $(date)"

Rscript - << 'REOF'
library(decontam)
library(phyloseq)

# ============================================================
# Config
# ============================================================
BASE_DIR <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"
KEEP_LIST <- file.path(BASE_DIR, "control_quality_viz/genuine_controls_keep.csv")
DECONTAM_DIR <- file.path(BASE_DIR, "decontam_results")
dir.create(DECONTAM_DIR, showWarnings=FALSE, recursive=TRUE)

# Load kept control sample list
keep <- read.csv(KEEP_LIST, stringsAsFactors=FALSE)
cat(sprintf("Kept control samples: %d across %d studies\n", nrow(keep), length(unique(keep$bioproject))))

# Build per-study control sample list
studies <- unique(keep$bioproject)
cat(sprintf("Processing %d studies with controls...\n", length(studies)))

all_results <- list()

for (study in studies) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Study: %s\n", study))

  # Check if full ASV table exists
  asv_file <- file.path("/beegfs/dongbiao/oral", study, "results", "ASVs_counts.tsv")
  fasta_file <- file.path("/beegfs/dongbiao/oral", study, "results", "ASVs.fa")

  if (!file.exists(asv_file) || !file.exists(fasta_file)) {
    cat(sprintf("  SKIP: ASV table or FASTA not found\n"))
    next
  }

  # Load ASV table
  asv_tab <- read.table(asv_file, header=TRUE, row.names=1, check.names=FALSE, sep="\t")
  # asv_tab: rows = ASVs, cols = samples

  # Load FASTA for sequences
  asv_seqs <- character()
  fasta_lines <- readLines(fasta_file)
  for (i in seq(1, length(fasta_lines), 2)) {
    header <- sub("^>", "", fasta_lines[i])
    seq <- fasta_lines[i+1]
    asv_seqs[[header]] <- seq
  }

  # Determine which samples are controls for this study
  study_controls <- keep[keep$bioproject == study, "sample"]
  # The counts table uses sample names like "ERRxxx.fastq.gz" or "ERRxxx"
  # Match by checking if any control sample name is a prefix of the column name

  sample_names <- colnames(asv_tab)
  is_control <- logical(length(sample_names))
  matched_controls <- 0

  for (j in seq_along(sample_names)) {
    sn <- sample_names[j]
    # Strip .fastq.gz, .R1.filtered.fastq.gz suffixes for matching
    sn_base <- gsub("\\.fastq\\.gz$", "", sn)
    sn_base <- gsub("\\.R1\\.filtered\\.fastq\\.gz$", "", sn_base)
    # Also try matching against study_controls entries (which may have various formats)
    for (ctrl in study_controls) {
      ctrl_base <- gsub("\\.fastq\\.gz$", "", ctrl)
      ctrl_base <- gsub("\\.R1\\.filtered\\.fastq\\.gz$", "", ctrl_base)
      if (sn_base == ctrl_base || grepl(ctrl_base, sn, fixed=TRUE) || grepl(sn_base, ctrl, fixed=TRUE)) {
        is_control[j] <- TRUE
        matched_controls <- matched_controls + 1
        break
      }
    }
  }

  cat(sprintf("  Total samples: %d\n", length(sample_names)))
  cat(sprintf("  Control samples matched: %d / %d expected\n", matched_controls, length(study_controls)))

  if (matched_controls == 0) {
    cat(sprintf("  SKIP: no controls matched in ASV table\n"))
    next
  }
  if (matched_controls == length(sample_names)) {
    cat(sprintf("  SKIP: ALL samples are controls (no oral samples to compare)\n"))
    next
  }

  # Convert to presence/absence
  asv_pa <- asv_tab > 0
  asv_pa <- as.matrix(asv_pa)
  mode(asv_pa) <- "integer"

  # Run decontam prevalence test
  cat(sprintf("  Running decontam prevalence test...\n"))

  # decontam expects: isContaminant(seqtab, conc=NULL, method="prevalence", neg=is_neg)
  # where neg=TRUE for negative control samples
  contam_result <- tryCatch({
    isContaminant(asv_pa, method="prevalence", neg=is_control, threshold=0.1)
  }, error = function(e) {
    cat(sprintf("  ERROR in decontam: %s\n", e$message))
    return(NULL)
  })

  if (is.null(contam_result)) next

  # Add ASV IDs and p-values
  contam_result$asv_id <- rownames(asv_tab)

  # Add sequences
  contam_result$sequence <- sapply(contam_result$asv_id, function(x) {
    if (x %in% names(asv_seqs)) asv_seqs[[x]] else ""
  })

  # Summary
  n_contam <- sum(contam_result$contaminant, na.rm=TRUE)
  n_total <- nrow(contam_result)
  cat(sprintf("  ASVs flagged as contaminants: %d / %d (%.1f%%)\n",
              n_contam, n_total, 100 * n_contam / n_total))

  # Show p-value distribution
  p_breaks <- c(0, 0.01, 0.05, 0.1, 0.5, 1.0)
  for (k in 2:length(p_breaks)) {
    n_in_range <- sum(contam_result$p >= p_breaks[k-1] & contam_result$p < p_breaks[k], na.rm=TRUE)
    cat(sprintf("    p in [%.2f, %.2f): %d\n", p_breaks[k-1], p_breaks[k], n_in_range))
  }

  # Export per-study results
  out_file <- file.path(DECONTAM_DIR, sprintf("%s_decontam.csv", study))
  write.csv(contam_result, out_file, row.names=FALSE)
  cat(sprintf("  -> %s\n", out_file))

  # Store for combined summary
  all_results[[study]] <- contam_result
}

# ============================================================
# Combined summary
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("COMBINED decontam SUMMARY\n"))
cat(sprintf("========================================\n"))

n_studies_processed <- length(all_results)
cat(sprintf("Studies processed: %d / %d\n", n_studies_processed, length(studies)))

# Collect all contaminant ASV sequences from across all studies
all_contam_sequences <- character()
all_contam_studies <- character()
all_contam_pvals <- numeric()

for (study in names(all_results)) {
  res <- all_results[[study]]
  contam_rows <- res[res$contaminant & !is.na(res$contaminant), ]
  if (nrow(contam_rows) > 0) {
    all_contam_sequences <- c(all_contam_sequences, contam_rows$sequence)
    all_contam_studies <- c(all_contam_studies, rep(study, nrow(contam_rows)))
    all_contam_pvals <- c(all_contam_pvals, contam_rows$p)
  }
}

cat(sprintf("\nTotal contaminant-flagged ASVs across all studies: %d\n", length(all_contam_sequences)))

# Deduplicate by sequence
unique_seqs <- unique(all_contam_sequences)
cat(sprintf("Unique contaminant sequences: %d\n", length(unique_seqs)))

# Count in how many studies each unique sequence appears
seq_to_studies <- list()
for (i in seq_along(all_contam_sequences)) {
  seq <- all_contam_sequences[i]
  study <- all_contam_studies[i]
  if (is.null(seq_to_studies[[seq]])) {
    seq_to_studies[[seq]] <- character()
  }
  seq_to_studies[[seq]] <- unique(c(seq_to_studies[[seq]], study))
}

n_studies_per_seq <- sapply(seq_to_studies, length)
cat(sprintf("  In 1 study:  %d\n", sum(n_studies_per_seq == 1)))
cat(sprintf("  In 2 studies: %d\n", sum(n_studies_per_seq == 2)))
cat(sprintf("  In 3+ studies: %d\n", sum(n_studies_per_seq >= 3)))

# Export combined contaminant FASTA with metadata
combined_fasta <- file.path(DECONTAM_DIR, "decontam_flagged_contaminants.fasta")
fh <- file(combined_fasta, "w")
for (i in seq_along(unique_seqs)) {
  seq <- unique_seqs[i]
  n_studies <- n_studies_per_seq[i]
  writeLines(sprintf(">decontam_%06d n_studies=%d studies=%s",
                     i, n_studies, paste(seq_to_studies[[seq]], collapse="|")), fh)
  writeLines(seq, fh)
}
close(fh)
cat(sprintf("\nCombined FASTA: %s (%d sequences)\n", combined_fasta, length(unique_seqs)))

cat(sprintf("\n=== decontam done at %s ===\n", Sys.time()))
REOF

echo ""
echo "=== Finished at $(date) ==="
ls -lh "${DECONTAM_DIR}/"
