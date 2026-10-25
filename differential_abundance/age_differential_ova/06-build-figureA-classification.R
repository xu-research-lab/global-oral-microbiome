# ==============================================================================
# 06-build-figureA-classification.R  (age OVA da_v2)
# Build figureA_classification_{2,3}method.tsv for the saliva age-stage
# consensus, with the SAME rule as the site / continent Figure A scripts:
#
#   robust (n_methods >= min) genus whose mean relative abundance is a UNIQUE
#   maximum (Enriched) / minimum (Depleted) across the six age stages, with a
#   significant MaAsLin2 coefficient of the matching sign supported by
#   >= min significant methods.
#
# No DA re-run; reads the stage-02 consensus tables only and reuses
# ova_extreme_counts() from differential_abundance/ova_plot_helpers.R.
# Group means are colMeans over the HC saliva samples of each age group
# (same scope as the Figure B heatmap of the age OVA pipeline).
#
# Usage (from project root):
#   Rscript differential_abundance/age_differential_ova/06-build-figureA-classification.R
# Output:
#   <out>/figures_saliva_target_top10_globalpeak_maaslin2coef/
#         figureA_classification_{2,3}method.tsv
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(data.table)
})

source("differential_abundance/age_differential_ova/da_helpers.R")  # AGE_ORDER, parse_ova_comparison, load_age_metadata, COL_*, RELABUND_TSV, MIN_SAMPLES_PER_GROUP, DA_ALPHA
source("differential_abundance/ova_plot_helpers.R")

root <- Sys.getenv("DA_RESULTS_ROOT")
if (!nzchar(root)) root <- "differential_abundance/result"
out_dir <- file.path(root, "age_hc_ova", "age_da_consensus")
figures_dir <- file.path(out_dir, "figures_saliva_target_top10_globalpeak_maaslin2coef")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

cat("==========================================\n")
cat("Age OVA figure-A classification (unique extremum + matching sign)\n")
cat("Start: ", format(Sys.time()), "\n", sep = "")
cat("==========================================\n\n")

enh <- fread(file.path(out_dir, "consensus_all_enhanced.tsv"), data.table = FALSE)
cons <- fread(file.path(out_dir, "consensus_all.tsv"), data.table = FALSE)
cat("    Enhanced rows:", nrow(enh), "  consensus rows:", nrow(cons), "\n")

enh$site <- vapply(lapply(enh$comparison, parse_ova_comparison),
                   `[[`, character(1), "site")
enh$target <- vapply(lapply(enh$comparison, parse_ova_comparison),
                     `[[`, character(1), "target")
if (any(is.na(enh$target))) stop("Unparsable OVA comparison in enhanced table")

full <- fread(file.path(out_dir, "full_table.tsv"), data.table = FALSE)
tested_targets <- unique(vapply(lapply(unique(full$comparison), parse_ova_comparison),
                                `[[`, character(1), "target"))
tested_targets <- AGE_ORDER[AGE_ORDER %in% tested_targets]
if (!length(tested_targets)) stop("No tested OVA targets found in full_table.tsv")
cat("    Tested targets: ", paste(tested_targets, collapse = ", "), "\n", sep = "")

# ---- per-age-stage mean relative abundance over the HC saliva scope ----
# Same scope as the age Figure B heatmap: load_age_metadata(hc_only=TRUE),
# saliva only, colMeans over all samples of each stage (zeros included).
robust_feats <- unique(c(enh$feature[enh$n_methods >= 2],
                         cons$feature[cons$n_methods >= 2]))
pats <- paste0(robust_feats, "\t")
cmd <- paste0("head -1 ", RELABUND_TSV, "; grep -F -e '",
              paste(pats, collapse = "' -e '"), "' ", RELABUND_TSV)
rel <- fread(cmd = cmd, sep = "\t", header = TRUE, check.names = FALSE,
             data.table = FALSE, quote = "")

meta <- load_age_metadata(hc_only = TRUE)
meta <- meta[meta[[COL_SITE]] == "saliva" & meta[[COL_AGE]] %in% AGE_ORDER, , drop = FALSE]
samples_all <- intersect(colnames(rel), rownames(meta))
sdf_all <- as.data.frame(t(rel[, samples_all, drop = FALSE]))
colnames(sdf_all) <- rel[[1]]
age_vec_all <- meta[rownames(sdf_all), COL_AGE]

mean_abundance_df <- do.call(rbind, lapply(AGE_ORDER, function(tgt) {
  idx <- which(age_vec_all == tgt)
  if (length(idx) == 0L) stop("No abundance samples for age group: ", tgt)
  means <- colMeans(sdf_all[idx, , drop = FALSE])
  data.frame(feature = names(means), target = tgt,
             target_mean_abundance = as.numeric(means))
}))
cat("    Features in mean table: ", nrow(mean_abundance_df) / length(AGE_ORDER),
    "  HC saliva samples: ", nrow(sdf_all), "\n", sep = "")

# ---- same Figure A classification rule as site / continent OVA ----
for (min_count_methods in c(2L, 3L)) {
  ver_count <- paste0(min_count_methods, "method")
  count_result <- ova_extreme_counts(enh, mean_abundance_df, AGE_ORDER,
                                     tested_targets, min_count_methods,
                                     alpha = DA_ALPHA)
  fwrite(count_result$detail,
         file.path(figures_dir, paste0("figureA_classification_", ver_count, ".tsv")),
         sep = "\t")
  ok <- count_result$counts %>% filter(n_genera > 0)
  cat("    ", ver_count, ": total passing = ",
      sum(count_result$counts$n_genera), "\n", sep = "")
  for (i in seq_len(nrow(ok)))
    cat("        ", ok$target[i], " / ", ok$direction[i],
        " : ", ok$n_genera[i], "\n", sep = "")
}
cat("\nDone: figureA_classification_{2,3}method.tsv written to ", figures_dir, "\n")