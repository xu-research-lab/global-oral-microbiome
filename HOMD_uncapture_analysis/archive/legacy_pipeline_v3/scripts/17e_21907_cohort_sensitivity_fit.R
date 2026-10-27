#!/usr/bin/env Rscript
# Cluster 21907: does the OSCC contrast survive the way the cohort is chosen?
#
# Two definitions of the cohort, both starting from the OSCC/HC cohort of
# 17b_plot_cluster_prevalence.R (all OSCC samples plus the HC samples of the
# same bioprojects):
#
#   S1  ever  -- the bioprojects in which the cluster is detected in at least
#                one sample.  Projects that never carry it say nothing about it
#                and only dilute the denominator; eight of the seventeen fall
#                away, leaving nine (three of which have HC samples).
#   S2  both  -- of those, only the bioprojects that have BOTH arms, so every
#                case is compared against controls from its own study.
#
# For each, two models are fitted with MaAsLin3, which is what the rest of the
# disease analyses use:
#
#   abundance   log relative abundance among samples where the cluster is
#               detected (non-zero)
#   prevalence  presence / absence over every sample of the cohort
#
# Covariates are the sampling site and a bioproject random intercept.  Age and
# sex are NOT used: they are missing for 578 of the 787 S1 samples, so a
# complete-case cut on them would silently throw away three quarters of the
# cohort and change what is being compared.  Site is complete.  In S2 site is
# constant inside every bioproject, so it is absorbed by the study term and
# dropped from the formula rather than fitted twice.
#
# Output
#   17e_21907_cohort_sensitivity_results.tsv    one row per scenario x model
#   17e_21907_cohort_sensitivity_samples.tsv    per-sample data for the figure
suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
REL  <- file.path(PIPE, "14_cluster_relabund_dedup.tsv")
META <- file.path(PIPE, "14_cluster_sample_metadata.tsv")
ROOT <- file.path(PIPE, "maaslin3_sensitivity")
CLUSTER <- "c21907"
STATUS  <- "OSCC"

OUT_RES <- file.path(PIPE, "17e_21907_cohort_sensitivity_results.tsv")
OUT_SMP <- file.path(PIPE, "17e_21907_cohort_sensitivity_samples.tsv")

rel <- as.matrix(read.delim(REL, sep = "\t", row.names = 1, check.names = FALSE))
storage.mode(rel) <- "double"
meta <- read.delim(META, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(sample_key = "character"))
stopifnot(identical(colnames(rel), meta$sample_key))
stopifnot(CLUSTER %in% rownames(rel))

st <- meta$status_cleaned
bp <- meta$bioproject_id
di <- which(st == STATUS)
hc <- which(st == "HC" & bp %in% bp[di])
coh <- c(di, hc)
detected <- rel[CLUSTER, ] > 0

# ---- the two cohort definitions --------------------------------------------
ever <- unique(bp[coh][detected[coh]])
both <- ever[vapply(ever, function(b) any(bp[di] == b) && any(bp[hc] == b),
                    logical(1))]
SCEN <- list(
  S1 = list(label = "All cohorts carrying the cluster",
            projects = sort(ever)),
  S2 = list(label = "Cohorts with both cases and controls",
            projects = sort(both)))

samples <- bind_rows(lapply(names(SCEN), function(k) {
  idx <- coh[bp[coh] %in% SCEN[[k]]$projects]
  data.frame(scenario = k,
             sample_key = meta$sample_key[idx],
             group = ifelse(st[idx] == "HC", "HC", "Disease"),
             bioproject_id = bp[idx],
             site = meta$site_cleaned[idx],
             rel_abund = rel[CLUSTER, idx],
             detected = detected[idx],
             stringsAsFactors = FALSE)
}))
write.table(samples, OUT_SMP, sep = "\t", quote = FALSE, row.names = FALSE)

# ---- MaAsLin3 fitting ------------------------------------------------------
dir.create(ROOT, showWarnings = FALSE, recursive = TRUE)
results <- list(); status <- list()

for (k in names(SCEN)) {
  d <- samples[samples$scenario == k, ]
  d$condition <- factor(d$group, levels = c("HC", "Disease"))
  md <- d[, c("condition", "site", "bioproject_id")]
  rownames(md) <- d$sample_key
  md$site <- factor(md$site, levels = names(sort(table(md$site),
                                                 decreasing = TRUE)))
  n_proj <- length(unique(md$bioproject_id))
  # site is study-level in S2: every project has one sampling site there, so a
  # site term and a project term would be the same variable twice
  site_varies <- any(tapply(md$site, md$bioproject_id,
                            function(z) length(unique(z))) > 1)
  form <- if (site_varies) "condition + site + (1|bioproject_id)" else
    "condition + (1|bioproject_id)"

  # samples x features; ZZ_other completes each row to the whole community so
  # TSS normalises by the true total rather than by this one cluster
  f <- data.frame(rel_abund = d$rel_abund)
  f$ZZ_other <- pmax(0, 1 - f$rel_abund)
  rownames(f) <- d$sample_key
  names(f)[1] <- CLUSTER

  outdir <- file.path(ROOT, k)
  unlink(outdir, recursive = TRUE)
  cat(sprintf("[17e] %s: %d cases / %d HC, %d projects, formula: %s\n",
              k, sum(d$group == "Disease"), sum(d$group == "HC"), n_proj, form))
  err <- tryCatch({
    maaslin3(input_data = f, input_metadata = md, output = outdir,
             formula = form, min_abundance = 0, min_prevalence = 0,
             max_significance = 0.05, normalization = "TSS", transform = "LOG",
             correction = "BH", standardize = TRUE,
             small_random_effects = n_proj > 1,
             plot_summary_plot = FALSE, plot_associations = FALSE,
             cores = 4, save_models = FALSE, verbosity = "WARN")
    NULL
  }, error = function(e) conditionMessage(e))

  if (!is.null(err)) {
    cat(sprintf("[17e] %s FAILED: %s\n", k, err))
    status[[k]] <- data.frame(scenario = k, model = NA, error = err)
    next
  }
  res <- read.delim(file.path(outdir, "all_results.tsv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  names(res)[1] <- "feature"
  res <- res[res$feature != "ZZ_other" & res$metadata == "condition", ]
  res$scenario <- k
  res$n_disease <- sum(d$group == "Disease")
  res$n_hc <- sum(d$group == "HC")
  res$n_projects <- n_proj
  res$formula <- form
  results[[k]] <- res
  cat(sprintf("[17e] %s: %d model rows\n", k, nrow(res)))
}

if (length(results)) {
  out <- bind_rows(results) %>%
    select(scenario, model, coef, stderr, pval_individual, qval_individual,
           pval_joint, qval_joint, N, N_not_zero, n_disease, n_hc,
           n_projects, formula, error)
  write.table(out, OUT_RES, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA")
  cat("[17e] wrote ", OUT_RES, "\n", sep = "")
} else {
  stop("no model fitted")
}
