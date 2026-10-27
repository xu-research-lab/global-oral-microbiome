#!/usr/bin/env Rscript
# MaAsLin3 differential abundance and prevalence, one disease cohort at a time.
#
# For each cohort selected in 15_cohorts/results/15_cohorts.tsv the samples are
# the disease arm plus the HC samples of the same bioprojects.  Both MaAsLin3
# models are fitted: the abundance model (log relative abundance among detected
# samples) and the prevalence model (presence/absence).  The disease contrast is
# adjusted for age, sex, sampling site, sequencing region and continent, with
# bioproject as a random intercept.
#
# The covariates are the cleaned metadata columns.  They are not complete:
# age is missing for about half the samples and sex for two in five, so each
# cohort is reduced to the samples complete on the covariates that vary within
# it.  A covariate with a single level in a cohort is dropped from that cohort's
# model (it cannot be estimated), as is the random effect when the cohort comes
# from a single bioproject.  Cohorts left with too few samples in either arm are
# skipped and the reason is recorded in 16_maaslin3_status.tsv.
#
# Sequencing region and continent are constant within every bioproject in this
# dataset, i.e. they are study-level variables.  They cannot be estimated on top
# of a study term -- with the random intercept the study variance goes to zero
# and the prevalence fit becomes singular, which makes its standard errors
# meaningless.  So in multi-project cohorts only covariates that vary within at
# least one project enter the model, and the study term carries the region and
# continent adjustment.  In single-project cohorts region and continent are
# constant and drop out anyway (site is kept where it varies).
#
# The abundance table fed in is the 486-cluster relative abundance plus an
# "other" column that makes each sample sum to one.  MaAsLin3's TSS step then
# divides by the whole-community total, so the values it models are the true
# relative abundances rather than proportions of the 486 clusters alone; the
# "other" feature is dropped from the combined results.
#
# Output: 16_maaslin3/results/maaslin3/<cohort>/ (MaAsLin3's own folder) and
#         16_maaslin3_all_results.tsv, 16_maaslin3_significant.tsv,
#         16_maaslin3_status.tsv
suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline"
REL  <- file.path(PIPE, "14_feature_profile/results/14_cluster_relabund.tsv")
META <- file.path(PIPE, "14_feature_profile/results/14_cluster_sample_metadata.tsv")
COHORT_TSV <- file.path(PIPE, "15_cohorts/results/15_cohorts.tsv")
ROOT <- file.path(PIPE, "16_maaslin3", "results", "maaslin3")
OUT_ALL <- file.path(PIPE, "16_maaslin3/results/16_maaslin3_all_results.tsv")
OUT_SIG <- file.path(PIPE, "16_maaslin3/results/16_maaslin3_significant.tsv")
OUT_STATUS <- file.path(PIPE, "16_maaslin3/results/16_maaslin3_status.tsv")

COV <- c("age_cleaned", "sex_cleaned", "site_cleaned",
         "sequencing_region_cleaned", "continent_group")
MIN_ARM <- 20          # samples needed in each arm after the complete-case cut
MIN_PREV <- 0.01       # a feature must be present in >1% of the cohort
CORES <- as.integer(Sys.getenv("MAASLIN_CORES", "8"))

only <- commandArgs(trailingOnly = TRUE)   # optional: run a single cohort
dir.create(ROOT, showWarnings = FALSE, recursive = TRUE)

meta <- read.delim(META, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(sample_key = "character"))
rel <- as.matrix(read.delim(REL, sep = "\t", row.names = 1, check.names = FALSE))
storage.mode(rel) <- "double"
stopifnot(identical(colnames(rel), meta$sample_key))

cohorts <- read.delim(COHORT_TSV, sep = "\t", stringsAsFactors = FALSE) %>%
  arrange(desc(n_disease + n_hc))
if (length(only)) cohorts <- cohorts[cohorts$cohort %in% only, ]
stopifnot(nrow(cohorts) > 0)

st <- meta$status_cleaned
results <- list(); status <- list()

for (k in seq_len(nrow(cohorts))) {
  d <- cohorts$cohort[k]
  tag <- gsub("[^A-Za-z0-9]+", "_", d)
  outdir <- file.path(ROOT, tag)
  di <- which(st == d)
  hi <- which(st == "HC" & meta$bioproject_id %in% meta$bioproject_id[di])
  idx <- c(di, hi)

  md <- meta[idx, ]
  for (v in COV) md[[v]][md[[v]] == ""] <- NA
  ok <- complete.cases(md[, COV])
  md <- md[ok, ]
  md$condition <- factor(ifelse(md$status_cleaned == "HC", "HC", "Disease"),
                         levels = c("HC", "Disease"))
  n_d <- sum(md$condition == "Disease"); n_h <- sum(md$condition == "HC")

  # Covariates that actually vary here, and -- when the cohort spans several
  # projects -- that vary *within* at least one project.  Sequencing region and
  # continent are constant within every project in this dataset: they are
  # project-level variables, so they cannot be estimated alongside a project
  # term (random or fixed) and are absorbed by it.  Keeping them would drive
  # the project variance to the boundary and leave the model singular.
  cov_use <- COV[vapply(COV, function(v) length(unique(md[[v]])) > 1, logical(1))]
  if (length(unique(md$bioproject_id)) > 1) {
    within <- vapply(cov_use, function(v) {
      any(tapply(md[[v]], md$bioproject_id, function(z) length(unique(z))) > 1)
    }, logical(1))
    cov_use <- cov_use[within]
  }
  # MaAsLin3 treats any term whose name contains "group" as its own group()
  # syntax, so continent_group is passed under the plain name; it is the same
  # column, unchanged
  if ("continent_group" %in% cov_use) {
    md$continent <- md$continent_group
    cov_use[cov_use == "continent_group"] <- "continent"
  }
  # categorical covariates as factors with the commonest level as reference, so
  # every coefficient reads against the dominant category
  for (v in cov_use) {
    if (is.character(md[[v]]) || is.factor(md[[v]])) {
      md[[v]] <- factor(md[[v]],
                        levels = names(sort(table(md[[v]]), decreasing = TRUE)))
    }
  }
  n_proj <- length(unique(md$bioproject_id))

  # a cohort can survive the complete-case cut in a single project that happens
  # to have no supported cluster at all; there is nothing to model there
  n_detected <- sum(rowMeans(rel[, md$sample_key, drop = FALSE] > 0) > 0)
  if (n_d >= MIN_ARM && n_h >= MIN_ARM && n_detected == 0) {
    status[[k]] <- data.frame(cohort = d, n_disease = n_d, n_hc = n_h,
                              n_disease_raw = length(di), n_hc_raw = length(hi),
                              n_projects = n_proj,
                              covariates = paste(cov_use, collapse = ","),
                              ranef = n_proj > 1,
                              status = "skipped (no cluster detected in the complete-case samples)")
    cat(sprintf("[16] SKIP %-32s no cluster detected in %d complete-case samples\n",
                d, length(md$sample_key)))
    next
  }

  if (n_d < MIN_ARM || n_h < MIN_ARM) {
    status[[k]] <- data.frame(cohort = d, n_disease = n_d, n_hc = n_h,
                              n_disease_raw = length(di), n_hc_raw = length(hi),
                              n_projects = n_proj,
                              covariates = paste(cov_use, collapse = ","),
                              ranef = n_proj > 1, status = "skipped (too few samples)")
    cat(sprintf("[16] SKIP %-32s %d disease / %d HC after complete cases\n", d, n_d, n_h))
    next
  }

  # with a formula, random effects have to be written into it; the reference
  # level of group is set by the factor levels above (HC first)
  form <- paste("condition", paste(cov_use, collapse = " + "), sep = " + ")
  ranef <- n_proj > 1
  if (ranef) form <- paste0(form, " + (1|bioproject_id)")

  # samples x features; "other" completes each row to the whole community
  f <- t(rel[, md$sample_key, drop = FALSE])
  f <- as.data.frame(f)
  f$ZZ_other <- pmax(0, 1 - rowSums(f))
  stopifnot(max(abs(rowSums(f) - 1)) < 1e-6)

  # MaAsLin3 matches samples between the two inputs by row name, so the
  # metadata rows have to carry the sample keys (MaAsLin3 expects samples as
  # rows in the metadata and as rows in the feature table here)
  md_in <- md[, c("condition", cov_use, "bioproject_id")]
  rownames(md_in) <- md$sample_key

  cat(sprintf("[16] %-32s %d disease / %d HC, %d projects, formula: %s%s\n",
              d, n_d, n_h, n_proj, form, ""))
  unlink(outdir, recursive = TRUE)   # so a rerun never reads old output
  t0 <- Sys.time()
  err <- tryCatch({
  maaslin3(
    input_data = f,
    input_metadata = md_in,
    output = outdir,
    formula = form,
    min_abundance = 0,
    min_prevalence = MIN_PREV,
    max_significance = 0.05,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = TRUE,
    # TRUE is MaAsLin3's default: the abundance model keeps bioproject as a
    # random intercept, and the prevalence model enters it as a fixed effect.
    # A logistic random intercept is not estimable here -- with a handful of
    # projects and sparse clusters it separates or returns boundary variances
    # (standard errors in the hundreds).  Treating project as a fixed effect
    # adjusts for it just the same and fits cleanly; the two models still share
    # one formula, so their coefficient sets agree.
    small_random_effects = ranef,
    plot_summary_plot = FALSE,
    plot_associations = FALSE,
    cores = CORES,
    save_models = FALSE,
    verbosity = "WARN"
  )
  NULL }, error = function(e) conditionMessage(e))
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (!is.null(err)) {
    status[[k]] <- data.frame(cohort = d, n_disease = n_d, n_hc = n_h,
                              n_disease_raw = length(di), n_hc_raw = length(hi),
                              n_projects = n_proj,
                              covariates = paste(cov_use, collapse = ","),
                              ranef = n_proj > 1,
                              status = paste("failed:", gsub("\\s+", " ", err)))
    cat(sprintf("[16] FAILED %-30s %s\n", d, err))
    next
  }

  res <- read.delim(file.path(outdir, "all_results.tsv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  names(res)[1] <- "feature"
  res <- res[res$feature != "ZZ_other" & res$metadata == "condition", ]
  # the error column can hold line breaks from lme4; left alone they would
  # break the unquoted table into junk rows when it is read back
  for (v in names(res)) if (is.character(res[[v]])) {
    res[[v]] <- gsub("[\t\r\n]+", " ", res[[v]])
  }
  res$cohort <- d
  results[[k]] <- res

  status[[k]] <- data.frame(cohort = d, n_disease = n_d, n_hc = n_h,
                            n_disease_raw = length(di), n_hc_raw = length(hi),
                            n_projects = n_proj,
                            covariates = paste(cov_use, collapse = ","),
                            ranef = n_proj > 1,
                            status = sprintf("fitted (%.1f min)", mins))
  cat(sprintf("[16] done %-30s %.1f min, %d rows\n", d, mins, nrow(res)))
}

all <- if (length(results)) bind_rows(results) else data.frame()
if (file.exists(OUT_ALL)) {
  old <- read.delim(OUT_ALL, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE, colClasses = c(feature = "character"))
  old <- old[!old$cohort %in% cohorts$cohort, ]
  # the annotation is joined again below, so the copy written by an earlier run
  # has to go; kept, the join would rename both sides to .x/.y and every partial
  # rerun would add another pair of columns
  old <- old[, !grepl("^(annotation|phylum)", names(old)), drop = FALSE]
  all <- bind_rows(old, all)
}
if (nrow(all)) {
  ann <- read.delim(file.path(PIPE, "14_feature_profile/results/14_supported_cluster_annotation.tsv"),
                    stringsAsFactors = FALSE,
                    colClasses = c(cluster_id = "character"))
  all <- left_join(all, ann[, c("cluster_id", "annotation", "phylum")],
                   by = c("feature" = "cluster_id"))
  write.table(all, OUT_ALL, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA")
}
# ---- the called set ---------------------------------------------------------
# Rows whose individual q value passes 0.05, one per cohort x phylotype x
# model.  The joint q value is carried alongside but is not the filter: it
# asks whether either model is associated, while this table lists each model
# that is.  Sorted by qval_individual, so the strongest calls come first.
if (nrow(all)) {
  sig <- all[!is.na(all$qval_individual) & all$qval_individual < 0.05, ]
  sig$direction <- ifelse(sig$coef > 0, "Disease", "HC")
  sig <- sig[order(sig$qval_individual, method = "radix"), ]
  write.table(sig, OUT_SIG, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA")
  cat(sprintf("[16] wrote %s (%d rows)\n", OUT_SIG, nrow(sig)))
}
if (length(status)) {
  st_out <- bind_rows(status)
  if (file.exists(OUT_STATUS)) {
    old_st <- read.delim(OUT_STATUS, sep = "\t", stringsAsFactors = FALSE)
    st_out <- bind_rows(old_st[!old_st$cohort %in% cohorts$cohort, ], st_out)
  }
  write.table(st_out, OUT_STATUS, sep = "\t", quote = FALSE, row.names = FALSE)
}
cat(sprintf("[16] wrote %s (%d rows across %d cohorts)\n", OUT_ALL, nrow(all),
            n_distinct(all$cohort)))
