#!/usr/bin/env Rscript
# The significant clusters of the head-and-neck and periodontal cohorts, re-fitted
# on the cohorts that have controls.
#
# 17_maaslin3_disease.R takes every sample of the disease status plus the HC
# samples of the same bioprojects.  A bioproject that contributes cases but no
# controls can then inflate the case arm without contributing anything to the
# within-study contrast; cluster 21907 showed what that does (17e).  This script
# re-fits the clusters that came out significant there on the restricted cohort:
# only the bioprojects that have BOTH arms, so every case is compared against
# controls from its own study.  Within that cohort a cluster is fitted on the
# bioprojects in which it was detected, because a study where the cluster never
# occurs cannot speak to its association with disease; the same samples are what
# 17f_plot_ctrl_only.R draws the prevalence bars on.
#
# Model: condition + site + (1|bioproject_id), dropping the random intercept for
# a cluster detected in a single study (with one study there is no
# between-study contrast to adjust, and the intercept is not estimable).  Age and
# sex are not used.  They are missing for most samples of exactly the cohorts in
# question -- with them a complete-case cut leaves 117 cases / 243 HC from a
# single HNC project, which cannot answer the question the restriction is meant
# to ask.  Site is complete
# and varies within projects, so it is fitted.  Because the covariates differ
# from the original model the p-values are not a like-for-like correction of the
# published ones; they are the same clusters tested on a cohort that no longer
# mixes studies with and without controls.
#
# Output: 17f_ctrl_only_results.tsv, 17f_ctrl_only_samples.tsv
suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
REL  <- file.path(PIPE, "14_cluster_relabund_dedup.tsv")
META <- file.path(PIPE, "14_cluster_sample_metadata.tsv")
RES  <- file.path(PIPE, "17_maaslin3_all_results.tsv")
ROOT <- file.path(PIPE, "maaslin3_ctrl_only")

DISEASES <- c("Head and Neck Cancer", "Periodontal Disease")
OUT_RES  <- file.path(PIPE, "17f_ctrl_only_results.tsv")
OUT_SMP  <- file.path(PIPE, "17f_ctrl_only_samples.tsv")

rel <- as.matrix(read.delim(REL, sep = "\t", row.names = 1, check.names = FALSE))
storage.mode(rel) <- "double"
meta <- read.delim(META, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(sample_key = "character"))
stopifnot(identical(colnames(rel), meta$sample_key))
res0 <- read.delim(RES, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(feature = "character"))

# the group numbering below is not stable between runs, so an old set directory
# would otherwise be left behind and read as if it belonged to this run
unlink(ROOT, recursive = TRUE)
dir.create(ROOT, showWarnings = FALSE, recursive = TRUE)
st <- meta$status_cleaned
bp <- meta$bioproject_id
all_res <- list(); all_smp <- list()

for (d in DISEASES) {
  di <- which(st == d)
  pr <- unique(bp[di])
  hi <- which(st == "HC" & bp %in% pr)
  both <- pr[vapply(pr, function(b) any(bp[hi] == b), logical(1))]
  idx <- c(di[bp[di] %in% both], hi[bp[hi] %in% both])

  md <- meta[idx, ]
  md$condition <- factor(ifelse(md$status_cleaned == "HC", "HC", "Disease"),
                         levels = c("HC", "Disease"))
  md$site <- factor(md$site_cleaned,
                    levels = names(sort(table(md$site_cleaned),
                                        decreasing = TRUE)))
  n_proj <- length(unique(md$bioproject_id))

  # the clusters the original cohort called significant in either model
  sig <- res0 %>% filter(cohort == d, !is.na(qval_individual),
                         qval_individual < 0.05) %>%
    distinct(feature) %>% pull(feature)
  sig <- intersect(sig, rownames(rel))
  if (!length(sig)) next
  cat(sprintf("[17f] %-22s %d clusters; %d cases / %d HC, %d projects\n",
              d, length(sig), sum(md$condition == "Disease"),
              sum(md$condition == "HC"), n_proj))

  # Each cluster is tested only on the samples of the bioprojects in which it
  # was detected.  A study in which the cluster never occurred cannot say
  # anything about the cluster's association with disease; left in, those
  # samples made the prevalence model compare one study against the others
  # rather than cases against controls, which is what a cluster found in a
  # single bioproject was being tested on.  The samples of a cluster's own
  # studies are also what the figure draws its bars on, so the percentages and
  # the p value are now computed on the same samples.
  #
  # maaslin3 fits every feature it is given on one sample table, so the clusters
  # are grouped by the set of bioprojects they were detected in and a call is
  # made per group.  Within a group the random intercept is estimable only when
  # the group spans more than one study, so a single-study group is fitted
  # without it.
  proj_of <- function(id) {
    v <- rel[id, md$sample_key] > 0
    sort(unique(md$bioproject_id[v]))
  }
  proj <- setNames(lapply(sig, proj_of), sig)
  sets <- split(sig, vapply(proj, paste, character(1), collapse = ","))
  cat(sprintf("[17f]   %d distinct sets of detecting bioprojects\n", length(sets)))

  rows <- list()
  for (j in seq_along(sets)) {
    ids <- sets[[j]]
    pr  <- proj[[ids[1]]]
    smp <- md[md$bioproject_id %in% pr, ]
    fml <- if (length(pr) > 1) "condition + site + (1|bioproject_id)" else
      "condition + site"

    f <- t(rel[ids, smp$sample_key, drop = FALSE])
    f <- as.data.frame(f)
    f$ZZ_other <- pmax(0, 1 - rowSums(f))
    stopifnot(max(abs(rowSums(f) - 1)) < 1e-6)
    rownames(f) <- smp$sample_key

    md_in <- smp[, c("condition", "site", "bioproject_id")]
    md_in$site <- droplevels(md_in$site)
    md_in$bioproject_id <- droplevels(factor(md_in$bioproject_id))
    rownames(md_in) <- smp$sample_key

    outdir <- file.path(ROOT, paste(gsub("[^A-Za-z0-9]+", "_", d),
                                    sprintf("set%02d", j), sep = "_"))
    unlink(outdir, recursive = TRUE)
    maaslin3(input_data = f, input_metadata = md_in, output = outdir,
             formula = fml, min_abundance = 0, min_prevalence = 0,
             max_significance = 0.05, normalization = "TSS", transform = "LOG",
             correction = "BH", standardize = TRUE,
             small_random_effects = length(pr) > 1,
             plot_summary_plot = FALSE, plot_associations = FALSE,
             cores = 4, save_models = FALSE, verbosity = "WARN")

    r <- read.delim(file.path(outdir, "all_results.tsv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
    names(r)[1] <- "feature"
    r <- r[r$feature != "ZZ_other" & r$metadata == "condition", ]
    r$feature_set <- if (length(pr) > 1) "multi" else "single"
    r$projects <- paste(pr, collapse = ",")
    r$n_projects_detected <- length(pr)
    r$n_samples_used <- nrow(smp)
    rows[[j]] <- r
    cat(sprintf("[17f]   set %d: %2d clusters, %2d samples from %d project(s), %s\n",
                j, length(ids), nrow(smp), length(pr), fml))
  }
  r <- bind_rows(rows)
  r$disease <- d
  r$n_disease <- sum(md$condition == "Disease")
  r$n_hc <- sum(md$condition == "HC")
  r$n_projects <- n_proj
  r$formula <- ifelse(r$feature_set == "multi",
                      "condition + site + (1|bioproject_id)",
                      "condition + site")
  all_res[[d]] <- r

  smp <- md[, c("sample_key", "bioproject_id", "site_cleaned", "condition")]
  smp$disease <- d
  all_smp[[d]] <- smp
  cat(sprintf("[17f] %-22s %d model rows\n", d, nrow(r)))
}

out <- bind_rows(all_res) %>%
  select(disease, feature, model, feature_set, coef, stderr, pval_individual,
         qval_individual, pval_joint, qval_joint, N, N_not_zero,
         n_disease, n_hc, n_projects, projects, n_projects_detected,
         n_samples_used, formula, error)
write.table(out, OUT_RES, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(bind_rows(all_smp), OUT_SMP, sep = "\t", quote = FALSE,
            row.names = FALSE)
cat("[17f] wrote ", OUT_RES, " and ", OUT_SMP, "\n", sep = "")
