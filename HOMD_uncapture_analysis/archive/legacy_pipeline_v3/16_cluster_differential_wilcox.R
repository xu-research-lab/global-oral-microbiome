#!/usr/bin/env Rscript
# Differential relative abundance of the supported clusters, disease vs HC.
#
# Cohorts are the same as in 15_cluster_cohort_detection: a disease status
# (status_cleaned) with more than 100 subjects, compared against the HC samples
# of the same bioprojects.  Here the 15 cohorts with the most samples (disease
# + HC) are used; ICU has no HC and so cannot be tested and is skipped.
#
# The test is a two-sided Wilcoxon rank-sum test on the sample-level relative
# abundance (14_cluster_relabund_dedup.tsv, cluster reads / all ASV reads in the
# sample), run for every cluster in every cohort, with Benjamini-Hochberg
# correction applied across the 486 clusters within each cohort.
#
# Samples are the unit of observation, as in the matrices; a subject that
# contributed several samples enters several times, so the p-values treat the
# samples as independent and are anti-conservative for such subjects.
#
# Output: 16_cluster_differential_wilcox.tsv
suppressPackageStartupMessages({
  library(dplyr)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
REL  <- file.path(PIPE, "14_cluster_relabund_dedup.tsv")
META <- file.path(PIPE, "14_cluster_sample_metadata.tsv")
OUT  <- file.path(PIPE, "16_cluster_differential_wilcox.tsv")
N_COHORT <- 15

rel <- as.matrix(read.delim(REL, sep = "\t", row.names = 1,
                            check.names = FALSE, stringsAsFactors = FALSE))
storage.mode(rel) <- "double"
meta <- read.delim(META, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(sample_key = "character"))
stopifnot(identical(colnames(rel), meta$sample_key))

st <- meta$status_cleaned
subj <- table(st)
# "" is the unannotated remainder (Python reads it as NaN and folds it into NA);
# "Other" is a catch-all rather than a disease, so neither can be a cohort
diseases <- names(subj)[nzchar(names(subj)) &
                        !names(subj) %in% c("HC", "NA", "Other") & subj > 100]

# cohort size = the samples the test is run on (disease arm + its HC)
arm <- lapply(diseases, function(d) {
  di <- which(st == d)
  hi <- which(st == "HC" & meta$bioproject_id %in% meta$bioproject_id[di])
  data.frame(cohort = d, n_dis = length(di), n_hc = length(hi),
             n_subj = length(unique(meta$subject_id[di])),
             n_total = length(di) + length(hi))
})
arm <- bind_rows(arm) %>% arrange(desc(n_total))
arm <- arm[arm$n_hc > 0, ]
arm <- head(arm, N_COHORT)
cat("[16] cohorts (samples, disease + HC):\n")
print(arm, row.names = FALSE)

res <- vector("list", nrow(arm))
for (k in seq_len(nrow(arm))) {
  d <- arm$cohort[k]
  di <- which(st == d)
  hi <- which(st == "HC" & meta$bioproject_id %in% meta$bioproject_id[di])
  x <- rel[, di, drop = FALSE]
  y <- rel[, hi, drop = FALSE]

  p <- numeric(nrow(rel)); w <- numeric(nrow(rel))
  for (i in seq_len(nrow(rel))) {
    xi <- x[i, ]; yi <- y[i, ]
    # a cluster absent from both arms has nothing to rank; p = 1 by convention
    if (all(c(xi, yi) == 0)) { p[i] <- 1; w[i] <- NA_real_; next }
    t <- suppressWarnings(wilcox.test(xi, yi, exact = FALSE))
    p[i] <- t$p.value
    w[i] <- unname(t$statistic)
  }

  md <- rowMeans(x); mh <- rowMeans(y)
  eps <- 1e-9
  res[[k]] <- data.frame(
    cohort = d,
    cluster_id = rownames(rel),
    n_disease = length(di),
    n_hc = length(hi),
    n_subjects_disease = length(unique(meta$subject_id[di])),
    n_studies = length(unique(meta$bioproject_id[di])),
    prev_disease = rowMeans(x > 0),
    prev_hc = rowMeans(y > 0),
    mean_disease = md,
    mean_hc = mh,
    median_disease = apply(x, 1, median),
    median_hc = apply(y, 1, median),
    log2fc = log2((md + eps) / (mh + eps)),
    W = w,
    p = p,
    p_adj = p.adjust(p, method = "BH"),
    stringsAsFactors = FALSE
  )
  cat(sprintf("[16] %2d/%d %-32s %d vs %d samples, %d clusters, %d with FDR<0.05\n",
              k, nrow(arm), d, length(di), length(hi), nrow(rel),
              sum(res[[k]]$p_adj < 0.05)))
}

out <- bind_rows(res)
write.table(out, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[16] wrote %s: %d rows\n", OUT, nrow(out)))
cat(sprintf("[16] FDR<0.05 per cohort: %s\n",
            paste(sprintf("%s=%d", arm$cohort,
                          sapply(res, function(z) sum(z$p_adj < 0.05))),
                  collapse = ", ")))
