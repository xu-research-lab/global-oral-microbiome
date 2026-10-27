#!/usr/bin/env Rscript
# Cohorts that enter the disease models.
#
# A cohort is a disease status (status_cleaned) with more than 100 samples,
# compared against the HC samples of the same BioProjects.  ICU has no HC and
# so cannot be compared; it is dropped by the n_hc > 0 rule.  The 15 cohorts
# with the most samples (disease + HC) are kept.
#
# This table is the only input step 16 needs from the disease side: it fixes
# which cohorts are modelled and their sample counts.  Sample membership of a
# cohort is recomputed there from the same metadata, so the counts here are
# descriptive, not a filter on the model.
#
# Output: 15_cohorts/results/15_cohorts.tsv
suppressPackageStartupMessages({
  library(dplyr)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline"
META <- file.path(PIPE, "14_feature_profile/results/14_cluster_sample_metadata.tsv")
OUT  <- file.path(PIPE, "15_cohorts/results/15_cohorts.tsv")
N_COHORT <- 15

meta <- read.delim(META, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c(sample_key = "character"))

st <- meta$status_cleaned
subj <- table(st)
# "" is the unannotated remainder (Python reads it as NaN and folds it into NA);
# "Other" is a catch-all rather than a disease, so neither can be a cohort
diseases <- names(subj)[nzchar(names(subj)) &
                        !names(subj) %in% c("HC", "NA", "Other") & subj > 100]

# cohort size = the samples the model is run on (disease arm + its HC)
arm <- lapply(diseases, function(d) {
  di <- which(st == d)
  hi <- which(st == "HC" & meta$bioproject_id %in% meta$bioproject_id[di])
  data.frame(cohort = d, n_disease = length(di), n_hc = length(hi),
             n_subject_disease = length(unique(meta$subject_id[di])),
             n_studies_disease = length(unique(meta$bioproject_id[di])),
             stringsAsFactors = FALSE)
})
arm <- bind_rows(arm) %>% arrange(desc(n_disease + n_hc))
arm <- arm[arm$n_hc > 0, ]
arm <- head(arm, N_COHORT)

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.table(arm, OUT, sep = "\t", quote = FALSE, row.names = FALSE)

cat("[15] cohorts (samples, disease + HC):\n")
print(arm, row.names = FALSE)
cat(sprintf("[15] wrote %s: %d cohorts\n", OUT, nrow(arm)))
