#!/usr/bin/env Rscript
# Phase 2: decontam (prevalence mode) within each blank-bearing BioProject.
#
# Gold-standard, blank-based contaminant identification for the ~23 BioProjects that carry
# negative controls (see 25_negative_control_audit.py). Runs decontam::isContaminant in
# prevalence mode per BioProject, then aggregates: a feature is a high-confidence
# cross-lab contaminant if flagged in multiple independent BioProjects. That aggregated
# catalog is later broadcast (as a FLAG, not an auto-delete) to blank-free studies by
# 26_flag_contaminants_low_homd.py.
#
# Inputs
#   --counts    wide TSV: rows = features, col 1 = feature_id, remaining cols = sample_ids
#   --metadata  TSV: sample_id, bioproject, is_blank (TRUE/FALSE)
#   --threshold decontam prevalence P* (default 0.5; stricter = more aggressive)
# Output
#   <outdir>/decontam_per_bioproject.tsv   feature x bioproject contaminant flags
#   <outdir>/decontam_contaminant_catalog.tsv  features flagged in >= min_projects BPs
#
# Cluster step (needs R + decontam). Not executed locally.

suppressMessages({
  library(optparse); library(decontam)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--outdir", type = "character", default = "novel_oral_taxa/03_decontam"),
  make_option("--threshold", type = "double", default = 0.5),
  make_option("--min-projects", type = "integer", default = 2)
)))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
counts <- read.delim(opt$counts, row.names = 1, check.names = FALSE)
md <- read.delim(opt$metadata, stringsAsFactors = FALSE)
md$is_blank <- as.logical(md$is_blank)

flags <- list()
for (bp in unique(md$bioproject)) {
  sm <- md[md$bioproject == bp, ]
  if (sum(sm$is_blank) < 1 || sum(!sm$is_blank) < 2) next
  ids <- intersect(sm$sample_id, colnames(counts))
  if (length(ids) < 3) next
  sub <- t(as.matrix(counts[, ids, drop = FALSE]))                 # samples x features
  neg <- sm$is_blank[match(ids, sm$sample_id)]
  keep <- colSums(sub) > 0
  sub <- sub[, keep, drop = FALSE]
  res <- isContaminant(sub, neg = neg, method = "prevalence", threshold = opt$threshold)
  flags[[bp]] <- data.frame(feature_id = rownames(res),
                            contaminant = res$contaminant, p = res$p, bioproject = bp)
}
if (length(flags) == 0) stop("no BioProject had both blanks and >=2 true samples")

all_flags <- do.call(rbind, flags)
write.table(all_flags, file.path(opt$outdir, "decontam_per_bioproject.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# aggregate: contaminant in >= min_projects independent BioProjects
agg <- aggregate(contaminant ~ feature_id, data = all_flags, FUN = function(x) sum(x, na.rm = TRUE))
colnames(agg)[2] <- "n_projects_flagged"
catalog <- agg[agg$n_projects_flagged >= opt$`min-projects`, ]
write.table(catalog[order(-catalog$n_projects_flagged), ],
            file.path(opt$outdir, "decontam_contaminant_catalog.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("BioProjects processed: %d; cross-lab contaminant features (>=%d projects): %d\n",
            length(flags), opt$`min-projects`, nrow(catalog)))
