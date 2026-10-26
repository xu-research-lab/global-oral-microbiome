#!/usr/bin/env Rscript
# sparcc_cohort.R — genus-level SparCC for one (group, cohort) unit of the
# HC/PD/DC/OSCC all-site networks (body site ignored). Parameters identical
# to the per-site / allnet / disease pipelines (SpiecEasi sparcc:
# iter=50, inner_iter=10, th=0.01).
#
# Usage:  Rscript sparcc_cohort.R <counts.tsv> <out_csv>
#   <counts.tsv> : per-unit genus x sample raw-count matrix (prepared by
#                  prepare_network.py)
#   <out_csv>    : output genus x genus SparCC correlation matrix
#
# Atomic write (tmp + rename) so an interrupted run never leaves a half-written
# matrix that a later run would mistake for a valid cache.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
infile  <- args[1]
outfile <- args[2]

library(SpiecEasi)

m <- as.matrix(read.delim(infile, row.names = 1, check.names = FALSE, quote = ""))
if (nrow(m) < 2) stop("fewer than 2 features after filtering: ", infile)
if (ncol(m) < 3) stop("too few samples (", ncol(m), "): ", infile)
if (any(colSums(m) == 0)) {
  stop("sample(s) with all-zero counts after feature filtering: ", infile)
}

sparcc_res <- sparcc(t(m), iter = 50, inner_iter = 10, th = 0.01)
cor_df <- as.data.frame(sparcc_res$Cor)
rownames(cor_df) <- rownames(m)
colnames(cor_df) <- rownames(m)   # Cor is feature x feature: use feature names

outdir <- dirname(outfile)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}
tmpfile <- paste0(outfile, ".tmp")
write.csv(cor_df, tmpfile, row.names = TRUE)
if (!file.rename(tmpfile, outfile)) {
  file.copy(tmpfile, outfile, overwrite = TRUE)
  unlink(tmpfile)
}
cat(sprintf("sparcc done: %d genera, %d samples -> %s\n",
            nrow(m), ncol(m), outfile))
