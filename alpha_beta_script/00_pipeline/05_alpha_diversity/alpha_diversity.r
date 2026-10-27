library(rhdf5)
library(ggplot2)
library(dplyr)

biom_path <- "/home/yemai/xulaoshi/oral20260528/alpha/table_gg2_id99_cleaned_62782_115345_all.biom"
meta_path <- "/home/yemai/xulaoshi/oral20260528/alpha/metadata-v2-merge-cleaned.csv"

# ── 1. Read sparse matrix components ──────────────────────────────────────────
cat("Reading biom file...\n")
sample_ids <- as.character(h5read(biom_path, "/sample/ids"))
data_vec   <- h5read(biom_path, "/sample/matrix/data")
indices    <- as.integer(h5read(biom_path, "/sample/matrix/indices"))
indptr     <- as.integer(h5read(biom_path, "/sample/matrix/indptr"))
n_samples  <- length(sample_ids)
cat("Samples:", n_samples, "\n")

# ── 2. Compute alpha diversity directly from sparse representation ────────────
# For each sample, extract its non-zero OTU counts and compute metrics.
cat("Computing alpha diversity per sample...\n")

observed_vec   <- numeric(n_samples)
shannon_vec    <- numeric(n_samples)
simpson_vec    <- numeric(n_samples)
inv_simpson_vec <- numeric(n_samples)
chao1_vec      <- numeric(n_samples)

for (j in seq_len(n_samples)) {
  start <- indptr[j] + 1       # R is 1-based, indptr is 0-based
  end   <- indptr[j + 1]       # 0-based exclusive → R inclusive
  if (start > end) {
    observed_vec[j]    <- 0
    shannon_vec[j]     <- 0
    simpson_vec[j]     <- 0
    inv_simpson_vec[j] <- 0
    chao1_vec[j]       <- 0
    next
  }
  counts <- data_vec[start:end]
  # Observed richness (number of non-zero OTUs)
  observed_vec[j] <- length(counts)
  # Total reads
  N <- sum(counts)
  if (N == 0) next
  p <- counts / N
  # Shannon
  shannon_vec[j] <- -sum(p * log(p))
  # Simpson (1 - D)
  simpson_vec[j] <- 1 - sum(p^2)
  # InvSimpson (1/D)
  inv_simpson_vec[j] <- 1 / sum(p^2)
  # Chao1
  f1 <- sum(counts == 1)
  f2 <- sum(counts == 2)
  if (f2 > 0) {
    chao1_vec[j] <- observed_vec[j] + f1^2 / (2 * f2)
  } else {
    chao1_vec[j] <- observed_vec[j] + f1 * (f1 - 1) / 2
  }
  if (j %% 10000 == 0) cat("  processed", j, "of", n_samples, "samples\n")
}
cat("Done computing alpha diversity.\n")

H5close()

# ── 3. Build data frame ──────────────────────────────────────────────────────
alpha_df <- data.frame(
  sample_id  = sample_ids,
  Observed    = observed_vec,
  Shannon     = shannon_vec,
  Simpson     = simpson_vec,
  InvSimpson  = inv_simpson_vec,
  Chao1       = chao1_vec,
  stringsAsFactors = FALSE
)

# ── 4. Read metadata & merge ──────────────────────────────────────────────────
cat("Reading metadata...\n")
meta <- read.csv(meta_path, stringsAsFactors = FALSE)

alpha_df <- merge(alpha_df, meta, by = "sample_id", all.x = TRUE)
alpha_df$bioproject <- as.factor(alpha_df$bioproject)

cat("Alpha diversity for", nrow(alpha_df), "samples across",
    length(levels(alpha_df$bioproject)), "datasets\n")

# ── 5. Output CSV ────────────────────────────────────────────────────────────
csv_path <- "/home/yemai/xulaoshi/oral20260528/alpha/alpha_diversity_data.csv"
write.csv(alpha_df, csv_path, row.names = FALSE)