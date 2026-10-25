# =============================================================================
# 01-run-three-method-da.R  (sex_differential_continent/pooled - da_v2 tsv)
# POOLED variant of ../01-run-three-method-da.R: ONE saliva F_vs_M
# comparison with ALL continents merged. Main effect = sex_cleaned;
# continent_group is added as a COVARIATE (no sex:continent interaction:
# Africa has 0 male saliva samples, an interaction would be rank deficient;
# prepare_covars() also folds F-only continents to "Other" automatically).
#
# Everything else mirrors the parent pipeline (da_v2, da_v4 aligned):
#   - counts tsv loaded once (no .qza / qiime2R OOM peak); reads>=1000 export
#   - Maaslin2 relabund derived in-memory per pair (TSS becomes a no-op)
#   - ANCOM-BC2: phyloseq, genus-fittability pre-filter, struc_zero=FALSE,
#     collinearity retry (drop offending covariates, then main-effect-only)
#   - conditional (1|bioproject_id) random effect when a BioProject spans
#     both sexes in the pooled saliva data
#   - age/bmi covariates use the pre-binned *_group columns
#
# Comparison (single):  saliva + pooled + F_vs_M   (ref = F, target = M =>
# log2FoldChange > 0 => enriched in MALES). Results land under
#   results root/sex_pooled/     (--scope=all)
#   results root/sex_pooled_hc/  (--scope=hc)
#   {sex_linda,sex_maaslin2,sex_ancombc2}_pairwise/saliva_pooled_F_vs_M/
#
# Extra flag vs the parent driver:  --smoke  subsamples a few BioProjects
# x capped samples per sex (<= 6 BioProjects, <= 40 samples per sex) for an
# end-to-end smoke run.
#
# Usage (from project root):
#   Rscript differential_abundance/sex_differential_continent/pooled/01-run-three-method-da.R --scope=hc
#   Rscript differential_abundance/sex_differential_continent/pooled/01-run-three-method-da.R --scope=hc --method=ancombc2 --smoke
#
# Idempotent: each method/pair skips if its analysis.done sentinel exists.
# =============================================================================

suppressPackageStartupMessages({
  library(MicrobiomeStat)
  library(Maaslin2)
  library(ANCOMBC)
  library(phyloseq)
  library(dplyr)
  library(data.table)
})

source("differential_abundance/sex_differential_continent/da_helpers.R")

set.seed(DA_SEED)

# ==========================================================================
# 0. Command-line arguments
# ==========================================================================
args <- commandArgs(trailingOnly = TRUE)
methods_to_run <- "all"
only_arg <- NULL
smoke <- "--smoke" %in% args
scope <- parse_scope_arg(args)
for (a in args) {
  if (grepl("^--method=", a)) methods_to_run <- sub("^--method=", "", a)
  else if (grepl("^--only=", a)) only_arg <- sub("^--only=", "", a)
}
ALL_METHODS <- c("linda", "maaslin2", "ancombc2")
if (methods_to_run == "all") methods_to_run <- ALL_METHODS
methods_to_run <- unlist(strsplit(methods_to_run, ",", fixed = TRUE))
bad <- setdiff(methods_to_run, ALL_METHODS)
if (length(bad) > 0)
  stop("Unknown --method value(s): ", paste(bad, collapse = ", "),
       " (valid: all | ", paste(ALL_METHODS, collapse = ","), ")")

# ==========================================================================
# 1. Configuration
# ==========================================================================
results_root <- if (scope == "hc") rd("sex_pooled_hc") else rd("sex_pooled")
group_kind_tag <- if (scope == "hc") "sex (pooled continents), saliva (HC only)" else "sex (pooled continents), saliva (all)"
main_col     <- COL_SEX
g1_col       <- "Sex1"
g2_col       <- "Sex2"
covar_candidates <- c(COL_CONTINENT,   # pooled: continent is a COVARIATE now
                      COL_AGE, COL_PLATFORM, COL_SEQ_REGION, COL_BMI)
              # sex = main effect; site fixed = saliva; continent = stratum
              # (continent_group excluded: no within-stratum variance)
n_cores <- 8L   # LinDA n.cores / ANCOMBC2 n_cl
main_results_filename <- "sex_results.tsv"

method_dirs <- list(
  linda    = file.path(results_root, "sex_linda_pairwise"),
  maaslin2 = file.path(results_root, "sex_maaslin2_pairwise"),
  ancombc2 = file.path(results_root, "sex_ancombc2_pairwise")
)
for (m in methods_to_run)
  if (!dir.exists(method_dirs[[m]])) dir.create(method_dirs[[m]], recursive = TRUE)

cat("==========================================\n")
cat("Three-method sex DA (pooled F vs M, saliva) — da_v2 tsv\n")
cat("Scope: ", group_kind_tag, "\n", sep = "")
cat("Methods: ", paste(methods_to_run, collapse = ", "), "\n", sep = "")
cat("Start: ", format(Sys.time()), "\n")
cat("==========================================\n\n")

# ==========================================================================
# 2. Load data (once, shared by all methods)
# ==========================================================================
cat("[1/2] Loading data...\n")
genus_counts <- load_genus_table("counts")
metadata     <- load_sex_metadata(hc_only = (scope == "hc"))
aligned      <- align_samples(metadata, genus_counts)
metadata     <- aligned$meta
genus_counts <- aligned$mat
rm(aligned); invisible(gc())

# keep only samples with a usable sex label
keep_sex <- metadata[[COL_SEX]] %in% c(SEX_REF, SEX_TGT)
metadata <- metadata[keep_sex, , drop = FALSE]
genus_counts <- genus_counts[, metadata[[COL_SAMPLE]], drop = FALSE]

# --smoke: subsample to a few BioProjects x capped samples per sex so the
# whole three-method chain finishes in minutes (<= 6 BioProjects, <= 40/sex).
if (smoke) {
  set.seed(DA_SEED)
  bps <- head(names(sort(table(metadata[[COL_PROJECT]]), decreasing = TRUE)), 6L)
  keep_idx <- unlist(lapply(c(SEX_REF, SEX_TGT), function(s) {
    cand <- which(metadata[[COL_SEX]] == s & metadata[[COL_PROJECT]] %in% bps)
    if (length(cand) > 40L) cand <- sample(cand, 40L)
    cand
  }), use.names = FALSE)
  keep_idx     <- sort(keep_idx)
  metadata     <- metadata[keep_idx, , drop = FALSE]
  genus_counts <- genus_counts[, metadata[[COL_SAMPLE]], drop = FALSE]
  cat("    [SMOKE] subsampled to ", nrow(metadata), " samples\n", sep = "")
}


# =========================================================================
# 2b. POOLED: a single saliva F_vs_M comparison across ALL continents
# =========================================================================
n_f <- sum(metadata[[COL_SEX]] == SEX_REF, na.rm = TRUE)
n_m <- sum(metadata[[COL_SEX]] == SEX_TGT, na.rm = TRUE)
if (n_f < MIN_SAMPLES_PER_SEX || n_m < MIN_SAMPLES_PER_SEX)
  stop("Pooled saliva needs both sexes >= ", MIN_SAMPLES_PER_SEX, " samples!")
pairs_list <- list(list(
  site = SITE_SALIVA, continent = "pooled",
  g1 = SEX_REF, g2 = SEX_TGT,
  comparison_name = paste0(SITE_SALIVA, "_pooled_", SEX_REF, "_vs_", SEX_TGT)
))
cat("    Continents in the pool (covariate; F-only folds to Other):\n")
print(table(factor(metadata[[COL_CONTINENT]]),
            factor(metadata[[COL_SEX]], levels = c(SEX_REF, SEX_TGT))))
cat("    Pooled saliva samples: F=", n_f, " / M=", n_m, "\n", sep = "")
cat("    Total samples: ", nrow(metadata), "\n", sep = "")
cat("    Genus features: ", nrow(genus_counts), "\n", sep = "")
cat("    Comparisons: ", length(pairs_list), "\n", sep = "")

# optional: run a SINGLE continent (for SLURM array parallelism)
only_idx <- NA_integer_
if (!is.null(only_arg)) {
  if (grepl("^[0-9]+$", only_arg)) {
    only_idx <- as.integer(only_arg)
  } else {
    nm <- vapply(pairs_list, function(p) p$comparison_name, "")
    m <- which(nm == only_arg)
    if (length(m) == 1) only_idx <- m
  }
  if (is.na(only_idx) || only_idx < 1 || only_idx > length(pairs_list))
    stop("--only out of range (1:", length(pairs_list), "): ", only_arg)
  cat("    [array mode] running only continent ", only_idx, ": ",
      pairs_list[[only_idx]]$comparison_name, "\n", sep = "")
  pairs_list <- pairs_list[only_idx]
}

# ==========================================================================
# 3. Shared per-pair preparation (covariates, formula, random effect)
# ==========================================================================
prep_pair <- function(continent, g1, g2) {
  comparison_name <- paste0(sanitize(SITE_SALIVA), "_", sanitize(continent),
                            "_", g1, "_vs_", g2)
  keep <- metadata[[COL_SITE]] == SITE_SALIVA &
          metadata[[COL_SEX]] %in% c(g1, g2)
  pair_meta  <- metadata[keep, , drop = FALSE]
  pair_counts <- genus_counts[, pair_meta[[COL_SAMPLE]], drop = FALSE]
  pair_meta[[main_col]] <- factor(as.character(pair_meta[[main_col]]),
                                  levels = c(g1, g2))
  stopifnot(identical(colnames(pair_counts), pair_meta[[COL_SAMPLE]]))

  prep <- prepare_covars(pair_meta, main_col, covar_candidates)
  pair_meta    <- prep$meta
  valid_covars <- prep$covars

  fixed_terms <- if (length(valid_covars) > 0)
    paste(c(main_col, valid_covars), collapse = " + ") else main_col
  # da_v4 style: random effect appended to the formula string (lme4 syntax)
  re <- add_random_effect(paste0("~ ", fixed_terms), pair_meta, main_col)
  pair_meta[[COL_PROJECT]] <- droplevels(factor(as.character(pair_meta[[COL_PROJECT]])))

  for (cov in c(main_col, valid_covars)) {
    if (is.factor(pair_meta[[cov]])) pair_meta[[cov]] <- droplevels(pair_meta[[cov]])
    if (length(unique(as.character(pair_meta[[cov]]))) < 2)
      stop("Factor '", cov, "' has <2 levels for ", comparison_name)
  }

  list(site = SITE_SALIVA, continent = continent, g1 = g1, g2 = g2,
       comparison_name = comparison_name,
       meta = pair_meta, counts = pair_counts,
       valid_covars = valid_covars, fixed_terms = fixed_terms, re = re)
}

# ---- per-pair summary/config writers (shared by all methods) ---------------
make_summary_row <- function(pp, n1, n2, formula_label) {
  data.frame(
    Sex1 = pp$g1, Sex2 = pp$g2, Site = pp$site, Continent = pp$continent,
    Comparison = pp$comparison_name,
    N1 = n1, N2 = n2, N_Significant = NA_integer_,
    N_Enriched_Sex1 = NA_integer_, N_Enriched_Sex2 = NA_integer_,
    Formula = formula_label, stringsAsFactors = FALSE
  )
}
fill_summary <- function(summary_row, output_subdir, res) {
  summary_row$N_Significant <- res$n_sig
  summary_row$N_Enriched_Sex1 <- res$n_up_g1
  summary_row$N_Enriched_Sex2 <- res$n_up_g2
  fwrite(summary_row, file.path(output_subdir, "summary.tsv"), sep = "\t", quote = FALSE)
}

# ==========================================================================
# 4a. LinDA
# ==========================================================================
run_pair_linda <- function(pair_index, pp) {
  output_subdir <- file.path(method_dirs$linda, pp$comparison_name)
  if (file.exists(file.path(output_subdir, "analysis.done"))) return("skip")

  pair_meta  <- pp$meta
  pair_genus <- pp$counts

  # lib>0 / rowSums>0 filter (da_v4 mirrors); reads>=1000 already in tsv
  lib <- colSums(pair_genus)
  keep_s <- lib > 0
  pair_genus <- pair_genus[, keep_s, drop = FALSE]
  pair_meta  <- pair_meta[keep_s, , drop = FALSE]
  keep_f <- rowSums(pair_genus) > 0
  pair_genus <- pair_genus[keep_f, , drop = FALSE]

  n1 <- sum(pair_meta[[main_col]] == pp$g1)
  n2 <- sum(pair_meta[[main_col]] == pp$g2)
  group_coef <- ancombc_group_coef(main_col, pp$g2)  # "<col><level>" convention
  formula_str <- pp$re$full_formula   # "~ base + (1|bioproject_id)" when shared

  cat("    Formula: ", formula_str, "\n", sep = "")
  cat("    Samples: ", pp$g1, " (", n1, ") vs ", pp$g2, " (", n2, ")\n", sep = "")
  cat("    LinDA: feature.dat.type=count, prev.filter=", PREV_CUT,
      ", is.winsor=FALSE, n.cores=", n_cores, "\n", sep = "")

  if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)

  config_df <- data.frame(
    pair_index = pair_index, comparison = pp$comparison_name,
    site = pp$site, continent = pp$continent, sex1 = pp$g1, sex2 = pp$g2,
    n1 = n1, n2 = n2,
    formula = formula_str, feature_dat_type = "count",
    prev_filter = PREV_CUT, is_winsor = FALSE, outlier_pct = 0.03,
    p_adj_method = "BH", alpha = DA_ALPHA, n_cores = n_cores,
    counts_tsv = COUNTS_TSV, metadata_file = METADATA_FILE,
    scope = scope, stringsAsFactors = FALSE
  )

  ok <- tryCatch({
    fit <- MicrobiomeStat::linda(
      feature.dat        = pair_genus,
      meta.dat           = pair_meta,
      formula            = formula_str,
      feature.dat.type   = "count",
      prev.filter        = PREV_CUT,
      mean.abund.filter  = 0,
      max.abund.filter   = 0,
      is.winsor          = FALSE,
      outlier.pct        = 0.03,
      p.adj.method       = "BH",
      alpha              = DA_ALPHA,
      n.cores            = n_cores,
      verbose            = TRUE
    )

    # bias (often unnamed -> align to $variables / output names)
    bias_vals  <- as.numeric(fit$bias)
    bias_names <- names(fit$bias)
    if (is.null(bias_names) || length(bias_names) != length(bias_vals)) {
      if (!is.null(fit$variables) && length(fit$variables) == length(bias_vals)) {
        bias_names <- fit$variables
      } else if (!is.null(names(fit$output)) &&
                 length(names(fit$output)) == length(bias_vals)) {
        bias_names <- names(fit$output)
      } else {
        bias_names <- paste0("V", seq_along(bias_vals))
      }
    }
    fwrite(data.frame(variable = bias_names, bias = bias_vals, stringsAsFactors = FALSE),
           file.path(output_subdir, "bias.tsv"), sep = "\t")

    common_df <- linda_to_common(fit, group_coef)
    all_res <- harmonize_results(common_df, group1 = pp$g1, group2 = pp$g2,
                                 comparison = pp$comparison_name,
                                 target_site = pp$site,
                                 group1_col = g1_col, group2_col = g2_col)

    summary_row <- make_summary_row(pp, n1, n2, formula_str)
    res <- write_comparison_outputs(
      all_res, output_subdir, pp$comparison_name,
      group1 = pp$g1, group2 = pp$g2,
      group1_col = g1_col, group2_col = g2_col,
      main_results_var_pattern = paste0("^", main_col),
      summary_row = summary_row, run_config_row = config_df,
      main_results_filename = main_results_filename
    )
    fill_summary(summary_row, output_subdir, res)

    cat("    SUCCESS!\n")
    TRUE
  }, error = function(e) { cat("    ERROR: ", e$message, "\n", sep = ""); FALSE })

  if (isTRUE(ok)) "ok" else "fail"
}

# ==========================================================================
# 4b. Maaslin2
# ==========================================================================
run_pair_maaslin2 <- function(pair_index, pp) {
  output_subdir <- file.path(method_dirs$maaslin2, pp$comparison_name)
  if (file.exists(file.path(output_subdir, "analysis.done"))) return("skip")

  pair_meta <- pp$meta
  # relabund derived in-memory per pair (cols sum to 1): avoids loading the
  # second ~2 GB relabund tsv; Maaslin2's TSS normalization becomes a no-op.
  pair_genus <- sweep(pp$counts, 2, colSums(pp$counts), "/")

  valid_covars <- pp$valid_covars
  fixed_effects  <- c(main_col, valid_covars)
  random_effects <- if (pp$re$has_random) COL_PROJECT else NULL

  n1 <- sum(pair_meta[[main_col]] == pp$g1)
  n2 <- sum(pair_meta[[main_col]] == pp$g2)
  cat("    Formula: ", pp$re$full_formula, "\n", sep = "")
  cat("    Samples: ", pp$g1, " (", n1, ") vs ", pp$g2, " (", n2, ")\n", sep = "")
  cat("    Fixed: ", paste(fixed_effects, collapse = ", "),
      " | Random: ", ifelse(pp$re$has_random, COL_PROJECT, "none"), "\n", sep = "")

  if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)

  config_df <- data.frame(
    pair_index = pair_index, comparison = pp$comparison_name,
    site = pp$site, continent = pp$continent, sex1 = pp$g1, sex2 = pp$g2,
    n1 = n1, n2 = n2,
    formula = pp$re$full_formula, fixed_effects = paste(fixed_effects, collapse = ","),
    random_effects = ifelse(pp$re$has_random, COL_PROJECT, "none"),
    normalization = "TSS", transform = "LOG", analysis_method = "LM",
    correction = "BH", min_abundance = MIN_ABUNDANCE, min_prevalence = PREV_CUT,
    max_significance = DA_ALPHA, standardize = FALSE,
    counts_tsv = COUNTS_TSV, relabund = "derived_from_counts",
    metadata_file = METADATA_FILE, scope = scope,
    stringsAsFactors = FALSE
  )

  ok <- tryCatch({
    fit <- Maaslin2::Maaslin2(
      input_data       = t(pair_genus),
      input_metadata   = pair_meta[, c(COL_SAMPLE, main_col, valid_covars,
                                       if (pp$re$has_random) COL_PROJECT else NULL), drop = FALSE],
      output           = output_subdir,
      fixed_effects    = fixed_effects,
      random_effects   = random_effects,
      normalization    = "TSS",
      transform        = "LOG",
      analysis_method  = "LM",
      correction       = "BH",
      min_abundance    = MIN_ABUNDANCE,
      min_prevalence   = PREV_CUT,
      max_significance = DA_ALPHA,
      standardize      = FALSE,
      plot_heatmap     = FALSE,
      plot_scatter     = FALSE,
      cores            = 1L
    )

    maaslin_all <- read.delim(file.path(output_subdir, "all_results.tsv"),
                              check.names = FALSE, stringsAsFactors = FALSE)
    common_df <- maaslin2_to_common(maaslin_all, main_col, pp$g2)
    all_res <- harmonize_results(common_df, group1 = pp$g1, group2 = pp$g2,
                                 comparison = pp$comparison_name,
                                 target_site = pp$site,
                                 group1_col = g1_col, group2_col = g2_col)

    summary_row <- make_summary_row(pp, n1, n2, pp$re$full_formula)
    res <- write_comparison_outputs(
      all_res, output_subdir, pp$comparison_name,
      group1 = pp$g1, group2 = pp$g2,
      group1_col = g1_col, group2_col = g2_col,
      main_results_var_pattern = paste0("^", main_col),
      summary_row = summary_row, run_config_row = config_df,
      main_results_filename = main_results_filename
    )
    fill_summary(summary_row, output_subdir, res)

    # remove Maaslin2's heavy intermediate artifacts to keep disk bounded
    for (d in c("features", "fits")) unlink(file.path(output_subdir, d), recursive = TRUE)
    unlink(file.path(output_subdir, c("maaslin2.log", "significant_results.tsv")))

    cat("    SUCCESS!\n")
    TRUE
  }, error = function(e) { cat("    ERROR: ", e$message, "\n", sep = ""); FALSE })

  if (isTRUE(ok)) "ok" else "fail"
}

# ==========================================================================
# 4c. ANCOMBC2 (with collinearity retry)
# ==========================================================================
run_pair_ancombc2 <- function(pair_index, pp) {
  output_subdir <- file.path(method_dirs$ancombc2, pp$comparison_name)
  if (file.exists(file.path(output_subdir, "analysis.done"))) return("skip")

  pair_meta  <- pp$meta
  # da_v4 genus-fittability pre-filter (>=2 BioProjects AND >=2 groups)
  pair_counts <- filter_fittable(pp$counts, pair_meta, main_col)

  valid_covars <- pp$valid_covars
  rand_formula <- if (pp$re$has_random) paste0("(1|", COL_PROJECT, ")") else NULL
  group_coef   <- ancombc_group_coef(main_col, pp$g2)

  n1 <- sum(pair_meta[[main_col]] == pp$g1)
  n2 <- sum(pair_meta[[main_col]] == pp$g2)
  cat("    Formula: ", pp$re$full_formula, "\n", sep = "")
  cat("    Samples: ", pp$g1, " (", n1, ") vs ", pp$g2, " (", n2, ")\n", sep = "")

  if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)

  # inner fit: runs ancombc2 with the supplied covariate set and writes outputs.
  # Returns TRUE on success; throws on error (caught by caller for retry).
  fit_once <- function(covars, msg_tag = "fit") {
    ff <- if (length(covars) > 0)
      paste0(main_col, " + ", paste(covars, collapse = " + ")) else main_col
    keep_cols <- c(main_col, covars, if (pp$re$has_random) COL_PROJECT else NULL)
    ps <- phyloseq(
      otu_table(pair_counts, taxa_are_rows = TRUE),
      sample_data(pair_meta[, keep_cols, drop = FALSE])
    )
    t0 <- Sys.time()
    fit <- ANCOMBC::ancombc2(
      data = ps, assay.type = "counts", fix_formula = ff,
      rand_formula = rand_formula, p_adj_method = "BH", prv_cut = PREV_CUT,
      lib_cut = LIB_CUT, group = NULL, struc_zero = FALSE, neg_lb = FALSE,
      alpha = DA_ALPHA, n_cl = n_cores, global = FALSE, pairwise = FALSE, verbose = FALSE
    )
    cat("    ANCOMBC2 elapsed [", msg_tag, "]: ",
        round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n", sep = "")
    res_wide <- as.data.frame(fit$res)
    common_df <- ancombc2_to_common(res_wide, group_coef)
    all_res <- harmonize_results(common_df, group1 = pp$g1, group2 = pp$g2,
                                 comparison = pp$comparison_name,
                                 target_site = pp$site,
                                 group1_col = g1_col, group2_col = g2_col)
    config_df <- data.frame(
      pair_index = pair_index, comparison = pp$comparison_name,
      site = pp$site, continent = pp$continent, sex1 = pp$g1, sex2 = pp$g2,
      n1 = n1, n2 = n2,
      fix_formula = ff, rand_formula = ifelse(pp$re$has_random, rand_formula, "none"),
      p_adj_method = "BH", prv_cut = PREV_CUT, lib_cut = LIB_CUT,
      group_coef = group_coef, struc_zero = FALSE, neg_lb = FALSE,
      alpha = DA_ALPHA, n_cl = n_cores,
      counts_tsv = COUNTS_TSV, metadata_file = METADATA_FILE, scope = scope,
      stringsAsFactors = FALSE
    )
    summary_row <- make_summary_row(
      pp, n1, n2,
      paste0(ff, if (pp$re$has_random) paste0(" + ", rand_formula) else "")
    )
    res <- write_comparison_outputs(
      all_res, output_subdir, pp$comparison_name,
      group1 = pp$g1, group2 = pp$g2,
      group1_col = g1_col, group2_col = g2_col,
      main_results_var_pattern = paste0("^", main_col),
      summary_row = summary_row, run_config_row = config_df,
      main_results_filename = main_results_filename
    )
    fill_summary(summary_row, output_subdir, res)
    cat("    SUCCESS!\n")
    TRUE
  }

  # First attempt; on collinearity error, drop the offending covariate(s) and
  # retry (up to a few rounds, since removing one can expose the next).
  ok <- FALSE
  cur_covars <- valid_covars
  for (attempt in 1:4) {
    err_obj <- tryCatch(
      { fit_once(cur_covars, msg_tag = paste0("attempt", attempt)); NULL },
      error = function(e) e
    )
    if (is.null(err_obj)) { ok <- TRUE; break }
    cat("    ERROR: ", conditionMessage(err_obj), "\n", sep = "")
    drop <- collinear_covars_from_error(conditionMessage(err_obj), cur_covars)
    if (length(drop) == 0L) { ok <- FALSE; break }
    cat("    Collinearity detected; dropping covariate(s): ",
        paste(drop, collapse = ", "), " and retrying.\n", sep = "")
    cur_covars <- setdiff(cur_covars, drop)
    ok <- FALSE
  }
  if (!isTRUE(ok)) {
    # final attempt: main effect only, no covariates (last resort)
    err_obj <- tryCatch(
      { cat("    Final attempt: main effect only (no covariates).\n")
        fit_once(character(0), msg_tag = "mainonly"); NULL },
      error = function(e) e
    )
    ok <- is.null(err_obj)
    if (!ok) cat("    ERROR: ", conditionMessage(err_obj), "\n", sep = "")
  }

  if (isTRUE(ok)) "ok" else "fail"
}

# =========================================================================
# 5. Run the single pooled comparison (method x 1 pair)
# =========================================================================
RUNNERS <- list(linda = run_pair_linda, maaslin2 = run_pair_maaslin2,
                ancombc2 = run_pair_ancombc2)

cat("\n[2/2] Running pooled saliva DA (F vs M)...\n")
any_fail <- FALSE
for (m in methods_to_run) {
  cat("\n############################################\n")
  cat("### METHOD: ", toupper(m), "\n", sep = "")
  cat("############################################\n")
  n_success <- 0L; n_skip <- 0L; n_fail <- 0L
  for (i in seq_along(pairs_list)) {
    pp <- pairs_list[[i]]
    comparison_name <- pp$comparison_name

    cat("\n------------------------------------------\n")
    cat("[", i, "/", length(pairs_list), "] ", comparison_name, "\n", sep = "")

    res <- RUNNERS[[m]](i, prep_pair(pp$continent, pp$g1, pp$g2))
    if (res == "ok") { n_success <- n_success + 1L; cat("    Output: ",
        file.path(method_dirs[[m]], comparison_name), "\n", sep = "") }
    else if (res == "skip") { cat("    Already completed. Skipping.\n"); n_skip <- n_skip + 1L }
    else n_fail <- n_fail + 1L
    invisible(gc(verbose = FALSE))
  }
  cat("\n", toupper(m), " pooled finished\n", sep = "")
  cat("  Success: ", n_success, "\n  Skipped: ", n_skip, "\n  Failed: ", n_fail, "\n", sep = "")
  if (n_fail > 0) any_fail <- TRUE
}

cat("\n==========================================\n")
cat("All requested methods finished (scope: ", scope, ")\n", sep = "")
cat("End: ", format(Sys.time()), "\n")
cat("==========================================\n")
if (any_fail) quit(save = "no", status = 1)
