# ==============================================================================
# 01-run-three-method-da.R  (continent OVA da_v2 — tsv based)
# Merged driver for the three-method differential abundance (DA) pipeline:
# LinDA + Maaslin2 + ANCOMBC2 for CONTINENT, SALIVA ONLY. Main effect =
# continent_group, run as OVA (one-vs-ALL): each valid continent (>=
# MIN_SAMPLES_PER_GROUP samples) is the TARGET, and the reference "Other"
# pools the samples of ALL OTHER VALID continents. No nested stratum.
# Covariates: age_group, sex, platform, seq_region, bmi_group (continent =
# main effect, excluded; site excluded — no variation under saliva-only
# filter). Conditional (1|bioproject_id) random effect when a BioProject
# spans both groups (target vs Other).
#
# Contrast convention: reference = "Other", target = the FOCAL continent
# -> log2FoldChange > 0  =>  enriched in the FOCAL continent.
# Comparison name "<continent>_enrich" (e.g. Asia_enrich).
#
# Scope (--scope=all|hc): "all" keeps every sample with a valid continent;
# "hc" restricts to status_cleaned == "HC". Results land under
#   <DA_RESULTS_ROOT>/continent_ova/   (all)  or  <DA_RESULTS_ROOT>/continent_hc_ova/ (hc)
#   {continent_linda,continent_maaslin2,continent_ancombc2}_ova/<continent>_enrich/
#
# SLURM array support: --task=N (1-based GLOBAL index across BOTH scopes).
# Task order = all-scope comparisons (in CONTINENT_ORDER order) followed by
# hc-scope comparisons. The mapping is computed at runtime from the
# tsv-aligned metadata, so the sbatch only needs to pass
# --task=$SLURM_ARRAY_TASK_ID.
#
# Usage (from project root):
#   Rscript differential_abundance/continent_differential_ova/01-run-three-method-da.R
#       # all three methods, all OVA comparisons, scope=all
#   Rscript differential_abundance/continent_differential_ova/01-run-three-method-da.R --scope=hc
#   Rscript differential_abundance/continent_differential_ova/01-run-three-method-da.R --method=linda
#   Rscript differential_abundance/continent_differential_ova/01-run-three-method-da.R \
#       --method=ancombc2 --task=1            # global comparison index (array)
#   Rscript differential_abundance/continent_differential_ova/01-run-three-method-da.R \
#       --method=ancombc2 --only=1            # single comparison by 1-based index WITHIN --scope
#   Rscript differential_abundance/continent_differential_ova/01-run-three-method-da.R \
#       --method=ancombc2 --only=Asia_enrich
#
# Idempotent: each method/comparison skips if its analysis.done sentinel exists.
# ==============================================================================

suppressPackageStartupMessages({
  library(MicrobiomeStat)
  library(Maaslin2)
  library(ANCOMBC)
  library(phyloseq)
  library(dplyr)
  library(data.table)
})

source("differential_abundance/continent_differential_ova/da_helpers.R")

set.seed(DA_SEED)

# ==========================================================================
# 0. Command-line arguments
# ==========================================================================
args <- commandArgs(trailingOnly = TRUE)
methods_to_run <- "all"
only_arg <- NULL
task_arg <- NULL
scope <- parse_scope_arg(args)
for (a in args) {
  if (grepl("^--method=", a)) methods_to_run <- sub("^--method=", "", a)
  else if (grepl("^--only=", a)) only_arg <- sub("^--only=", "", a)
  else if (grepl("^--task=", a)) task_arg <- sub("^--task=", "", a)
}
ALL_METHODS <- c("linda", "maaslin2", "ancombc2")
if (methods_to_run == "all") methods_to_run <- ALL_METHODS
methods_to_run <- strsplit(methods_to_run, ",", fixed = TRUE)[[1]]
bad <- setdiff(methods_to_run, ALL_METHODS)
if (length(bad) > 0)
  stop("Unknown --method value(s): ", paste(bad, collapse = ", "),
       " (valid: all | ", paste(ALL_METHODS, collapse = ","), ")")

# ==========================================================================
# 1. Load counts ONCE (shared by scopes and methods)
# ==========================================================================
cat("[1/3] Loading genus counts tsv...\n")
genus_counts <- load_genus_table("counts")
cat("    Genus features: ", nrow(genus_counts),
    "; samples in table: ", ncol(genus_counts), "\n", sep = "")

# Build the aligned metadata + OVA comparison list for one scope (counts
# matrix reused).
setup_scope <- function(sc) {
  meta <- load_age_metadata(hc_only = (sc == "hc"))
  aligned <- align_samples(meta, genus_counts)
  meta <- aligned$meta
  ova <- make_age_ova_comparisons(meta, SITE_ORDER)
  list(meta = meta, ova = ova)
}

# ==========================================================================
# 2. Determine scope + comparisons (single scope, or global --task across both)
# ==========================================================================
cat("[2/3] Building OVA comparison list...\n")
if (!is.null(task_arg)) {
  if (!grepl("^[0-9]+$", task_arg)) stop("--task must be a positive integer: ", task_arg)
  task_idx <- as.integer(task_arg)
  setup_all <- setup_scope("all")
  setup_hc  <- setup_scope("hc")
  n_all <- length(setup_all$ova)
  n_hc  <- length(setup_hc$ova)
  cat("    Global task table: all-scope comparisons 1-", n_all,
      ", hc-scope comparisons ", n_all + 1, "-", n_all + n_hc, "\n", sep = "")
  if (task_idx < 1 || task_idx > n_all + n_hc) {
    # Oversized array: extra tasks exit cleanly so SLURM marks them done.
    cat("--task ", task_idx, " out of range (1:", n_all + n_hc,
        ") — nothing to do, exiting cleanly.\n", sep = "")
    quit(save = "no", status = 0)
  }
  if (task_idx <= n_all) {
    scope <- "all"
    scope_setup <- setup_all
    pairs_list <- setup_all$ova[task_idx]
  } else {
    scope <- "hc"
    scope_setup <- setup_hc
    pairs_list <- setup_hc$ova[task_idx - n_all]
  }
  cat("    [array mode] task ", task_idx, " -> scope=", scope, ", comparison index=",
      if (task_idx <= n_all) task_idx else task_idx - n_all, "\n", sep = "")
} else {
  scope_setup <- setup_scope(scope)
  pairs_list <- scope_setup$ova
}

metadata <- scope_setup$meta

if (length(pairs_list) < 1) stop("No continent has >= ",
                                 MIN_SAMPLES_PER_GROUP, " samples!")

cat("    Scope: ", if (scope == "hc") "HC only" else "all samples", "\n", sep = "")
cat("    OVA comparisons selected: ", length(pairs_list), " of a possible ",
    length(make_age_ova_comparisons(metadata, SITE_ORDER)), " in this scope\n", sep = "")
if (length(pairs_list) <= 6 || !is.null(task_arg) || !is.null(only_arg)) {
  cat("    Comparisons:\n")
  for (p in pairs_list) {
    n1 <- sum(metadata[[COL_CONTINENT]] %in% p$other_groups)
    n2 <- sum(metadata[[COL_CONTINENT]] == p$g2)
    cat("      ", p$comparison_name, "  (n_Other=", n1, ", n_", p$g2, "=", n2, ")\n", sep = "")
  }
}

# optional: run a SINGLE comparison within the chosen scope (--only, manual runs)
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
  cat("    [single mode] running only comparison ", only_idx, ": ",
      pairs_list[[only_idx]]$comparison_name, "\n", sep = "")
  pairs_list <- pairs_list[only_idx]
}

# ==========================================================================
# 3. Configuration (per selected scope)
# ==========================================================================
results_root <- continent_root(scope)
scope_tag <- if (scope == "hc") "continent OVA (HC only)" else "continent OVA (all)"
main_col     <- COL_CONTINENT
g1_col       <- "Continent1"   # "Other" (pooled other valid continents, reference)
g2_col       <- "Continent2"   # focal continent (target)
covar_candidates <- c(COL_AGE, COL_SEX, COL_PLATFORM, COL_SEQ_REGION, COL_BMI)
                      # continent is the main effect; site excluded (no variation
                      # under the saliva-only filter)
n_cores <- 8L   # LinDA n.cores / ANCOMBC2 n_cl
main_results_filename <- "continent_results.tsv"

method_dirs <- list(
  linda    = file.path(results_root, "continent_linda_ova"),
  maaslin2 = file.path(results_root, "continent_maaslin2_ova"),
  ancombc2 = file.path(results_root, "continent_ancombc2_ova")
)
for (m in methods_to_run)
  if (!dir.exists(method_dirs[[m]])) dir.create(method_dirs[[m]], recursive = TRUE)

cat("==========================================\n")
cat("Three-method continent DA (OVA: each continent vs Other, saliva only) — da_v2 tsv\n")
cat("Scope: ", scope_tag, "\n", sep = "")
cat("Methods: ", paste(methods_to_run, collapse = ", "), "\n", sep = "")
cat("Start: ", format(Sys.time()), "\n")
cat("==========================================\n\n")

# ==========================================================================
# 4. Shared per-comparison preparation (covariates, formula, random effect)
# ==========================================================================
prep_ova <- function(site, target, other_groups) {
  comparison_name <- paste0(sanitize(target), "_enrich")
  # samples of the target site + all other VALID sites (no nesting: site IS
  # the main grouping variable)
  keep <- metadata[[main_col]] %in% c(target, other_groups)
  pair_meta  <- metadata[keep, , drop = FALSE]
  pair_counts <- genus_counts[, pair_meta[[COL_SAMPLE]], drop = FALSE]
  # recode the main effect: reference = "Other" (pooled other valid groups),
  # target = the focal age group -> coef of target level > 0 = enriched there
  pair_meta[[main_col]] <- factor(
    ifelse(as.character(pair_meta[[main_col]]) == target, target, "Other"),
    levels = c("Other", target))
  stopifnot(identical(colnames(pair_counts), pair_meta[[COL_SAMPLE]]))

  prep <- prepare_covars(pair_meta, main_col, covar_candidates)
  pair_meta    <- prep$meta
  valid_covars <- prep$covars

  fixed_terms <- if (length(valid_covars) > 0)
    paste(c(main_col, valid_covars), collapse = " + ") else main_col
  re <- add_random_effect(paste0("~ ", fixed_terms), pair_meta, main_col)
  pair_meta[[COL_PROJECT]] <- droplevels(factor(as.character(pair_meta[[COL_PROJECT]])))

  for (cov in c(main_col, valid_covars)) {
    if (is.factor(pair_meta[[cov]])) pair_meta[[cov]] <- droplevels(pair_meta[[cov]])
    if (length(unique(as.character(pair_meta[[cov]]))) < 2)
      stop("Factor '", cov, "' has <2 levels for ", comparison_name)
  }

  # keep the fixed design full rank (LinDA crashes on aliased columns)
  valid_covars <- drop_aliased_covars(pair_meta, main_col, valid_covars)
  fixed_terms <- if (length(valid_covars) > 0)
    paste(c(main_col, valid_covars), collapse = " + ") else main_col
  re <- add_random_effect(paste0("~ ", fixed_terms), pair_meta, main_col)

  list(site = site, g1 = "Other", g2 = target,
       other_groups = other_groups,
       comparison_name = comparison_name,
       meta = pair_meta, counts = pair_counts,
       valid_covars = valid_covars, fixed_terms = fixed_terms, re = re)
}

# ---- per-pair summary/config writers (shared by all methods) ---------------
make_summary_row <- function(pp, n1, n2, formula_label) {
  data.frame(
    Continent1 = pp$g1, Continent2 = pp$g2, Comparison = pp$comparison_name,
    N1 = n1, N2 = n2, N_Significant = NA_integer_,
    N_Enriched_Continent1 = NA_integer_, N_Enriched_Continent2 = NA_integer_,
    Formula = formula_label, stringsAsFactors = FALSE
  )
}
fill_summary <- function(summary_row, output_subdir, res) {
  summary_row$N_Significant <- res$n_sig
  summary_row$N_Enriched_Continent1 <- res$n_up_g1
  summary_row$N_Enriched_Continent2 <- res$n_up_g2
  fwrite(summary_row, file.path(output_subdir, "summary.tsv"), sep = "\t", quote = FALSE)
}

# ==========================================================================
# 5a. LinDA
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
  formula_str <- pp$re$full_formula

  cat("    Formula: ", formula_str, "\n", sep = "")
  cat("    Samples: ", pp$g1, " (", n1, ") vs ", pp$g2, " (", n2, ")\n", sep = "")
  cat("    LinDA: feature.dat.type=count, prev.filter=", PREV_CUT,
      ", is.winsor=FALSE, n.cores=", n_cores, "\n", sep = "")

  if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)

  config_df <- data.frame(
    pair_index = pair_index, comparison = pp$comparison_name,
    target = pp$g2, n1 = n1, n2 = n2,
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

    group_coef <- resolve_group_coef(names(fit$output), main_col, pp$g2)
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
# 5b. Maaslin2
# ==========================================================================
run_pair_maaslin2 <- function(pair_index, pp) {
  output_subdir <- file.path(method_dirs$maaslin2, pp$comparison_name)
  if (file.exists(file.path(output_subdir, "analysis.done"))) return("skip")

  pair_meta <- pp$meta
  # relabund derived in-memory per pair (cols sum to 1); Maaslin2's TSS
  # normalization becomes a no-op.
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
    target = pp$g2, n1 = n1, n2 = n2,
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

    for (d in c("features", "fits")) unlink(file.path(output_subdir, d), recursive = TRUE)
    unlink(file.path(output_subdir, c("maaslin2.log", "significant_results.tsv")))

    cat("    SUCCESS!\n")
    TRUE
  }, error = function(e) { cat("    ERROR: ", e$message, "\n", sep = ""); FALSE })

  if (isTRUE(ok)) "ok" else "fail"
}

# ==========================================================================
# 5c. ANCOMBC2 (with collinearity retry)
# ==========================================================================
run_pair_ancombc2 <- function(pair_index, pp) {
  output_subdir <- file.path(method_dirs$ancombc2, pp$comparison_name)
  if (file.exists(file.path(output_subdir, "analysis.done"))) return("skip")

  pair_meta  <- pp$meta
  # da_v4 genus-fittability pre-filter (>=2 BioProjects AND >=2 groups)
  pair_counts <- filter_fittable(pp$counts, pair_meta, main_col)

  valid_covars <- pp$valid_covars
  rand_formula <- if (pp$re$has_random) paste0("(1|", COL_PROJECT, ")") else NULL

  n1 <- sum(pair_meta[[main_col]] == pp$g1)
  n2 <- sum(pair_meta[[main_col]] == pp$g2)
  cat("    Formula: ", pp$re$full_formula, "\n", sep = "")
  cat("    Samples: ", pp$g1, " (", n1, ") vs ", pp$g2, " (", n2, ")\n", sep = "")

  if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)

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
    group_coef <- resolve_group_coef(names(res_wide), main_col, pp$g2)
    common_df <- ancombc2_to_common(res_wide, group_coef)
    all_res <- harmonize_results(common_df, group1 = pp$g1, group2 = pp$g2,
                                 comparison = pp$comparison_name,
                                 target_site = pp$site,
                                 group1_col = g1_col, group2_col = g2_col)
    config_df <- data.frame(
      pair_index = pair_index, comparison = pp$comparison_name,
      target = pp$g2, n1 = n1, n2 = n2,
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

# ==========================================================================
# 6. Run all (method x pair)
# ==========================================================================
RUNNERS <- list(linda = run_pair_linda, maaslin2 = run_pair_maaslin2,
                ancombc2 = run_pair_ancombc2)

cat("\n[3/3] Running DA for all selected OVA comparisons...\n")
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

    res <- RUNNERS[[m]](i, prep_ova(pp$site, pp$g2, pp$other_groups))
    if (res == "ok") { n_success <- n_success + 1L; cat("    Output: ",
        file.path(method_dirs[[m]], comparison_name), "\n", sep = "") }
    else if (res == "skip") { cat("    Already completed. Skipping.\n"); n_skip <- n_skip + 1L }
    else n_fail <- n_fail + 1L
    invisible(gc(verbose = FALSE))
  }
  cat("\n", toupper(m), " OVA finished\n", sep = "")
  cat("  Success: ", n_success, "\n  Skipped: ", n_skip, "\n  Failed: ", n_fail, "\n", sep = "")
  if (n_fail > 0) any_fail <- TRUE
}

cat("\n==========================================\n")
cat("All requested methods finished (scope: ", scope, ")\n", sep = "")
cat("End: ", format(Sys.time()), "\n")
cat("==========================================\n")
if (any_fail) quit(save = "no", status = 1)