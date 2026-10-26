# Comparison preparation and shared LinDA, MaAsLin2 and ANCOM-BC2 fits.

prep_ova <- function(site, target, other_groups) {
  comparison_name <- if (analysis == "age") {
    paste0(sanitize(site), "_enrich_", sanitize_age_level(target))
  } else {
    paste0(sanitize(target), "_enrich")
  }
  keep <- metadata[[main_col]] %in% c(target, other_groups)
  if (analysis == "age") keep <- keep & metadata[[COL_SITE]] == site
  pair_meta <- metadata[keep, , drop = FALSE]
  pair_counts <- genus_counts[, pair_meta[[COL_SAMPLE]], drop = FALSE]
  pair_meta[[main_col]] <- factor(ifelse(as.character(pair_meta[[main_col]]) == target, target, "Other"),
    levels = c("Other", target)
  )
  stopifnot(identical(colnames(pair_counts), pair_meta[[COL_SAMPLE]]))
  prep <- prepare_covars(pair_meta, main_col, covar_candidates)
  pair_meta <- prep$meta
  valid_covars <- prep$covars
  fixed_terms <- if (length(valid_covars) > 0) {
    paste(c(main_col, valid_covars), collapse = " + ")
  } else {
    main_col
  }
  re <- add_random_effect(paste0("~ ", fixed_terms), pair_meta, main_col)
  pair_meta[[COL_PROJECT]] <- droplevels(factor(as.character(pair_meta[[COL_PROJECT]])))
  for (cov in c(main_col, valid_covars)) {
    if (is.factor(pair_meta[[cov]])) {
      pair_meta[[cov]] <- droplevels(pair_meta[[cov]])
    }
    if (length(unique(as.character(pair_meta[[cov]]))) < 2) {
      stop("Factor '", cov, "' has <2 levels for ", comparison_name)
    }
  }
  valid_covars <- drop_aliased_covars(pair_meta, main_col, valid_covars)
  fixed_terms <- if (length(valid_covars) > 0) {
    paste(c(main_col, valid_covars), collapse = " + ")
  } else {
    main_col
  }
  re <- add_random_effect(paste0("~ ", fixed_terms), pair_meta, main_col)
  list(
    site = site, g1 = "Other", g2 = target, other_groups = other_groups, comparison_name = comparison_name,
    meta = pair_meta, counts = pair_counts, valid_covars = valid_covars, fixed_terms = fixed_terms,
    re = re
  )
}

prep_pair <- function(continent, g1, g2) {
  comparison_name <- paste0(sanitize(SITE_SALIVA), "_", sanitize(continent), "_", g1, "_vs_", g2)
  keep <- metadata[[COL_SITE]] == SITE_SALIVA & metadata[[COL_SEX]] %in% c(g1, g2)
  pair_meta <- metadata[keep, , drop = FALSE]
  pair_counts <- genus_counts[, pair_meta[[COL_SAMPLE]], drop = FALSE]
  pair_meta[[main_col]] <- factor(as.character(pair_meta[[main_col]]), levels = c(g1, g2))
  stopifnot(identical(colnames(pair_counts), pair_meta[[COL_SAMPLE]]))
  prep <- prepare_covars(pair_meta, main_col, covar_candidates)
  pair_meta <- prep$meta
  valid_covars <- prep$covars
  fixed_terms <- if (length(valid_covars) > 0) {
    paste(c(main_col, valid_covars), collapse = " + ")
  } else {
    main_col
  }
  re <- add_random_effect(paste0("~ ", fixed_terms), pair_meta, main_col)
  pair_meta[[COL_PROJECT]] <- droplevels(factor(as.character(pair_meta[[COL_PROJECT]])))
  for (cov in c(main_col, valid_covars)) {
    if (is.factor(pair_meta[[cov]])) {
      pair_meta[[cov]] <- droplevels(pair_meta[[cov]])
    }
    if (length(unique(as.character(pair_meta[[cov]]))) < 2) {
      stop("Factor '", cov, "' has <2 levels for ", comparison_name)
    }
  }
  list(
    site = SITE_SALIVA, continent = continent, g1 = g1, g2 = g2, comparison_name = comparison_name,
    meta = pair_meta, counts = pair_counts, valid_covars = valid_covars, fixed_terms = fixed_terms,
    re = re
  )
}

comparison_identity <- function(pp) {
  switch(analysis,
    age = list(site = pp$site, age1 = pp$g1, age2 = pp$g2),
    sex = list(site = pp$site, continent = pp$continent, sex1 = pp$g1, sex2 = pp$g2),
    list(target = pp$g2)
  )
}

make_summary_row <- function(pp, n1, n2, formula_label) {
  identity <- setNames(list(pp$g1, pp$g2), c(g1_col, g2_col))
  if (analysis %in% c("age", "sex")) identity$Site <- pp$site
  if (analysis == "sex") identity$Continent <- pp$continent
  data.frame(identity,
    Comparison = pp$comparison_name, N1 = n1, N2 = n2,
    N_Significant = NA_integer_,
    setNames(list(NA_integer_, NA_integer_), paste0("N_Enriched_", c(g1_col, g2_col))),
    Formula = formula_label, stringsAsFactors = FALSE
  )
}
fill_summary <- function(summary_row, output_subdir, res) {
  summary_row$N_Significant <- res$n_sig
  summary_row[[paste0("N_Enriched_", g1_col)]] <- res$n_up_g1
  summary_row[[paste0("N_Enriched_", g2_col)]] <- res$n_up_g2
  fwrite(summary_row, file.path(output_subdir, "summary.tsv"), sep = "\t", quote = FALSE)
}

resolve_fit_coef <- function(available, group_col, target) {
  if (analysis == "sex") {
    ancombc_group_coef(group_col, target)
  } else {
    resolve_group_coef(available, group_col, target)
  }
}

run_pair_linda <- function(pair_index, pp) {
  output_subdir <- file.path(method_dirs$linda, pp$comparison_name)
  if (file.exists(file.path(output_subdir, "analysis.done"))) {
    return("skip")
  }
  pair_meta <- pp$meta
  pair_genus <- pp$counts
  lib <- colSums(pair_genus)
  keep_s <- lib > 0
  pair_genus <- pair_genus[, keep_s, drop = FALSE]
  pair_meta <- pair_meta[keep_s, , drop = FALSE]
  keep_f <- rowSums(pair_genus) > 0
  pair_genus <- pair_genus[keep_f, , drop = FALSE]
  n1 <- sum(pair_meta[[main_col]] == pp$g1)
  n2 <- sum(pair_meta[[main_col]] == pp$g2)
  formula_str <- pp$re$full_formula
  cat("    Formula: ", formula_str, "\n", sep = "")
  cat("    Samples: ", pp$g1, " (", n1, ") vs ", pp$g2, " (", n2, ")\n", sep = "")
  cat("    LinDA: feature.dat.type=count, prev.filter=", PREV_CUT, ", is.winsor=FALSE, n.cores=", n_cores,
    "\n",
    sep = ""
  )
  if (!dir.exists(output_subdir)) {
    dir.create(output_subdir, recursive = TRUE)
  }
  config_df <- data.frame(
    pair_index = pair_index, comparison = pp$comparison_name, comparison_identity(pp),
    n1 = n1, n2 = n2, formula = formula_str, feature_dat_type = "count", prev_filter = PREV_CUT,
    is_winsor = FALSE, outlier_pct = 0.03, p_adj_method = "BH", alpha = DA_ALPHA, n_cores = n_cores,
    counts_tsv = COUNTS_TSV, metadata_file = METADATA_FILE, scope = scope, stringsAsFactors = FALSE
  )
  ok <- tryCatch(
    {
      fit <- MicrobiomeStat::linda(
        feature.dat = pair_genus, meta.dat = pair_meta, formula = formula_str,
        feature.dat.type = "count", prev.filter = PREV_CUT, mean.abund.filter = 0, max.abund.filter = 0,
        is.winsor = FALSE, outlier.pct = 0.03, p.adj.method = "BH", alpha = DA_ALPHA, n.cores = n_cores,
        verbose = TRUE
      )
      bias_vals <- as.numeric(fit$bias)
      bias_names <- names(fit$bias)
      if (is.null(bias_names) || length(bias_names) != length(bias_vals)) {
        if (!is.null(fit$variables) && length(fit$variables) == length(bias_vals)) {
          bias_names <- fit$variables
        } else if (!is.null(names(fit$output)) && length(names(fit$output)) == length(bias_vals)) {
          bias_names <- names(fit$output)
        } else {
          bias_names <- paste0("V", seq_along(bias_vals))
        }
      }
      fwrite(data.frame(variable = bias_names, bias = bias_vals, stringsAsFactors = FALSE), file.path(
        output_subdir,
        "bias.tsv"
      ), sep = "\t")
      group_coef <- resolve_fit_coef(names(fit$output), main_col, pp$g2)
      common_df <- linda_to_common(fit, group_coef)
      all_res <- harmonize_results(common_df,
        group1 = pp$g1, group2 = pp$g2, comparison = pp$comparison_name,
        target_site = pp$site, group1_col = g1_col, group2_col = g2_col
      )
      summary_row <- make_summary_row(pp, n1, n2, formula_str)
      res <- write_comparison_outputs(all_res, output_subdir, pp$comparison_name,
        group1 = pp$g1, group2 = pp$g2,
        group1_col = g1_col, group2_col = g2_col, main_results_var_pattern = paste0("^", main_col),
        summary_row = summary_row, run_config_row = config_df, main_results_filename = main_results_filename
      )
      fill_summary(summary_row, output_subdir, res)
      cat("    SUCCESS!\n")
      TRUE
    },
    error = function(e) {
      cat("    ERROR: ", e$message, "\n", sep = "")
      FALSE
    }
  )
  if (isTRUE(ok)) {
    "ok"
  } else {
    "fail"
  }
}

run_pair_maaslin2 <- function(pair_index, pp) {
  output_subdir <- file.path(method_dirs$maaslin2, pp$comparison_name)
  if (file.exists(file.path(output_subdir, "analysis.done"))) {
    return("skip")
  }
  pair_meta <- pp$meta
  pair_genus <- sweep(pp$counts, 2, colSums(pp$counts), "/")
  valid_covars <- pp$valid_covars
  fixed_effects <- c(main_col, valid_covars)
  random_effects <- if (pp$re$has_random) {
    COL_PROJECT
  } else {
    NULL
  }
  n1 <- sum(pair_meta[[main_col]] == pp$g1)
  n2 <- sum(pair_meta[[main_col]] == pp$g2)
  cat("    Formula: ", pp$re$full_formula, "\n", sep = "")
  cat("    Samples: ", pp$g1, " (", n1, ") vs ", pp$g2, " (", n2, ")\n", sep = "")
  cat("    Fixed: ", paste(fixed_effects, collapse = ", "), " | Random: ", ifelse(pp$re$has_random,
    COL_PROJECT, "none"
  ), "\n", sep = "")
  if (!dir.exists(output_subdir)) {
    dir.create(output_subdir, recursive = TRUE)
  }
  config_df <- data.frame(
    pair_index = pair_index, comparison = pp$comparison_name, comparison_identity(pp),
    n1 = n1, n2 = n2, formula = pp$re$full_formula, fixed_effects = paste(fixed_effects, collapse = ","),
    random_effects = ifelse(pp$re$has_random, COL_PROJECT, "none"), normalization = "TSS", transform = "LOG",
    analysis_method = "LM", correction = "BH", min_abundance = MIN_ABUNDANCE, min_prevalence = PREV_CUT,
    max_significance = DA_ALPHA, standardize = FALSE, counts_tsv = COUNTS_TSV, relabund = "derived_from_counts",
    metadata_file = METADATA_FILE, scope = scope, stringsAsFactors = FALSE
  )
  ok <- tryCatch(
    {
      fit <- Maaslin2::Maaslin2(
        input_data = t(pair_genus), input_metadata = pair_meta[, c(
          COL_SAMPLE,
          main_col, valid_covars, if (pp$re$has_random) COL_PROJECT else NULL
        ), drop = FALSE], output = output_subdir,
        fixed_effects = fixed_effects, random_effects = random_effects, normalization = "TSS", transform = "LOG",
        analysis_method = "LM", correction = "BH", min_abundance = MIN_ABUNDANCE, min_prevalence = PREV_CUT,
        max_significance = DA_ALPHA, standardize = FALSE, plot_heatmap = FALSE, plot_scatter = FALSE,
        cores = 1L
      )
      maaslin_all <- read.delim(file.path(output_subdir, "all_results.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
      common_df <- maaslin2_to_common(maaslin_all, main_col, pp$g2)
      all_res <- harmonize_results(common_df,
        group1 = pp$g1, group2 = pp$g2, comparison = pp$comparison_name,
        target_site = pp$site, group1_col = g1_col, group2_col = g2_col
      )
      summary_row <- make_summary_row(pp, n1, n2, pp$re$full_formula)
      res <- write_comparison_outputs(all_res, output_subdir, pp$comparison_name,
        group1 = pp$g1, group2 = pp$g2,
        group1_col = g1_col, group2_col = g2_col, main_results_var_pattern = paste0("^", main_col),
        summary_row = summary_row, run_config_row = config_df, main_results_filename = main_results_filename
      )
      fill_summary(summary_row, output_subdir, res)
      for (d in c("features", "fits")) unlink(file.path(output_subdir, d), recursive = TRUE)
      unlink(file.path(output_subdir, c("maaslin2.log", "significant_results.tsv", "all_results.tsv")))
      cat("    SUCCESS!\n")
      TRUE
    },
    error = function(e) {
      cat("    ERROR: ", e$message, "\n", sep = "")
      FALSE
    }
  )
  if (isTRUE(ok)) {
    "ok"
  } else {
    "fail"
  }
}

run_pair_ancombc2 <- function(pair_index, pp) {
  output_subdir <- file.path(method_dirs$ancombc2, pp$comparison_name)
  if (file.exists(file.path(output_subdir, "analysis.done"))) {
    return("skip")
  }
  pair_meta <- pp$meta
  pair_counts <- filter_fittable(pp$counts, pair_meta, main_col)
  valid_covars <- pp$valid_covars
  rand_formula <- if (pp$re$has_random) {
    paste0("(1|", COL_PROJECT, ")")
  } else {
    NULL
  }
  n1 <- sum(pair_meta[[main_col]] == pp$g1)
  n2 <- sum(pair_meta[[main_col]] == pp$g2)
  cat("    Formula: ", pp$re$full_formula, "\n", sep = "")
  cat("    Samples: ", pp$g1, " (", n1, ") vs ", pp$g2, " (", n2, ")\n", sep = "")
  if (!dir.exists(output_subdir)) {
    dir.create(output_subdir, recursive = TRUE)
  }
  fit_once <- function(covars, msg_tag = "fit") {
    ff <- if (length(covars) > 0) {
      paste0(main_col, " + ", paste(covars, collapse = " + "))
    } else {
      main_col
    }
    keep_cols <- c(main_col, covars, if (pp$re$has_random) COL_PROJECT else NULL)
    ps <- phyloseq(otu_table(pair_counts, taxa_are_rows = TRUE), sample_data(pair_meta[, keep_cols,
      drop = FALSE
    ]))
    t0 <- Sys.time()
    fit <- ANCOMBC::ancombc2(
      data = ps, assay.type = "counts", fix_formula = ff, rand_formula = rand_formula,
      p_adj_method = "BH", prv_cut = PREV_CUT, lib_cut = LIB_CUT, group = NULL, struc_zero = FALSE,
      neg_lb = FALSE, alpha = DA_ALPHA, n_cl = n_cores, global = FALSE, pairwise = FALSE, verbose = FALSE
    )
    cat("    ANCOMBC2 elapsed [", msg_tag, "]: ", round(
      as.numeric(Sys.time() - t0, units = "secs"),
      1
    ), "s\n", sep = "")
    res_wide <- as.data.frame(fit$res)
    group_coef <- resolve_fit_coef(names(res_wide), main_col, pp$g2)
    common_df <- ancombc2_to_common(res_wide, group_coef)
    all_res <- harmonize_results(common_df,
      group1 = pp$g1, group2 = pp$g2, comparison = pp$comparison_name,
      target_site = pp$site, group1_col = g1_col, group2_col = g2_col
    )
    config_df <- data.frame(
      pair_index = pair_index, comparison = pp$comparison_name, comparison_identity(pp),
      n1 = n1, n2 = n2, fix_formula = ff, rand_formula = ifelse(pp$re$has_random, rand_formula,
        "none"
      ), p_adj_method = "BH", prv_cut = PREV_CUT, lib_cut = LIB_CUT, group_coef = group_coef,
      struc_zero = FALSE, neg_lb = FALSE, alpha = DA_ALPHA, n_cl = n_cores, counts_tsv = COUNTS_TSV,
      metadata_file = METADATA_FILE, scope = scope, stringsAsFactors = FALSE
    )
    summary_row <- make_summary_row(pp, n1, n2, paste0(ff, if (pp$re$has_random) {
      paste0(" + ", rand_formula)
    } else {
      ""
    }))
    res <- write_comparison_outputs(all_res, output_subdir, pp$comparison_name,
      group1 = pp$g1, group2 = pp$g2,
      group1_col = g1_col, group2_col = g2_col, main_results_var_pattern = paste0("^", main_col),
      summary_row = summary_row, run_config_row = config_df, main_results_filename = main_results_filename
    )
    fill_summary(summary_row, output_subdir, res)
    cat("    SUCCESS!\n")
    TRUE
  }
  ok <- FALSE
  cur_covars <- valid_covars
  for (attempt in 1:4) {
    err_obj <- tryCatch(
      {
        fit_once(cur_covars, msg_tag = paste0("attempt", attempt))
        NULL
      },
      error = function(e) e
    )
    if (is.null(err_obj)) {
      ok <- TRUE
      break
    }
    cat("    ERROR: ", conditionMessage(err_obj), "\n", sep = "")
    drop <- collinear_covars_from_error(conditionMessage(err_obj), cur_covars)
    if (length(drop) == 0L) {
      ok <- FALSE
      break
    }
    cat("    Collinearity detected; dropping covariate(s): ", paste(drop, collapse = ", "), " and retrying.\n",
      sep = ""
    )
    cur_covars <- setdiff(cur_covars, drop)
    ok <- FALSE
  }
  if (!isTRUE(ok)) {
    err_obj <- tryCatch(
      {
        cat("    Final attempt: main effect only (no covariates).\n")
        fit_once(character(0), msg_tag = "mainonly")
        NULL
      },
      error = function(e) e
    )
    ok <- is.null(err_obj)
    if (!ok) {
      cat("    ERROR: ", conditionMessage(err_obj), "\n", sep = "")
    }
  }
  if (isTRUE(ok)) {
    "ok"
  } else {
    "fail"
  }
}
