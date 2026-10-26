# Load each analysis into its own environment so settings never leak across runs.
da_profile <- function(analysis, project_root,
                       result_root = Sys.getenv("DA_RESULTS_ROOT", "")) {
  stopifnot(analysis %in% c("age", "continent", "site", "sex"))
  project_root <- normalizePath(project_root, mustWork = TRUE)
  script_dir <- file.path(project_root, "differential_abundance/scripts")
  env <- new.env(parent = globalenv())
  sys.source(file.path(script_dir, "common.R"), envir = env)
  if (analysis != "sex") sys.source(file.path(script_dir, "ova_helpers.R"), envir = env)
  sys.source(file.path(script_dir, "profiles", paste0(analysis, ".R")), envir = env)
  env$analysis <- analysis
  env$script_dir <- script_dir
  env$scope <- if (analysis %in% c("age", "sex")) "hc" else "all"
  env$DATA_DIR <- file.path(project_root, "data/metadata-and-biom-v2-qc")
  env$COUNTS_TSV <- file.path(env$DATA_DIR, "genus_table_counts.tsv")
  env$RELABUND_TSV <- file.path(env$DATA_DIR, "genus_table_relabund.tsv")
  env$METADATA_FILE <- file.path(env$DATA_DIR, "metadata_oral_qc_dedup.tsv")
  if (!nzchar(result_root)) result_root <- file.path(project_root, "differential_abundance/result")
  env$result_root <- result_root
  branch <- switch(analysis,
    age = "age_hc_ova",
    continent = "continent_ova",
    site = "site_ova",
    sex = "sex_pooled_hc"
  )
  env$results_root <- file.path(result_root, branch)
  prefix <- if (analysis == "sex") "sex_pooled" else analysis
  env$out_dir <- file.path(env$results_root, paste0(prefix, "_da_consensus"))
  figure_subdir <- switch(analysis,
    age = "figures_saliva_target_top10_globalpeak_maaslin2coef",
    sex = "figures",
    "figures_target_top10_globalpeak_maaslin2coef"
  )
  env$figures_dir <- file.path(env$out_dir, figure_subdir)
  env$figure_stem <- switch(analysis,
    age = "figureB_abundance_saliva_2method_combined_zscore",
    sex = "figure_lefse_maaslin2_coef",
    "figureB_abundance_2method_combined"
  )
  env$main_col <- env$group_col <- switch(analysis,
    age = env$COL_AGE,
    continent = env$COL_CONTINENT,
    site = env$COL_SITE,
    sex = env$COL_SEX
  )
  group_label <- switch(analysis,
    age = "Age",
    continent = "Continent",
    site = "Site",
    sex = "Sex"
  )
  env$g1_col <- paste0(group_label, "1")
  env$g2_col <- paste0(group_label, "2")
  env$main_results_filename <- env$main_res_file <- paste0(analysis, "_results.tsv")
  suffix <- if (analysis == "sex") "pairwise" else "ova"
  methods <- c("linda", "maaslin2", "ancombc2")
  env$method_dirs <- setNames(lapply(methods, function(m) {
    file.path(env$results_root, paste(analysis, m, suffix, sep = "_"))
  }), methods)
  env$covar_candidates <- switch(analysis,
    age = c(env$COL_SEX, env$COL_PLATFORM, env$COL_SEQ_REGION, env$COL_BMI, env$COL_CONTINENT),
    continent = c(env$COL_AGE, env$COL_SEX, env$COL_PLATFORM, env$COL_SEQ_REGION, env$COL_BMI),
    site = c(env$COL_AGE, env$COL_SEX, env$COL_PLATFORM, env$COL_SEQ_REGION, env$COL_BMI, env$COL_CONTINENT),
    sex = c(env$COL_CONTINENT, env$COL_AGE, env$COL_PLATFORM, env$COL_SEQ_REGION, env$COL_BMI)
  )
  env$n_cores <- 8L
  env$alpha <- env$DA_ALPHA
  env$scope_tag <- paste(analysis, env$scope)
  env$group_order <- switch(analysis,
    age = env$AGE_ORDER,
    continent = env$CONTINENT_ORDER,
    site = env$SITE_ORDER,
    sex = c("F", "M")
  )
  env$group_name <- switch(analysis,
    age = "Age group",
    continent = "Continent",
    site = "Site",
    sex = "Sex"
  )
  env
}
