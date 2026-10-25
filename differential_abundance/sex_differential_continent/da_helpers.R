# ==============================================================================
# da_helpers.R  (sex-within-continent da_v2 — adapted from
#                sex_differential/da_helpers.R)
# Shared preprocessing + harmonization helpers for the three-method
# differential abundance (DA) analysis of SEX: LinDA + Maaslin2 + ANCOMBC2.
# Main effect = sex_cleaned (F vs M), run PER CONTINENT within the SALIVA site
# only (one F_vs_M comparison per valid continent). Covariates exclude sex
# (main), site (fixed = saliva) AND continent (stratum, no within-stratum
# variance):
#   age_cleaned_group, platform_cleaned, sequencing_region_cleaned,
#   bmi_cleaned_group.
#
# DESIGN (identical to the sex/site da_v2 scripts except where noted):#   - Data loaded from plain .tsv via read.table (ONE table per driver), NOT
#     .qza via qiime2R (sparse->dense coercion OOM risk).
#   - reads>=1000 sample filter is applied at TSV EXPORT time (Python,
#     prepare_genus_tables.py); R drivers only do lib>0 / rowSums>0 here.
#   - 1% feature prevalence stays INSIDE each tool (prev_cut / min_prevalence /
#     prv_cut = 0.01), per comparison.
#   - age / BMI covariates: pre-binned *_group metadata columns.
#   - Conditional (1|bioproject_id) random effect when a BioProject spans
#     BOTH sexes within a continent (saliva) — the key guard against
#     study/batch confounding.
#   - SCOPE: --scope=all (all samples with sex) or --scope=hc (status_cleaned
#     == "HC" only). Results root: rd("sex_continent") vs rd("sex_continent_hc").
#
# Locked constants match the geography/site da_v2 scripts verbatim:
#   seed 20260615; alpha 0.05; prev_cut 0.01; min_samples_per_sex = 10
#   LinDA (MicrobiomeStat::linda): feature.dat.type="count", is.winsor=FALSE,
#          prev.filter=0.01, n.cores=8, BH, alpha=0.05
#   Maaslin2: TSS+LOG+LM, BH, min_prevalence=0.01, standardize=FALSE, cores=1
#   ANCOMBC2: fix_formula binary contrast, (1|bioproject_id), prv_cut=0.01,
#             lib_cut=0, struc_zero=FALSE, neg_lb=FALSE, n_cl=8
#
# Contrast convention (used consistently in outputs):
#   reference = F, target = M  ->  log2FoldChange > 0  =>  enriched in MALES
#   comparison name = "saliva_<continent>_F_vs_M" (e.g. saliva_North_America_F_vs_M)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

# ---- Locked constants --------------------------------------------------------
DA_SEED   <- 20260615
DA_ALPHA  <- 0.05

# ---- INPUT FILTERS ----------------------------------------------------------
PREV_CUT                 <- 0.01   # feature prevalence (per-comparison): LinDA prev.filter, Maaslin2 min_prevalence, ANCOMBC2 prv_cut
MIN_ABUNDANCE            <- 0      # feature abundance floor (Maaslin2 min_abundance); 0 = no floor
LIB_CUT                  <- 0      # library-size floor (ANCOMBC2 lib_cut, LinDA lib cut); 0 = no floor
MIN_SAMPLES_PER_SEX      <- 10L    # min samples PER SEX within a site to keep that site's F_vs_M comparison
DROP_BIOPROJECTS         <- NULL   # optional bioproject_id vector to exclude (NULL = off)
DROP_SITES               <- NULL   # optional sites to drop (NULL = off)
# <<< end input filters

COL_SAMPLE     <- "sample_id"
COL_CONTINENT  <- "continent_group"
COL_SITE       <- "site_cleaned"
COL_SITE_GROUP <- "site_group_cleaned"
COL_AGE        <- "age_cleaned_group"   # pre-binned (da_v4 aligned)
COL_SEX        <- "sex_cleaned"
COL_STATUS     <- "status_cleaned"
COL_BMI        <- "bmi_cleaned_group"   # pre-binned (da_v4 aligned)
COL_PROJECT    <- "bioproject_id"
COL_PLATFORM   <- "platform_cleaned"
COL_SEQ_REGION <- "sequencing_region_cleaned"

SEX_REF  <- "F"   # reference level of the contrast
SEX_TGT  <- "M"   # target level of the contrast (coef > 0 => enriched in M)

SITE_ORDER <- c(
  "saliva", "tooth", "oral mucosa", "oral swab", "tongue", "gingiva"
)

# Fixed site for this whole pipeline (sex within continent, saliva only).
SITE_SALIVA <- "saliva"

# Continents eligible for the sex-within-continent contrast (saliva, HC).
# Africa is EXCLUDED: it has 0 male saliva-HC samples, so the both-sexes gate
# in valid_sex_continents() would drop it anyway; omitting it here keeps
# logs/plots clean. Order matches the geography da_v2 CONTINENT_ORDER minus
# Africa.
CONTINENT_ORDER <- c(
  "Asia", "Europe", "North America", "Oceania", "South America"
)

# data paths (project-root relative)
DATA_DIR     <- "data/metadata-and-biom-v2-qc"
COUNTS_TSV   <- file.path(DATA_DIR, "genus_table_counts.tsv")
RELABUND_TSV <- file.path(DATA_DIR, "genus_table_relabund.tsv")
METADATA_FILE <- file.path(DATA_DIR, "metadata_oral_qc_dedup.tsv")

# harmonized result schema (must match LinDA all_results.tsv columns exactly)
COMMON_COLS <- c(
  "feature", "genus_short", "variable",
  "baseMean", "log2FoldChange", "lfcSE", "stat",
  "pvalue", "padj", "reject", "df",
  "direction", "enriched_in"
)

sanitize <- function(x) gsub(" ", "_", x)

# Resolve a results subdirectory under the optional DA_RESULTS_ROOT.
# (Same contract as the geography/site da_helpers.R.)
rd <- function(subpath) {
  root <- Sys.getenv("DA_RESULTS_ROOT", "")
  if (nchar(root) > 0) file.path(root, subpath) else file.path("differential_abundance", "result", subpath)
}

# Scope -> results subpath: "all" -> rd("sex/..."), "hc" -> rd("sex_hc/...")
sex_root <- function(scope) {
  stopifnot(scope %in% c("all", "hc"))
  file.path("results", "da_v2", if (scope == "hc") "sex_hc" else "sex")
}

# Scope -> results subpath for the sex-within-continent pipeline:
#   "all" -> rd("sex_continent/..."), "hc" -> rd("sex_continent_hc/...")
continent_root <- function(scope) {
  stopifnot(scope %in% c("all", "hc"))
  file.path("results", "da_v2", if (scope == "hc") "sex_continent_hc" else "sex_continent")
}

# Parse the --scope CLI arg (default "all"); stops on unknown values.
parse_scope_arg <- function(args) {
  scope <- "all"
  for (a in args) if (grepl("^--scope=", a)) scope <- sub("^--scope=", "", a)
  if (!scope %in% c("all", "hc"))
    stop("Unknown --scope value: ", scope, " (valid: all | hc)")
  scope
}

# ---- 1. Data loading (tsv, one table per driver) ----------------------------
# Identical to the geography/site da_v2 loader: read ONE genus table from tsv
# as a dense double matrix. The reads>=1000 sample filter was already applied
# at export (prepare_genus_tables.py).
# `which`: "counts" -> integer counts (LinDA / ANCOMBC2), "relabund" ->
#   relative abundance cols sum to 1 (Maaslin2 / abundance panels).
load_genus_table <- function(which = c("counts", "relabund")) {
  which <- match.arg(which)
  path <- if (which == "counts") COUNTS_TSV else RELABUND_TSV
  if (!file.exists(path)) stop("Genus tsv not found: ", path,
                               " (run prepare_genus_tables.py / sbatch first)")
  df <- read.table(path, header = TRUE, sep = "\t", row.names = 1,
                   check.names = FALSE, comment.char = "",
                   stringsAsFactors = FALSE)
  m <- as.matrix(df)
  storage.mode(m) <- "double"
  m
}

# ---- 2. Metadata cleaning ---------------------------------------------------
# Sex-within-continent scope: fixed to the SALIVA site only (Oral group), no
# main-col filter here — each per-continent F_vs_M comparison subsets its own
# continent and the two sex levels. hc_only = TRUE restricts to
# status_cleaned == "HC". Uses the pre-binned age_cleaned_group /
# bmi_cleaned_group columns directly.
load_sex_metadata <- function(hc_only = FALSE) {
  if (!file.exists(METADATA_FILE)) stop("Metadata not found: ", METADATA_FILE)
  metadata_raw <- fread(
    METADATA_FILE,
    sep = "\t", header = TRUE,
    quote = "", fill = TRUE, check.names = FALSE,
    data.table = FALSE
  )

  required_cols <- c(COL_SAMPLE, COL_SITE, COL_SEX, COL_PROJECT)
  missing_cols <- setdiff(required_cols, names(metadata_raw))
  if (length(missing_cols) > 0) {
    stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # optional covariate columns: create if absent so pipeline does not break
  for (opt_col in c(COL_AGE, COL_BMI, COL_PLATFORM, COL_SEQ_REGION,
                    COL_CONTINENT, COL_SITE_GROUP, COL_STATUS)) {
    if (!opt_col %in% names(metadata_raw)) {
      metadata_raw[[opt_col]] <- NA
    }
  }

  metadata_clean <- metadata_raw %>%
    filter(!is.na(.data[[COL_SAMPLE]]) & .data[[COL_SAMPLE]] != "") %>%
    distinct(.data[[COL_SAMPLE]], .keep_all = TRUE) %>%
    mutate(
      !!COL_SITE        := trimws(as.character(.data[[COL_SITE]])),
      !!COL_SEX         := trimws(as.character(.data[[COL_SEX]])),
      !!COL_STATUS      := trimws(as.character(.data[[COL_STATUS]])),
      !!COL_AGE         := trimws(as.character(.data[[COL_AGE]])),
      !!COL_BMI         := trimws(as.character(.data[[COL_BMI]])),
      !!COL_CONTINENT   := trimws(as.character(.data[[COL_CONTINENT]])),
      !!COL_PLATFORM    := trimws(as.character(.data[[COL_PLATFORM]])),
      !!COL_SEQ_REGION  := trimws(as.character(.data[[COL_SEQ_REGION]])),
      !!COL_PROJECT     := trimws(as.character(.data[[COL_PROJECT]])),
      !!COL_SITE_GROUP  := trimws(as.character(.data[[COL_SITE_GROUP]]))
    ) %>%
    mutate(
      !!COL_AGE         := ifelse(.data[[COL_AGE]] %in% c("", "NA", "na"), "Unknown", .data[[COL_AGE]]),
      !!COL_BMI         := ifelse(.data[[COL_BMI]] %in% c("", "NA", "na"), "Unknown", .data[[COL_BMI]]),
      !!COL_PLATFORM    := ifelse(.data[[COL_PLATFORM]] %in% c("", "NA", "na"), "Unknown", .data[[COL_PLATFORM]]),
      !!COL_SEQ_REGION  := ifelse(.data[[COL_SEQ_REGION]] %in% c("", "NA", "na"), "Unknown", .data[[COL_SEQ_REGION]]),
      !!COL_CONTINENT   := ifelse(.data[[COL_CONTINENT]] %in% c("", "NA", "na"), "Unknown", .data[[COL_CONTINENT]])
    )

  # valid oral sites only (Oral group)
  metadata_clean <- metadata_clean %>%
    filter(
      .data[[COL_SITE]] %in% SITE_ORDER &
        !is.na(.data[[COL_SITE]]) & .data[[COL_SITE]] != ""
    )
  sg <- trimws(as.character(metadata_clean[[COL_SITE_GROUP]]))
  keep_sg <- is.na(sg) | sg %in% c("", "NA", "Oral")
  metadata_clean <- metadata_clean[keep_sg, , drop = FALSE]

  # fix to the SALIVA site (this pipeline is saliva-only)
  n_pre <- nrow(metadata_clean)
  metadata_clean <- metadata_clean[
    !is.na(metadata_clean[[COL_SITE]]) &
      metadata_clean[[COL_SITE]] == SITE_SALIVA, , drop = FALSE]
  cat("    Saliva-only: kept ", nrow(metadata_clean), " / ", n_pre,
      " oral samples (site == ", SITE_SALIVA, ")\n", sep = "")

  # HC scope
  if (hc_only) {
    n0 <- nrow(metadata_clean)
    metadata_clean <- metadata_clean[
      !is.na(metadata_clean[[COL_STATUS]]) &
        metadata_clean[[COL_STATUS]] == "HC", , drop = FALSE]
    cat("    HC scope: kept ", nrow(metadata_clean), " / ", n0,
        " samples (status_cleaned == HC)\n", sep = "")
  }

  metadata_clean <- as.data.frame(metadata_clean)
  rownames(metadata_clean) <- metadata_clean[[COL_SAMPLE]]
  metadata_clean
}

# Intersect metadata (row-named by sample_id) with a genus table matrix and
# return the common-samples subset. No reads filter here (done at tsv export).
align_samples <- function(meta, genus_mat) {
  # --- optional pre-filters (gated by INPUT FILTERS block; default off) ---
  if (!is.null(DROP_BIOPROJECTS)) {
    n0 <- nrow(meta)
    meta <- meta[!(as.character(meta[[COL_PROJECT]]) %in% DROP_BIOPROJECTS), , drop = FALSE]
    if (nrow(meta) < n0) cat("    Dropped ", n0 - nrow(meta),
                             " samples by DROP_BIOPROJECTS\n", sep = "")
  }
  if (COL_SITE %in% names(meta) && !is.null(DROP_SITES)) {
    n0 <- nrow(meta)
    meta <- meta[!(as.character(meta[[COL_SITE]]) %in% DROP_SITES), , drop = FALSE]
    if (nrow(meta) < n0) cat("    Dropped ", n0 - nrow(meta),
                             " samples by DROP_SITES\n", sep = "")
  }

  common <- intersect(rownames(meta), colnames(genus_mat))
  if (length(common) == 0) stop("No common samples between metadata and genus table!")
  meta <- meta[common, , drop = FALSE]
  mat  <- genus_mat[, common, drop = FALSE]
  list(meta = meta, mat = mat)
}

# Valid continents for the sex contrast (within saliva, already fixed by the
# loader): BOTH sexes have >= min_n samples in the continent (preserving
# CONTINENT_ORDER). Continents absent from CONTINENT_ORDER (e.g. Africa, which
# has 0 males) are dropped automatically.
valid_sex_continents <- function(meta, min_n = MIN_SAMPLES_PER_SEX) {
  keep <- meta[[COL_SEX]] %in% c(SEX_REF, SEX_TGT)
  counts <- table(factor(meta[[COL_CONTINENT]][keep], levels = CONTINENT_ORDER),
                  factor(meta[[COL_SEX]][keep], levels = c(SEX_REF, SEX_TGT)))
  ok <- rowSums(counts >= min_n) == 2
  setNames(names(ok)[ok], NULL)
}

# ---- 3. Pair generation -----------------------------------------------------
# One F_vs_M comparison per valid continent (within saliva). Each element:
# list(site = SITE_SALIVA, continent, g1, g2, comparison_name).
make_sex_continent_pairs <- function(valid_continents) {
  lapply(valid_continents, function(continent) {
    list(
      site = SITE_SALIVA,
      continent = continent,
      g1 = SEX_REF, g2 = SEX_TGT,
      comparison_name = paste0(sanitize(SITE_SALIVA), "_",
                               sanitize(continent), "_",
                               SEX_REF, "_vs_", SEX_TGT)
    )
  })
}

# ---- 4. Covariate preparation (collapse exclusive levels to "Other") -------
# Identical to the geography/site da_v2 version.
prepare_covars <- function(meta, group_col, covars, bmi_col = COL_BMI,
                           age_col = COL_AGE) {
  valid_covars <- character(0)
  meta_out <- meta

  for (cov in covars) {
    vals <- as.character(meta_out[[cov]])
    vals[is.na(vals) | vals %in% c("", "NA", "na")] <- "Unknown"

    cross_tab <- table(meta_out[[group_col]], vals)
    # levels present in both groups
    shared_levels <- colnames(cross_tab)[colSums(cross_tab > 0) == 2]
    exclusive_levels <- setdiff(colnames(cross_tab), shared_levels)

    if (length(exclusive_levels) > 0) {
      vals[vals %in% exclusive_levels] <- "Other"
      cat("    NOTE: ", cov, " exclusive levels collapsed to Other: ",
          paste(exclusive_levels, collapse = ", "), "\n", sep = "")
    }

    vals <- droplevels(factor(vals))
    tab <- sort(table(vals), decreasing = TRUE)
    if (length(tab) < 2) {
      cat("    WARNING: ", cov, " has only 1 usable level, dropping.\n", sep = "")
      next
    }

    cross2 <- table(meta_out[[group_col]], vals)
    if (sum(colSums(cross2 > 0) == 2) < 1) {
      cat("    WARNING: ", cov,
          " has no level shared by the two sexes after collapse, dropping.\n", sep = "")
      next
    }

    ref_level <- names(tab)[1]
    if (cov == bmi_col && "Normal" %in% names(tab)) ref_level <- "Normal"
    if (cov == age_col && "Young_Adult" %in% names(tab)) ref_level <- "Young_Adult"
    if (ref_level %in% c("Other", "Unknown")) {
      prefer <- setdiff(names(tab), c("Other", "Unknown"))
      if (length(prefer) > 0) ref_level <- prefer[1]
    }

    other_levels <- setdiff(names(tab), ref_level)
    meta_out[[cov]] <- factor(as.character(vals), levels = c(ref_level, other_levels))
    valid_covars <- c(valid_covars, cov)
    cat("    [OK] Kept ", cov, " (Ref: ", ref_level,
        "; levels: ", paste(levels(meta_out[[cov]]), collapse = ", "), ")\n", sep = "")
  }

  list(meta = meta_out, covars = valid_covars)
}

# ---- 5. Conditional random effect for BioProject ----------------------------
# Returns list(fixed_formula, rand_formula, full_formula, has_random, n_shared)
add_random_effect <- function(base_formula, meta, group_col,
                              project_col = COL_PROJECT) {
  meta[[project_col]] <- droplevels(factor(as.character(meta[[project_col]])))
  project_tab <- table(meta[[group_col]], meta[[project_col]])
  shared_projects <- sum(colSums(project_tab > 0) == 2)

  if (shared_projects > 0) {
    full_formula <- paste0(base_formula, " + (1|", project_col, ")")
    rand_formula <- paste0("(1|", project_col, ")")
    cat("    Using random effect for BioProject (",
        shared_projects, " projects in both groups).\n", sep = "")
    list(
      fixed_formula = base_formula,
      rand_formula  = rand_formula,
      full_formula = full_formula,
      has_random  = TRUE,
      n_shared    = shared_projects
    )
  } else {
    cat("    No BioProjects spanning both groups; fixed-effects only.\n")
    list(
      fixed_formula = base_formula,
      rand_formula  = NULL,
      full_formula  = base_formula,
      has_random   = FALSE,
      n_shared     = 0L
    )
  }
}

# ---- 6. Result harmonization ------------------------------------------------
# genus_short: strip to genus label after last 'g__'
make_genus_short <- function(feature) {
  g <- sub(".*g__", "", feature)
  sub(";.*", "", g)
}

# Given a per-feature data.frame with at least:
#   feature, coef (->log2FoldChange), se (->lfcSE), stat, pvalue, padj
# and the contrast labels (group1=reference=F, group2=target=M so coef>0 =>
# enriched in group2), produce the harmonized common-schema data.frame.
harmonize_results <- function(df, group1, group2, comparison,
                              target_site = NA_character_,
                              group1_col = "Group1", group2_col = "Group2") {
  if (nrow(df) == 0) {
    return(data.frame(
      feature = character(0), genus_short = character(0), variable = character(0),
      baseMean = numeric(0), log2FoldChange = numeric(0), lfcSE = numeric(0),
      stat = numeric(0), pvalue = numeric(0), padj = numeric(0), reject = logical(0),
      df = numeric(0), direction = character(0), enriched_in = character(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- data.frame(
    feature        = df$feature,
    genus_short    = make_genus_short(df$feature),
    variable       = if (!is.null(df$variable)) df$variable else NA_character_,
    baseMean       = if (!is.null(df$baseMean)) df$baseMean else NA_real_,
    log2FoldChange = df$log2FoldChange,
    lfcSE          = if (!is.null(df$lfcSE)) df$lfcSE else NA_real_,
    stat           = if (!is.null(df$stat)) df$stat else NA_real_,
    pvalue         = df$pvalue,
    padj           = df$padj,
    reject         = if (!is.null(df$reject)) df$reject else (df$padj < DA_ALPHA),
    df             = if (!is.null(df$df)) df$df else NA_real_,
    stringsAsFactors = FALSE
  )

  # coef sign: group2 is the contrast level -> >0 enriched in group2
  out$direction <- ifelse(
    is.na(out$log2FoldChange), NA_character_,
    ifelse(out$log2FoldChange > 0, "up_in_target",
           ifelse(out$log2FoldChange < 0, "up_in_reference", "zero"))
  )
  out$enriched_in <- ifelse(
    is.na(out$log2FoldChange), NA_character_,
    ifelse(out$log2FoldChange > 0, group2,
           ifelse(out$log2FoldChange < 0, group1, "none"))
  )

  out[[group1_col]] <- group1
  out[[group2_col]] <- group2
  out$Comparison <- comparison
  if (!is.na(target_site)) out$Target_Site <- target_site

  keep <- c(COMMON_COLS, group1_col, group2_col, "Comparison",
            if (!is.na(target_site)) "Target_Site" else NULL)
  out[, keep, drop = FALSE]
}

# ---- 6b. Method-specific converters to intermediate common format ------------

# MicrobiomeStat::linda returns $output = named list of per-coefficient
# data.frames; group_coef is the coefficient for the main contrast
# ("<group_col><group2>", e.g. "sex_cleanedM").
linda_to_common <- function(linda_fit, group_coef) {
  if (is.null(linda_fit$output) || !group_coef %in% names(linda_fit$output)) {
    stop("LinDA coefficient not in output: ", group_coef,
         " (have: ", paste(names(linda_fit$output), collapse=", "), ")")
  }
  res <- linda_fit$output[[group_coef]]
  data.frame(
    feature        = rownames(res),
    variable       = group_coef,
    baseMean       = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    lfcSE          = res$lfcSE,
    stat           = res$stat,
    pvalue         = res$pvalue,
    padj           = res$padj,
    reject         = res$reject,
    df             = res$df,
    stringsAsFactors = FALSE
  )
}

# Maaslin2 all_results.tsv rows for the group contrast: metadata == group_col
# and value == group2. Restore original taxonomy feature ids after Maaslin2's
# make.names() mangling (cached reverse map from the counts tsv).
.maalind2_rev_map_cache <- NULL
restore_feature_names <- function(feature_vec) {
  if (is.null(.maalind2_rev_map_cache)) {
    feats <- fread(COUNTS_TSV, sep = "\t", header = TRUE, select = "V1",
                   showProgress = FALSE)$V1
    mk <- make.names(feats)
    if (anyDuplicated(mk))
      stop("make.names not 1:1 on feature ids; cannot restore safely.")
    .maalind2_rev_map_cache <<- setNames(feats, mk)
  }
  m <- .maalind2_rev_map_cache
  known <- feature_vec %in% names(m)
  ifelse(known, m[feature_vec], feature_vec)
}

maaslin2_to_common <- function(maaslin_all, group_col, group2) {
  grp_rows <- maaslin_all[
    maaslin_all$metadata == group_col &
    as.character(maaslin_all$value) == group2, , drop = FALSE
  ]
  if (nrow(grp_rows) == 0) {
    grp_rows <- maaslin_all[maaslin_all$metadata == group_col, , drop = FALSE]
  }
  feat_fixed <- restore_feature_names(as.character(grp_rows$feature))
  data.frame(
    feature        = feat_fixed,
    variable       = paste0(group_col, group2),
    baseMean       = NA_real_,
    log2FoldChange = grp_rows$coef,
    lfcSE          = grp_rows$stderr,
    stat           = grp_rows$coef / pmax(grp_rows$stderr, .Machine$double.xmin),
    pvalue         = grp_rows$pval,
    padj           = grp_rows$qval,
    reject         = !is.na(grp_rows$qval) & grp_rows$qval < DA_ALPHA,
    df             = NA_real_,
    stringsAsFactors = FALSE
  )
}

# ANCOMBC2 fit$res wide columns lfc_<coef>/se_/W_/p_/q_ ; group_coef e.g.
# "sex_cleanedM".
ancombc2_to_common <- function(ancombc_res, group_coef) {
  lfc_col <- paste0("lfc_", group_coef)
  se_col  <- paste0("se_", group_coef)
  W_col   <- paste0("W_", group_coef)
  p_col   <- paste0("p_", group_coef)
  q_col   <- paste0("q_", group_coef)

  data.frame(
    feature        = ancombc_res$taxon,
    variable       = group_coef,
    baseMean       = NA_real_,
    log2FoldChange = ancombc_res[[lfc_col]],
    lfcSE          = ancombc_res[[se_col]],
    stat           = ancombc_res[[W_col]],
    pvalue         = ancombc_res[[p_col]],
    padj           = ancombc_res[[q_col]],
    reject         = !is.na(ancombc_res[[q_col]]) & ancombc_res[[q_col]] < DA_ALPHA,
    df             = NA_real_,
    stringsAsFactors = FALSE
  )
}

# Build the ANCOMBC2 coefficient string for the group main effect:
# "<group_col><group2>" (no separator; matches lmer-style naming).
ancombc_group_coef <- function(group_col, group2) {
  paste0(group_col, group2)
}

# Parse an ANCOMBC2 "Estimation failed" error message and return the subset of
# `valid_covars` that should be DROPPED to fix the collinearity.
collinear_covars_from_error <- function(error_msg, valid_covars) {
  if (length(valid_covars) == 0) return(character(0))
  lines <- strsplit(error_msg, "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines)
  hdr <- grep("Estimation failed", lines)
  fin <- grep("Please ensure|collinearity|missing values", lines)
  if (length(hdr) == 0) return(character(0))
  end <- if (length(fin)) min(fin[fin > hdr[1]]) - 1L else length(lines)
  tokens <- lines[(hdr[1] + 1L):end]
  tokens <- tokens[nzchar(tokens)]
  drop <- character(0)
  for (tk in tokens) {
    cand <- valid_covars[vapply(valid_covars, function(v)
      substr(tk, 1L, nchar(v)) == v, logical(1))]
    if (length(cand) >= 1L) drop <- unique(c(drop, cand[which.max(nchar(cand))]))
  }
  drop
}

# ---- 7. Per-comparison output writer (shared) -------------------------------
write_comparison_outputs <- function(all_res, output_subdir, comparison,
                                     group1, group2,
                                     group1_col, group2_col,
                                     main_results_var_pattern,
                                     summary_row, run_config_row,
                                     main_results_filename = "main_results.tsv",
                                     sanitize_names = TRUE) {
  if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)

  fwrite(all_res, file.path(output_subdir, "all_results.tsv"), sep = "\t", quote = FALSE)

  main_res <- all_res[grepl(main_results_var_pattern, all_res$variable), , drop = FALSE]
  main_res <- main_res[order(main_res$padj, -abs(main_res$log2FoldChange)), , drop = FALSE]
  fwrite(main_res, file.path(output_subdir, main_results_filename), sep = "\t", quote = FALSE)
  if (main_results_filename != "main_results.tsv") {
    fwrite(main_res, file.path(output_subdir, "main_results.tsv"), sep = "\t", quote = FALSE)
  }

  sig <- main_res[!is.na(main_res$padj) & main_res$padj < DA_ALPHA, , drop = FALSE]
  if (nrow(sig) > 0) {
    fwrite(sig, file.path(output_subdir, "significant_results.tsv"), sep = "\t", quote = FALSE)
    g1_name <- if (sanitize_names) sanitize(group1) else group1
    g2_name <- if (sanitize_names) sanitize(group2) else group2
    sig_g1 <- sig[sig$enriched_in == group1, , drop = FALSE]
    sig_g2 <- sig[sig$enriched_in == group2, , drop = FALSE]
    if (nrow(sig_g1) > 0) {
      fwrite(sig_g1, file.path(output_subdir, paste0("significant_enriched_in_", g1_name, ".tsv")),
             sep = "\t", quote = FALSE)
    }
    if (nrow(sig_g2) > 0) {
      fwrite(sig_g2, file.path(output_subdir, paste0("significant_enriched_in_", g2_name, ".tsv")),
             sep = "\t", quote = FALSE)
    }
  }

  fwrite(summary_row, file.path(output_subdir, "summary.tsv"), sep = "\t", quote = FALSE)
  fwrite(run_config_row, file.path(output_subdir, "run_config.tsv"), sep = "\t", quote = FALSE)
  writeLines("done", file.path(output_subdir, "analysis.done"))

  n_sig <- nrow(sig)
  n_up_g1 <- if (n_sig > 0) sum(sig$enriched_in == group1, na.rm = TRUE) else 0L
  n_up_g2 <- if (n_sig > 0) sum(sig$enriched_in == group2, na.rm = TRUE) else 0L
  cat("    Significant taxa (padj < ", DA_ALPHA, "): ", n_sig, "\n", sep = "")
  cat("      enriched in ", group1, ": ", n_up_g1, "\n", sep = "")
  cat("      enriched in ", group2, ": ", n_up_g2, "\n", sep = "")

  if (n_sig > 0) {
    top_n <- min(10L, n_sig)
    cat("    Top ", top_n, " by padj:\n", sep = "")
    print(sig[seq_len(top_n), c("genus_short", "log2FoldChange", "padj", "enriched_in")],
          row.names = FALSE)
  }
  invisible(list(n_sig = n_sig, n_up_g1 = n_up_g1, n_up_g2 = n_up_g2))
}

# ---- 8. ANCOMBC2 genus-fittability pre-filter (da_v4 aligned) ---------------
filter_fittable <- function(counts, meta, group_col, project_col = COL_PROJECT) {
  detected <- counts > 0
  bp_ind  <- model.matrix(~ 0 + meta[[project_col]])
  grp_ind <- model.matrix(~ 0 + meta[[group_col]])
  keep_fit <- rowSums((detected %*% bp_ind)  > 0) >= 2 &
              rowSums((detected %*% grp_ind) > 0) >= 2
  cat("    Unfittable genera dropped: ", sum(!keep_fit), " / ", length(keep_fit),
      "\n", sep = "")
  counts[keep_fit, , drop = FALSE]
}