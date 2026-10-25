# ==============================================================================
# da_helpers.R  (age OVA da_v2 — adapted from age_differential da_helpers.R)
# Shared preprocessing + harmonization helpers for the three-method
# differential abundance (DA) analysis of AGE: LinDA + Maaslin2 + ANCOMBC2.
# Main effect = age_cleaned_group (pre-binned life stages), run as OVA
# (one-vs-ALL) comparisons WITHIN saliva (the only site in this pipeline):
# each valid age group (>= MIN_SAMPLES_PER_GROUP samples in the site) is the
# TARGET, and the reference "Other" pools the samples of ALL OTHER VALID age
# groups in the same site (invalid <10-sample groups are excluded entirely,
# consistent with the OVO pipeline). Covariates exclude age (main) and site
# (within-site stratum):
#   sex_cleaned, platform_cleaned, sequencing_region_cleaned,
#   bmi_cleaned_group, continent_group.
#
# Why binned age (not continuous): the cohort spans 0-6 months to >60 years and
# the microbiome-age relationship is strongly nonlinear (infant succession ->
# stable adulthood -> elderly drift); a single linear slope would be
# uninterpretable. Bins absorb the nonlinearity.
#
# DESIGN (identical to the geography/site/sex da_v2 scripts except where noted):
#   - Data loaded from plain .tsv via read.table (ONE table per driver).
#   - reads>=1000 sample filter applied at TSV export (prepare_genus_tables.py).
#   - 1% feature prevalence stays INSIDE each tool (prev_cut = 0.01).
#   - bmi covariate: pre-binned bmi_cleaned_group column.
#   - Conditional (1|bioproject_id) random effect when a BioProject spans both
#     groups (target vs Other) of a comparison.
#   - SCOPE: --scope=all or --scope=hc (status_cleaned == "HC" only).
#     Results root: rd("age_ova") vs rd("age_hc_ova").
#
# Locked constants match the da_v2 scripts verbatim:
#   seed 20260615; alpha 0.05; prev_cut 0.01; min_samples_per_group = 10
#   LinDA (MicrobiomeStat::linda): feature.dat.type="count", is.winsor=FALSE,
#          prev.filter=0.01, n.cores=8, BH, alpha=0.05
#   Maaslin2: TSS+LOG+LM, BH, min_prevalence=0.01, standardize=FALSE, cores=1
#   ANCOMBC2: fix_formula binary contrast, (1|bioproject_id), prv_cut=0.01,
#             lib_cut=0, struc_zero=FALSE, neg_lb=FALSE, n_cl=8
#
# Contrast convention (used consistently in outputs):
#   reference = "Other" (all other valid age groups in the site),
#   target    = the focal age group
#   ->  log2FoldChange > 0  =>  enriched in the FOCAL age group
#   comparison name = "<site>_enrich_<g>" with sanitized level labels
#   (e.g. saliva_enrich_18-60_years; full-width ＞ -> "gt").
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

# ---- Locked constants --------------------------------------------------------
DA_SEED   <- 20260615
DA_ALPHA  <- 0.05

# ---- INPUT FILTERS ----------------------------------------------------------
PREV_CUT                 <- 0.01   # feature prevalence (per-comparison)
MIN_ABUNDANCE            <- 0      # feature abundance floor (Maaslin2); 0 = none
LIB_CUT                  <- 0      # library-size floor (ANCOMBC2/LinDA); 0 = none
MIN_SAMPLES_PER_GROUP    <- 10L    # min samples PER AGE GROUP within the site to keep that group
                                   # (both as OVA targets AND inside the "Other" reference pool)
DROP_BIOPROJECTS         <- NULL   # optional bioproject_id vector to exclude (NULL = off)
DROP_SITES               <- NULL   # optional sites to drop (NULL = off)
# <<< end input filters

COL_SAMPLE     <- "sample_id"
COL_CONTINENT  <- "continent_group"
COL_SITE       <- "site_cleaned"
COL_SITE_GROUP <- "site_group_cleaned"
COL_AGE        <- "age_cleaned_group"
COL_SEX        <- "sex_cleaned"
COL_STATUS     <- "status_cleaned"
COL_BMI        <- "bmi_cleaned_group"
COL_PROJECT    <- "bioproject_id"
COL_PLATFORM   <- "platform_cleaned"
COL_SEQ_REGION <- "sequencing_region_cleaned"

# Canonical life-stage order (NOTE: the metadata uses the FULL-WIDTH "＞"
# character in "＞60 years" — do not "fix" it to ASCII here; it must match the
# metadata values exactly).
AGE_ORDER <- c(
  "0-6 months", "6-36 months", "3-12 years", "12-18 years",
  "18-60 years", "＞60 years"
)

# OVA pipeline is SALIVA-ONLY (user request 2026-09-28): restricting
# SITE_ORDER to saliva filters the metadata at load time, so all downstream
# stages (comparisons, merge, figures) automatically cover saliva only.
SITE_ORDER <- c("saliva")

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

# Sanitize an age-group label for use in comparison names / directory names:
# spaces -> "_", full-width ＞ -> "gt" (ASCII-safe).
sanitize_age_level <- function(x) {
  x <- gsub(" ", "_", x)
  gsub("＞", "gt", x, fixed = TRUE)
}

# Resolve a results subdirectory under the optional DA_RESULTS_ROOT.
rd <- function(subpath) {
  root <- Sys.getenv("DA_RESULTS_ROOT", "")
  if (nchar(root) > 0) file.path(root, subpath) else file.path("differential_abundance", "result", subpath)
}

# Scope -> results subpath: "all" -> rd("age_ova/..."), "hc" -> rd("age_hc_ova/...")
age_root <- function(scope) {
  stopifnot(scope %in% c("all", "hc"))
  rd(if (scope == "hc") "age_hc_ova" else "age_ova")
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
# Age scope: keep all valid oral sites (Oral group); no main-col filter here —
# each per-site comparison subsets its own site and its two age groups.
# hc_only = TRUE restricts to status_cleaned == "HC".
load_age_metadata <- function(hc_only = FALSE) {
  if (!file.exists(METADATA_FILE)) stop("Metadata not found: ", METADATA_FILE)
  metadata_raw <- fread(
    METADATA_FILE,
    sep = "\t", header = TRUE,
    quote = "", fill = TRUE, check.names = FALSE,
    data.table = FALSE
  )

  required_cols <- c(COL_SAMPLE, COL_SITE, COL_AGE, COL_PROJECT)
  missing_cols <- setdiff(required_cols, names(metadata_raw))
  if (length(missing_cols) > 0) {
    stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  for (opt_col in c(COL_SEX, COL_BMI, COL_PLATFORM, COL_SEQ_REGION,
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
      !!COL_AGE         := trimws(as.character(.data[[COL_AGE]])),
      !!COL_STATUS      := trimws(as.character(.data[[COL_STATUS]])),
      !!COL_SEX         := trimws(as.character(.data[[COL_SEX]])),
      !!COL_BMI         := trimws(as.character(.data[[COL_BMI]])),
      !!COL_CONTINENT   := trimws(as.character(.data[[COL_CONTINENT]])),
      !!COL_PLATFORM    := trimws(as.character(.data[[COL_PLATFORM]])),
      !!COL_SEQ_REGION  := trimws(as.character(.data[[COL_SEQ_REGION]])),
      !!COL_PROJECT     := trimws(as.character(.data[[COL_PROJECT]])),
      !!COL_SITE_GROUP  := trimws(as.character(.data[[COL_SITE_GROUP]]))
    ) %>%
    mutate(
      !!COL_AGE         := ifelse(.data[[COL_AGE]] %in% c("", "NA", "na"), "Unknown", .data[[COL_AGE]]),
      !!COL_SEX         := ifelse(.data[[COL_SEX]] %in% c("", "NA", "na"), "Unknown", .data[[COL_SEX]]),
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
# return the common-samples subset.
align_samples <- function(meta, genus_mat) {
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

# ---- 3. OVA comparison generation -------------------------------------------
# Valid age groups within a site: >= min_n samples in that site, in AGE_ORDER.
valid_age_groups <- function(meta, site, min_n = MIN_SAMPLES_PER_GROUP) {
  keep <- meta[[COL_SITE]] == site & meta[[COL_AGE]] %in% AGE_ORDER
  counts <- table(factor(meta[[COL_AGE]][keep], levels = AGE_ORDER))
  names(counts)[counts >= min_n]
}

# One OVA comparison per valid age group in the site: TARGET = the focal
# group g, REFERENCE "Other" = all samples of the OTHER VALID age groups in
# the same site (invalid <min_n groups excluded entirely). Each entry keeps
# `other_groups` so prep_ova can pool the reference samples.
make_site_age_ova <- function(site, valid_groups) {
  lapply(valid_groups, function(g) {
    list(
      site = site,
      g1 = "Other",                 # reference (pooled other valid groups)
      g2 = g,                       # target (focal age group)
      other_groups = setdiff(valid_groups, g),
      comparison_name = paste0(sanitize(site), "_enrich_",
                               sanitize_age_level(g))
    )
  })
}

# All OVA comparisons across a set of sites (site order = SITE_ORDER).
make_age_ova_comparisons <- function(meta, sites,
                                     min_n = MIN_SAMPLES_PER_GROUP) {
  out <- list()
  for (site in sites) {
    vg <- valid_age_groups(meta, site, min_n)
    out <- c(out, make_site_age_ova(site, vg))
  }
  out
}

# Parse an OVA comparison name "<site>_enrich_<g>" -> list(site, target),
# matched against SITE_ORDER x AGE_ORDER sanitized labels.
parse_ova_comparison <- function(comp) {
  for (site in SITE_ORDER) {
    pre <- paste0(sanitize(site), "_enrich_")
    if (startsWith(comp, pre)) {
      rest <- substr(comp, nchar(pre) + 1L, nchar(comp))
      for (g in AGE_ORDER) {
        if (identical(sanitize_age_level(g), rest))
          return(list(site = site, target = g))
      }
    }
  }
  list(site = NA_character_, target = NA_character_)
}

# ---- 4. Covariate preparation (collapse exclusive levels to "Other") -------
prepare_covars <- function(meta, group_col, covars, bmi_col = COL_BMI) {
  valid_covars <- character(0)
  meta_out <- meta

  for (cov in covars) {
    vals <- as.character(meta_out[[cov]])
    vals[is.na(vals) | vals %in% c("", "NA", "na")] <- "Unknown"

    cross_tab <- table(meta_out[[group_col]], vals)
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
          " has no level shared by the two age groups after collapse, dropping.\n", sep = "")
      next
    }

    ref_level <- names(tab)[1]
    if (cov == bmi_col && "Normal" %in% names(tab)) ref_level <- "Normal"
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

# Drop covariates whose dummy columns make the fixed design rank-deficient.
# LinDA fits ONE multivariate lm across all features; with an aliased column
# its coef summary has fewer rows than the covariance blocks and it crashes
# ("length of 'dimnames' [1] not equal to array extent"). Dropping the
# offending covariates up front keeps the design full rank for all methods.
drop_aliased_covars <- function(meta, main_col, valid_covars, project_col = COL_PROJECT) {
  # random effect is handled separately by add_random_effect; rank check on
  # the fixed design only.
  repeat {
    terms <- c(main_col, valid_covars)
    ff <- paste0("~ ", paste(terms, collapse = " + "))
    X <- tryCatch(model.matrix(as.formula(ff), data = meta),
                  error = function(e) NULL)
    if (is.null(X) || qr(X)$rank == ncol(X)) break
    q <- qr(X)
    aliased <- setdiff(colnames(X), colnames(X)[q$pivot[seq_len(q$rank)]])
    drop_cov <- unique(vapply(aliased, function(a) {
      hit <- valid_covars[startsWith(a, valid_covars)]
      if (length(hit) > 0) hit[which.max(nchar(hit))] else NA_character_
    }, character(1)))
    drop_cov <- drop_cov[!is.na(drop_cov) & drop_cov %in% valid_covars]
    if (length(drop_cov) == 0) break
    cat("    NOTE: dropping rank-deficient covariate(s): ",
        paste(drop_cov, collapse = ", "), "\n", sep = "")
    valid_covars <- setdiff(valid_covars, drop_cov)
    if (length(valid_covars) == 0) break
  }
  valid_covars
}

# ---- 5. Conditional random effect for BioProject ----------------------------
add_random_effect <- function(base_formula, meta, group_col,
                              project_col = COL_PROJECT) {
  meta[[project_col]] <- droplevels(factor(as.character(meta[[project_col]])))
  n_levels <- nlevels(meta[[project_col]])
  project_tab <- table(meta[[group_col]], meta[[project_col]])
  shared_projects <- sum(colSums(project_tab > 0) == 2)

  # A single-level factor (all samples from one BioProject) makes the random
  # effect degenerate and crashes model.matrix ("contrasts can be applied
  # only to factors with 2 or more levels"); fall back to fixed effects.
  if (shared_projects > 0 && n_levels >= 2) {
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
make_genus_short <- function(feature) {
  g <- sub(".*g__", "", feature)
  sub(";.*", "", g)
}

# group1 = reference ("Other" pooled valid groups), group2 = target (focal
# age group), so coef > 0 => enriched in the FOCAL age group.
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

# Resolve the design-matrix coefficient name for the main contrast across
# tools/levels: exact "<col><level>" first, then a longest matching prefix
# (levels contain spaces / full-width chars that tools may reformat).
resolve_group_coef <- function(avail_names, group_col, group2) {
  exact <- paste0(group_col, group2)
  if (exact %in% avail_names) return(exact)
  pref <- paste0(group_col, substr(group2, 1, 3))
  cand <- avail_names[startsWith(avail_names, pref)]
  if (length(cand) >= 1) return(cand[which.max(nchar(cand))])
  exact  # let the caller fail with a clear message
}

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
    g1_name <- if (sanitize_names) sanitize_age_level(group1) else group1
    g2_name <- if (sanitize_names) sanitize_age_level(group2) else group2
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
  # model.matrix(~0+single-level factor) crashes ("contrasts can be applied
  # only to factors with 2 or more levels"); with one BioProject the
  # >=2-projects criterion is trivially untestable -> use a group-only filter.
  bp_levels <- nlevels(droplevels(factor(as.character(meta[[project_col]]))))
  grp_ind <- model.matrix(~ 0 + meta[[group_col]])
  keep_fit <- if (bp_levels >= 2) {
    bp_ind <- model.matrix(~ 0 + meta[[project_col]])
    rowSums((detected %*% bp_ind)  > 0) >= 2 &
    rowSums((detected %*% grp_ind) > 0) >= 2
  } else {
    rowSums((detected %*% grp_ind) > 0) >= 2
  }
  cat("    Unfittable genera dropped: ", sum(!keep_fit), " / ", length(keep_fit),
      "\n", sep = "")
  counts[keep_fit, , drop = FALSE]
}