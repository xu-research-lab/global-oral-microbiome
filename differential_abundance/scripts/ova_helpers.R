# Metadata cleaning and model preparation shared by OVA analyses.

AGE_ORDER <- c("0-6 months", "6-36 months", "3-12 years", "12-18 years", "18-60 years", "＞60 years")

sanitize_age_level <- function(x) {
  x <- gsub(" ", "_", x)
  gsub("＞", "gt", x, fixed = TRUE)
}

load_age_metadata <- function(hc_only = FALSE) {
  if (!file.exists(METADATA_FILE)) {
    stop("Metadata not found: ", METADATA_FILE)
  }
  metadata_raw <- fread(METADATA_FILE,
    sep = "\t", header = TRUE, quote = "", fill = TRUE, check.names = FALSE,
    data.table = FALSE
  )
  required_cols <- c(COL_SAMPLE, COL_SITE, COL_AGE, COL_PROJECT)
  missing_cols <- setdiff(required_cols, names(metadata_raw))
  if (length(missing_cols) > 0) {
    stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  for (opt_col in c(
    COL_SEX, COL_BMI, COL_PLATFORM, COL_SEQ_REGION, COL_CONTINENT, COL_SITE_GROUP,
    COL_STATUS
  )) {
    if (!opt_col %in% names(metadata_raw)) {
      metadata_raw[[opt_col]] <- NA
    }
  }
  metadata_clean <- metadata_raw %>%
    filter(!is.na(.data[[COL_SAMPLE]]) & .data[[COL_SAMPLE]] != "") %>%
    distinct(.data[[COL_SAMPLE]], .keep_all = TRUE) %>%
    mutate(
      `:=`(!!COL_SITE, trimws(as.character(.data[[COL_SITE]]))),
      `:=`(!!COL_AGE, trimws(as.character(.data[[COL_AGE]]))), `:=`(!!COL_STATUS, trimws(as.character(.data[[COL_STATUS]]))),
      `:=`(!!COL_SEX, trimws(as.character(.data[[COL_SEX]]))), `:=`(!!COL_BMI, trimws(as.character(.data[[COL_BMI]]))),
      `:=`(!!COL_CONTINENT, trimws(as.character(.data[[COL_CONTINENT]]))), `:=`(!!COL_PLATFORM, trimws(as.character(.data[[COL_PLATFORM]]))),
      `:=`(!!COL_SEQ_REGION, trimws(as.character(.data[[COL_SEQ_REGION]]))), `:=`(!!COL_PROJECT, trimws(as.character(.data[[COL_PROJECT]]))),
      `:=`(!!COL_SITE_GROUP, trimws(as.character(.data[[COL_SITE_GROUP]])))
    ) %>%
    mutate(
      `:=`(
        !!COL_AGE,
        ifelse(.data[[COL_AGE]] %in% c("", "NA", "na"), "Unknown", .data[[COL_AGE]])
      ), `:=`(
        !!COL_SEX,
        ifelse(.data[[COL_SEX]] %in% c("", "NA", "na"), "Unknown", .data[[COL_SEX]])
      ), `:=`(
        !!COL_BMI,
        ifelse(.data[[COL_BMI]] %in% c("", "NA", "na"), "Unknown", .data[[COL_BMI]])
      ), `:=`(
        !!COL_PLATFORM,
        ifelse(.data[[COL_PLATFORM]] %in% c("", "NA", "na"), "Unknown", .data[[COL_PLATFORM]])
      ), `:=`(
        !!COL_SEQ_REGION,
        ifelse(.data[[COL_SEQ_REGION]] %in% c("", "NA", "na"), "Unknown", .data[[COL_SEQ_REGION]])
      ),
      `:=`(!!COL_CONTINENT, ifelse(.data[[COL_CONTINENT]] %in% c("", "NA", "na"), "Unknown", .data[[COL_CONTINENT]]))
    )
  metadata_clean <- metadata_clean %>% filter(.data[[COL_SITE]] %in% SITE_ORDER & !is.na(.data[[COL_SITE]]) &
    .data[[COL_SITE]] != "")
  sg <- trimws(as.character(metadata_clean[[COL_SITE_GROUP]]))
  keep_sg <- is.na(sg) | sg %in% c("", "NA", "Oral")
  metadata_clean <- metadata_clean[keep_sg, , drop = FALSE]
  if (hc_only) {
    n0 <- nrow(metadata_clean)
    metadata_clean <- metadata_clean[!is.na(metadata_clean[[COL_STATUS]]) & metadata_clean[[COL_STATUS]] ==
      "HC", , drop = FALSE]
    cat("    HC scope: kept ", nrow(metadata_clean), " / ", n0, " samples (status_cleaned == HC)\n",
      sep = ""
    )
  }
  metadata_clean <- as.data.frame(metadata_clean)
  rownames(metadata_clean) <- metadata_clean[[COL_SAMPLE]]
  metadata_clean
}

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
      cat("    NOTE: ", cov, " exclusive levels collapsed to Other: ", paste(exclusive_levels,
        collapse = ", "
      ), "\n", sep = "")
    }
    vals <- droplevels(factor(vals))
    tab <- sort(table(vals), decreasing = TRUE)
    if (length(tab) < 2) {
      cat("    WARNING: ", cov, " has only 1 usable level, dropping.\n", sep = "")
      next
    }
    cross2 <- table(meta_out[[group_col]], vals)
    if (sum(colSums(cross2 > 0) == 2) < 1) {
      cat("    WARNING: ", cov, " has no level shared by the two age groups after collapse, dropping.\n",
        sep = ""
      )
      next
    }
    ref_level <- names(tab)[1]
    if (cov == bmi_col && "Normal" %in% names(tab)) {
      ref_level <- "Normal"
    }
    if (ref_level %in% c("Other", "Unknown")) {
      prefer <- setdiff(names(tab), c("Other", "Unknown"))
      if (length(prefer) > 0) {
        ref_level <- prefer[1]
      }
    }
    other_levels <- setdiff(names(tab), ref_level)
    meta_out[[cov]] <- factor(as.character(vals), levels = c(ref_level, other_levels))
    valid_covars <- c(valid_covars, cov)
    cat("    [OK] Kept ", cov, " (Ref: ", ref_level, "; levels: ", paste(levels(meta_out[[cov]]),
      collapse = ", "
    ), ")\n", sep = "")
  }
  list(meta = meta_out, covars = valid_covars)
}

drop_aliased_covars <- function(meta, main_col, valid_covars, project_col = COL_PROJECT) {
  repeat {
    terms <- c(main_col, valid_covars)
    ff <- paste0("~ ", paste(terms, collapse = " + "))
    X <- tryCatch(model.matrix(as.formula(ff), data = meta), error = function(e) NULL)
    if (is.null(X) || qr(X)$rank == ncol(X)) {
      break
    }
    q <- qr(X)
    aliased <- setdiff(colnames(X), colnames(X)[q$pivot[seq_len(q$rank)]])
    drop_cov <- unique(vapply(aliased, function(a) {
      hit <- valid_covars[startsWith(a, valid_covars)]
      if (length(hit) > 0) {
        hit[which.max(nchar(hit))]
      } else {
        NA_character_
      }
    }, character(1)))
    drop_cov <- drop_cov[!is.na(drop_cov) & drop_cov %in% valid_covars]
    if (length(drop_cov) == 0) {
      break
    }
    cat("    NOTE: dropping rank-deficient covariate(s): ", paste(drop_cov, collapse = ", "), "\n",
      sep = ""
    )
    valid_covars <- setdiff(valid_covars, drop_cov)
    if (length(valid_covars) == 0) {
      break
    }
  }
  valid_covars
}

add_random_effect <- function(base_formula, meta, group_col, project_col = COL_PROJECT) {
  meta[[project_col]] <- droplevels(factor(as.character(meta[[project_col]])))
  n_levels <- nlevels(meta[[project_col]])
  project_tab <- table(meta[[group_col]], meta[[project_col]])
  shared_projects <- sum(colSums(project_tab > 0) == 2)
  if (shared_projects > 0 && n_levels >= 2) {
    full_formula <- paste0(base_formula, " + (1|", project_col, ")")
    rand_formula <- paste0("(1|", project_col, ")")
    cat("    Using random effect for BioProject (", shared_projects, " projects in both groups).\n",
      sep = ""
    )
    list(
      fixed_formula = base_formula, rand_formula = rand_formula, full_formula = full_formula,
      has_random = TRUE, n_shared = shared_projects
    )
  } else {
    cat("    No BioProjects spanning both groups; fixed-effects only.\n")
    list(
      fixed_formula = base_formula, rand_formula = NULL, full_formula = base_formula, has_random = FALSE,
      n_shared = 0L
    )
  }
}

resolve_group_coef <- function(avail_names, group_col, group2) {
  exact <- paste0(group_col, group2)
  if (exact %in% avail_names) {
    return(exact)
  }
  pref <- paste0(group_col, substr(group2, 1, 3))
  cand <- avail_names[startsWith(avail_names, pref)]
  if (length(cand) >= 1) {
    return(cand[which.max(nchar(cand))])
  }
  exact
}

filter_fittable <- function(counts, meta, group_col, project_col = COL_PROJECT) {
  detected <- counts > 0
  bp_levels <- nlevels(droplevels(factor(as.character(meta[[project_col]]))))
  grp_ind <- model.matrix(~ 0 + meta[[group_col]])
  keep_fit <- if (bp_levels >= 2) {
    bp_ind <- model.matrix(~ 0 + meta[[project_col]])
    rowSums((detected %*% bp_ind) > 0) >= 2 & rowSums((detected %*% grp_ind) > 0) >= 2
  } else {
    rowSums((detected %*% grp_ind) > 0) >= 2
  }
  cat("    Unfittable genera dropped: ", sum(!keep_fit), " / ", length(keep_fit), "\n", sep = "")
  counts[keep_fit, , drop = FALSE]
}
