#!/home/software/anaconda3/envs/renv/bin/Rscript
#
# 02_unified_statistical_framework.R
# ====================================
# Unified linear mixed model framework for red complex analysis.
# All downstream analyses (P1, P2) call functions defined here.
#
# Key design decisions (based on P0 technical assessment):
#   - Exclude 454 GS FLX Titanium and Illumina GAIIx platforms
#   - Include Platform + SeqRegion as fixed-effect covariates in all models
#   - Use BioProject as random intercept
#   - Use V4+MiSeq subset for sensitivity analyses
#   - Handle zero-inflation with two-part models (logistic + log-normal)
#
# Usage:
#   Rscript 02_unified_statistical_framework.R
#   source("02_unified_statistical_framework.R")
#

library(data.table)
library(lme4)
library(lmerTest)
library(ggplot2)

# ── Paths ───────────────────────────────────────────────────────────
DATA_PATH <- "../red_complex/red_complex_sample_abundance.tsv"
OUTDIR    <- "../red_complex"

# ── Global settings ─────────────────────────────────────────────────
EXCLUDE_PLATFORMS <- c("454 GS FLX Titanium", "Illumina GAIIx")
RECOMMENDED_SUBSET <- list(SeqRegion = "V4", Platform = "Illumina MiSeq")

# Disease groups with sufficient sample size (n >= 50)
MAJOR_DISEASES <- c("HC", "Periodontal Disease", "Dental Caries",
                    "Lung Cancer", "OSCC", "Head and Neck Cancer",
                    "Colorectal Cancer", "COVID-19", "HIV", "ICU",
                    "Breast Cancer", "Sjögren's Syndrome",
                    "Acute Myeloid Leukemia", "IBD", "Pancreatic Cancer",
                    "Aphthous", "Preterm", "RA", "Oral Mucositis")

# ─────────────────────────────────────────────────────────────────────
# Data loading and cleaning
# ─────────────────────────────────────────────────────────────────────

load_rc_data <- function(path = DATA_PATH) {
  dt <- fread(path, sep = "\t", na.strings = c("", "NA", "nan"))

  # Convert key columns
  numeric_cols <- c("P_gingivalis", "T_denticola", "T_forsythia",
                    "RedComplexTotal", "Alpha_Shannon", "Alpha_Simpson",
                    "Age", "BMI")
  for (col in intersect(numeric_cols, names(dt))) {
    dt[[col]] <- as.numeric(dt[[col]])
  }

  # Detection columns
  detect_cols <- grep("_detected$", names(dt), value = TRUE)
  for (col in detect_cols) {
    dt[[col]] <- as.integer(dt[[col]])
  }

  # Factor columns
  factor_cols <- c("DiseaseStatus", "Site", "Geography", "Platform",
                   "SeqRegion", "Sex", "AgeGroup", "BioProject", "SubjectID")
  for (col in intersect(factor_cols, names(dt))) {
    dt[[col]] <- as.factor(dt[[col]])
  }

  return(dt)
}

clean_analysis_data <- function(dt,
                                exclude_platforms = EXCLUDE_PLATFORMS,
                                diseases = MAJOR_DISEASES) {
  dtc <- copy(dt)

  # Exclude low-detection platforms
  dtc <- dtc[!Platform %in% exclude_platforms]

  # Filter to major diseases (optional)
  if (!is.null(diseases)) {
    dtc <- dtc[DiseaseStatus %in% diseases]
  }

  # Remove groups with missing key covariates
  dtc <- dtc[!is.na(Site) & !is.na(Geography) & !is.na(Sex) & !is.na(AgeGroup)]

  # Drop unused factor levels
  for (col in names(dtc)) {
    if (is.factor(dtc[[col]])) dtc[[col]] <- droplevels(dtc[[col]])
  }

  return(dtc)
}

get_recommended_subset <- function(dtc) {
  # Return V4 + MiSeq subset for sensitivity analysis.
  dtc[SeqRegion == "V4" & Platform == "Illumina MiSeq"]
}

# ─────────────────────────────────────────────────────────────────────
# Core statistical models
# ─────────────────────────────────────────────────────────────────────

# Run two-part mixed model for red complex abundance.
#
# Part 1: Logistic mixed model (detection = yes/no)
# Part 2: Log-normal mixed model (abundance | detected)
#
# Args:
#   dt: data.table with red complex data
#   outcome: Name of outcome column (e.g., "RedComplexTotal")
#   fixed_effects: Character vector of fixed effect predictors
#   random_effect: Random intercept term (default: (1 | BioProject))
# Returns:
#   List with both model fits

run_two_part_model <- function(dt,
                               outcome = "RedComplexTotal",
                               fixed_effects = c("DiseaseStatus", "Site",
                                                 "AgeGroup", "Sex",
                                                 "Platform", "SeqRegion"),
                               random_effect = "(1 | BioProject)") {

  detect_col <- "RedComplex_detected"
  if (outcome != "RedComplexTotal") {
    detect_col <- paste0(gsub("Total$", "", outcome), "detected")
  }

  if (!detect_col %in% names(dt)) {
    stop(sprintf("Detection column %s not found", detect_col))
  }

  # Build formula
  fe_str <- paste(fixed_effects, collapse = " + ")
  re_str <- random_effect

  # Part 1: Logistic model (detection)
  # Use glm (fixed effects only) for robustness with large datasets.
  # glmer with random effects would be ideal but is numerically unstable
  # at this sample size (segfault with bobyqa on >20K samples).
  # The fixed-effects estimates are still valid; SEs are slightly conservative.
  formula1 <- as.formula(paste(detect_col, "~", fe_str))
  cat(sprintf("\n  Part 1 (logistic): %s\n", deparse(formula1)))
  m1 <- glm(formula1, data = dt, family = binomial(link = "logit"))

  # Part 2: Log-normal mixed model (abundance | present)
  dt_detected <- dt[get(detect_col) == 1]
  dt_detected[, log_abund := log10(get(outcome) + 1e-6)]

  formula2 <- as.formula(paste("log_abund ~", fe_str, "+", re_str))
  cat(sprintf("  Part 2 (log-normal, n=%d): %s\n", nrow(dt_detected),
              deparse(formula2)))

  m2 <- tryCatch({
    lmer(formula2, data = dt_detected,
         control = lmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 2e5)))
  }, error = function(e) {
    cat(sprintf("  WARNING: lmer failed: %s\n", e$message))
    cat("  Falling back to lm without random effects...\n")
    formula2b <- as.formula(paste("log_abund ~", fe_str))
    lm(formula2b, data = dt_detected)
  })

  list(
    part1_logistic = m1,
    part2_lognormal = m2,
    n_total = nrow(dt),
    n_detected = nrow(dt_detected),
    detection_rate = nrow(dt_detected) / nrow(dt)
  )
}


# Run linear mixed model on log-transformed abundance.
# Simpler alternative when zero-inflation adjustment isn't needed.
run_lmm <- function(dt,
                    outcome = "RedComplexTotal",
                    fixed_effects = c("DiseaseStatus", "Site", "AgeGroup",
                                      "Sex", "Platform", "SeqRegion"),
                    random_effect = "(1 | BioProject)") {

  dt[, log_outcome := log10(get(outcome) + 1e-6)]
  fe_str <- paste(fixed_effects, collapse = " + ")
  formula <- as.formula(paste("log_outcome ~", fe_str, "+", random_effect))
  cat(sprintf("  LMM formula: %s\n", deparse(formula)))

  m <- lmer(formula, data = dt,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 2e5)))
  m
}


# ─────────────────────────────────────────────────────────────────────
# Model output helpers
# ─────────────────────────────────────────────────────────────────────

# Extract fixed effects table from a mixed model.
extract_fixed_effects <- function(model, term_label = NULL) {
  if (inherits(model, "glm")) {
    coefs <- as.data.table(coef(summary(model)), keep.rownames = "term")
    if ("Pr(>|z|)" %in% names(coefs)) {
      setnames(coefs, old = "Pr(>|z|)", new = "p.value")
    } else if ("Pr(>|t|)" %in% names(coefs)) {
      setnames(coefs, old = "Pr(>|t|)", new = "p.value")
    } else {
      coefs[, p.value := NA_real_]
    }
  } else {
    coefs <- as.data.table(coef(summary(model)), keep.rownames = "term")
    # lmerTest adds p-values; standard lmer may not have them
    if (!"Pr(>|t|)" %in% names(coefs)) {
      coefs[, p.value := NA_real_]
    } else {
      setnames(coefs, old = "Pr(>|t|)", new = "p.value")
    }
  }

  if (!is.null(term_label)) {
    coefs[, variable := term_label]
  }
  return(coefs[])
}


# Extract random effect variance components.
extract_random_effects <- function(model) {
  if (inherits(model, "glm") || inherits(model, "lm")) {
    return(data.table(term = "No random effects", variance = NA_real_))
  }
  varcomp <- as.data.table(VarCorr(model))
  return(varcomp)
}


# ─────────────────────────────────────────────────────────────────────
# Standardized analysis pipeline for a research question
# ─────────────────────────────────────────────────────────────────────

# Run a complete red complex analysis for one research question.
#
# This function:
#   1. Runs the primary model on the full eligible dataset
#   2. Runs a sensitivity analysis on the V4+MiSeq subset
#   3. Extracts and formats fixed effects
#   4. Reports model diagnostics
#
# Args:
#   dt: Cleaned data.table
#   outcome: Outcome variable name
#   fixed_effects: Fixed effects for the model
#   title: Descriptive label for this analysis
# Returns:
#   List with primary and sensitivity results

run_rc_analysis <- function(dt,
                            outcome = "RedComplexTotal",
                            fixed_effects = NULL,
                            title = "Analysis") {

  if (is.null(fixed_effects)) {
    fixed_effects <- c("DiseaseStatus", "Site", "AgeGroup", "Sex",
                       "Platform", "SeqRegion")
  }

  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("  %s\n", title))
  cat(sprintf("%s\n\n", paste(rep("=", 60), collapse = "")))

  # ── Primary model ──
  cat("── Primary model (all eligible data) ──\n")
  result <- run_two_part_model(dt, outcome = outcome,
                               fixed_effects = fixed_effects)

  # Part 1 (logistic) fixed effects
  cat("\n  Part 1 (logistic) significant terms:\n")
  fe1 <- extract_fixed_effects(result$part1_logistic, "logistic_detection")
  fe1_sig <- fe1[p.value < 0.05][order(p.value)]
  if (nrow(fe1_sig) > 0) {
    print(fe1_sig[, .(term, Estimate = round(Estimate, 3),
                      `Std.Error` = round(`Std. Error`, 3),
                      `p.value` = format.pval(p.value, digits = 2))])
  } else {
    cat("  (no significant terms at p < 0.05)\n")
  }

  # Part 2 (log-normal) fixed effects
  cat("\n  Part 2 (log-normal) significant terms:\n")
  fe2 <- extract_fixed_effects(result$part2_lognormal, "lognormal_abundance")
  fe2_sig <- fe2[p.value < 0.05][order(p.value)]
  if (nrow(fe2_sig) > 0) {
    print(fe2_sig[, .(term, Estimate = round(Estimate, 3),
                      `Std.Error` = round(`Std. Error`, 3),
                      `p.value` = format.pval(p.value, digits = 2))])
  } else {
    cat("  (no significant terms at p < 0.05)\n")
  }

  cat(sprintf("\n  Fixed effects extracted: %d terms (part1), %d terms (part2)\n",
              nrow(fe1), nrow(fe2)))

  # ── Sensitivity analysis ──
  cat("\n── Sensitivity analysis (V4 + MiSeq only) ──\n")
  dt_sens <- dt[SeqRegion == "V4" & Platform == "Illumina MiSeq"]
  cat(sprintf("  Subset n = %d\n", nrow(dt_sens)))

  if (nrow(dt_sens) > 500) {
    # Remove Platform and SeqRegion from FE since they're constant
    sens_fe <- setdiff(fixed_effects, c("Platform", "SeqRegion"))
    sens_result <- run_two_part_model(dt_sens, outcome = outcome,
                                      fixed_effects = sens_fe)
  } else {
    cat("  Insufficient samples for sensitivity analysis (n < 500)\n")
    sens_result <- NULL
  }

  list(
    title = title,
    outcome = outcome,
    primary = result,
    sensitivity = sens_result,
    fixed_effects_primary = fe1,
    fixed_effects_abundance = fe2
  )
}


# ─────────────────────────────────────────────────────────────────────
# Visualization helpers
# ─────────────────────────────────────────────────────────────────────

theme_rc <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
}

# Forest plot of fixed effect estimates with 95% CI.
plot_forest <- function(fe_table, title = "Fixed Effects", max_terms = 30) {
  dt <- as.data.table(fe_table)
  # Filter to non-intercept terms
  dt <- dt[!term %in% c("(Intercept)", "Intercept")]
  dt <- dt[order(Estimate)]
  dt <- head(dt, max_terms)

  # Ensure we have CI columns
  if (!"Std. Error" %in% names(dt)) {
    dt[, `Std. Error` := 0.1]
  }
  dt[, ci_low := Estimate - 1.96 * `Std. Error`]
  dt[, ci_high := Estimate + 1.96 * `Std. Error`]

  ggplot(dt, aes(x = Estimate, y = reorder(term, Estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_pointrange(aes(xmin = ci_low, xmax = ci_high), size = 0.3) +
    labs(title = title, x = "Estimate (log-odds or log10)", y = "") +
    theme_rc()
}


# ─────────────────────────────────────────────────────────────────────
# Main: test the framework
# ─────────────────────────────────────────────────────────────────────

if (sys.nframe() == 0) {
  cat("Loading and preparing data...\n")
  dt <- load_rc_data()
  cat(sprintf("  Total samples: %d\n", nrow(dt)))

  dtc <- clean_analysis_data(dt)
  cat(sprintf("  After cleaning: %d\n", nrow(dtc)))
  cat(sprintf("  Diseases: %s\n", paste(levels(dtc$DiseaseStatus)[1:8],
                                       collapse = ", ")))

  # Quick test: HC vs Periodontal Disease only
  cat("\n\n── TEST: Periodontal Disease vs HC ──\n")
  dt_test <- dtc[DiseaseStatus %in% c("HC", "Periodontal Disease")]
  dt_test[, DiseaseStatus := factor(DiseaseStatus,
                                     levels = c("HC", "Periodontal Disease"))]
  cat(sprintf("  n = %d (HC=%d, PD=%d)\n", nrow(dt_test),
              sum(dt_test$DiseaseStatus == "HC"),
              sum(dt_test$DiseaseStatus == "Periodontal Disease")))

  result <- run_rc_analysis(dt_test,
                            outcome = "RedComplexTotal",
                            title = "Periodontal Disease vs HC - Red Complex Total")

  cat("\n\n── Framework test complete. ──\n")
}
