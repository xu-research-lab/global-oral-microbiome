#!/home/software/anaconda3/envs/renv/bin/Rscript
#
# 12_species_level_analysis.R (NEW)
# ==================================
# Addresses the "any member vs. all three" validity concern.
#
# Core questions:
#   1. Is the disease-RC association driven by one dominant species?
#   2. Does "any detected" OR differ from "all three detected" OR?
#   3. Does the number of co-occurring members (0→1→2→3) scale with disease?
#
# Methods:
#   A. Species-specific binary logistic (Pg / Td / Tf + AllThree as outcomes)
#   B. Proportional odds model (RC_member_count 0-3 as ordinal outcome)
#   C. Proportional odds assumption test (Brant test)
#   D. Side-by-side OR comparison across all 5 outcome definitions
#

source("02_unified_statistical_framework.R")
library(data.table)
library(ggplot2)

# For proportional odds model — use MASS::polr (available in base R)
# or ordinal::clm. We'll implement manually to avoid dependency issues.
OUTDIR <- "../red_complex"

cat("Loading data...\n")
dt <- load_rc_data()
dtc <- clean_analysis_data(dt)
cat(sprintf("Cleaned: %d samples\n", nrow(dtc)))

# ── Helper: run binary logistic for a given detection outcome ───
run_binary_logistic <- function(sub, disease, outcome_col = "RedComplex_detected") {
  sub[, y := get(outcome_col)]
  m <- tryCatch({
    glm(y ~ DiseaseStatus + Site + AgeGroup + Sex + Platform + SeqRegion,
        data = sub, family = binomial())
  }, error = function(e) NULL)
  if (is.null(m)) return(NULL)

  coef_name <- paste0("DiseaseStatus", disease)
  s <- summary(m)$coefficients[coef_name, , drop = FALSE]
  if (nrow(s) == 0) return(NULL)
  ci <- tryCatch(confint.default(m)[coef_name, ], error = function(e) c(NA, NA))

  data.table(
    disease = disease,
    n_disease = sum(sub$DiseaseStatus == disease),
    n_hc = sum(sub$DiseaseStatus == "HC"),
    outcome = outcome_col,
    OR = exp(s[1, "Estimate"]),
    OR_lo = exp(ci[1]),
    OR_hi = exp(ci[2]),
    pvalue = s[1, "Pr(>|z|)"]
  )
}

# ── A. Species-specific disease associations ────────────────────
cat("\n══════════════════════════════════════════════════════\n")
cat("A. Species-Specific Detection OR by Disease\n")
cat("══════════════════════════════════════════════════════\n\n")

outcomes <- c("RedComplex_detected", "P_gingivalis_detected",
              "T_denticola_detected", "T_forsythia_detected", "AllThree_detected")
outcome_labels <- c("Any RC", "P. gingivalis", "T. denticola", "T. forsythia", "All Three")

disease_counts <- dtc[, .N, by = DiseaseStatus][N >= 50 & DiseaseStatus != "HC" & DiseaseStatus != ""]
valid_diseases <- as.character(disease_counts$DiseaseStatus)

# Run all models
species_results <- rbindlist(lapply(valid_diseases, function(d) {
  sub <- dtc[DiseaseStatus %in% c("HC", d)]
  sub[, DiseaseStatus := factor(DiseaseStatus, levels = c("HC", d))]
  rbindlist(lapply(outcomes, function(oc) {
    run_binary_logistic(sub, d, oc)
  }), fill = TRUE)
}))

# FDR correction within each outcome
for (oc in outcomes) {
  idx <- species_results$outcome == oc
  if (sum(idx) > 1) {
    species_results[idx, p_fdr := p.adjust(pvalue, method = "BH")]
    species_results[idx, sig_fdr := p_fdr < 0.05]
  }
}

# Print comparison table for key diseases
cat("Side-by-side OR comparison (all outcomes FDR-corrected):\n\n")
for (d in valid_diseases) {
  sub <- species_results[disease == d]
  if (nrow(sub) < 2) next
  max_or <- sub[outcome == "RedComplex_detected"]$OR
  if (length(max_or) == 0 || is.na(max_or) || max_or < 1.5) next  # only strong associations

  cat(sprintf("  %s (any-RC OR=%.2f):\n", d, max_or))
  print(sub[, .(
    Outcome = outcome_labels[match(outcome, outcomes)],
    OR = round(OR, 2),
    CI = sprintf("%.2f–%.2f", OR_lo, OR_hi),
    p_FDR = format.pval(p_fdr, digits = 2),
    Sig = ifelse(sig_fdr, "*", "")
  )])
  cat("\n")
}

fwrite(species_results, file.path(OUTDIR, "species_level_or_comparison.csv"))

# ── B. Proportional odds model (ordinal logistic) ──────────────
cat("══════════════════════════════════════════════════════\n")
cat("B. Proportional Odds Model (RC_member_count 0→3)\n")
cat("══════════════════════════════════════════════════════\n\n")

cat("What this model does:\n")
cat("  logit(P(Y ≥ k)) = α_k + β·Disease + covariates, k = 1,2,3\n")
cat("  exp(β) = cumulative OR for being at a HIGHER member count\n")
cat("  Assumption: β is the SAME across all k thresholds\n\n")

# RC_member_count distribution per disease
cat("RC_member_count distribution:\n")
rc_dist <- dtc[, .N, by = .(DiseaseStatus, RC_member_count)][order(DiseaseStatus, RC_member_count)]
rc_dist[, pct := round(100 * N / sum(N), 1), by = DiseaseStatus]
dist_wide <- dcast(rc_dist, DiseaseStatus ~ RC_member_count, value.var = "pct", sep = "_")
cat("  (row%):\n")
print(dist_wide)
fwrite(rc_dist, file.path(OUTDIR, "species_rc_member_count_distribution.csv"))

# Proportional odds model using polr from MASS (base R package)
library(MASS)

ordinal_results <- rbindlist(lapply(valid_diseases, function(d) {
  sub <- dtc[DiseaseStatus %in% c("HC", d)]
  sub[, DiseaseStatus := factor(DiseaseStatus, levels = c("HC", d))]
  n_d <- sum(sub$DiseaseStatus == d)
  if (n_d < 30) return(NULL)
  # Skip if the disease has no variation in RC_member_count (all same level)
  if (uniqueN(sub$RC_member_count) < 2) return(NULL)

  # Proportional odds model on RC_member_count
  m <- tryCatch({
    polr(as.ordered(RC_member_count) ~ DiseaseStatus + Site + AgeGroup +
           Sex + Platform + SeqRegion,
         data = sub, method = "logistic", Hess = TRUE)
  }, error = function(e) NULL)
  if (is.null(m)) return(NULL)

  coef_name <- paste0("DiseaseStatus", d)
  coefs <- as.data.table(coef(summary(m)), keep.rownames = "term")
  if (!coef_name %in% coefs$term) return(NULL)

  s_row <- coefs[term == coef_name]
  ci <- tryCatch(confint.default(m)[coef_name, ], error = function(e) c(NA, NA))
  # polr gives t-value, not p-value — compute via normal approximation
  pval <- 2 * pnorm(-abs(s_row$`t value`))

  data.table(
    disease = d,
    n_disease = n_d,
    n_hc = sum(sub$DiseaseStatus == "HC"),
    or_val = exp(s_row$Value),
    or_ci_low = exp(ci[1]),
    or_ci_high = exp(ci[2]),
    pvalue = pval,
    n_total = nrow(sub)
  )
}))

if (nrow(ordinal_results) > 0) {
  ordinal_results[, p_fdr := p.adjust(pvalue, method = "BH")]
  ordinal_results[, sig_fdr := p_fdr < 0.05]
  setorder(ordinal_results, -or_val)

  cat("\nProportional odds model results:\n")
  cat("  OR = cumulative odds for higher RC_member_count in disease vs HC\n\n")
  print(ordinal_results[, .(
    disease,
    N_disease = n_disease,
    Cum_OR = round(or_val, 2),
    `95% CI` = sprintf("%.2f–%.2f", or_ci_low, or_ci_high),
    p = format.pval(pvalue, digits = 2),
    p_FDR = format.pval(p_fdr, digits = 2),
    Sig = ifelse(sig_fdr, "*", "")
  )])

  fwrite(ordinal_results, file.path(OUTDIR, "species_proportional_odds_results.csv"))

  # ── C. Test proportional odds assumption ────────────────────
  cat("\n\nC. Testing the Proportional Odds Assumption\n")
  cat("  (Does the same β apply across all thresholds?)\n\n")

  # Brant test approach: fit separate binary models for each threshold
  # and compare coefficients
  brant_test <- rbindlist(lapply(valid_diseases[1:min(6, length(valid_diseases))], function(d) {
    sub <- dtc[DiseaseStatus %in% c("HC", d)]
    sub[, DiseaseStatus := factor(DiseaseStatus, levels = c("HC", d))]
    n_d <- sum(sub$DiseaseStatus == d)
    if (n_d < 30) return(NULL)

    # Three binary models: y≥1, y≥2, y=3
    thresholds <- c(">=1", ">=2", ">=3")
    th_inds <- list(1:3, 2:3, 3:3)
    coefs <- rbindlist(lapply(seq_along(thresholds), function(i) {
      sub[, y := as.integer(RC_member_count %in% th_inds[[i]])]
      m <- tryCatch({
        glm(y ~ DiseaseStatus + Site + AgeGroup + Sex + Platform + SeqRegion,
            data = sub, family = binomial())
      }, error = function(e) NULL)
      if (is.null(m)) return(NULL)
      s <- coef(summary(m))
      cn <- paste0("DiseaseStatus", d)
      if (!cn %in% rownames(s)) return(NULL)
      data.table(threshold = thresholds[i], beta = s[cn, "Estimate"], se = s[cn, "Std. Error"])
    }))
    if (nrow(coefs) < 3) return(NULL)
    # Test for heterogeneity of β across thresholds
    betas <- coefs$beta
    ses <- coefs$se
    # Cochran's Q test for heterogeneity
    w <- 1 / ses^2
    beta_w <- sum(w * betas) / sum(w)
    Q <- sum(w * (betas - beta_w)^2)
    brant_p <- pchisq(Q, 2, lower.tail = FALSE)

    data.table(
      disease = d,
      beta_ge1 = betas[1], beta_ge2 = betas[2], beta_eq3 = betas[3],
      se_ge1 = ses[1], se_ge2 = ses[2], se_eq3 = ses[3],
      brant_Q = Q, brant_p = brant_p,
      po_assumption = ifelse(brant_p > 0.05,
        "OK (same β across thresholds)",
        "VIOLATED (β differs across thresholds)"
      )
    )
  }))

  if (nrow(brant_test) > 0) {
    cat("  Brant-like heterogeneity test:\n")
    cat("  (Compares β for Y≥1, Y≥2, Y=3 — null: they are equal)\n\n")
    brant_test[, `:=`(
      beta_ge1 = round(beta_ge1, 3),
      beta_ge2 = round(beta_ge2, 3),
      beta_eq3 = round(beta_eq3, 3),
      brant_Q = round(brant_Q, 2),
      brant_p = format.pval(brant_p, digits = 2)
    )]
    print(brant_test[, .(disease, beta_ge1, beta_ge2, beta_eq3, brant_Q, brant_p, po_assumption)])
    fwrite(brant_test, file.path(OUTDIR, "species_proportional_odds_brant_test.csv"))
  }
}

# ── D. Side-by-side OR comparison: 5 outcomes + proportional odds ─
cat("\n\n══════════════════════════════════════════════════════\n")
cat("D. Integrated Comparison Across All Outcome Definitions\n")
cat("══════════════════════════════════════════════════════\n\n")

# Merge all results
binary_wide <- dcast(species_results[, .(disease, outcome, OR, OR_lo, OR_hi, sig_fdr)],
                     disease ~ outcome, value.var = c("OR", "OR_lo", "OR_hi", "sig_fdr"),
                     sep = "__")
if (nrow(ordinal_results) > 0) {
  ordinal_wide <- ordinal_results[, .(disease, polr_OR = or_val, polr_lo = or_ci_low,
                                       polr_hi = or_ci_high, polr_sig = sig_fdr)]
  combined <- merge(binary_wide, ordinal_wide, by = "disease", all = TRUE)
} else {
  combined <- binary_wide
}

# For top diseases, print a clean comparison
cat("  OR Comparison Table (top diseases with any-RC OR > 1.5):\n\n")
comp_diseases <- species_results[outcome == "RedComplex_detected" & OR > 1.5]$disease
for (d in unique(species_results$disease)) {
  sub <- species_results[disease == d]
  polr_row <- ordinal_results[disease == d]
  any_or <- sub[outcome == "RedComplex_detected"]$OR
  if (length(any_or) == 0 || is.na(any_or)) next
  if (any_or < 1.5 && any_or > 0.67) next

  # Format the summary line
  parts <- c()
  for (i in seq_along(outcomes)) {
    v <- sub[outcome == outcomes[i]]
    if (nrow(v) > 0 && !is.na(v$OR)) {
      parts <- c(parts, sprintf("%s: %.2f%s",
        outcome_labels[i], v$OR, ifelse(v$sig_fdr, "*", "")))
    }
  }
  if (nrow(polr_row) > 0 && !is.na(polr_row$or_val)) {
    parts <- c(parts, sprintf("Ordinal: %.2f%s", polr_row$or_val, ifelse(polr_row$sig_fdr, "*", "")))
  }

  cat(sprintf("  %-25s → %s\n", d, paste(parts, collapse = " | ")))
}

cat("\n\n──────────────────────────────────────────────────\n")
cat("INTERPRETATION GUIDE\n")
cat("──────────────────────────────────────────────────\n")
cat("  If any-RC OR ≈ all-three OR:\n")
cat("    → Disease association is DRIVEN by true co-occurrence\n")
cat("  If any-RC OR >> all-three OR:\n")
cat("    → Association is driven by SINGLE-MEMBER presence\n")
cat("  If proportional odds model OR ≈ any-RC OR:\n")
cat("    → Each additional member increases disease odds equally\n")
cat("  If proportional odds OR < any-RC OR:\n")
cat("    → Most of the effect is at the 0→1 step (marginal colonization)\n")
cat("──────────────────────────────────────────────────\n")

# ── E. Visualization ──────────────────────────────────────────
cat("\n\nE. Visualizations\n")

# Forest plot: compare ORs across outcomes for top diseases
top_diseases <- unique(species_results[outcome == "RedComplex_detected"][order(-OR)]$disease)
top_diseases <- top_diseases[1:min(10, length(top_diseases))]

plot_data <- species_results[disease %in% top_diseases]
plot_data[, outcome_label := factor(
  outcome_labels[match(outcome, outcomes)],
  levels = outcome_labels)]
plot_data[, disease := factor(disease, levels = top_diseases)]

# Add proportional odds
if (nrow(ordinal_results) > 0 && any(ordinal_results$disease %in% top_diseases)) {
  polr_plot <- ordinal_results[disease %in% top_diseases]
  polr_plot[, outcome_label := "Ordinal (cumulative)"]
  polr_plot[, y_OR := or_val]
  polr_plot[, y_lo := or_ci_low]
  polr_plot[, y_hi := or_ci_high]
  polr_plot[, sig_fdr := sig_fdr]
  polr_plot[, disease := factor(disease, levels = top_diseases)]
  # Add to plot_data
  plot_data <- rbind(
    plot_data[, .(disease, outcome_label, OR, OR_lo, OR_hi, sig_fdr)],
    polr_plot[, .(disease, outcome_label, OR = or_val, OR_lo = or_ci_low,
                   OR_hi = or_ci_high, sig_fdr)]
  )
}

p1 <- ggplot(plot_data, aes(x = OR, y = outcome_label, color = sig_fdr)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = OR_lo, xmax = OR_hi), size = 0.4, position = position_dodge(width = 0.3)) +
  scale_x_log10() +
  facet_wrap(~ disease, scales = "free_x", ncol = 2) +
  scale_color_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#95a5a6", "NA" = "#95a5a6")) +
  labs(title = "Disease-Red Complex Association: By Outcome Definition",
       subtitle = "Each panel = one disease; colors = FDR < 0.05",
       x = "Odds Ratio (log scale)", y = "Outcome definition", color = "FDR<0.05") +
  theme_rc() +
  theme(axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 8))

ggsave(file.path(OUTDIR, "species_level_forest_by_outcome.png"), p1,
       width = 12, height = max(5, length(top_diseases) * 1.5), dpi = 150)
cat("  Saved: species_level_forest_by_outcome.png\n")

# RC_member_count distribution by disease (stacked bar)
dist_plot <- rc_dist
dist_plot[, DiseaseStatus := factor(DiseaseStatus,
  levels = names(sort(xtabs(N ~ DiseaseStatus, rc_dist), decreasing = TRUE)))]
dist_plot[, RC_member_count := factor(RC_member_count)]

p2 <- ggplot(dist_plot[DiseaseStatus %in% names(sort(table(dist_plot$DiseaseStatus),
                                                       decreasing = TRUE))[1:15]],
             aes(x = DiseaseStatus, y = N, fill = RC_member_count)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("0" = "#ecf0f1", "1" = "#fdd49e", "2" = "#fd8d3c", "3" = "#c0392b")) +
  labs(title = "Red Complex Member Count Distribution by Disease",
       subtitle = "0 = none, 1 = single, 2 = double, 3 = all three",
       x = "", y = "Proportion", fill = "Member count") +
  theme_rc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTDIR, "species_member_count_distribution.png"), p2,
       width = 10, height = 6, dpi = 150)
cat("  Saved: species_member_count_distribution.png\n")

cat("\nDone.\n")
