#!/home/software/anaconda3/envs/renv/bin/Rscript
#
# 03_multi_disease_landscape.R (P1.1) — RIGOR-ENHANCED v2
# ========================================================
# New in v2:
#   + Benjamini-Hochberg FDR correction across all disease tests
#   + Site-stratified analysis (saliva-only, tooth-only)
#   + Separation detection & flagging (Preterm, etc.)
#   + Forest plot uses FDR significance (not raw p)
#

source("02_unified_statistical_framework.R")

library(data.table)
library(ggplot2)

OUTDIR <- "../red_complex"

cat("\nLoading data...\n")
dt <- load_rc_data()
dtc <- clean_analysis_data(dt)
cat(sprintf("  Cleaned: %d samples\n", nrow(dtc)))

# ── Helper: run one disease vs HC ───────────────────────────────────

run_disease_model <- function(sub, disease_name, adjust_for_site = TRUE) {
  n_disease <- sum(sub$DiseaseStatus == disease_name)
  n_hc <- sum(sub$DiseaseStatus == "HC")
  results <- list()

  # Build covariate string
  covars <- "AgeGroup + Sex + Platform + SeqRegion"
  if (adjust_for_site) covars <- paste("Site +", covars)

  # Logistic model
  logit_ok <- tryCatch({
    formula_str <- paste("RedComplex_detected ~ DiseaseStatus +", covars)
    m <- glm(as.formula(formula_str), data = sub, family = binomial())
    coef_name <- paste0("DiseaseStatus", disease_name)
    if (coef_name %in% names(coef(m))) {
      s <- summary(m)$coefficients[coef_name, ]
      ci <- suppressMessages(confint.default(m)[coef_name, ])
      sep_flag <- (abs(s["Estimate"]) > 10) || (s["Std. Error"] > 50)
      data.table(
        disease = disease_name, n_disease = n_disease, n_hc = n_hc,
        model = "logistic",
        coef = s["Estimate"], se = s["Std. Error"],
        or_val = exp(s["Estimate"]),
        or_ci_low = exp(ci[1]), or_ci_high = exp(ci[2]),
        pvalue = s["Pr(>|z|)"],
        separation = sep_flag
      )
    }
  }, error = function(e) NULL)
  if (!is.null(logit_ok)) results$logistic <- logit_ok

  # Linear model (among detected)
  sub_det <- sub[RedComplex_detected == 1]
  if (nrow(sub_det) >= 50 && sum(sub_det$DiseaseStatus == disease_name) >= 20) {
    lm_ok <- tryCatch({
      sub_det[, log10_rc := log10(RedComplexTotal + 1e-6)]
      formula_str <- paste("log10_rc ~ DiseaseStatus +", covars)
      m <- lm(as.formula(formula_str), data = sub_det)
      coef_name <- paste0("DiseaseStatus", disease_name)
      if (coef_name %in% names(coef(m))) {
        s <- summary(m)$coefficients[coef_name, ]
        ci <- suppressMessages(confint.default(m)[coef_name, ])
        data.table(
          disease = disease_name,
          n_detected_disease = sum(sub_det$DiseaseStatus == disease_name),
          n_detected_hc = sum(sub_det$DiseaseStatus == "HC"),
          model = "lm_log10",
          coef = s["Estimate"], se = s["Std. Error"],
          fc = 10^s["Estimate"],
          fc_ci_low = 10^ci[1], fc_ci_high = 10^ci[2],
          pvalue = s["Pr(>|t|)"]
        )
      }
    }, error = function(e) NULL)
    if (!is.null(lm_ok)) results$lm <- lm_ok
  }

  return(results)
}

# ── Multi-disease (full dataset) ────────────────────────────────────

cat("\n── Adjusted Effect Sizes (full dataset) ──\n")

disease_counts <- dtc[, .N, by = DiseaseStatus]
valid_diseases <- disease_counts[N >= 50 & DiseaseStatus != "HC" & DiseaseStatus != ""]

logistic_results <- list()
lm_results <- list()

for (disease in as.character(valid_diseases$DiseaseStatus)) {
  sub <- dtc[DiseaseStatus %in% c("HC", disease)]
  sub[, DiseaseStatus := droplevels(DiseaseStatus)]
  sub[, DiseaseStatus := factor(DiseaseStatus, levels = c("HC", disease))]
  n_disease <- sum(sub$DiseaseStatus == disease)

  res <- run_disease_model(sub, disease)

  if (!is.null(res$logistic)) {
    cat(sprintf("  %-30s n=%-5d  OR=%5.2f [%5.2f-%5.2f]  p=%s",
                disease, n_disease,
                res$logistic$or_val, res$logistic$or_ci_low,
                res$logistic$or_ci_high,
                format.pval(res$logistic$pvalue, digits=2)))
    if (isTRUE(res$logistic$separation)) cat("  ⚠️ SEPARATION")
    cat("\n")
    logistic_results[[disease]] <- res$logistic
  }
  if (!is.null(res$lm)) {
    cat(sprintf("    %-28s (detected)  FC=%5.2f [%5.2f-%5.2f]  p=%s\n",
                "", res$lm$fc, res$lm$fc_ci_low, res$lm$fc_ci_high,
                format.pval(res$lm$pvalue, digits=2)))
    lm_results[[disease]] <- res$lm
  }
}

logistic_dt <- rbindlist(logistic_results, fill = TRUE)
lm_dt <- rbindlist(lm_results, fill = TRUE)

# ── FDR Correction ──────────────────────────────────────────────────

if (nrow(logistic_dt) > 0) {
  logistic_dt[, p_fdr := p.adjust(pvalue, method = "BH")]
  logistic_dt[, sig_fdr := p_fdr < 0.05]
}
if (nrow(lm_dt) > 0) {
  lm_dt[, p_fdr := p.adjust(pvalue, method = "BH")]
  lm_dt[, sig_fdr := p_fdr < 0.05]
}

cat(sprintf("\n  --- FDR correction (Benjamini-Hochberg) ---\n"))
cat(sprintf("  Logistic: %d/%d significant after FDR\n",
            sum(logistic_dt$sig_fdr, na.rm = TRUE), nrow(logistic_dt)))
cat(sprintf("  Abundance: %d/%d significant after FDR\n",
            sum(lm_dt$sig_fdr, na.rm = TRUE), nrow(lm_dt)))

# Flag diseases that lost significance after FDR
if (nrow(logistic_dt) > 0) {
  lost_fdr <- logistic_dt[pvalue < 0.05 & sig_fdr == FALSE]
  if (nrow(lost_fdr) > 0) {
    cat(sprintf("  ⚠️  Lost significance after FDR: %s\n",
                paste(lost_fdr$disease, collapse = ", ")))
  }
}

fwrite(logistic_dt, file.path(OUTDIR, "multi_disease_logistic_effects.csv"))
fwrite(lm_dt, file.path(OUTDIR, "multi_disease_abundance_effects.csv"))

# ── Site-stratified analysis ────────────────────────────────────────

cat("\n── Site-Stratified Analysis ──\n")

site_stratified <- list()
for (target_site in c("saliva", "tooth")) {
  cat(sprintf("\n  === %s-only ===\n", target_site))
  dt_site <- dtc[Site == target_site]
  site_diseases <- dt_site[, .N, by = DiseaseStatus][N >= 30 & DiseaseStatus != "HC"]

  site_logistic <- list()
  for (disease in as.character(site_diseases$DiseaseStatus)) {
    sub <- dt_site[DiseaseStatus %in% c("HC", disease)]
    sub[, DiseaseStatus := droplevels(DiseaseStatus)]
    sub[, DiseaseStatus := factor(DiseaseStatus, levels = c("HC", disease))]
    n_disease <- sum(sub$DiseaseStatus == disease)

    if (n_disease < 20) next

    res <- run_disease_model(sub, disease, adjust_for_site = FALSE)
    if (!is.null(res$logistic)) {
      res$logistic[, site_stratum := target_site]
      site_logistic[[paste(disease, target_site)]] <- res$logistic
      cat(sprintf("    %-28s n=%-4d OR=%5.2f [%5.2f-%5.2f] p=%s\n",
                  disease, n_disease,
                  res$logistic$or_val, res$logistic$or_ci_low,
                  res$logistic$or_ci_high,
                  format.pval(res$logistic$pvalue, digits=2)))
    }
  }

  if (length(site_logistic) > 0) {
    site_dt <- rbindlist(site_logistic, fill = TRUE)
    site_dt[, p_fdr := p.adjust(pvalue, method = "BH")]
    site_dt[, sig_fdr := p_fdr < 0.05]
    fname <- sprintf("multi_disease_logistic_%s.csv", target_site)
    fwrite(site_dt, file.path(OUTDIR, fname))
    cat(sprintf("    → %d/%d significant after FDR (saved: %s)\n",
                sum(site_dt$sig_fdr, na.rm = TRUE), nrow(site_dt), fname))
    site_stratified[[target_site]] <- site_dt
  }
}

# ── Forest Plot (FDR-corrected) ─────────────────────────────────────

cat("\n── Forest Plot ──\n")

# Filter to clean data with valid finite CIs (exclude separation cases)
plot_logistic <- logistic_dt[order(or_val)]
plot_logistic <- plot_logistic[is.finite(or_ci_low) & is.finite(or_ci_high) & or_ci_low > 1e-6 & or_ci_high < 1e6]
plot_logistic[, disease := factor(disease, levels = unique(disease))]
plot_logistic[, sig_label := ifelse(sig_fdr, "FDR<0.05", "FDR>=0.05")]

# Compute sensible axis limits
or_xmin <- max(0.005, min(plot_logistic$or_ci_low[is.finite(plot_logistic$or_ci_low)], na.rm = TRUE) * 0.5)
or_xmax <- min(200, max(plot_logistic$or_ci_high[is.finite(plot_logistic$or_ci_high)], na.rm = TRUE) * 2)

p1 <- ggplot(plot_logistic, aes(x = or_val, y = disease, color = sig_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = or_ci_low, xmax = or_ci_high), size = 0.4) +
  scale_x_log10(limits = c(or_xmin, or_xmax),
                breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16),
                labels = c("0.05", "0.1", "0.25", "0.5", "1", "2", "4", "8", "16")) +
  scale_color_manual(values = c("FDR<0.05" = "#c0392b", "FDR>=0.05" = "#95a5a6")) +
  labs(title = "Detection Odds Ratio vs HC (adjusted + FDR corrected)",
       subtitle = sprintf("%d diseases (separation cases excluded); red = FDR < 0.05", nrow(plot_logistic)),
       x = "Odds Ratio (log scale)", y = "", color = "") +
  theme_rc()

ggsave(file.path(OUTDIR, "multi_disease_forest_or.png"), p1,
       width = 10, height = max(5, nrow(plot_logistic) * 0.35), dpi = 150)

# FC plot — filter extreme CIs and set reasonable x-axis limits
plot_lm <- lm_dt[order(fc)]
# Remove rows with extreme or missing CIs (e.g., separation cases)
plot_lm <- plot_lm[is.finite(fc_ci_low) & is.finite(fc_ci_high) & fc_ci_low > 0.01 & fc_ci_high < 100]
plot_lm[, disease := factor(disease, levels = unique(disease))]
plot_lm[, sig_label := ifelse(sig_fdr, "FDR<0.05", "FDR>=0.05")]

# Determine axis limits with padding
x_min <- max(0.02, min(plot_lm$fc_ci_low, na.rm = TRUE) * 0.5)
x_max <- min(50, max(plot_lm$fc_ci_high, na.rm = TRUE) * 1.5)

p2 <- ggplot(plot_lm, aes(x = fc, y = disease, color = sig_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = fc_ci_low, xmax = fc_ci_high), size = 0.4) +
  scale_x_log10(limits = c(x_min, x_max),
                breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
                labels = c("0.1×", "0.25×", "0.5×", "1×", "2×", "4×", "8×")) +
  scale_color_manual(values = c("FDR<0.05" = "#2980b9", "FDR>=0.05" = "#95a5a6")) +
  labs(title = "Abundance Fold Change vs HC (among detected, FDR corrected)",
       subtitle = sprintf("n = %d diseases with valid CIs; extreme CI values clipped", nrow(plot_lm)),
       x = "Fold Change (log scale)", y = "", color = "") +
  theme_rc()

ggsave(file.path(OUTDIR, "multi_disease_forest_fc.png"), p2,
       width = 10, height = max(5, nrow(plot_lm) * 0.35), dpi = 150)

# ── Detection bar chart ─────────────────────────────────────────────

rates <- dtc[, .(
  n = .N,
  rc_detection_rate = mean(RedComplex_detected),
  all3_rate = mean(AllThree_detected)
), by = DiseaseStatus][order(-rc_detection_rate)]

rates[, `:=`(
  none_rate = 1 - rc_detection_rate,
  partial_rate = rc_detection_rate - all3_rate,
  disease_label = sprintf("%s (n=%d)", DiseaseStatus, n)
)]
rates[, disease_label := factor(disease_label,
        levels = rates[order(rc_detection_rate)]$disease_label)]

bar_data <- melt(rates, measure.vars = c("none_rate", "partial_rate", "all3_rate"),
                 variable.name = "detection_type", value.name = "proportion")
bar_data[, detection_type := factor(detection_type,
  levels = c("all3_rate", "partial_rate", "none_rate"),
  labels = c("All three", "1-2 members", "None detected"))]

p3 <- ggplot(bar_data, aes(x = proportion, y = disease_label, fill = detection_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("None detected" = "#ecf0f1",
                               "1-2 members" = "#f39c12",
                               "All three" = "#c0392b")) +
  labs(title = "Red Complex Detection by Disease Status",
       x = "Proportion of samples", y = "", fill = "") +
  theme_rc() + theme(axis.text.y = element_text(size = 8))

ggsave(file.path(OUTDIR, "multi_disease_detection_bar.png"), p3,
       width = 12, height = 8, dpi = 150)

# ── Summary ─────────────────────────────────────────────────────────

cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")
cat("TOP FINDINGS (FDR-corrected)\n")
cat(paste(rep("=", 60), collapse = "") %+% "\n")

cat("\n  Significant after FDR — Detection (logistic):\n")
sig_log <- logistic_dt[sig_fdr == TRUE][order(p_fdr)]
if (nrow(sig_log) > 0) {
  print(sig_log[, .(disease, OR = round(or_val, 2),
                    CI = sprintf("%.2f-%.2f", or_ci_low, or_ci_high),
                    p_FDR = format.pval(p_fdr, digits = 2))])
} else cat("    (none — all nominal findings may be false positives)\n")

cat("\n  Significant after FDR — Abundance (lm):\n")
sig_lm <- lm_dt[sig_fdr == TRUE][order(p_fdr)]
if (nrow(sig_lm) > 0) {
  print(sig_lm[, .(disease, FC = round(fc, 2),
                   CI = sprintf("%.2f-%.2f", fc_ci_low, fc_ci_high),
                   p_FDR = format.pval(p_fdr, digits = 2))])
} else cat("    (none)\n")

# Separation warnings
sep_dt <- logistic_dt[separation == TRUE]
if (nrow(sep_dt) > 0) {
  cat(sprintf("\n  ⚠️  Separation detected: %s (OR estimates unreliable, excluded from FDR ranking)\n",
              paste(sep_dt$disease, collapse = ", ")))
}

cat("\nDone.\n")
