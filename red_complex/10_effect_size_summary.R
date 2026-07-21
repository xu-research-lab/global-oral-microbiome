#!/home/software/anaconda3/envs/renv/bin/Rscript
#
# 10_effect_size_summary.R (NEW, optimized)
# ==========================================
# Consolidate effect sizes from existing output CSVs
# Rank by magnitude not p-value, flag large-N significance trap
#

library(data.table); library(ggplot2)
source("02_unified_statistical_framework.R")  # for theme_rc
OUTDIR <- "../red_complex"

# ── Helper: classify effect size ────────────────────────────────
grade_or <- function(or) fcase(
  is.na(or), NA_character_, or >= 3 | or <= 0.33, "Large",
  or >= 2 | or <= 0.5, "Moderate", or >= 1.5 | or <= 0.67, "Small", default = "Negligible")
grade_icc <- function(icc) fcase(
  is.na(icc), NA_character_, icc >= 0.75, "Large", icc >= 0.5, "Moderate",
  icc >= 0.25, "Small", default = "Negligible")
grade_auc <- function(auc) fcase(
  is.na(auc), NA_character_, auc >= 0.9, "Large", auc >= 0.8, "Moderate",
  auc >= 0.7, "Small", default = "Negligible")
grade_phi <- function(phi) fcase(
  is.na(phi), NA_character_, phi >= 0.5, "Large", phi >= 0.3, "Moderate",
  phi >= 0.1, "Small", default = "Negligible")

all_effects <- list()

# ── 1. Disease OR from Fix 1 CSV ────────────────────────────────
f1 <- file.path(OUTDIR, "multi_disease_logistic_effects.csv")
if (file.exists(f1)) {
  d1 <- fread(f1)
  effects1 <- d1[, .(
    module = "Disease Assoc.", finding = paste("Detection OR:", disease),
    metric = "OR", value = or_val, lo = or_ci_low, hi = or_ci_high,
    pvalue = pvalue, p_fdr = p_fdr, magnitude = grade_or(or_val))]
  all_effects[["disease_logistic"]] <- effects1[!is.na(value)]
}

f1b <- file.path(OUTDIR, "multi_disease_abundance_effects.csv")
if (file.exists(f1b)) {
  d1b <- fread(f1b)
  effects1b <- d1b[, .(
    module = "Disease Assoc.", finding = paste("Abundance FC:", disease),
    metric = "FC", value = fc, lo = fc_ci_low, hi = fc_ci_high,
    pvalue = pvalue, p_fdr = p_fdr, magnitude = grade_or(fc))]
  all_effects[["disease_abundance"]] <- effects1b[!is.na(value)]
}

# ── 2. Longitudinal ICC from Fix 2 CSV ─────────────────────────
f2 <- file.path(OUTDIR, "longitudinal_icc_mixed_model.csv")
if (file.exists(f2)) {
  d2 <- fread(f2)
  effects2 <- d2[, .(
    module = "Longitudinal", finding = paste("ICC:", measure),
    metric = "ICC", value = ICC, lo = NA_real_, hi = NA_real_,
    pvalue = NA_real_, p_fdr = NA_real_, magnitude = grade_icc(ICC))]
  all_effects[["icc"]] <- effects2[!is.na(value)]
}

# ── 3. Diagnostic AUC/AUPRC from Fix 3 CSV ─────────────────────
f3 <- file.path(OUTDIR, "roc_site_stratified_full.csv")
if (file.exists(f3)) {
  d3 <- fread(f3)
  effects3a <- d3[, .(
    module = "Diagnostic", finding = paste("PD AUC:", Site, "-", Predictor),
    metric = "AUC", value = AUC, lo = AUC_low, hi = AUC_high,
    pvalue = NA_real_, p_fdr = NA_real_, magnitude = grade_auc(AUC))]
  effects3b <- d3[!is.na(AUPRC), .(
    module = "Diagnostic", finding = paste("PD AUPRC:", Site, "-", Predictor),
    metric = "AUPRC", value = AUPRC, lo = NA_real_, hi = NA_real_,
    pvalue = NA_real_, p_fdr = NA_real_, magnitude = fcase(
      AUPRC >= 0.5, "Large", AUPRC >= 0.3, "Moderate", AUPRC >= 0.1, "Small", default = "Negligible"))]
  all_effects[["auc"]] <- rbind(effects3a, effects3b, fill = TRUE)[!is.na(value)]
}

# ── 4. Co-occurrence Phi from Fix 5 CSV ────────────────────────
f4 <- file.path(OUTDIR, "cooccurrence_global_pairs.csv")
if (file.exists(f4)) {
  d4 <- fread(f4)
  effects4 <- d4[, .(
    module = "Co-occurrence", finding = paste("Phi:", species1, "-", species2),
    metric = "Phi", value = phi, lo = NA_real_, hi = NA_real_,
    pvalue = p_fisher, p_fdr = NA_real_, magnitude = grade_phi(phi))]
  all_effects[["phi"]] <- effects4[!is.na(value)]
}

# CMH-adjusted OR
f4b <- file.path(OUTDIR, "cooccurrence_mh_adjusted.csv")
if (file.exists(f4b)) {
  d4b <- fread(f4b)
  effects4b <- d4b[, .(
    module = "Co-occurrence", finding = paste("CMH OR:", species1, "-", species2),
    metric = "OR_cmh", value = cmh_OR, lo = NA_real_, hi = NA_real_,
    pvalue = cmh_p, p_fdr = NA_real_, magnitude = grade_or(cmh_OR))]
  all_effects[["cmh"]] <- effects4b[!is.na(value)]
}

# ── 5. Geography OR from Fix 4 CSV ─────────────────────────────
f5 <- file.path(OUTDIR, "geography_adjusted_or.csv")
if (file.exists(f5)) {
  d5 <- fread(f5)
  d5_geo <- d5[!grepl("ref", term)]
  if (nrow(d5_geo) > 0) {
    effects5 <- d5_geo[, .(
      module = "Geography", finding = paste("OR:", gsub("Geography", "", term)),
      metric = "OR", value = OR, lo = OR_ci_low, hi = OR_ci_high,
      pvalue = pvalue, p_fdr = NA_real_, magnitude = grade_or(OR))]
    all_effects[["geo"]] <- effects5[!is.na(value)]
  }
}

# ── 6. Age effect from Fix 6 data ──────────────────────────────
f6 <- file.path(OUTDIR, "age_site_interaction_test.csv")
if (file.exists(f6)) {
  d6 <- fread(f6)
  effects6 <- d6[, .(
    module = "Age Trajectory", finding = paste("Age×Site interaction"),
    metric = "p-value", value = pvalue, lo = NA_real_, hi = NA_real_,
    pvalue = pvalue, p_fdr = NA_real_,
    magnitude = fcase(pvalue < 0.001, "Large", pvalue < 0.01, "Moderate", pvalue < 0.05, "Small", default = "NS"))]
  all_effects[["age"]] <- effects6[!is.na(value)]
}

# ── 7. Combine and rank ─────────────────────────────────────────
effects_dt <- rbindlist(all_effects, fill = TRUE)
effects_dt <- effects_dt[!is.na(value)]

# Magnitude score
effects_dt[, mag_score := fcase(
  magnitude == "Large", 3, magnitude == "Moderate", 2,
  magnitude == "Small", 1, default = 0)]

# Evidence grade
effects_dt[, evidence_grade := fcase(
  magnitude == "Large" & (!is.na(pvalue) & pvalue < 0.001), "A+",
  magnitude == "Large" & (!is.na(pvalue) & pvalue < 0.05), "A",
  magnitude == "Large" & is.na(pvalue), "A (no p)",
  magnitude == "Moderate" & (!is.na(pvalue) & pvalue < 0.001), "A-",
  magnitude == "Moderate" & (!is.na(pvalue) & pvalue < 0.05), "B+",
  magnitude == "Moderate" & is.na(pvalue), "B (no p)",
  magnitude == "Small" & (!is.na(pvalue) & pvalue < 0.05), "C+",
  magnitude == "Small", "C", default = "D")]

setorder(effects_dt, -mag_score, pvalue, na.last = TRUE)

cat(sprintf("\n  Total findings: %d across %d modules\n", nrow(effects_dt), uniqueN(effects_dt$module)))

cat("\n  TOP 25 FINDINGS BY EFFECT SIZE:\n")
cat(paste(rep("-", 80), collapse = "") %+% "\n")
top25 <- head(effects_dt, 25)
print(top25[, .(
  Grade = evidence_grade, Module = module, Finding = finding,
  Metric = metric, Value = round(value, 3), Magnitude = magnitude,
  p = ifelse(is.na(pvalue), "-", format.pval(pvalue, digits = 2))
)])

fwrite(effects_dt, file.path(OUTDIR, "effect_size_summary_all.csv"))

# ── 8. Large-N trap ────────────────────────────────────────────
trap <- effects_dt[!is.na(pvalue) & pvalue < 0.05 &
                     magnitude %in% c("Negligible", "Small") &
                     value > 0.5 & value < 2.0]
cat(sprintf("\n  Findings p<0.05 but small/negligible: %d\n", nrow(trap)))
if (nrow(trap) > 0) print(trap[, .(finding, value = round(value, 3), magnitude, p = format.pval(pvalue, digits = 2))])

# ── 9. Grade distribution ──────────────────────────────────────
grade_s <- effects_dt[, .N, by = evidence_grade][order(evidence_grade)]
cat("\n  Evidence grade distribution:\n"); print(grade_s)

module_s <- effects_dt[, .(N = .N, best = evidence_grade[1]), by = module][order(module)]
cat("\n  Per-module summary:\n"); print(module_s)

fwrite(grade_s, file.path(OUTDIR, "effect_size_evidence_grades.csv"))
fwrite(module_s, file.path(OUTDIR, "effect_size_module_summary.csv"))

# ── 10. Plots ──────────────────────────────────────────────────
plot_eff <- head(effects_dt, 30)
plot_eff[, fshort := substr(finding, 1, 55)]

p1 <- ggplot(plot_eff, aes(x = value, y = reorder(fshort, mag_score))) +
  geom_point(aes(color = evidence_grade), size = 3) +
  geom_text(aes(label = sprintf("%.2f", value)), hjust = -0.3, size = 2.5) +
  scale_color_manual(values = c("A+"="#c0392b","A"="#e74c3c","A (no p)"="#e67e22",
    "A-"="#e67e22","B+"="#f39c12","B (no p)"="#2ecc71",
    "C+"="#95a5a6","C"="#bdc3c7","D"="#ecf0f1")) +
  labs(title="Top Findings by Effect Size (Not p-value)", x="Effect Size", y="", color="Grade") +
  theme_rc() + theme(axis.text.y = element_text(size = 6))
ggsave(file.path(OUTDIR, "effect_size_ranking.png"), p1, width = 12, height = 9, dpi = 150)

p2 <- ggplot(grade_s, aes(x = evidence_grade, y = N, fill = evidence_grade)) +
  geom_bar(stat = "identity") +
  labs(title = "Evidence Grade Distribution",
       subtitle = sprintf("%d findings, %d modules", nrow(effects_dt), uniqueN(effects_dt$module)),
       x = "Grade", y = "Count") + theme_rc() + theme(legend.position = "none")
ggsave(file.path(OUTDIR, "effect_size_grade_distribution.png"), p2, width = 8, height = 5, dpi = 150)

cat(sprintf("\n  Saved: effect_size_ranking.png, effect_size_grade_distribution.png\nDone.\n"))
