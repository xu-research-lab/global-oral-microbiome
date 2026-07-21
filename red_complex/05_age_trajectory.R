#!/home/software/anaconda3/envs/renv/bin/Rscript
#
# 05_age_trajectory.R (P1.3) — RIGOR-ENHANCED v2
# ================================================
# New in v2:
#   + Formal age × site interaction test (GAM factor-smooth interaction)
#   + ANOVA model comparison (with vs without interaction)
#   + Fixed HC/PD comparison: shared reference levels, consistent formula
#   + Age×disease interaction test
#

source("02_unified_statistical_framework.R")

library(data.table)
library(ggplot2)
library(mgcv)

OUTDIR <- "../red_complex"

cat("Loading data...\n")
dt <- load_rc_data()
dtc <- clean_analysis_data(dt)
dtc[, Age := as.numeric(Age)]
dtc <- dtc[!is.na(Age) & Age > 0 & Age < 100]
cat(sprintf("  With valid age: %d samples\n", nrow(dtc)))

# ── 1. Age distribution summary (HC) ────────────────────────────
cat("\n── 1. Age Groups (HC only) ──\n")
hc <- dtc[DiseaseStatus == "HC"]
hc_age <- hc[, .(n = .N, rc_detection = mean(RedComplex_detected),
  all3_detection = mean(AllThree_detected), rc_mean = mean(RedComplexTotal)),
  by = AgeGroup][order(AgeGroup)]
print(hc_age)
fwrite(hc_age, file.path(OUTDIR, "age_trajectory_hc_summary.csv"))

# ── 2. GAM smooth: age vs RC in HC ──────────────────────────────
cat("\n── 2. GAM: Age vs Red Complex in HC ──\n")
cat("  Fitting GAM for detection probability...\n")
gam_detect <- gam(RedComplex_detected ~ s(Age, bs = "tp") + Site + Platform + SeqRegion + Sex,
                  data = hc, family = binomial(), method = "REML")
cat(sprintf("  Detection GAM deviance explained: %.1f%%\n", summary(gam_detect)$dev.expl * 100))

hc_det <- hc[RedComplex_detected == 1]
cat("  Fitting GAM for log10 abundance (detected only)...\n")
gam_abund <- gam(log10(RedComplexTotal + 1e-6) ~ s(Age, bs = "tp") + Site + Platform + SeqRegion + Sex,
                 data = hc_det, method = "REML")
cat(sprintf("  Abundance GAM deviance explained: %.1f%%\n", summary(gam_abund)$dev.expl * 100))

# Predictions
age_grid <- data.table(Age = seq(0, 95, by = 0.5))
age_grid[, `:=`(Site = factor("saliva", levels = levels(hc$Site)),
  Platform = factor("Illumina MiSeq", levels = levels(hc$Platform)),
  SeqRegion = factor("V4", levels = levels(hc$SeqRegion)),
  Sex = factor("F", levels = levels(hc$Sex)))]
pred_detect <- predict(gam_detect, newdata = age_grid, type = "link", se.fit = TRUE)
age_grid[, `:=`(detect_logit = pred_detect$fit, detect_se = pred_detect$se.fit,
  detect_prob = plogis(pred_detect$fit),
  detect_low = plogis(pred_detect$fit - 1.96 * pred_detect$se.fit),
  detect_high = plogis(pred_detect$fit + 1.96 * pred_detect$se.fit))]
pred_abund <- predict(gam_abund, newdata = age_grid, type = "link", se.fit = TRUE)
age_grid[, `:=`(abund_log10 = pred_abund$fit, abund_se = pred_abund$se.fit,
  abund_low = 10^(pred_abund$fit - 1.96 * pred_abund$se.fit),
  abund_high = 10^(pred_abund$fit + 1.96 * pred_abund$se.fit))]
fwrite(age_grid, file.path(OUTDIR, "age_trajectory_gam_predictions.csv"))

# ── 3. AGE × SITE INTERACTION (NEW) ─────────────────────────────
cat("\n── 3. Age × Site Interaction Test ──\n")

site_age_n <- hc[, .(n = .N, age_range = max(Age) - min(Age)), by = Site][n >= 200 & age_range >= 50]
major_sites <- as.character(site_age_n$Site)
cat(sprintf("  Sites for interaction test: %s\n", paste(major_sites, collapse = ", ")))

hc_interact <- hc[Site %in% major_sites]
hc_interact[, Site := droplevels(Site)]

# Model 1: No interaction (smooth age effect, same across sites)
cat("\n  Fitting model WITHOUT interaction...\n")
m_no_int <- gam(RedComplex_detected ~ s(Age, bs = "tp") + Site + Platform + SeqRegion + Sex,
                data = hc_interact, family = binomial(), method = "REML")
cat(sprintf("  No-interaction model: deviance explained = %.1f%%\n", summary(m_no_int)$dev.expl * 100))

# Model 2: Factor-smooth interaction
cat("  Fitting model WITH age×site factor-smooth interaction...\n")
m_with_int <- tryCatch({
  gam(RedComplex_detected ~ s(Age, bs = "tp", by = Site) + Site + Platform + SeqRegion + Sex,
      data = hc_interact, family = binomial(), method = "REML")
}, error = function(e) {
  cat(sprintf("  Factor-smooth interaction failed: %s\n", e$message))
  cat("  Trying ordered factor approach...\n")
  hc_interact[, Site_ordered := factor(Site, ordered = TRUE)]
  gam(RedComplex_detected ~ s(Age, bs = "tp") + s(Age, bs = "tp", by = Site_ordered) +
      Site + Platform + SeqRegion + Sex,
      data = hc_interact, family = binomial(), method = "REML")
})

aov_p <- NA
if (!is.null(m_with_int)) {
  cat(sprintf("  With-interaction model: deviance explained = %.1f%%\n",
              summary(m_with_int)$dev.expl * 100))
  anova_result <- anova(m_no_int, m_with_int, test = "Chisq")
  cat("\n  ANOVA: model comparison (no interaction vs with interaction):\n")
  print(anova_result)
  aov_p <- anova_result[2, "Pr(>Chi)"]
  cat(sprintf("\n  Age×Site interaction significance: p = %s\n", format.pval(aov_p, digits = 4)))
  if (aov_p < 0.05) {
    cat("  → SIGNIFICANT: age trajectory differs by oral site\n")
  } else {
    cat("  → NOT significant: age trajectory is consistent across oral sites\n")
  }

  int_test <- data.table(
    test = "Age x Site interaction",
    dev_no_int = anova_result[1, "Deviance"], dev_with_int = anova_result[2, "Deviance"],
    resid_df_no_int = anova_result[1, "Resid. Df"], resid_df_with_int = anova_result[2, "Resid. Df"],
    dev_diff = anova_result[2, "Deviance"] - anova_result[1, "Deviance"],
    df_diff = anova_result[2, "Resid. Df"] - anova_result[1, "Resid. Df"],
    pvalue = anova_result[2, "Pr(>Chi)"],
    n_samples = nrow(hc_interact)
  )
  fwrite(int_test, file.path(OUTDIR, "age_site_interaction_test.csv"))
}

# Site-specific smooth predictions
cat("\n  Generating site-specific predictions...\n")
site_smooth <- rbindlist(lapply(major_sites, function(s) {
  sub <- hc[Site == s]
  if (nrow(sub) < 200) return(NULL)
  tryCatch({
    covars <- "s(Age, bs = 'tp')"
    for (cv in c("Platform", "SeqRegion", "Sex")) {
      if (cv %in% names(sub) && nlevels(sub[[cv]]) >= 2) covars <- paste(covars, "+", cv)
    }
    formula_str <- paste("RedComplex_detected ~", covars)
    gam_site <- gam(as.formula(formula_str), data = sub, family = binomial(), method = "REML")
    age_site <- data.table(Site = s, Age = seq(max(0.5, min(sub$Age)), min(95, max(sub$Age)), by = 0.5))
    for (col in c("Platform", "SeqRegion")) {
      if (col %in% names(sub) && nlevels(sub[[col]]) >= 2)
        age_site[[col]] <- factor(as.character(sub[[col]][1]), levels = levels(sub[[col]]))
    }
    pred <- predict(gam_site, newdata = age_site, type = "link", se.fit = TRUE)
    age_site[, `:=`(detect_prob = plogis(pred$fit),
      detect_low = plogis(pred$fit - 1.96 * pred$se.fit),
      detect_high = plogis(pred$fit + 1.96 * pred$se.fit), n_samples = nrow(sub))]
  }, error = function(e) {
    cat(sprintf("  GAM failed for site %s: %s\n", s, e$message)); NULL
  })
}))
fwrite(site_smooth, file.path(OUTDIR, "age_trajectory_by_site.csv"))

# ── 4. HC vs PD: FIXED COMPARISON (NEW) ─────────────────────────
cat("\n── 4. Disease vs HC: Age-adjusted RC ──\n")

pd <- dtc[DiseaseStatus == "Periodontal Disease" & !is.na(Age) & Age > 0 & Age < 100]

# Use consistent covariates (Site + Sex) for fair comparison
shared_formula <- RedComplex_detected ~ s(Age, bs = "tp") + Site + Sex

gam_hc_cmp <- gam(shared_formula, data = hc, family = binomial(), method = "REML")
gam_pd_cmp <- gam(shared_formula, data = pd, family = binomial(), method = "REML")

cat(sprintf("  HC GAM deviance: %.1f%%, PD GAM deviance: %.1f%%\n",
            summary(gam_hc_cmp)$dev.expl * 100, summary(gam_pd_cmp)$dev.expl * 100))

# Test age×disease interaction
cat("\n  Testing age × disease interaction...\n")
hc_pd_combined <- rbind(
  hc[, .(Age, Site, Sex, RedComplex_detected, DiseaseStatus = "HC")],
  pd[, .(Age, Site, Sex, RedComplex_detected, DiseaseStatus = "Periodontal Disease")]
)
hc_pd_combined[, DiseaseStatus := factor(DiseaseStatus, levels = c("HC", "Periodontal Disease"))]
hc_pd_combined[, Site := droplevels(Site)]
hc_pd_combined[, Sex := droplevels(Sex)]

m_disease_no_int <- gam(RedComplex_detected ~ s(Age, bs = "tp") + DiseaseStatus + Site + Sex,
                        data = hc_pd_combined, family = binomial(), method = "REML")

disease_int_p <- NA
m_disease_with_int <- tryCatch({
  gam(RedComplex_detected ~ s(Age, bs = "tp", by = DiseaseStatus) + DiseaseStatus + Site + Sex,
      data = hc_pd_combined, family = binomial(), method = "REML")
}, error = function(e) { NULL })

if (!is.null(m_disease_with_int)) {
  aov_disease <- anova(m_disease_no_int, m_disease_with_int, test = "Chisq")
  cat("  ANOVA: age trajectory × disease status:\n")
  print(aov_disease)
  disease_int_p <- aov_disease[2, "Pr(>Chi)"]
  cat(sprintf("\n  Age×Disease interaction: p = %s\n", format.pval(disease_int_p, digits = 4)))
  if (disease_int_p < 0.05) {
    cat("  → SIGNIFICANT: RC age trajectory differs between HC and PD\n")
  } else {
    cat("  → NOT significant: same age trajectory shape, just different intercept\n")
  }
}

# Generate comparable predictions with shared reference levels
ref_site <- if ("saliva" %in% levels(hc$Site)) "saliva" else levels(hc$Site)[1]
ref_sex <- if ("F" %in% levels(hc$Sex)) "F" else levels(hc$Sex)[1]

age_hc_pred <- data.table(Age = seq(0.5, 95, by = 0.5))
age_hc_pred[, `:=`(Site = factor(ref_site, levels = levels(hc$Site)),
                    Sex = factor(ref_sex, levels = levels(hc$Sex)))]
pred_hc_cmp <- predict(gam_hc_cmp, newdata = age_hc_pred, type = "link", se.fit = TRUE)
age_hc_pred[, `:=`(DiseaseStatus = "HC",
  detect_prob = plogis(pred_hc_cmp$fit),
  detect_low = plogis(pred_hc_cmp$fit - 1.96 * pred_hc_cmp$se.fit),
  detect_high = plogis(pred_hc_cmp$fit + 1.96 * pred_hc_cmp$se.fit))]

age_pd_pred <- data.table(Age = seq(10, 90, by = 0.5))
age_pd_pred[, `:=`(Site = factor(ref_site, levels = levels(pd$Site)),
                    Sex = factor(ref_sex, levels = levels(pd$Sex)))]
pred_pd_cmp <- predict(gam_pd_cmp, newdata = age_pd_pred, type = "link", se.fit = TRUE)
age_pd_pred[, `:=`(DiseaseStatus = "Periodontal Disease",
  detect_prob = plogis(pred_pd_cmp$fit),
  detect_low = plogis(pred_pd_cmp$fit - 1.96 * pred_pd_cmp$se.fit),
  detect_high = plogis(pred_pd_cmp$fit + 1.96 * pred_pd_cmp$se.fit))]

disease_age <- rbind(
  age_hc_pred[, .(Age, DiseaseStatus, detect_prob, detect_low, detect_high)],
  age_pd_pred[, .(Age, DiseaseStatus, detect_prob, detect_low, detect_high)]
)
fwrite(disease_age, file.path(OUTDIR, "age_trajectory_hc_vs_pd.csv"))

# ── 5. Visualization ────────────────────────────────────────────
cat("\n── 5. Visualizations ──\n")

p1 <- ggplot(age_grid, aes(x = Age, y = detect_prob)) +
  geom_ribbon(aes(ymin = detect_low, ymax = detect_high), fill = "#3498db", alpha = 0.2) +
  geom_line(color = "#2980b9", linewidth = 1.2) +
  labs(title = "Red Complex Detection Probability by Age (Healthy Controls)",
       subtitle = sprintf("n = %d, GAM smooth with 95%% CI", nrow(hc)),
       x = "Age (years)", y = "Detection Probability") +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) + theme_rc()
ggsave(file.path(OUTDIR, "age_trajectory_hc_detection.png"), p1, width = 8, height = 5, dpi = 150)

p2 <- ggplot(hc[Site %in% major_sites & Age <= 95], aes(x = Age, y = RedComplex_detected, color = Site, fill = Site)) +
  geom_smooth(method = "loess", span = 0.5, alpha = 0.15) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") +
  labs(title = "Red Complex Detection by Age and Oral Site (HC)",
       subtitle = sprintf("LOESS smooth, n = %d", nrow(hc)), x = "Age (years)",
       y = "Detection Probability", color = "Site", fill = "Site") + theme_rc()
ggsave(file.path(OUTDIR, "age_trajectory_by_site.png"), p2, width = 9, height = 5, dpi = 150)

p3 <- ggplot(disease_age, aes(x = Age, y = detect_prob, color = DiseaseStatus, fill = DiseaseStatus)) +
  geom_ribbon(aes(ymin = detect_low, ymax = detect_high), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = c("HC" = "#3498db", "Periodontal Disease" = "#e74c3c")) +
  scale_fill_manual(values = c("HC" = "#3498db", "Periodontal Disease" = "#e74c3c")) +
  labs(title = "Red Complex Detection: HC vs Periodontal Disease",
       subtitle = "Consistent covariates (Site + Sex), with 95% CI",
       x = "Age (years)", y = "Detection Probability") + theme_rc()
ggsave(file.path(OUTDIR, "age_trajectory_hc_vs_pd.png"), p3, width = 8, height = 5, dpi = 150)

p4 <- ggplot(hc_det[Age <= 95], aes(x = Age, y = RedComplexTotal + 1e-6)) +
  geom_point(alpha = 0.05, size = 0.5, color = "#3498db") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp"), color = "#c0392b", fill = "#c0392b", alpha = 0.2) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "Red Complex Abundance by Age (HC, detected only)",
       subtitle = sprintf("n = %d, log scale", nrow(hc_det)), x = "Age (years)",
       y = "Relative Abundance (log10 scale)") + theme_rc()
ggsave(file.path(OUTDIR, "age_trajectory_hc_abundance.png"), p4, width = 8, height = 5, dpi = 150)
cat("  All plots saved.\n")

# ── 6. Key findings ─────────────────────────────────────────────
cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\nKEY FINDINGS\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")

early <- hc_age[AgeGroup %in% c("0-6 months", "6-36 months")]
cat(sprintf("\n  1. Early colonization: 0-6mo=%.1f%%, 6-36mo=%.1f%%\n",
            early[AgeGroup == "0-6 months"]$rc_detection * 100,
            early[AgeGroup == "6-36 months"]$rc_detection * 100))

peak_age <- age_grid[which.max(detect_prob)]
cat(sprintf("\n  2. Peak detection: %.0f%% at age %.0f\n", peak_age$detect_prob * 100, peak_age$Age))

cat("\n  3. Age × Site interaction:\n")
if (!is.na(aov_p)) {
  cat(sprintf("     GAM interaction test: p = %s\n", format.pval(aov_p, digits = 4)))
} else {
  cat("     Interaction model did not converge\n")
}

cat("\n  4. Age × Disease interaction:\n")
if (!is.na(disease_int_p)) {
  cat(sprintf("     HC vs PD age trajectory: p = %s\n", format.pval(disease_int_p, digits = 4)))
} else {
  cat("     Interaction model did not converge\n")
}

cat("\n  5. Site-specific early colonization (LOESS at age 5):\n")
for (s in major_sites[1:min(5, length(major_sites))]) {
  sub_site <- hc[Site == s & Age <= 95 & Age >= 0.5]
  if (nrow(sub_site) > 50) {
    lo <- loess(RedComplex_detected ~ Age, data = sub_site, span = 0.5)
    pred5 <- predict(lo, newdata = data.frame(Age = 5))
    cat(sprintf("     %s: %.1f%%\n", s, pmax(0, pred5) * 100))
  }
}

cat("\nDone.\n")
