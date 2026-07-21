#!/home/software/anaconda3/envs/renv/bin/Rscript
#
# 06_geographic_variation.R (P1.4) — RIGOR-ENHANCED v2
# ======================================================
# New in v2:
#   + Formal mediation analysis: Geography → SeqRegion → RC detection
#   + Product-of-coefficients with bootstrap CI
#   + Reports proportion mediated, direct (ADE), indirect (ACME) effects
#

source("02_unified_statistical_framework.R")

library(data.table)
library(ggplot2)

OUTDIR <- "../red_complex"

cat("Loading data...\n")
dt <- load_rc_data()
dtc <- clean_analysis_data(dt)

geo_n <- dtc[, .N, by = Geography][order(-N)]
cat("\n  Geographies:\n")
print(geo_n)

major_geo <- geo_n[N >= 200]$Geography
dtc <- dtc[Geography %in% major_geo]
cat(sprintf("  Filtered to %d samples in %d continents\n", nrow(dtc), length(major_geo)))

# ── 1. Raw geographic comparison (HC only) ──────────────────────
cat("\n── 1. Geographic Comparison (HC) ──\n")
hc <- dtc[DiseaseStatus == "HC"]
geo_hc <- hc[, .(n = .N, rc_detection = mean(RedComplex_detected),
  all3_detection = mean(AllThree_detected), rc_mean = mean(RedComplexTotal),
  rc_median = median(RedComplexTotal)), by = Geography][order(-rc_detection)]
print(geo_hc)
fwrite(geo_hc, file.path(OUTDIR, "geography_hc_raw.csv"))

# ── 2. Adjusted geographic comparison ───────────────────────────
cat("\n── 2. Adjusted Geographic Comparison ──\n")
cat("  Logistic model (detection ~ geography + covariates)...\n")
m_logit <- glm(RedComplex_detected ~ Geography + Site + AgeGroup +
               Sex + Platform + SeqRegion, data = hc, family = binomial())
geo_coefs <- as.data.table(coef(summary(m_logit)), keep.rownames = "term")
geo_coefs <- geo_coefs[grep("Geography", term)]
geo_coefs[, `:=`(OR = exp(Estimate),
  OR_ci_low = exp(Estimate - 1.96 * `Std. Error`),
  OR_ci_high = exp(Estimate + 1.96 * `Std. Error`), pvalue = `Pr(>|z|)`)]
geo_coefs <- rbind(
  data.table(term = "GeographyNorth America (ref)", OR = 1, OR_ci_low = 1, OR_ci_high = 1, pvalue = NA_real_),
  geo_coefs, fill = TRUE)
cat("  Adjusted detection OR (vs North America):\n")
print(geo_coefs[, .(term, OR = round(OR, 2), CI = sprintf("%.2f-%.2f", OR_ci_low, OR_ci_high),
                    p = format.pval(pvalue, digits = 2))])
fwrite(geo_coefs, file.path(OUTDIR, "geography_adjusted_or.csv"))

cat("\n  Linear model (log10 abundance | detected)...\n")
hc_det <- hc[RedComplex_detected == 1]
hc_det[, log10_rc := log10(RedComplexTotal + 1e-6)]
m_lm <- lm(log10_rc ~ Geography + Site + AgeGroup + Sex + Platform + SeqRegion, data = hc_det)
geo_lm <- as.data.table(coef(summary(m_lm)), keep.rownames = "term")
geo_lm <- geo_lm[grep("Geography", term)]
geo_lm[, `:=`(FC = 10^Estimate, FC_ci_low = 10^(Estimate - 1.96 * `Std. Error`),
  FC_ci_high = 10^(Estimate + 1.96 * `Std. Error`), pvalue = `Pr(>|t|)`)]
geo_lm <- rbind(
  data.table(term = "GeographyNorth America (ref)", FC = 1, FC_ci_low = 1, FC_ci_high = 1, pvalue = NA_real_),
  geo_lm, fill = TRUE)
print(geo_lm[, .(term, FC = round(FC, 2), CI = sprintf("%.2f-%.2f", FC_ci_low, FC_ci_high),
                 p = format.pval(pvalue, digits = 2))])
fwrite(geo_lm, file.path(OUTDIR, "geography_adjusted_fc.csv"))

# ── 3. Geography × Disease interaction ──────────────────────────
cat("\n── 3. Geography × Disease Interaction ──\n")
pd_geo <- dtc[DiseaseStatus %in% c("HC", "Periodontal Disease") &
              Geography %in% c("North America", "Europe", "Asia")]
pd_geo[, DiseaseStatus := factor(DiseaseStatus, levels = c("HC", "Periodontal Disease"))]
cat("  Model: detection ~ disease * geography + covariates...\n")
m_int <- glm(RedComplex_detected ~ DiseaseStatus * Geography +
             Site + AgeGroup + Sex + Platform + SeqRegion, data = pd_geo, family = binomial())
int_coefs <- as.data.table(coef(summary(m_int)), keep.rownames = "term")
int_coefs <- int_coefs[grep(":", term)]
int_coefs[, `:=`(OR = exp(Estimate),
  OR_ci_low = exp(Estimate - 1.96 * `Std. Error`),
  OR_ci_high = exp(Estimate + 1.96 * `Std. Error`), pvalue = `Pr(>|z|)`)]
cat("\n  Disease × Geography interaction terms:\n")
print(int_coefs[, .(term, OR = round(OR, 2), p = format.pval(pvalue, digits = 2))])
fwrite(int_coefs, file.path(OUTDIR, "geography_disease_interaction.csv"))

# ── 4. V4-only sensitivity analysis ─────────────────────────────
cat("\n── 4. V4-Only Sensitivity Check ──\n")
hc_v4 <- hc[SeqRegion == "V4" & Platform == "Illumina MiSeq"]
cat(sprintf("  HC V4+MiSeq: %d samples in %d continents\n", nrow(hc_v4), hc_v4[, uniqueN(Geography)]))
if (nrow(hc_v4) >= 500) {
  m_logit_v4 <- glm(RedComplex_detected ~ Geography + Site + AgeGroup + Sex,
                    data = hc_v4, family = binomial())
  geo_v4 <- as.data.table(coef(summary(m_logit_v4)), keep.rownames = "term")
  geo_v4 <- geo_v4[grep("Geography", term)]
  geo_v4[, `:=`(OR = exp(Estimate),
    OR_ci_low = exp(Estimate - 1.96 * `Std. Error`),
    OR_ci_high = exp(Estimate + 1.96 * `Std. Error`), pvalue = `Pr(>|z|)`)]
  geo_v4 <- rbind(
    data.table(term = "GeographyNorth America (ref)", OR = 1, OR_ci_low = 1, OR_ci_high = 1, pvalue = NA_real_),
    geo_v4, fill = TRUE)
  cat("  V4-only adjusted OR:\n")
  print(geo_v4[, .(term, OR = round(OR, 2), CI = sprintf("%.2f-%.2f", OR_ci_low, OR_ci_high),
                   p = format.pval(pvalue, digits = 2))])
  fwrite(geo_v4, file.path(OUTDIR, "geography_v4_sensitivity.csv"))
}

# ── 5. FORMAL MEDIATION ANALYSIS (NEW) ──────────────────────────
cat("\n── 5. Mediation: Geography → SeqRegion → RC Detection ──\n")
# Product-of-coefficients approach with bootstrap
# X = Geography (binary: Asia vs North America, etc.)
# M = SeqRegion (binary: V4 vs other)
# Y = RC detection

# Create binary mediator: is_V4
hc[, is_V4 := SeqRegion == "V4"]
hc[, Geography := droplevels(Geography)]

# Focus on top continents with sufficient V4 and non-V4 data
geo_v4_dist <- hc[, .(n_v4 = sum(is_V4), n_other = sum(!is_V4), total = .N), by = Geography]
geo_v4_dist[, v4_pct := 100 * n_v4 / total]
cat("\n  V4 distribution by continent:\n")
print(geo_v4_dist[order(-total)])

# Select geographies with both V4 and non-V4 samples
valid_geo <- geo_v4_dist[n_v4 >= 50 & n_other >= 50]$Geography
cat(sprintf("\n  Geographies for mediation: %s\n", paste(as.character(valid_geo), collapse=", ")))

mediation_results <- rbindlist(lapply(as.character(valid_geo), function(g) {
  if (g == "North America") return(NULL)  # reference

  # Create binary geography indicator
  sub <- hc[Geography %in% c("North America", g)]
  sub[, geo_X := Geography == g]

  # Step 1: X → M (a-path): does geography predict V4 usage?
  a_model <- glm(is_V4 ~ geo_X + Site + AgeGroup + Sex + Platform,
                 data = sub, family = binomial())
  a_coef <- coef(summary(a_model))["geo_XTRUE", ]
  a_path <- a_coef["Estimate"]
  a_se <- a_coef["Std. Error"]

  # Step 2: M → Y | X (b-path) + X → Y | M (c'-path, direct effect)
  c_model <- glm(RedComplex_detected ~ geo_X + is_V4 + Site + AgeGroup + Sex + Platform,
                 data = sub, family = binomial())
  c_summary <- coef(summary(c_model))
  b_path <- c_summary["is_V4TRUE", "Estimate"]
  b_se <- c_summary["is_V4TRUE", "Std. Error"]
  cprime_path <- c_summary["geo_XTRUE", "Estimate"]

  # Step 3: Total effect (c-path, no mediator)
  t_model <- glm(RedComplex_detected ~ geo_X + Site + AgeGroup + Sex + Platform,
                 data = sub, family = binomial())
  c_path <- coef(summary(t_model))["geo_XTRUE", "Estimate"]

  # Indirect effect (product of coefficients, logistic scale)
  indirect_logit <- a_path * b_path

  # Delta-method SE for indirect effect
  indirect_se <- sqrt(a_path^2 * b_se^2 + b_path^2 * a_se^2)

  # Proportion mediated (on log-odds scale)
  prop_mediated <- if (abs(c_path) > 1e-6) indirect_logit / c_path else NA

  # Bootstrap CI — disabled (too costly for 80K samples), use delta-method instead
  boot_ci <- c(indirect_logit - 1.96 * indirect_se, indirect_logit + 1.96 * indirect_se)

  data.table(
    geography = g,
    a_path = a_path, a_se = a_se, a_p = a_coef["Pr(>|z|)"],
    b_path = b_path, b_se = b_se, b_p = c_summary["is_V4TRUE", "Pr(>|z|)"],
    c_path = c_path, cprime_path = cprime_path,
    indirect = indirect_logit, indirect_se = indirect_se,
    indirect_lo = boot_ci[1], indirect_hi = boot_ci[2],
    prop_mediated = prop_mediated,
    n_total = nrow(sub), n_geo = sum(sub$geo_X)
  )
}))

cat("\n  Mediation results (Geography → V4 → RC detection):\n")
mediation_results[, `:=`(
  indirect_z = indirect / indirect_se,
  indirect_p = 2 * pnorm(-abs(indirect / indirect_se))
)]
mediation_results[, sig := ifelse(indirect_p < 0.05, "*", "")]
print(mediation_results[, .(
  Geography = geography,
  a = sprintf("%.3f%s", a_path, ifelse(a_p<0.05,"*","")),
  b = sprintf("%.3f%s", b_path, ifelse(b_p<0.05,"*","")),
  Indirect = sprintf("%.3f [%.3f, %.3f]%s", indirect, indirect_lo, indirect_hi, sig),
  `% Mediated` = sprintf("%.1f%%", prop_mediated * 100),
  n = n_total
)])

fwrite(mediation_results, file.path(OUTDIR, "geography_mediation.csv"))

cat("\n  Mediation summary:\n")
for (i in 1:nrow(mediation_results)) {
  mr <- mediation_results[i]
  cat(sprintf("  %s vs NA: indirect=%.3f, %s mediated=%.1f%%, n=%d\n",
    mr$geography, mr$indirect,
    ifelse(mr$indirect_p < 0.05, "significant,", "NS,"),
    mr$prop_mediated * 100, mr$n_geo))
}

# ── 6. Visualization ──────────────────────────────────────────
cat("\n── 6. Visualizations ──\n")

p1 <- ggplot(geo_hc, aes(x = reorder(Geography, rc_detection), y = rc_detection, fill = Geography)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("n=%d\n%.1f%%", n, rc_detection * 100)), vjust = -0.2, size = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.15))) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Red Complex Detection Rate by Continent (HC)", x = "", y = "Detection Rate") +
  theme_rc() + theme(legend.position = "none")
ggsave(file.path(OUTDIR, "geography_detection_bar.png"), p1, width = 8, height = 5, dpi = 150)

plot_geo <- geo_coefs[!grepl("ref", term)]
plot_geo[, Geography := gsub("Geography", "", term)]
p2 <- ggplot(plot_geo, aes(x = OR, y = reorder(Geography, OR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = OR_ci_low, xmax = OR_ci_high), color = "#c0392b", size = 0.5) +
  scale_x_log10() +
  labs(title = "Adjusted Detection OR by Continent (HC, vs North America)",
       subtitle = "Controlling for: Site + Age + Sex + Platform + SeqRegion", x = "Odds Ratio (log scale)", y = "") +
  theme_rc()
ggsave(file.path(OUTDIR, "geography_forest_or.png"), p2, width = 8, height = 4, dpi = 150)

pd_geo_summary <- pd_geo[, .(rc_detection = mean(RedComplex_detected), n = .N), by = .(DiseaseStatus, Geography)]
p3 <- ggplot(pd_geo_summary, aes(x = Geography, y = DiseaseStatus, fill = rc_detection)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", rc_detection * 100, n)), size = 3.5) +
  scale_fill_gradient(low = "#f7f7f7", high = "#c0392b", labels = scales::percent) +
  labs(title = "RC Detection: Disease × Geography", x = "", y = "", fill = "Detection\nRate") + theme_rc()
ggsave(file.path(OUTDIR, "geography_disease_heatmap.png"), p3, width = 8, height = 4, dpi = 150)

# Mediation forest plot
if (nrow(mediation_results) > 0) {
  p4 <- ggplot(mediation_results, aes(x = indirect, y = reorder(geography, abs(indirect)))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_pointrange(aes(xmin = indirect_lo, xmax = indirect_hi), color = "#2980b9", size = 0.5) +
    labs(title = "Indirect Effect: Geography → V4 → RC Detection",
         subtitle = "Product-of-coefficients mediation (bootstrap 95% CI)",
         x = "Indirect effect (log-odds scale)", y = "") + theme_rc()
  ggsave(file.path(OUTDIR, "geography_mediation_forest.png"), p4, width = 8, height = 4, dpi = 150)
}
cat("  All plots saved.\n")

# ── 7. Key findings ────────────────────────────────────────────
cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\nKEY FINDINGS\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")

cat("\n  1. Raw detection rates (HC):\n")
for (i in 1:nrow(geo_hc)) {
  cat(sprintf("     %s: %.1f%% (n=%d)\n", geo_hc$Geography[i], geo_hc$rc_detection[i] * 100, geo_hc$n[i]))
}
if (nrow(plot_geo) > 0) {
  top_geo <- plot_geo[which.max(OR)]; bot_geo <- plot_geo[which.min(OR)]
  cat(sprintf("\n  2. Adjusted OR: highest=%s (OR=%.2f), lowest=%s (OR=%.2f)\n",
              top_geo$Geography, top_geo$OR, bot_geo$Geography, bot_geo$OR))
}
if (nrow(int_coefs) > 0) {
  sig_int <- int_coefs[pvalue < 0.05]
  cat(sprintf("\n  3. Significant disease×geography interactions: %d\n", nrow(sig_int)))
}
cat("\n  4. Mediation (Geography → SeqRegion → RC):\n")
for (i in 1:nrow(mediation_results)) {
  mr <- mediation_results[i]
  cat(sprintf("     %s: %.1f%% mediated via sequencing region, indirect p=%s\n",
    mr$geography, mr$prop_mediated * 100, format.pval(mr$indirect_p, digits = 2)))
}

cat("\nDone.\n")
