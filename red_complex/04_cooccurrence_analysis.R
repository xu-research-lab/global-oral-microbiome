#!/home/software/anaconda3/envs/renv/bin/Rscript
#
# 04_cooccurrence_analysis.R (P1.2) — RIGOR-ENHANCED v2
# ======================================================
# New in v2:
#   + Cochran-Mantel-Haenszel (CMH) test stratified by Site
#   + Site-adjusted odds ratios vs raw OR comparison
#   + Breslow-Day test for homogeneity of OR across strata
#   + Partial phi controlling for site confounding
#

source("02_unified_statistical_framework.R")

library(data.table)
library(ggplot2)
library(scales)

OUTDIR <- "../red_complex"

cat("Loading data...\n")
dt <- load_rc_data()
dtc <- clean_analysis_data(dt)
cat(sprintf("  Analysis set: %d samples\n", nrow(dtc)))

# ── 1. Global pairwise associations ─────────────────────────────
cat("\n── 1. Pairwise Co-occurrence (Fisher + Phi) ──\n")

pairwise_cooccurrence <- function(dt, sp1, sp2) {
  tab <- table(dt[[paste0(sp1, "_detected")]], dt[[paste0(sp2, "_detected")]])
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
  ft <- fisher.test(tab)
  n <- sum(tab)
  chi2 <- suppressWarnings(chisq.test(tab)$statistic)
  phi <- sqrt(chi2 / n)
  sp1_given_sp2 <- tab[2, 2] / sum(tab[, 2])
  sp2_given_sp1 <- tab[2, 2] / sum(tab[2, ])
  data.table(species1 = sp1, species2 = sp2, n_total = n,
    n_both = tab[2, 2], n_neither = tab[1, 1], n_sp1_only = tab[2, 1], n_sp2_only = tab[1, 2],
    or_val = as.numeric(ft$estimate), or_ci_low = ft$conf.int[1], or_ci_high = ft$conf.int[2],
    phi = phi, p_fisher = ft$p.value, p_sp1_given_sp2 = sp1_given_sp2, p_sp2_given_sp1 = sp2_given_sp1)
}

pairs <- list(c("P_gingivalis", "T_denticola"), c("P_gingivalis", "T_forsythia"), c("T_denticola", "T_forsythia"))

global_pairs <- rbindlist(lapply(pairs, function(p) pairwise_cooccurrence(dtc, p[1], p[2])))
cat("\n  Global pairwise associations:\n")
print(global_pairs[, .(species1, species2, n_total, n_both, OR = round(or_val, 1),
  Phi = round(phi, 3), p = format.pval(p_fisher, digits = 2))])
fwrite(global_pairs, file.path(OUTDIR, "cooccurrence_global_pairs.csv"))

# ── 2. Co-occurrence by disease status ──────────────────────────
cat("\n── 2. Co-occurrence by Disease Status ──\n")
major_diseases <- dtc[, .N, by = DiseaseStatus][N >= 200]$DiseaseStatus
disease_pairs <- rbindlist(lapply(as.character(major_diseases), function(d) {
  sub <- dtc[DiseaseStatus == d]
  if (nrow(sub) < 100) return(NULL)
  rbindlist(lapply(pairs, function(p) {
    res <- pairwise_cooccurrence(sub, p[1], p[2])
    if (is.null(res)) return(NULL)
    res[, DiseaseStatus := d]; res
  }))
}))
disease_pairs[, pair := paste0(gsub("_", ".", substr(species1, 1, 2)), "-",
                                gsub("_", ".", substr(species2, 1, 2)))]
pg_tf <- disease_pairs[species1 == "P_gingivalis" & species2 == "T_forsythia"][order(-phi)]
cat("\n  Pg-Tf phi by disease:\n")
print(pg_tf[, .(DiseaseStatus, n_total, Phi = round(phi, 3), OR = round(or_val, 1))])
fwrite(disease_pairs, file.path(OUTDIR, "cooccurrence_by_disease.csv"))

# ── 3. Co-occurrence by oral site ───────────────────────────────
cat("\n── 3. Co-occurrence by Oral Site ──\n")
major_sites <- dtc[, .N, by = Site][N >= 200]$Site
site_pairs <- rbindlist(lapply(as.character(major_sites), function(s) {
  sub <- dtc[Site == s]
  if (nrow(sub) < 100) return(NULL)
  rbindlist(lapply(pairs, function(p) {
    res <- pairwise_cooccurrence(sub, p[1], p[2])
    if (is.null(res)) return(NULL)
    res[, Site := s]; res
  }))
}))
fwrite(site_pairs, file.path(OUTDIR, "cooccurrence_by_site.csv"))

# ── 3b. MANTEL-HAENSZEL ADJUSTED CO-OCCURRENCE (NEW) ───────────
cat("\n── 3b. Mantel-Haenszel Stratified Analysis ──\n")
cat("  Stratifying by Site to control for site-specific co-occurrence...\n")

# For each species pair, compute CMH OR and test
mh_results <- rbindlist(lapply(pairs, function(p) {
  sp1_col <- paste0(p[1], "_detected"); sp2_col <- paste0(p[2], "_detected")

  # Build 2x2 tables per site stratum
  site_strata <- dtc[, .N, by = Site][N >= 100]$Site

  # Collect per-stratum data
  stratum_tabs <- lapply(as.character(site_strata), function(s) {
    sub <- dtc[Site == s]
    tab <- table(sub[[sp1_col]], sub[[sp2_col]])
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
    tab
  })
  stratum_tabs <- stratum_tabs[!sapply(stratum_tabs, is.null)]

  if (length(stratum_tabs) < 2) return(NULL)

  # Mantel-Haenszel common OR
  # Manual CMH computation
  n11 <- as.numeric(sapply(stratum_tabs, function(x) x[2, 2]))
  n10 <- as.numeric(sapply(stratum_tabs, function(x) x[2, 1]))
  n01 <- as.numeric(sapply(stratum_tabs, function(x) x[1, 2]))
  n00 <- as.numeric(sapply(stratum_tabs, function(x) x[1, 1]))
  Nk <- n11 + n10 + n01 + n00

  # CMH OR = sum(n11*n00/Nk) / sum(n10*n01/Nk)
  cmh_num <- sum(n11 * n00 / Nk, na.rm = TRUE)
  cmh_denom <- sum(n10 * n01 / Nk, na.rm = TRUE)
  cmh_or <- if (cmh_denom > 0) cmh_num / cmh_denom else NA

  # CMH test statistic (use as.numeric to avoid integer overflow)
  exp_n11 <- (n11 + n10) * (n11 + n01) / Nk
  var_n11 <- (n11 + n10) * (n01 + n00) * (n10 + n00) * (n11 + n01) / (Nk^2 * (Nk - 1))
  cmh_stat <- (sum(n11 - exp_n11))^2 / sum(var_n11)
  cmh_p <- pchisq(cmh_stat, 1, lower.tail = FALSE)

  # Breslow-Day test for homogeneity of OR across strata
  # OR_k for each stratum
  or_k <- (n11 * n00) / (n10 * n01)
  # Woolf method for homogeneity (simpler than BD iteration)
  log_or_k <- log(or_k)
  log_or_k <- log_or_k[is.finite(log_or_k)]
  var_log_or_k <- 1/n11 + 1/n10 + 1/n01 + 1/n00
  var_log_or_k <- var_log_or_k[is.finite(log_or_k)]

  if (length(log_or_k) >= 2) {
    # Weighted mean of log OR
    w <- 1 / var_log_or_k[1:length(log_or_k)]
    log_or_pooled <- sum(w * log_or_k) / sum(w)
    # Woolf test for homogeneity (approximates Breslow-Day)
    woolf_stat <- sum(w * (log_or_k - log_or_pooled)^2)
    woolf_df <- length(log_or_k) - 1
    woolf_p <- pchisq(woolf_stat, woolf_df, lower.tail = FALSE)
  } else {
    woolf_stat <- NA; woolf_df <- NA; woolf_p <- NA
  }

  # Unadjusted (raw) OR for comparison
  raw_tab <- table(dtc[[sp1_col]], dtc[[sp2_col]])
  raw_or <- (raw_tab[2,2] * raw_tab[1,1]) / (raw_tab[2,1] * raw_tab[1,2])

  data.table(
    species1 = p[1], species2 = p[2],
    raw_OR = raw_or,
    cmh_OR = cmh_or,
    cmh_chi2 = cmh_stat,
    cmh_p = cmh_p,
    or_change_pct = 100 * (cmh_or - raw_or) / raw_or,
    woolf_chi2 = woolf_stat,
    woolf_df = woolf_df,
    woolf_p = woolf_p,
    n_strata = length(stratum_tabs)
  )
}))

cat("\n  Mantel-Haenszel adjusted vs raw OR:\n")
mh_results[, homogeneity := ifelse(!is.na(woolf_p) & woolf_p > 0.05, "Homogeneous", "Heterogeneous")]
print(mh_results[, .(
  Pair = paste0(substr(species1,1,2), "-", substr(species2,1,2)),
  Raw_OR = round(raw_OR, 1),
  CMH_OR = round(cmh_OR, 1),
  Change = sprintf("%+.1f%%", or_change_pct),
  CMH_p = format.pval(cmh_p, digits = 2),
  Homogeneity = homogeneity,
  Strata = n_strata
)])

fwrite(mh_results, file.path(OUTDIR, "cooccurrence_mh_adjusted.csv"))

# ── 3c. SITE-STRATIFIED CMH BY DISEASE (NEW) ────────────────────
cat("\n── 3c. CMH Stratified by Site, Within Disease Groups ──\n")

disease_cmh <- rbindlist(lapply(as.character(major_diseases), function(d) {
  sub_d <- dtc[DiseaseStatus == d]
  if (nrow(sub_d) < 200) return(NULL)

  rbindlist(lapply(pairs, function(p) {
    sp1_col <- paste0(p[1], "_detected"); sp2_col <- paste0(p[2], "_detected")
    site_strata <- sub_d[, .N, by = Site][N >= 30]$Site
    if (length(site_strata) < 2) return(NULL)

    tabs <- lapply(as.character(site_strata), function(s) {
      sub <- sub_d[Site == s]
      tab <- table(sub[[sp1_col]], sub[[sp2_col]])
      if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
      tab
    })
    tabs <- tabs[!sapply(tabs, is.null)]
    if (length(tabs) < 2) return(NULL)

    n11 <- as.numeric(sapply(tabs, function(x) x[2,2])); n10 <- as.numeric(sapply(tabs, function(x) x[2,1]))
    n01 <- as.numeric(sapply(tabs, function(x) x[1,2])); n00 <- as.numeric(sapply(tabs, function(x) x[1,1]))
    Nk <- n11 + n10 + n01 + n00

    cmh_num <- sum(n11 * n00 / Nk, na.rm = TRUE)
    cmh_denom <- sum(n10 * n01 / Nk, na.rm = TRUE)
    cmh_or <- if (cmh_denom > 0) cmh_num / cmh_denom else NA

    exp_n11 <- (n11 + n10) * (n11 + n01) / Nk
    var_n11 <- (n11 + n10) * (n01 + n00) * (n10 + n00) * (n11 + n01) / (Nk^2 * (Nk - 1))
    cmh_stat <- (sum(n11 - exp_n11))^2 / sum(var_n11)
    cmh_p <- pchisq(cmh_stat, 1, lower.tail = FALSE)

    # Raw OR within this disease
    raw_tab <- table(sub_d[[sp1_col]], sub_d[[sp2_col]])
    raw_or <- (raw_tab[2,2] * raw_tab[1,1]) / (raw_tab[2,1] * raw_tab[1,2])

    data.table(DiseaseStatus = d, species1 = p[1], species2 = p[2],
      n_total = nrow(sub_d), n_strata = length(tabs),
      raw_OR = raw_or, cmh_OR = cmh_or,
      or_change_pct = 100 * (cmh_or - raw_or) / raw_or,
      cmh_p = cmh_p)
  }))
}))

cat("\n  Selected disease CMH comparisons (HC vs PD):\n")
disease_cmh_show <- disease_cmh[DiseaseStatus %in% c("HC", "Periodontal Disease")]
disease_cmh_show[, pair := paste0(substr(species1,1,2), "-", substr(species2,1,2))]
print(disease_cmh_show[, .(DiseaseStatus, pair,
  Raw_OR = round(raw_OR, 1), CMH_OR = round(cmh_OR, 1),
  Change = sprintf("%+.1f%%", or_change_pct),
  CMH_p = format.pval(cmh_p, digits = 2))])

fwrite(disease_cmh, file.path(OUTDIR, "cooccurrence_cmb_by_disease.csv"))

# ── 4. Combination distribution ─────────────────────────────────
cat("\n── 4. 8-Combination Distribution ──\n")
combo_overall <- dtc[, .N, by = RC_combination][order(-N)]
combo_overall[, pct := 100 * N / sum(N)]
cat("\n  Overall combination distribution:\n"); print(combo_overall)
fwrite(combo_overall, file.path(OUTDIR, "cooccurrence_combinations_overall.csv"))

combo_by_disease <- dtc[, .N, by = .(DiseaseStatus, RC_combination)]
combo_by_disease[, pct := 100 * N / sum(N), by = DiseaseStatus]
fwrite(combo_by_disease, file.path(OUTDIR, "cooccurrence_combinations_by_disease.csv"))

cat("\n  Periodontal Disease vs HC combination comparison:\n")
pd_vs_hc <- combo_by_disease[DiseaseStatus %in% c("HC", "Periodontal Disease")]
pd_wide <- dcast(pd_vs_hc, RC_combination ~ DiseaseStatus, value.var = "pct")
pd_wide[, diff := `Periodontal Disease` - HC]
print(pd_wide[order(-diff)])

# ── 5. Conditional probability analysis ─────────────────────────
cat("\n── 5. Conditional Probabilities ──\n")
cond_probs <- function(sub, label) {
  pg_det <- sub[P_gingivalis_detected == 1]; td_det <- sub[T_denticola_detected == 1]; tf_det <- sub[T_forsythia_detected == 1]
  data.table(group = label, n = nrow(sub),
    pg_rate = sub[, mean(P_gingivalis_detected)], td_rate = sub[, mean(T_denticola_detected)], tf_rate = sub[, mean(T_forsythia_detected)],
    td_given_pg = pg_det[, mean(T_denticola_detected)], tf_given_pg = pg_det[, mean(T_forsythia_detected)],
    pg_given_td = td_det[, mean(P_gingivalis_detected)], tf_given_td = td_det[, mean(T_forsythia_detected)],
    pg_given_tf = tf_det[, mean(P_gingivalis_detected)], td_given_tf = tf_det[, mean(T_denticola_detected)],
    all3_given_pg_td = sub[P_gingivalis_detected == 1 & T_denticola_detected == 1, mean(T_forsythia_detected)],
    all3_given_pg_tf = sub[P_gingivalis_detected == 1 & T_forsythia_detected == 1, mean(T_denticola_detected)],
    all3_given_td_tf = sub[T_denticola_detected == 1 & T_forsythia_detected == 1, mean(P_gingivalis_detected)])
}
cond_combined <- rbind(cond_probs(dtc, "All"), cond_probs(dtc[DiseaseStatus == "HC"], "HC"),
                        cond_probs(dtc[DiseaseStatus == "Periodontal Disease"], "PD"))
cat("\n  Conditional detection probabilities:\n")
print(cond_combined[, .(group, n, Pg = round(pg_rate, 3), Td = round(td_rate, 3), Tf = round(tf_rate, 3),
  `Td|Pg` = round(td_given_pg, 3), `Tf|Pg` = round(tf_given_pg, 3), `All3|Pg+Td` = round(all3_given_pg_td, 3))])
fwrite(cond_combined, file.path(OUTDIR, "cooccurrence_conditional_probs.csv"))

# ── 6. Visualization ────────────────────────────────────────────
cat("\n── 6. Visualizations ──\n")

# ── Phi heatmap by disease (top diseases with n>=200, with numeric labels) ──
# Filter to diseases with decent sample size
major_dis <- dtc[, .N, by = DiseaseStatus][N >= 200]$DiseaseStatus
plot_phi <- disease_pairs[DiseaseStatus %in% major_dis]
# Short pair labels
plot_phi[, pair_label := paste0(substring(species1, 1, 2), "-", substring(species2, 1, 2))]
plot_phi_wide <- dcast(plot_phi, DiseaseStatus ~ pair_label, value.var = "phi")
phi_mat <- as.matrix(plot_phi_wide[, -"DiseaseStatus"])
rownames(phi_mat) <- plot_phi_wide$DiseaseStatus
phi_mat <- phi_mat[order(rowMeans(phi_mat, na.rm = TRUE)), , drop = FALSE]
phi_mat <- phi_mat[, order(colMeans(phi_mat, na.rm = TRUE)), drop = FALSE]

plot_phi_long <- melt(as.data.table(phi_mat, keep.rownames = "Disease"),
                       id.vars = "Disease", variable.name = "Pair", value.name = "Phi")
plot_phi_long[, Disease := factor(Disease, levels = rownames(phi_mat))]
plot_phi_long[, Pair := factor(Pair, levels = colnames(phi_mat))]
# Add N per disease
dis_n <- dtc[DiseaseStatus %in% rownames(phi_mat), .N, by = DiseaseStatus]
plot_phi_long <- merge(plot_phi_long, dis_n, by.x = "Disease", by.y = "DiseaseStatus", all.x = TRUE)
plot_phi_long[, DiseaseLabel := sprintf("%s (n=%d)", Disease, N)]

phi_max <- max(plot_phi_long$Phi, na.rm = TRUE)

p_phi <- ggplot(plot_phi_long, aes(x = Pair, y = DiseaseLabel, fill = Phi)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Phi)), size = 4.5, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#f7f7f7", "#fdd49e", "#fdbb84", "#fc8d59", "#e34a33", "#b2182b"),
    limits = c(0, phi_max + 0.05),
    guide = guide_colorbar(barwidth = 0.8, barheight = 6)) +
  labs(title = "Red Complex Pairwise Co-occurrence (Phi coefficient)",
       subtitle = sprintf("%d diseases (n≥200); numeric values = φ; higher = stronger co-occurrence", nrow(phi_mat)),
       x = "", y = "", fill = expression(phi)) +
  theme_rc() +
  theme(axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

ggsave(file.path(OUTDIR, "cooccurrence_phi_heatmap.png"), p_phi,
       width = 8, height = max(4, nrow(phi_mat) * 0.55), dpi = 150)

# HC vs PD phi comparison
phi_compare <- disease_pairs[DiseaseStatus %in% c("HC", "Periodontal Disease")]
p_phi_comp <- ggplot(phi_compare, aes(x = pair, y = phi, fill = DiseaseStatus)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("HC" = "#3498db", "Periodontal Disease" = "#e74c3c")) +
  labs(title = "Co-occurrence Strength: HC vs Periodontal Disease", x = "Species Pair", y = "Phi Coefficient") +
  ylim(0, 1) + theme_rc()
ggsave(file.path(OUTDIR, "cooccurrence_hc_vs_pd.png"), p_phi_comp, width = 8, height = 5, dpi = 150)

# CMH vs Raw OR comparison
if (nrow(mh_results) > 0) {
  mh_plot <- melt(mh_results, id.vars = c("species1", "species2"),
                   measure.vars = c("raw_OR", "cmh_OR"),
                   variable.name = "Method", value.name = "OR")
  mh_plot[, Method := ifelse(Method == "raw_OR", "Raw (unadjusted)", "CMH (site-adjusted)")]
  mh_plot[, pair := paste0(substr(species1,1,2), "-", substr(species2,1,2))]
  p_mh <- ggplot(mh_plot, aes(x = pair, y = OR, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    scale_fill_manual(values = c("Raw (unadjusted)" = "#95a5a6", "CMH (site-adjusted)" = "#2980b9")) +
    labs(title = "Co-occurrence OR: Raw vs Site-Adjusted (CMH)", x = "Species Pair", y = "Odds Ratio") +
    theme_rc()
  ggsave(file.path(OUTDIR, "cooccurrence_cmb_vs_raw.png"), p_mh, width = 8, height = 5, dpi = 150)
  cat("  Saved: cooccurrence_cmb_vs_raw.png\n")
}

# Combination distribution by disease
combo_top <- combo_by_disease[RC_combination != "Pg-_Td-_Tf-"]
combo_top <- combo_top[DiseaseStatus %in% names(sort(table(combo_top$DiseaseStatus), decreasing = TRUE))[1:12]]
p_combo <- ggplot(combo_top, aes(x = DiseaseStatus, y = pct, fill = RC_combination)) +
  geom_bar(stat = "identity", position = "fill") + scale_fill_brewer(palette = "Set2") +
  labs(title = "Red Complex Combination Distribution by Disease", x = "", y = "Proportion", fill = "Combination") +
  theme_rc() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTDIR, "cooccurrence_combinations_by_disease.png"), p_combo, width = 12, height = 6, dpi = 150)
cat("  All plots saved.\n")

# ── 7. Key findings ────────────────────────────────────────────
cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\nKEY FINDINGS\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")

strongest <- global_pairs[which.max(phi)]
cat(sprintf("\n  1. Strongest global pair: %s-%s (phi=%.3f, OR=%.1f)\n",
            strongest$species1, strongest$species2, strongest$phi, strongest$or_val))

cat("\n  2. Site-adjusted (CMH) co-occurrence:\n")
for (i in 1:nrow(mh_results)) {
  mr <- mh_results[i]
  cat(sprintf("     %s-%s: Raw OR=%.1f → CMH OR=%.1f (%+.1f%%), CMH p=%s, %s across sites\n",
    substr(mr$species1,1,2), substr(mr$species2,1,2),
    mr$raw_OR, mr$cmh_OR, mr$or_change_pct,
    format.pval(mr$cmh_p, digits=2), mr$homogeneity))
}

cat("\n  3. Does co-occurrence weaken in PD after site adjustment?\n")
for (p in pairs) {
  hc_row <- disease_cmh[DiseaseStatus == "HC" & species1 == p[1] & species2 == p[2]]
  pd_row <- disease_cmh[DiseaseStatus == "Periodontal Disease" & species1 == p[1] & species2 == p[2]]
  if (nrow(hc_row) > 0 && nrow(pd_row) > 0) {
    cat(sprintf("     %s-%s: HC CMH-OR=%.1f, PD CMH-OR=%.1f\n",
      substr(p[1],1,2), substr(p[2],1,2), hc_row$cmh_OR, pd_row$cmh_OR))
  }
}

cat("\nDone.\n")
