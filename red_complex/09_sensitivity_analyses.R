#!/home/software/anaconda3/envs/renv/bin/Rscript
#
# 09_sensitivity_analyses.R (NEW, optimized)
# ==========================================
# Cross-subset sensitivity: run key disease OR in 3 subsets
# Focus on 7 core diseases, logistic only, 3 subsets
#

source("02_unified_statistical_framework.R")
library(data.table); library(ggplot2)
OUTDIR <- "../red_complex"

cat("Loading data...\n")
dt <- load_rc_data()
dtc <- clean_analysis_data(dt)
cat(sprintf("  Full: %d samples\n", nrow(dtc)))

dt_v4 <- dtc[SeqRegion == "V4" & Platform == "Illumina MiSeq"]
dt_saliva <- dtc[Site == "saliva"]

subsets <- list("Full" = dtc, "V4+MiSeq" = dt_v4, "Saliva-only" = dt_saliva)

# 7 core diseases only (skip with-light-N diseases)
key_diseases <- c("Periodontal Disease","OSCC","Sjögren's Syndrome",
  "Breast Cancer","Lung Cancer","Colorectal Cancer","Aphthous","COVID-19","ICU","Dental Caries")

cat("\nRunning disease OR in subsets...\n")
sens_res <- rbindlist(lapply(names(subsets), function(sn) {
  cat(sprintf("  %s...\n", sn))
  sub <- subsets[[sn]]
  rbindlist(lapply(key_diseases, function(d) {
    if (!d %in% unique(sub$DiseaseStatus)) return(NULL)
    subd <- sub[DiseaseStatus %in% c("HC", d)]
    subd[, DiseaseStatus := factor(DiseaseStatus, levels = c("HC", d))]
    n_d <- sum(subd$DiseaseStatus == d)
    if (n_d < 20) return(NULL)
    res <- tryCatch({
      m <- glm(RedComplex_detected ~ DiseaseStatus + AgeGroup + Sex,
               data = subd, family = binomial())
      cn <- paste0("DiseaseStatus", d)
      s <- summary(m)$coefficients[cn, ]
      ci <- suppressMessages(confint.default(m)[cn, ])
      data.table(subset = sn, disease = d, n_disease = n_d, n_hc = sum(subd$DiseaseStatus == "HC"),
        OR = exp(s["Estimate"]), OR_lo = exp(ci[1]), OR_hi = exp(ci[2]), pvalue = s["Pr(>|z|)"])
    }, error = function(e) NULL)
    res
  }), fill = TRUE)
}))

# FDR per subset
sens_res[, p_fdr := p.adjust(pvalue, method = "BH"), by = subset]
sens_res[, sig_fdr := p_fdr < 0.05]

cat("\n  Cross-subset OR comparison:\n")
for (d in unique(sens_res$disease)) {
  sub <- sens_res[disease == d]
  if (nrow(sub) >= 2) {
    cat(sprintf("\n  %s:\n", d))
    print(sub[, .(subset, n_disease, OR = round(OR, 2), p_FDR = format.pval(p_fdr, digits = 2), sig = sig_fdr)])
  }
}

fwrite(sens_res, file.path(OUTDIR, "sensitivity_cross_subset.csv"))

# Robustness score
robust <- sens_res[, .(n_subsets = .N, n_sig = sum(sig_fdr, na.rm = TRUE),
  mean_OR = mean(OR, na.rm = TRUE), min_OR = min(OR, na.rm = TRUE), max_OR = max(OR, na.rm = TRUE),
  OR_cv = sd(OR, na.rm = TRUE) / abs(mean(OR, na.rm = TRUE))), by = disease]
robust[, robustness := fcase(
  n_subsets >= 2 & n_sig == n_subsets & OR_cv < 0.5, "A: Consistent",
  n_subsets >= 2 & n_sig >= n_subsets - 1 & OR_cv < 1.0, "B: Robust",
  n_sig >= 1, "C: Variable", default = "D: Unstable")]
setorder(robust, robustness, -mean_OR)
cat("\n  Robustness scores:\n")
print(robust[, .(disease, n_subsets, n_sig, mean_OR = round(mean_OR, 2), CV = round(OR_cv, 2), robustness)])
fwrite(robust, file.path(OUTDIR, "sensitivity_robustness_scores.csv"))

# Heatmap
p1 <- ggplot(sens_res, aes(x = subset, y = disease, fill = log2(OR) * ifelse(sig_fdr, 1, 0.3))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#2980b9", mid = "white", high = "#c0392b", midpoint = 0) +
  geom_text(aes(label = sprintf("%.1f%s", OR, ifelse(sig_fdr, "*", ""))), size = 3) +
  labs(title = "Sensitivity: Disease Detection OR Across Subsets",
       subtitle = "* = FDR < 0.05; color = log2(OR) direction", x = "", y = "", fill = "log2(OR)") + theme_rc()
ggsave(file.path(OUTDIR, "sensitivity_heatmap.png"), p1, width = 8, height = 8, dpi = 150)

cat("\nDone.\n")
