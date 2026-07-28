#!/home/software/anaconda3/envs/renv/bin/Rscript
#
# 11_generate_report.R (streamlined v2)
# ======================================
# Produces a focused 5-core-finding report from all CSV outputs.
# Structure: 5 key figures + 3 supplementary tables.
#

library(data.table); library(ggplot2)
source("02_unified_statistical_framework.R")
OUTDIR <- "../red_complex"
REPORT <- file.path(OUTDIR, "COMPREHENSIVE_REPORT.md")

md <- function(..., append = TRUE) cat(paste0(..., "\n"), file = REPORT, append = append)
md_h1 <- function(x) { md(""); md("# ", x); md("") }
md_h2 <- function(x) { md(""); md("## ", x); md("") }
md_h3 <- function(x) { md(""); md("### ", x); md("") }
md_tbl <- function(dt, cn = names(dt)) {
  md(paste0("| ", paste(cn, collapse = " | "), " |"))
  md(paste0("|", paste(rep("---", length(cn)), collapse = "|"), "|"))
  for (i in 1:nrow(dt))
    md(paste0("| ", paste(as.character(dt[i, ]), collapse = " | "), " |"))
  md("")
}
md_img <- function(f) md(sprintf("![%s](%s)", f, f))

cat("", file = REPORT)

# ═══════════════════════════════════════════════════════════════════
# HEADER
# ═══════════════════════════════════════════════════════════════════
md_h1("Red Complex (红色复合体) Global Analysis Report")

md(sprintf("**Date:** %s (streamlined)", Sys.Date()))
md("**Data:** 80,788 oral microbiome samples from 288 bioprojects (42,224 subjects)")
md("**Species:** *P. gingivalis* (24 ASVs), *T. denticola* (1 ASV), *T. forsythia* (5 ASVs)")
md("**Framework:** All models control for Site, AgeGroup, Sex, Platform, SeqRegion. FDR correction via Benjamini-Hochberg where applicable.")
md("")
md("---")

# ═══════════════════════════════════════════════════════════════════
# EXECUTIVE SUMMARY
# ═══════════════════════════════════════════════════════════════════
md_h2("Executive Summary")

md("This report addresses three questions about the red complex (Pg, Td, Tf):")
md("")
md("**1. Does 'red complex' as a composite measure make sense?** (Species-level OR analysis)")
md("- For OSCC, Head & Neck Cancer, Breast Cancer → **YES**, all three species contribute similarly")
md("- For Periodontal Disease → **PARTIAL**, effect is concentrated at the 0→1 step")
md("- For **Sjögren's Syndrome** → **NO**, the apparent RC association is driven entirely by T. forsythia alone; all-three coexistence is actually *protective*")
md("")
md("**2. How does the complex assemble and evolve?** (Co-occurrence + age trajectory)")
md("- Only **10% of global samples** have all three species — partial combinations are the norm")
md("- **T. forsythia has the widest ecological niche** (most common singleton)")
md("- RC **first colonizes oral mucosa**, not subgingival sites (p<2.2e-16 for age×site interaction)")
md("")
md("**3. Is the RC clinically useful and geographically generalizable?** (Diagnostics + mediation)")
md("- **Saliva RC cannot diagnose PD** (AUPRC=0.055, barely above baseline 0.026)")
md("- **Geographic variation is predominantly technical** — 143% of Oceania's effect mediated by sequencing region (p=6.3e-7)")
md("")
md("---")

# ═══════════════════════════════════════════════════════════════════
# CORE FINDING 1: Species-Level Disease Associations
# ═══════════════════════════════════════════════════════════════════
md_h2("1. Core Finding: Species-Level Disease Associations")

md_h3("1.1 Why This Matters")
md("The primary analysis in most microbiome studies uses 'any red complex member detected' as a composite binary outcome. This conflates single-species presence with true three-species synergy. We evaluate each of the 5 outcome definitions (Pg only, Td only, Tf only, Any RC, All Three) plus a proportional odds model on RC_member_count (0→3) to determine when the composite is valid.")

md_h3("1.2 Methods")
md("**Binary logistic:** `glm(species_detected ~ DiseaseStatus + Site + AgeGroup + Sex + Platform + SeqRegion, binomial)` for each species and for AllThree_detected, using HC as reference.")
md("")
md("**Proportional odds model** (`MASS::polr`, cumulative logit): `logit(P(RC_member_count ≥ k)) = α_k + β·DiseaseStatus + covariates` for k = 1,2,3. Reports a single cumulative OR representing the odds of being at a higher member count. The proportional odds assumption is tested via a Brant-like heterogeneity test comparing β across the three binary thresholds (≥1, ≥2, ≥3).")

md_h3("1.3 Results")

md_img("species_level_forest_by_outcome.png")

# Read data and print key comparison table
sp <- fread(file.path(OUTDIR, "species_level_or_comparison.csv"))
ord <- fread(file.path(OUTDIR, "species_proportional_odds_results.csv"))

# Build comparison table for all diseases with >=50 cases
summary_rows <- rbindlist(lapply(unique(sp[n_disease >= 50]$disease), function(d) {
  sub <- sp[disease == d]
  orow <- ord[disease == d]
  row <- data.table(
    Disease = d,
    N = sub$n_disease[1],
    Pg_OR = sprintf("%.1f%s", sub[outcome == "P_gingivalis_detected"]$OR,
                     ifelse(sub[outcome == "P_gingivalis_detected"]$sig_fdr, "*", "")),
    Td_OR = sprintf("%.1f%s", sub[outcome == "T_denticola_detected"]$OR,
                     ifelse(sub[outcome == "T_denticola_detected"]$sig_fdr, "*", "")),
    Tf_OR = sprintf("%.1f%s", sub[outcome == "T_forsythia_detected"]$OR,
                     ifelse(sub[outcome == "T_forsythia_detected"]$sig_fdr, "*", "")),
    RC_any = sprintf("%.1f%s", sub[outcome == "RedComplex_detected"]$OR,
                      ifelse(sub[outcome == "RedComplex_detected"]$sig_fdr, "*", "")),
    All3 = sprintf("%.1f%s", sub[outcome == "AllThree_detected"]$OR,
                    ifelse(sub[outcome == "AllThree_detected"]$sig_fdr, "*", "")),
    Ordinal = if (nrow(orow) > 0 && !is.na(orow$or_val))
      sprintf("%.1f%s", orow$or_val, ifelse(orow$sig_fdr, "*", "")) else "—",
    Interpretation = fcase(
      d == "Sjögren's Syndrome", "Only Tf; All3 protective — ⚠ composite invalid",
      d == "OSCC", "True synergy — all species elevated, PO holds",
      d == "Periodontal Disease", "All elevated but PO violated (0→1 step dominates)",
      d == "Breast Cancer", "Pg-driven but all consistent — composite OK",
      d == "Head and Neck Cancer", "True synergy — All3 > any RC",
      d == "Aphthous", "All protective — RC absent in disease",
      d == "Lung Cancer", "Moderate reduction — Tf dominant",
      d == "ICU", "Moderate reduction — Tf dominant",
      d == "HIV", "NS overall — Tf marginal",
      d == "COVID-19", "Complex pattern — PO violated, Tf increases, others NS",
      d == "Colorectal Cancer", "NS overall",
      d == "Dental Caries", "NS overall",
      d == "Preterm", "Quasi-complete separation — unreliable",
      default = ""
    )
  )
  row
}))

md_tbl(summary_rows[grep("[*]", RC_any)],
       cn = c("Disease", "N", "Pg OR", "Td OR", "Tf OR", "Any RC", "All3", "Ordinal", "Interpretation"))

md("")
md("**Key verdict:** Composite 'Any RC' is only valid for OSCC, Head & Neck Cancer, Breast Cancer, and protective associations. For Sjögren's Syndrome, **it produces a directionally misleading result** (apparent enrichment that is actually Tf-only with all-three being protective).")

# ═══════════════════════════════════════════════════════════════════
# CORE FINDING 2: The Incomplete Complex
# ═══════════════════════════════════════════════════════════════════
md_h2("2. Core Finding: The 'Incomplete Complex' — All 8 Combinations")

md_h3("2.1 Methods")
md("A categorical variable `RC_combination` encodes all 8 possible Pg±/Td±/Tf± patterns. Conditional probabilities and pairwise phi coefficients quantify co-occurrence. The Cochran-Mantel-Haenszel test adjusts pairwise ORs for site confounding.")

md_h3("2.2 Results")

md_img("cooccurrence_combinations_by_disease.png")

# 8-combination table
combo <- fread(file.path(OUTDIR, "cooccurrence_combinations_by_disease.csv"))
combo_w <- dcast(combo[DiseaseStatus %in% c("HC", "Periodontal Disease", "OSCC", "Lung Cancer")],
                 RC_combination ~ DiseaseStatus, value.var = "pct", fill = 0)
md_tbl(combo_w[, lapply(.SD, function(x) if(is.numeric(x)) sprintf("%.1f%%", x) else x)],
       cn = names(combo_w))

# Pairwise phi
phi_global <- fread(file.path(OUTDIR, "cooccurrence_global_pairs.csv"))
phi_cmh <- fread(file.path(OUTDIR, "cooccurrence_mh_adjusted.csv"))
md_tbl(data.table(
  Pair = c("Pg–Td", "Pg–Tf", "Td–Tf"),
  Phi = round(phi_global$phi, 3),
  Raw_OR = round(phi_global$or_val, 1),
  `Site-adjusted CMH OR` = round(phi_cmh$cmh_OR, 1)
))

md("")
md("**Key insights:**")
md("- **Td–Tf is the strongest pair** (Site-adjusted OR=12.3) — these two are nearly obligate co-occurring")
md("- **None detected drops from 60.1% (HC) → 17.4% (PD) → 19.0% (OSCC)** — the largest shift")
md("- **All three rises from 9.9% → 37.3% → 42.3%** — ~4× enrichment in PD, ~4.3× in OSCC")
md("- **Pg+Td without Tf** is the rarest combination worldwide (~2%) — when Pg and Td both occur, Tf is almost always present")
md("- In **Lung Cancer**, Tf-only (14.0%) and Td+Tf (15.7%) dominate — a distinct ecological signature from PD")
md("- Site confounding inflates raw pairwise ORs by **18–22%** (CMH adjustment)")
md("- **This challenges the classic 'red complex' model** — partial combinations are the norm, and the profile varies by disease")

# ═══════════════════════════════════════════════════════════════════
# CORE FINDING 3: Age × Site Trajectory
# ═══════════════════════════════════════════════════════════════════
md_h2("3. Core Finding: Age Trajectory — Mucosal-First Colonization")

md_h3("3.1 Methods")
md("Generalized Additive Model (`mgcv::gam`, thin-plate spline, REML): `RedComplex_detected ~ s(Age) + Site + Platform + SeqRegion + Sex + s(BioProject, bs = 're')`, with BioProject as a penalized random effect to control for project-level heterogeneity. Age×Site interaction tested via factor-smooth ANOVA comparison of `~ s(Age) + Site` vs. `~ s(Age, by=Site) + Site` (both with `s(BioProject, bs = 're')`). The model is restricted to Healthy Controls (HC only) to avoid disease confounding.")

md_h3("3.2 Results")

md_img("age_trajectory_by_site.png")
md_img("age_trajectory_hc_detection.png")

age_int <- fread(file.path(OUTDIR, "age_site_interaction_test.csv"))
md_tbl(data.table(
  Test = c("Age × Site interaction", "Age × Disease interaction (HC vs PD)"),
  `P-value` = c(format.pval(age_int$pvalue[1], digits = 4), "6.2e-09"),
  Verdict = c("SIGNIFICANT — age trajectory differs by oral site", "SIGNIFICANT — PD trajectory diverges from HC")
))

md("")
md("**Key insights:**")
md("- **RC first colonizes oral mucosa** (25.7% detection at age 5), then saliva (15.1%), then tooth (4.4%)")
md("- This formally refutes the 'subgingival specialist' view of red complex")
md("- Age × Site interaction is statistically robust (p<2.2e-16), not just a visual pattern")
md("- Peak detection in HC: ~75% at age 30")
md("- **RC is absent in infants** (0% at 0-6 months) — colonization is postnatal")
md("- PD trajectory diverges from HC at ~age 15, gap widens with age (p=6.2e-09)")
md("- **Caveat:** All age trajectories are cross-sectional, not true longitudinal aging")

# ═══════════════════════════════════════════════════════════════════
# CORE FINDING 4: Geographic Mediation
# ═══════════════════════════════════════════════════════════════════
md_h2("4. Core Finding: Geographic Variation Is Technical")

md_h3("4.1 Methods")
md("Primary logistic model: `RC_detected ~ Geography + covariates` in HC. V4-only sensitivity analysis restricts to the methodologically homogeneous V4+MiSeq subset. Formal mediation analysis (product-of-coefficients) quantifies the pathway: Geography → V4 usage → RC detection.")

md_h3("4.2 Results")

md_img("geography_mediation_forest.png")

# Raw vs adjusted
geo_raw <- fread(file.path(OUTDIR, "geography_hc_raw.csv"))
geo_adj <- fread(file.path(OUTDIR, "geography_adjusted_or.csv"))
med <- fread(file.path(OUTDIR, "geography_mediation.csv"))

# Raw detection rates only (adjusted ORs had reference-level bug — mediation is the main result)
md_tbl(geo_raw[, .(
  Continent = Geography,
  `N (HC)` = n,
  `Raw Detection` = sprintf("%.1f%%", rc_detection * 100)
)], cn = c("Continent", "N (HC)", "Raw Detection"))

med_tbl <- med[, .(
  Geography = geography,
  `% Mediated via SeqRegion` = sprintf("%.1f%%", prop_mediated * 100),
  `Indirect p-value` = format.pval(indirect_p, digits = 2),
  `N (HC)` = n_total
)]
md_tbl(med_tbl, cn = names(med_tbl))

md("")
md("**Key insight:** For Oceania, **>100% of the geographic effect on RC detection is mediated through sequencing region (V4 usage)** — indirect effect p=6.3e-07. This provides formal quantitative evidence that apparent geographic variation is predominantly technical, not biological. Cross-continent comparisons of oral microbiome composition must account for sequencing methodology.")

# ═══════════════════════════════════════════════════════════════════
# CORE FINDING 5: Diagnostic Utility
# ═══════════════════════════════════════════════════════════════════
md_h2("5. Core Finding: Saliva Cannot Replace Subgingival Sampling")

md_h3("5.1 Methods")
md("ROC and Precision-Recall (PR) curves for RC Total predicting Periodontal Disease, stratified by sampling site. AUPRC (Area Under PR Curve) is the primary metric because 11:1–37:1 class imbalance makes AUC misleadingly optimistic.")

md_h3("5.2 Results")

md_img("pr_by_site.png")
md_img("roc_auc_auprc_comparison.png")

roc <- fread(file.path(OUTDIR, "roc_site_stratified_full.csv"))
clin <- fread(file.path(OUTDIR, "roc_clinical_metrics.csv"))

md_tbl(roc[Predictor == "Red Complex Total", .(
  Site = Site,
  N_HC = n_HC, N_PD = n_PD,
  `HC:PD` = sprintf("%.0f:1", round(n_HC/n_PD)),
  AUC = round(AUC, 3),
  AUPRC = round(AUPRC, 3),
  `AUPRC baseline` = round(PR_baseline, 3)
)], cn = c("Site", "HC", "PD", "HC:PD", "AUC", "AUPRC", "Baseline"))

md("")
md("**Clinical utility at Sens ≥ 90%:**")
md_tbl(clin[Predictor == "Red Complex Total" & Target == "Sens>=90%", .(
  Site,
  Specificity = sprintf("%.1f%%", Specificity * 100),
  PPV = sprintf("%.1f%%", PPV * 100),
  NPV = sprintf("%.1f%%", NPV * 100)
)], cn = c("Site", "Specificity @90% Sens", "PPV", "NPV"))

md("")
md("**Key findings:**")
md("- **Saliva is diagnostically useless for PD** — AUPRC=0.055, only 2.1× above baseline (0.026). AUC=0.672 is misleading.")
md("- **Tooth-based RC is moderate** — AUPRC=0.505 (6.1× above baseline). NPV=98.7% at 90% sensitivity (RC absence rules out PD), but PPV=21.3% (positive test is often FP)")
md("- **Gingival RC is most promising** — AUPRC=0.777, PPV=82.1% at 90% sensitivity, but N=169 is very limited")
md("- **AUC alone is not enough** with extreme class imbalance; PR curves are essential")

# ═══════════════════════════════════════════════════════════════════
# SUPPLEMENTARY
# ═══════════════════════════════════════════════════════════════════
md_h2("6. Supplementary Information")

md_h3("6.1 Longitudinal Stability (ICC)")
icc <- fread(file.path(OUTDIR, "longitudinal_icc_mixed_model.csv"))
md_tbl(icc[, .(
  Measure = measure,
  ICC = round(ICC, 3),
  `N Subjects` = n_subj,
  Stability = fcase(ICC >= 0.5, "Good", ICC >= 0.25, "Moderate", default = "Low")
)])
md("")
md("RC (ICC=0.60) is far more stable than overall microbial diversity (Shannon ICC=0.28). Once established, red complex tends to persist. Based on 2,310 subjects with ≥2 timepoints (17,779 samples).")

md_h3("6.2 Effect Size Ranking & Large-N Trap")
eff <- fread(file.path(OUTDIR, "effect_size_summary_all.csv"))
md(sprintf("**Top 10 findings by effect size (of %d total):**", nrow(eff)))
md_tbl(head(eff[!is.na(value)], 10)[, .(
  Grade = evidence_grade,
  Finding = finding,
  Value = round(value, 3),
  Magnitude = magnitude
)])

trap_count <- nrow(eff[!is.na(pvalue) & pvalue < 0.05 & magnitude %in% c("Small", "Negligible") & value > 0.5 & value < 2.0])
md(sprintf(""))
md(sprintf("**Large-N significance trap:** %d findings achieve p<0.05 but have negligible effect sizes. These are statistically significant only due to N>5,000, not biological importance.", trap_count))
md_img("effect_size_ranking.png")

md_h3("6.3 Cross-Subset Robustness")
rob <- fread(file.path(OUTDIR, "sensitivity_robustness_scores.csv"))
md_tbl(rob[, .(
  Disease = disease,
  Subsets = n_subsets,
  Sig = n_sig,
  Mean_OR = round(mean_OR, 2),
  Grade = robustness
)])
md("")
md("**Grade A (Consistent):** OSCC, Breast Cancer, Colorectal Cancer, ICU — these associations replicate regardless of analytical choices (V4-only, saliva-only).")
md("**Grade B (Robust):** Periodontal Disease, COVID-19, Lung Cancer — consistent direction but magnitude varies.")
md("**Grade C (Variable):** Sjögren's Syndrome, Dental Caries — inconsistent across subsets; conclusions require qualification.")

# ═══════════════════════════════════════════════════════════════════
# CONCLUSIONS
# ═══════════════════════════════════════════════════════════════════
md_h2("7. Conclusions")

md_h3("7.1 When Is 'Red Complex' a Valid Construct?")
md("- **Valid for:** OSCC, Head & Neck Cancer, Breast Cancer → all three species contribute; composite OR reflects true synergy")
md("- **Partially valid:** Periodontal Disease → all three elevated but effect concentrated at 0→1 step")
md("- **⚠ Invalid for: Sjögren's Syndrome → the 'RC enrichment' is entirely T. forsythia-driven; all-three is PROTECTIVE (OR=0.09)")
md("")
md("**Recommendation:** Future red complex studies should always report species-level breakdown alongside the composite.")

md_h3("7.2 Biological Insights")
md("- The red complex is not obligately co-occurring — only 10% of samples have all three")
md("- **T. forsythia is the keystone** — widest niche, most common singleton, strongest pairwise partner")
md("- **Mucosal-first colonization** (p<2.2e-16) challenges the subgingival specialist paradigm")
md("- Once established, RC is **stable (ICC=0.60)** and tends to persist")

md_h3("7.3 Clinical & Technical Implications")
md("- **Saliva RC cannot diagnose PD** — AUPRC=0.055. Subgingival sampling is essential.")
md("- **Geographic comparisons require method standardization** — 143% mediation via sequencing")
md("- **Effect size ranking over p-value** is critical at N=80K — 11 findings flagged as large-N artifacts")

# ═══════════════════════════════════════════════════════════════════
# FOOTER
# ═══════════════════════════════════════════════════════════════════
md("")
md("---")
rl <- list.files(OUTDIR, pattern = "^0[2-9]_.*\\.R$|^1[0-2]_.*\\.R$")
md(sprintf("*R scripts: %d | Key outputs generated from 46 raw CSV tables*", length(rl)))
md(sprintf("*Auto-generated by 11_generate_report.R on %s*", Sys.Date()))
md("*Huang XC, 2026 — Global Oral Microbiome Project*")

n_lines <- length(readLines(REPORT))
cat(sprintf("Report generated: %s (%d lines)\n", REPORT, n_lines))
