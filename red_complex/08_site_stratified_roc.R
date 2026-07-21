#!/home/software/anaconda3/envs/renv/bin/Rscript
#
# 08_site_stratified_roc.R — RIGOR-ENHANCED v2
# ==============================================
# New in v2:
#   + Precision-Recall (PR) curves with AUPRC
#   + Calibration assessment (Brier, Hosmer-Lemeshow)
#   + PPV/NPV at Youden, Sens>=90%, Sens>=95% thresholds
#   + Class imbalance diagnostics
#

source("02_unified_statistical_framework.R")
library(data.table)
library(ggplot2)
OUTDIR <- "../red_complex"

cat("Loading data...\n")
dt <- load_rc_data()
dtc <- clean_analysis_data(dt)
dt_roc <- dtc[DiseaseStatus %in% c("HC", "Periodontal Disease")]
dt_roc[, DiseaseStatus := factor(DiseaseStatus, levels = c("HC", "Periodontal Disease"))]
dt_roc[, site_char := as.character(Site)]

site_n <- dt_roc[, .N, by = .(site_char, DiseaseStatus)]
cat("\n  Site x Disease sample sizes:\n")
print(dcast(site_n, site_char ~ DiseaseStatus, value.var = "N"))

valid_sites <- c("tooth", "saliva", "gingiva")
dt_roc <- dt_roc[site_char %in% valid_sites]
dt_roc[, is_pd := DiseaseStatus == "Periodontal Disease"]
cat(sprintf("  Analysis set: %d HC + %d PD across %d sites\n",
            sum(!dt_roc$is_pd), sum(dt_roc$is_pd), length(valid_sites)))

cat("\n── Class Imbalance ──\n")
for (st in valid_sites) {
  sub <- dt_roc[site_char == st]
  n_hc <- sum(!sub$is_pd); n_pd <- sum(sub$is_pd)
  cat(sprintf("  %-10s HC=%d PD=%d  prevalence=%.3f  imbalance=%.1f:1\n",
              st, n_hc, n_pd, n_pd/(n_hc+n_pd), n_hc/n_pd))
}

calc_all_metrics <- function(pred, truth) {
  truth <- as.integer(as.logical(truth))
  n_pos <- sum(truth == 1, na.rm = TRUE)
  n_neg <- sum(truth == 0, na.rm = TRUE)
  if (n_pos < 5 || n_neg < 5) {
    return(list(auc = NA, auc_se = NA, auc_low = NA, auc_high = NA,
                auprc = NA, pr_baseline = NA, brier = NA, brier_scaled = NA,
                tpr = c(0, 1), fpr = c(0, 1), precision = c(1, 0), recall = c(0, 1),
                n_pos = n_pos, n_neg = n_neg))
  }
  ord <- order(pred, decreasing = TRUE)
  truth_ord <- truth[ord]
  tpr <- cumsum(truth_ord) / n_pos
  fpr <- cumsum(1 - truth_ord) / n_neg
  auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
  q1 <- auc / (2 - auc); q2 <- 2 * auc^2 / (1 + auc)
  auc_se <- sqrt((auc*(1-auc) + (n_pos-1)*(q1-auc^2) + (n_neg-1)*(q2-auc^2)) / (n_pos*n_neg))
  tp <- cumsum(truth_ord)
  precision <- tp / seq_len(length(tp))
  recall <- tp / n_pos
  pr_baseline <- n_pos / (n_pos + n_neg)
  rc_ord <- order(recall); rc_u <- recall[rc_ord]; pr_u <- precision[rc_ord]
  keep <- !duplicated(rc_u); rc_u <- rc_u[keep]; pr_u <- pr_u[keep]
  auprc <- if (length(rc_u) >= 2) sum(diff(rc_u)*(pr_u[-1]+pr_u[-length(pr_u)])/2, na.rm=TRUE) else NA
  pc <- pmax(1e-10, pmin(1-1e-10, pred))
  brier <- mean((pc - truth)^2)
  brier_scaled <- 1 - brier/mean((mean(truth)-truth)^2)
  list(auc=auc, auc_se=auc_se, auc_low=max(0,auc-1.96*auc_se), auc_high=min(1,auc+1.96*auc_se),
       auprc=auprc, pr_baseline=pr_baseline, brier=brier, brier_scaled=brier_scaled,
       tpr=tpr, fpr=fpr, precision=precision, recall=recall, n_pos=n_pos, n_neg=n_neg)
}

predictors <- list("Red Complex Total"="RedComplexTotal","P. gingivalis"="P_gingivalis",
                   "T. denticola"="T_denticola","T. forsythia"="T_forsythia",
                   "All Three (0/1)"="AllThree_detected")

cat("\n── 2. ROC + PR + Calibration ──\n")
all_results <- list(); all_curves <- list()
for (st in valid_sites) {
  sub <- dt_roc[site_char == st]
  for (pname in names(predictors)) {
    pcol <- predictors[[pname]]
    pred_raw <- sub[[pcol]]
    if (pcol == "AllThree_detected") {
      pred_prob <- as.numeric(pred_raw)
    } else {
      pred_prob <- (rank(pred_raw, ties.method="average") - 0.5) / length(pred_raw)
    }
    res <- calc_all_metrics(pred_prob, sub$is_pd)
    all_results[[paste(st,pname)]] <- data.table(
      Site=st, Predictor=pname, AUC=res$auc, AUC_SE=res$auc_se,
      AUC_low=res$auc_low, AUC_high=res$auc_high, AUPRC=res$auprc,
      PR_baseline=res$pr_baseline, Brier=res$brier, Brier_scaled=res$brier_scaled,
      n_HC=res$n_neg, n_PD=res$n_pos, prevalence=res$n_pos/(res$n_pos+res$n_neg))
    n_pts <- min(150, length(res$fpr))
    idx <- floor(seq(1, length(res$fpr), length.out=n_pts))
    all_curves[[paste(st,pname,"roc")]] <- data.table(
      Site=st, Predictor=pname, Curve="ROC", X=res$fpr[idx], Y=res$tpr[idx], AUC=res$auc)
    ok <- is.finite(res$precision) & !is.na(res$precision)
    if (sum(ok) > 1) {
      n_pr <- min(150, sum(ok)); pr_idx <- floor(seq(1, sum(ok), length.out=n_pr))
      all_curves[[paste(st,pname,"pr")]] <- data.table(
        Site=st, Predictor=pname, Curve="PR", X=res$recall[ok][pr_idx],
        Y=res$precision[ok][pr_idx], AUC=res$auprc)
    }
  }
}
roc_results <- rbindlist(all_results, fill=TRUE)
roc_curves <- rbindlist(all_curves, fill=TRUE)
roc_results[, AUC_CI := sprintf("%.3f [%.3f-%.3f]", AUC, AUC_low, AUC_high)]

cat("\n  Full metrics:\n")
for (st in valid_sites) {
  cat(sprintf("\n  --- %s ---\n", st))
  sub <- roc_results[Site==st][order(-AUC)]
  print(sub[,.(Predictor, AUC_CI, AUPRC=round(AUPRC,3), Brier=round(Brier,4),
               Brier_scaled=round(Brier_scaled,3), Prevalence=round(prevalence,3))])
}
fwrite(roc_results, file.path(OUTDIR, "roc_site_stratified_full.csv"))
fwrite(roc_curves, file.path(OUTDIR, "roc_pr_curve_data.csv"))

cat("\n── 3. Clinical Metrics ──\n")
clinical_metrics <- rbindlist(lapply(valid_sites, function(st) {
  sub <- dt_roc[site_char == st]
  rbindlist(lapply(names(predictors), function(pname) {
    pcol <- predictors[[pname]]; pred <- sub[[pcol]]; truth <- as.integer(sub$is_pd)
    n_pos <- sum(truth==1); n_neg <- sum(truth==0)
    if (n_pos<5 || n_neg<5) return(data.table())
    ord <- order(pred, decreasing=TRUE); t_ord <- truth[ord]; p_ord <- pred[ord]
    tpr <- cumsum(t_ord)/n_pos; fpr <- cumsum(1-t_ord)/n_neg
    yi <- which.max(tpr - fpr)
    get_m <- function(idx) {
      tp_at <- sum(t_ord[1:idx]); fp_at <- idx - tp_at
      fn_at <- n_pos - tp_at; tn_at <- n_neg - fp_at
      data.table(Sensitivity=tpr[idx], Specificity=1-fpr[idx],
        PPV=if(idx>0) tp_at/idx else NA,
        NPV=if(idx<length(t_ord)) tn_at/(length(t_ord)-idx) else NA,
        Threshold=if(idx<=length(p_ord)) p_ord[idx] else NA)
    }
    i90 <- which(tpr>=0.90)[1]; i95 <- which(tpr>=0.95)[1]
    if(is.na(i90)) i90 <- length(t_ord); if(is.na(i95)) i95 <- length(t_ord)
    rbind(
      cbind(data.table(Site=st,Predictor=pname,Target="Youden"), get_m(yi)),
      cbind(data.table(Site=st,Predictor=pname,Target="Sens>=90%"), get_m(i90)),
      cbind(data.table(Site=st,Predictor=pname,Target="Sens>=95%"), get_m(i95)))
  }))
}))
cat("\n  RC Total clinical performance:\n")
for (st in valid_sites) {
  cat(sprintf("\n  --- %s ---\n", st))
  sub <- clinical_metrics[Site==st & Predictor=="Red Complex Total"]
  print(sub[,.(Target, Sensitivity=round(Sensitivity,3), Specificity=round(Specificity,3),
               PPV=round(PPV,3), NPV=round(NPV,3), Threshold=signif(Threshold,4))])
}
fwrite(clinical_metrics, file.path(OUTDIR, "roc_clinical_metrics.csv"))

cat("\n── 4. Calibration (Hosmer-Lemeshow) ──\n")
calib_data <- rbindlist(lapply(valid_sites, function(st) {
  sub <- dt_roc[site_char == st]
  rbindlist(lapply(names(predictors), function(pname) {
    pcol <- predictors[[pname]]; pred_raw <- sub[[pcol]]
    pred_prob <- if(pcol=="AllThree_detected") as.numeric(pred_raw) else
      (rank(pred_raw, ties.method="average")-0.5)/length(pred_raw)
    truth <- as.integer(sub$is_pd); g <- 10
    breaks <- quantile(pred_prob, probs=seq(0,1,length.out=g+1), type=2)
    breaks[1] <- -Inf; breaks[g+1] <- Inf
    bin_id <- findInterval(pred_prob, breaks, rightmost.closed=TRUE)
    dt_bins <- data.table(Site=st, Predictor=pname, Bin=1:g,
      pred_mean=as.numeric(tapply(pred_prob, bin_id, mean)),
      obs_rate=as.numeric(tapply(truth, bin_id, mean)),
      n_bin=as.integer(tapply(truth, bin_id, length)),
      exp_pos=as.numeric(tapply(pred_prob, bin_id, sum)),
      obs_pos=as.integer(tapply(truth, bin_id, sum)))
    dt_bins[, exp_neg := n_bin - exp_pos]; dt_bins[, obs_neg := n_bin - obs_pos]
    vb <- dt_bins[exp_pos>1 & exp_neg>1]
    if(nrow(vb)>=2) {
      hl_stat <- sum((vb$obs_pos-vb$exp_pos)^2/vb$exp_pos + (vb$obs_neg-vb$exp_neg)^2/vb$exp_neg)
      hl_df <- nrow(vb)-2; hl_p <- pchisq(hl_stat, hl_df, lower.tail=FALSE)
    } else { hl_stat <- NA; hl_df <- NA; hl_p <- NA }
    dt_bins[, `:=`(HL_stat=hl_stat, HL_df=hl_df, HL_p=hl_p)]; dt_bins
  }))
}))
hl_s <- unique(calib_data[,.(Site,Predictor,HL_stat,HL_df,HL_p)])
hl_s[, calibration := ifelse(!is.na(HL_p)&HL_p>0.05,"OK","POOR")]
cat("\n  Hosmer-Lemeshow:\n")
print(hl_s[,.(Site,Predictor,HL_stat=round(HL_stat,1),HL_df,
              HL_p=format.pval(HL_p,digits=2),calibration)])
fwrite(calib_data, file.path(OUTDIR, "roc_calibration_data.csv"))

cat("\n── 5. Plots ──\n")
p1 <- ggplot(roc_curves[Curve=="ROC"], aes(x=X, y=Y, color=Predictor)) +
  geom_line(linewidth=0.8) + geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  facet_wrap(~Site) + scale_color_brewer(palette="Set1") +
  labs(title="ROC Curves: Red Complex for Periodontal Disease", x="FPR", y="TPR", color="") +
  coord_equal() + theme_rc()
ggsave(file.path(OUTDIR,"roc_by_site.png"), p1, width=14, height=5, dpi=150)

pr_baselines <- roc_results[,.(baseline=PR_baseline[1]),by=Site]
p2 <- ggplot(roc_curves[Curve=="PR"], aes(x=X, y=Y, color=Predictor)) +
  geom_line(linewidth=0.8) + facet_wrap(~Site) + scale_color_brewer(palette="Set1") +
  geom_hline(data=pr_baselines, aes(yintercept=baseline), linetype="dashed", color="grey50") +
  labs(title="Precision-Recall Curves: Red Complex for PD", x="Recall", y="Precision", color="") +
  xlim(0,1) + ylim(0,1) + theme_rc() + theme(legend.position="bottom")
ggsave(file.path(OUTDIR,"pr_by_site.png"), p2, width=14, height=5, dpi=150)

mm <- melt(roc_results, id.vars=c("Site","Predictor"), measure.vars=c("AUC","AUPRC"),
           variable.name="Metric", value.name="Value")
p3 <- ggplot(mm, aes(x=Predictor, y=Value, fill=Site)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  facet_wrap(~Metric, scales="free_y") + scale_fill_brewer(palette="Set2") +
  geom_text(aes(label=sprintf("%.3f",Value)), position=position_dodge(0.7), vjust=-0.5, size=2.8) +
  labs(title="AUC and AUPRC: Red Complex for Periodontal Disease", x="", y="") + theme_rc()
ggsave(file.path(OUTDIR,"roc_auc_auprc_comparison.png"), p3, width=12, height=6, dpi=150)

cp <- calib_data[!is.na(pred_mean)]
p4 <- ggplot(cp, aes(x=pred_mean, y=obs_rate, color=Predictor)) +
  geom_point(aes(size=n_bin), alpha=0.7) + geom_line(linewidth=0.6) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  facet_wrap(~Site) + scale_color_brewer(palette="Set1") + scale_size_continuous(range=c(1,5)) +
  labs(title="Calibration: Predicted vs Observed PD Rate", x="Predicted probability",
       y="Observed rate", color="", size="N/bin") +
  xlim(0,1) + ylim(0,1) + coord_equal() + theme_rc() +
  theme(legend.position="bottom", legend.text=element_text(size=6))
ggsave(file.path(OUTDIR,"roc_calibration.png"), p4, width=14, height=6, dpi=150)

rcp <- rbindlist(lapply(valid_sites, function(st) {
  sub <- dt_roc[site_char==st]; pred <- sub[["RedComplexTotal"]]; truth <- as.integer(sub$is_pd)
  ord <- order(pred, decreasing=TRUE); t_ord <- truth[ord]
  tpr <- cumsum(t_ord)/sum(t_ord); ppv <- cumsum(t_ord)/seq_len(length(t_ord))
  idx <- floor(seq(1,length(tpr),length.out=min(500,length(tpr))))
  data.table(Site=st, TPR=tpr[idx], PPV=ppv[idx])
}))
p5 <- ggplot(rcp, aes(x=TPR, y=PPV, color=Site)) + geom_line(linewidth=0.8) +
  geom_hline(data=roc_results[Predictor=="Red Complex Total"],
             aes(yintercept=prevalence, color=Site), linetype="dashed") +
  labs(title="PPV vs Sensitivity Trade-off (RC Total)", x="Sensitivity", y="PPV", color="Site") +
  xlim(0,1) + ylim(0,1) + theme_rc()
ggsave(file.path(OUTDIR,"roc_ppv_vs_sensitivity.png"), p5, width=8, height=5, dpi=150)
cat("  All plots saved.\n")

cat("\n" %+% paste(rep("=",60),collapse="") %+% "\nKEY FINDINGS\n" %+% paste(rep("=",60),collapse="") %+% "\n")
best_auc <- roc_results[order(-AUC)][1]
best_auprc <- roc_results[order(-AUPRC)][1]
cat(sprintf("\n  1. Best AUC: %s in %s (AUC=%.3f, AUPRC=%.3f)\n",
            best_auc$Predictor, best_auc$Site, best_auc$AUC, best_auc$AUPRC))
cat("\n  2. Class imbalance impact (RC Total):\n")
for (st in valid_sites) {
  s <- roc_results[Site==st & Predictor=="Red Complex Total"]
  cat(sprintf("     %s: %.0f:1, AUC=%.3f, AUPRC=%.3f (baseline %.3f)\n",
              st, s$n_HC/s$n_PD, s$AUC, s$AUPRC, s$PR_baseline))
}
cat("\n  3. Calibration (RC Total):\n")
for (st in valid_sites) {
  s <- calib_data[Site==st & Predictor=="Red Complex Total" & Bin==1]
  br <- roc_results[Site==st & Predictor=="Red Complex Total"]
  if(nrow(s)>0 && nrow(br)>0) {
    cat(sprintf("     %s: HL p=%s, Brier=%.4f, scaled=%.3f\n",
                st, format.pval(s$HL_p,digits=2), br$Brier, br$Brier_scaled))
  }
}
cat("\n  4. Clinical utility (RC Total, Sens>=90%):\n")
for (st in valid_sites) {
  s <- clinical_metrics[Site==st & Predictor=="Red Complex Total" & Target=="Sens>=90%"]
  if(nrow(s)>0) {
    cat(sprintf("     %s: Spec=%.1f%%, PPV=%.1f%%, NPV=%.1f%%\n",
                st, s$Specificity*100, s$PPV*100, s$NPV*100))
  }
}
cat("\nDone.\n")
