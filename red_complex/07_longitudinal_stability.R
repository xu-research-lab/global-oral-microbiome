#!/home/software/anaconda3/envs/renv/bin/Rscript
source("02_unified_statistical_framework.R")
library(data.table); library(ggplot2); library(lme4)
OUTDIR <- "../red_complex"
dt <- load_rc_data()
dtl <- dt[IsTimeSeries == "yes"][!Platform %in% EXCLUDE_PLATFORMS][!is.na(Site) & !is.na(AgeGroup) & !is.na(DiseaseStatus)]
cat(sprintf("Longitudinal: %d samples, %d subjects\n", nrow(dtl), uniqueN(dtl$SubjectID)))

subj_summary <- dtl[, .(n_tp=.N, rc_det=mean(RedComplex_detected),
  rc_mean=mean(RedComplexTotal), rc_sd=sd(RedComplexTotal),
  rc_cv=if(.N>=3&&mean(RedComplexTotal)>0) sd(RedComplexTotal)/mean(RedComplexTotal) else NA_real_),
  by=.(BioProject,SubjectID)]
tp_dist <- subj_summary[, .N, by=n_tp][order(n_tp)]
cat(sprintf("2tp:%d 3tp:%d 4-10tp:%d 11+:%d\n",
  sum(tp_dist[n_tp==2]$N,na.rm=TRUE), sum(tp_dist[n_tp==3]$N,na.rm=TRUE),
  sum(tp_dist[n_tp>=4&n_tp<=10]$N,na.rm=TRUE), sum(tp_dist[n_tp>=11]$N,na.rm=TRUE)))
fwrite(subj_summary, file.path(OUTDIR,"longitudinal_subject_summary.csv"))

cat("\n--- Mixed-Model ICC ---\n")
dtl[, log10_rc := log10(RedComplexTotal + 1e-6)]
calc_icc <- function(col) {
  dtl[, y := log10(get(col) + 1e-6)]
  m <- tryCatch({lmer(y~1+(1|SubjectID),data=dtl,control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))},error=function(e)NULL)
  if(is.null(m)) return(data.table(measure=col,ICC=NA_real_,n_subj=NA_integer_,n_obs=NA_integer_))
  vc <- as.data.table(VarCorr(m))
  data.table(measure=col,ICC=vc[grp=="SubjectID"]$vcov/(vc[grp=="SubjectID"]$vcov+vc[grp=="Residual"]$vcov),n_subj=uniqueN(dtl$SubjectID),n_obs=nrow(dtl))
}
icc_mm <- rbindlist(lapply(c("RedComplexTotal","P_gingivalis","T_denticola","T_forsythia","Alpha_Shannon"),calc_icc))
print(icc_mm[,.(measure,ICC=round(ICC,3),n_subj,n_obs)])
fwrite(icc_mm, file.path(OUTDIR,"longitudinal_icc_mixed_model.csv"))

cat("\n--- Detection Transitions ---\n")
dtl[, tp_order := seq_len(.N), by=.(BioProject,SubjectID)]
pairs_all <- dtl[tp_order<=2,.(rc1=RedComplex_detected[1],rc2=RedComplex_detected[2],rc_v1=RedComplexTotal[1],rc_v2=RedComplexTotal[2]),by=.(BioProject,SubjectID)]
pairs_all <- pairs_all[!is.na(rc1)&!is.na(rc2)]
tt <- table(pairs_all$rc1,pairs_all$rc2); print(tt)
sp <- cor(pairs_all$rc_v1,pairs_all$rc_v2,method="spearman")
cat(sprintf("Spearman r (all %d pairs): %.3f\n",nrow(pairs_all),sp))
cat(sprintf("Persist detected: %.1f%%  Persist undetected: %.1f%%\n",
  100*tt[2,2]/sum(tt[2,]),100*tt[1,1]/sum(tt[1,])))

cat("\n--- Variance Decomposition ---\n")
m_null <- lmer(log10_rc~1+(1|SubjectID),data=dtl,control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vc_null <- as.data.table(VarCorr(m_null))[,pct:=100*vcov/sum(vcov)]; print(vc_null[,.(grp,vcov,pct)])
m_cov <- lmer(log10_rc~Site+DiseaseStatus+AgeGroup+(1|SubjectID),data=dtl,control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vc_cov <- as.data.table(VarCorr(m_cov))[,pct:=100*vcov/sum(vcov)]; print(vc_cov[,.(grp,vcov,pct)])
fwrite(vc_cov,file.path(OUTDIR,"longitudinal_variance_components.csv"))

cat("\n--- CV (mean-adjusted) ---\n")
cv_data <- subj_summary[n_tp>=3&!is.na(rc_cv)&is.finite(rc_cv)]
cv_data <- merge(cv_data,unique(dtl[,.(SubjectID,DiseaseStatus)]),by="SubjectID",all.x=TRUE)
cv_data[,`:=`(log_mean=log10(rc_mean+1e-8),disease_group=ifelse(DiseaseStatus=="HC","HC","Disease"))]
cv_data[,cv_adj:=residuals(lm(rc_cv~log_mean,data=cv_data))]
hc_cv<-cv_data[disease_group=="HC"]$cv_adj; dis_cv<-cv_data[disease_group=="Disease"]$cv_adj
if(length(hc_cv)>10&&length(dis_cv)>10){wt<-wilcox.test(hc_cv,dis_cv);cat(sprintf("Adj CV: HC median=%.3f Disease median=%.3f p=%s\n",median(hc_cv,na.rm=TRUE),median(dis_cv,na.rm=TRUE),format.pval(wt$pvalue,digits=3)))}

cat("\n--- Plots ---\n")
# 1. ICC bar
p1<-ggplot(icc_mm,aes(x=reorder(measure,ICC),y=ICC))+geom_bar(stat="identity",fill="#3498db",width=0.6)+geom_hline(yintercept=c(0.3,0.5,0.7),linetype=c("dashed","dotted","dashed"),color=c("grey50","grey30","grey50"))+labs(title="Mixed-Model ICC (All Longitudinal Subjects)",subtitle=sprintf("n=%d subjects, %d samples",uniqueN(dtl$SubjectID),nrow(dtl)),x="",y="ICC")+ylim(0,1)+theme_rc()
ggsave(file.path(OUTDIR,"longitudinal_icc_bar.png"),p1,width=8,height=5,dpi=150)

# 2. Spearman stability bar (all measures on consecutive pairs)
# Compute Spearman r for each species individually
sp_pg <- cor(pairs_all$rc_v1, pairs_all$rc_v2, method="spearman")
# Compute per-species Spearman
dtl_pg <- dtl[,.(v1=P_gingivalis[tp_order==1],v2=P_gingivalis[tp_order==2]),by=.(BioProject,SubjectID)]
dtl_pg <- dtl_pg[!is.na(v1)&!is.na(v2)]
dtl_td <- dtl[,.(v1=T_denticola[tp_order==1],v2=T_denticola[tp_order==2]),by=.(BioProject,SubjectID)]
dtl_td <- dtl_td[!is.na(v1)&!is.na(v2)]
dtl_tf <- dtl[,.(v1=T_forsythia[tp_order==1],v2=T_forsythia[tp_order==2]),by=.(BioProject,SubjectID)]
dtl_tf <- dtl_tf[!is.na(v1)&!is.na(v2)]
dtl_sh <- dtl[,.(v1=Alpha_Shannon[tp_order==1],v2=Alpha_Shannon[tp_order==2]),by=.(BioProject,SubjectID)]
dtl_sh <- dtl_sh[!is.na(v1)&!is.na(v2)]

stability_dt <- data.table(
  measure = c("P. gingivalis","Red Complex Total","T. forsythia","Shannon Diversity","T. denticola"),
  spearman_r = c(
    cor(dtl_pg$v1, dtl_pg$v2, method="spearman"),
    sp,
    cor(dtl_tf$v1, dtl_tf$v2, method="spearman"),
    cor(dtl_sh$v1, dtl_sh$v2, method="spearman"),
    cor(dtl_td$v1, dtl_td$v2, method="spearman")
  ),
  n_pairs = c(nrow(dtl_pg), nrow(pairs_all), nrow(dtl_tf), nrow(dtl_sh), nrow(dtl_td))
)
stability_dt <- stability_dt[order(-spearman_r)]

p2 <- ggplot(stability_dt, aes(x=reorder(measure, spearman_r), y=spearman_r)) +
  geom_bar(stat="identity", fill="#e74c3c", width=0.6) +
  geom_text(aes(label=sprintf("%.3f", spearman_r)), hjust=-0.1, size=3.5) +
  labs(title="Test-Retest Stability (Consecutive Timepoints)",
       subtitle=sprintf("Spearman r, n=%d consecutive pairs from %d subjects", nrow(pairs_all), uniqueN(pairs_all$SubjectID)),
       x="", y="Spearman r") + ylim(0, 1) + coord_flip() + theme_rc()
ggsave(file.path(OUTDIR,"longitudinal_stability_bar.png"),p2,width=8,height=5,dpi=150)

# 3. Detection transition bar chart
trans_dt <- data.table(
  transition = c("Stay undetected\n(0→0)", "Persist detected\n(1→1)", "Acquisition\n(0→1)", "Loss\n(1→0)"),
  count = as.integer(c(tt[1,1], tt[2,2], tt[1,2], tt[2,1])),
  pct = c(100*tt[1,1]/sum(tt), 100*tt[2,2]/sum(tt), 100*tt[1,2]/sum(tt), 100*tt[2,1]/sum(tt))
)
trans_dt[, transition := factor(transition, levels=transition[order(-pct)])]

p3 <- ggplot(trans_dt, aes(x=transition, y=count, fill=transition)) +
  geom_bar(stat="identity", width=0.6) +
  geom_text(aes(label=sprintf("%d\n(%.1f%%)", count, pct)), vjust=-0.2, size=3.5) +
  scale_fill_manual(values=c("Stay undetected\n(0→0)"="#ecf0f1","Persist detected\n(1→1)"="#c0392b",
                              "Acquisition\n(0→1)"="#f39c12","Loss\n(1→0)"="#3498db")) +
  labs(title="Red Complex Detection Transitions",
       subtitle=sprintf("Consecutive timepoint pairs, n=%d pairs from %d subjects", nrow(pairs_all), uniqueN(pairs_all$SubjectID)),
       x="", y="Number of pairs") + theme_rc() + theme(legend.position="none")
ggsave(file.path(OUTDIR,"longitudinal_transitions.png"),p3,width=7,height=5,dpi=150)

# Save transition data for report
fwrite(trans_dt, file.path(OUTDIR,"longitudinal_transitions.csv"))
cat("Done.\n")
