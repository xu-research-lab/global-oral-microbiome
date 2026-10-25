suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(MuMIn)
})

get_fcNS <- function() {
  if (!requireNamespace("lme4", quietly = TRUE)) return(NULL)
  ok <- tryCatch(exists("formatVC", envir = asNamespace("lme4")),
                 error = function(e) FALSE)
  if (!ok) return(NULL)
  environment(get("formatVC", envir = asNamespace("lme4")))
}
fcNS <- get_fcNS()
if (!exists("%||%", envir = baseenv(), inherits = FALSE) && !is.null(fcNS)) {
  for (nm in ls(fcNS, all.names = TRUE, sorted = FALSE)) {
    obj <- get0(nm, envir = fcNS, inherits = FALSE)
    if (!is.function(obj)) next
    if (any(grepl("%\\||\\|%", deparse(body(obj))))) {
      menv <- new.env(parent = environment(obj))
      menv[["%||%"]] <- function(x, y) if (is.null(x)) y else x
      environment(obj) <- menv
      unlockBinding(nm, fcNS)
      assign(nm, obj, envir = fcNS)
      lockBinding(nm, fcNS)
    }
  }
}

find_script_dir <- function() {
  if (interactive()) return(normalizePath(getwd()))
  args <- commandArgs(trailingOnly = FALSE)
  f <- args[grep("^--file=", args)]
  if (length(f) == 1) return(dirname(normalizePath(sub("^--file=", "", f))))
  normalizePath(getwd())
}
script_dir <- find_script_dir()

PROJECT_ROOT <- normalizePath(file.path(script_dir, ".."))
DATA_DIR     <- file.path(PROJECT_ROOT, "data", "metadata-and-biom-v2-qc")

out_dir <- file.path(script_dir, "results", "Fig3", "alpha_fig3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)
cat("Output dir:", out_dir, "\n")
cat("========================================\n")
cat("Step 1: Load and prepare data\n")
cat("========================================\n")

metadata_path <- file.path(DATA_DIR, "metadata_oral_qc_dedup.tsv")
cat("Metadata file:", metadata_path, "\n")
metadata <- read.table(metadata_path,
                       header = TRUE, sep = "\t",
                       fill = TRUE, comment.char = "", stringsAsFactors = FALSE, quote = "")

site_map <- c("saliva" = "Saliva", "tongue" = "Tongue", "tooth" = "Tooth",
              "oral swab" = "Oral_swab", "oral mucosa" = "Oral_Mucosa",
              "gingiva" = "Gingiva", "other" = "Others")
metadata$site_cleaned <- ifelse(metadata$site_cleaned %in% names(site_map),
                                site_map[metadata$site_cleaned], metadata$site_cleaned)

metadata$age_cleaned_group <- gsub("＞", ">", metadata$age_cleaned_group)

cat(sprintf("Total samples: %d\n", nrow(metadata)))

needed_cols <- c("Shannon", "sex_cleaned", "age_cleaned_group", "site_cleaned",
                 "sequencing_region_cleaned", "continent_group", "bioproject_id")
meta_sub <- metadata[, needed_cols]

meta_sub <- meta_sub[complete.cases(meta_sub), ]
cat(sprintf("Complete cases: %d\n", nrow(meta_sub)))

meta_sub <- meta_sub[meta_sub$sex_cleaned %in% c("M", "F"), ]
meta_sub <- meta_sub[!meta_sub$site_cleaned %in% c("Others", "Tonsils", "Other"), ]
meta_sub <- meta_sub[!is.na(meta_sub$Shannon), ]
cat(sprintf("After filtering: %d samples\n", nrow(meta_sub)))

meta_sub$age_cleaned_group <- factor(meta_sub$age_cleaned_group)
meta_sub$age_cleaned_group <- relevel(meta_sub$age_cleaned_group, ref = "18-60 years")

meta_sub$site_cleaned <- factor(meta_sub$site_cleaned)
meta_sub$site_cleaned <- relevel(meta_sub$site_cleaned, ref = "Saliva")

meta_sub$sequencing_region_cleaned <- factor(meta_sub$sequencing_region_cleaned)
meta_sub$sequencing_region_cleaned <- relevel(meta_sub$sequencing_region_cleaned, ref = "V4")

meta_sub$continent_group <- factor(meta_sub$continent_group)
meta_sub$continent_group <- relevel(meta_sub$continent_group, ref = "Asia")

meta_sub$sex_cleaned <- factor(meta_sub$sex_cleaned)
meta_sub$sex_cleaned <- relevel(meta_sub$sex_cleaned, ref = "F")

meta_sub$bioproject_id <- factor(meta_sub$bioproject_id)

cat("\nSample distribution:\n")
cat(sprintf("  Sex: F=%d M=%d\n", sum(meta_sub$sex_cleaned=="F"), sum(meta_sub$sex_cleaned=="M")))
cat(sprintf("  BioProjects: %d\n", length(unique(meta_sub$bioproject_id))))
cat(sprintf("  Age: ")); print(table(meta_sub$age_cleaned_group))
cat(sprintf("  Site: ")); print(table(meta_sub$site_cleaned))
cat(sprintf("  Region: ")); print(table(meta_sub$sequencing_region_cleaned))
cat(sprintf("  Continent: ")); print(table(meta_sub$continent_group))

write.table(meta_sub, "mixed_model_input_data.tsv", sep = "\t", quote = FALSE)

cat("\n========================================\n")
cat("Step 2: Fit mixed effects model\n")
cat("========================================\n")

model <- lmer(Shannon ~ sex_cleaned + age_cleaned_group + site_cleaned +
              sequencing_region_cleaned + continent_group + (1 | bioproject_id),
              data = meta_sub, REML = TRUE)

saveRDS(model, "mixed_model_results.rds")

cat("\nModel Summary:\n")
print(summary(model))

anova_results <- anova(model)
print(anova_results)
write.csv(anova_results, "model_anova.csv")

cat("\n========================================\n")
cat("Step 3: Fixed effects\n")
cat("========================================\n")

fe <- as.data.frame(coef(summary(model)))
fe$Variable <- rownames(fe)
names(fe)[1:5] <- c("Estimate", "Std.Error", "df", "t_value", "p_value")

ci <- confint(model, parm = "beta_", level = 0.95)
fe$CI_lower <- ci[, 1]
fe$CI_upper <- ci[, 2]
fe$FDR <- p.adjust(fe$p_value, method = "BH")

plot_effects <- fe[fe$Variable != "(Intercept)", ]
plot_effects$Factor <- NA
plot_effects$Factor[grepl("^sex_cleaned", plot_effects$Variable)] <- "Sex"
plot_effects$Factor[grepl("^age_cleaned_group", plot_effects$Variable)] <- "Age Group"
plot_effects$Factor[grepl("^site_cleaned", plot_effects$Variable)] <- "Site"
plot_effects$Factor[grepl("^sequencing_region_cleaned", plot_effects$Variable)] <- "16S region"
plot_effects$Factor[grepl("^continent_group", plot_effects$Variable)] <- "Continent"

plot_effects$Label <- gsub("^(sex_cleaned|age_cleaned_group|site_cleaned|sequencing_region_cleaned|continent_group)",
                            "", plot_effects$Variable)

plot_effects$Sig <- ifelse(plot_effects$FDR < 0.001, "***",
                    ifelse(plot_effects$FDR < 0.01, "**",
                    ifelse(plot_effects$FDR < 0.05, "*", "NS")))

print(plot_effects[, c("Variable", "Estimate", "p_value", "FDR", "Sig")])
write.csv(plot_effects, "fixed_effects_results.csv", row.names = FALSE)

cat("\n========================================\n")
cat("Step 4: Model fit statistics\n")
cat("========================================\n")
r2 <- tryCatch(
  r.squaredGLMM(model),
  error = function(e) {
    cat(sprintf("  NOTE: r.squaredGLMM failed (%s); computing R2 manually\n", conditionMessage(e)))
    tryCatch({
      vc <- VarCorr(model)
      vf <- var(as.vector(model.matrix(model) %*% fixef(model)))
      vr <- sum(unlist(lapply(vc, attr, "stddev"))^2)
      ve <- as.numeric(attr(vc, "sc"))^2
      c(Marginal.R2 = vf/(vf + vr + ve), Conditional.R2 = (vf + vr)/(vf + vr + ve))
    }, error = function(e2) {
      cat(sprintf("  NOTE: manual R2 also failed (%s); setting R2 to NA\n", conditionMessage(e2)))
      c(Marginal.R2 = NA_real_, Conditional.R2 = NA_real_)
    })
  })
cat(sprintf("Marginal R2: %.4f\n", r2[1]))
cat(sprintf("Conditional R2: %.4f\n", r2[2]))
aic_val <- tryCatch(AIC(model), error = function(e) NA)
cat(sprintf("AIC: %.1f\n", aic_val))
cat("\n========================================\n")
cat("Step 5: Forest plot\n")
cat("========================================\n")

main_factors <- c("Sex", "Age Group", "Site", "16S region", "Continent")

get_label <- function(lab, fac) {
  if (fac == "Sex" && lab == "M") return("Male")
  if (fac == "Age Group") {
    mapping <- c("0-6 months" = "0-6 mo",
                 "6-36 months" = "6-36 mo",
                 "3-12 years" = "3-12 yr",
                 "12-18 years" = "12-18 yr",
                 ">60 years" = ">60 yr")
    return(mapping[lab])
  }
  if (fac == "Site") {
    smap <- c("Oral_swab" = "Oral Swab", "Oral_Mucosa" = "Oral Mucosa")
    if (lab %in% names(smap)) return(smap[lab])
  }
  if (fac == "Continent") {
    cmap <- c("Asia" = "AS", "Europe" = "EU", "Africa" = "AF",
              "North America" = "NA", "South America" = "SA", "Oceania" = "OC")
    if (lab %in% names(cmap)) return(cmap[lab])
  }
  return(lab)
}

plot_effects$ILab <- mapply(get_label, plot_effects$Label, plot_effects$Factor)
plot_effects <- plot_effects[!is.na(plot_effects$ILab) & plot_effects$ILab != "", ]

pf <- plot_effects[plot_effects$Factor %in% main_factors, ]
pf <- pf[!(pf$Factor == "16S region" & pf$Label == ""), ]

region_seq <- c("V1-V2", "V6-V8", "V3-V4", "V1-V3", "V4-V5", "V3", "V2", "V1")
pf <- pf[!(pf$Factor == "16S region" & !pf$ILab %in% region_seq), ]
pf$yorder <- pf$Estimate
is16 <- pf$Factor == "16S region"
pf$yorder[is16] <- length(region_seq) - match(pf$ILab[is16], region_seq)

sex_idx <- pf$Factor == "Sex"
if (sum(sex_idx) == 1) {
  base_row <- pf[sex_idx, , drop = FALSE]
  pad <- base_row[rep(1, 2), ]
  pad$ILab <- c("__pad1", "__pad2")
  pad$Estimate <- NA_real_
  pad$CI_lower <- NA_real_
  pad$CI_upper <- NA_real_
  pad$Sig <- "NS"
  pad$yorder <- base_row$yorder[1] + c(1, -1)
  pf <- rbind(pf, pad)
}
y_breaks <- unique(pf$ILab[!grepl("^__pad", pf$ILab)])

colors <- c("***" = "#E41A1C", "**" = "#FF7F00", "*" = "#4DAF4A", "NS" = "gray50")

p <- ggplot(pf, aes(x = Estimate, y = reorder(ILab, yorder))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Sig), height = 0.3, size = 1) +
  geom_point(aes(color = Sig, size = Sig)) +
  scale_color_manual(values = colors) +
  scale_size_manual(values = c("***" = 1.8, "**" = 1.6, "*" = 1.4, "NS" = 1.1)) +
  scale_y_discrete(breaks = y_breaks) +
  facet_grid(rows = vars(Factor), scales = "free_y", space = "free_y") +
  labs(x = "Effect size (coefficient)", y = NULL,
       title = "Shannon diversity",
       subtitle = paste(strwrap("Reference: Female / 18-60 yr / Saliva / V4 / AS", 36), collapse = "\n")) +
  theme_bw(base_size = 18) +
  theme(strip.text.y = element_text(face = "bold", size = 14, angle = -90, hjust = 0.5,
                                    margin = margin(1, 3, 1, 3)),
        strip.background = element_rect(fill = "gray90"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        legend.key.size = unit(0.8, "lines"),
        legend.position = "bottom")

ggsave("shannon_forest_plot.pdf", p, width = 4.5, height = 8.5)
ggsave("shannon_forest_plot.png", p, width = 4.5, height = 8.5, dpi = 600)
cat("Saved: shannon_forest_plot.pdf/png\n")

cat("\n========================================\n")
cat("Step 6: Overview boxplots\n")
cat("========================================\n")

p_age <- ggplot(meta_sub, aes(x = age_cleaned_group, y = Shannon, fill = age_cleaned_group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon diversity by age group", x = NULL, y = "Shannon index") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(colour = "black", size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"))
ggsave("shannon_age.pdf", p_age, width = 4.2, height = 3.2)
cat("Saved: shannon_age.pdf\n")

p_site <- ggplot(meta_sub, aes(x = site_cleaned, y = Shannon, fill = site_cleaned)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon diversity by oral site", x = NULL, y = "Shannon index") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(colour = "black", size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"))
ggsave("shannon_site.pdf", p_site, width = 4.2, height = 3.2)
cat("Saved: shannon_site.pdf\n")

p_geo <- ggplot(meta_sub, aes(x = continent_group, y = Shannon, fill = continent_group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon diversity by continent", x = NULL, y = "Shannon index") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(colour = "black", size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"))
ggsave("shannon_continent.pdf", p_geo, width = 4.2, height = 3.2)
cat("Saved: shannon_continent.pdf\n")

p_region <- ggplot(meta_sub, aes(x = sequencing_region_cleaned, y = Shannon,
                                  fill = sequencing_region_cleaned)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon diversity by 16S region", x = NULL, y = "Shannon index") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(colour = "black", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"))
ggsave("shannon_amplicon.pdf", p_region, width = 5, height = 3.2)
cat("Saved: shannon_amplicon.pdf\n")

cat("\n=== Alpha Diversity Analysis Complete ===\n")
