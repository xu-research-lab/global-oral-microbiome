#!/usr/bin/env Rscript
# ============================================================================
# Alpha Diversity v2 — Mixed Effects Model (updated metadata)
# Formula: Shannon ~ sex_cleaned + age_cleaned_group + site_cleaned +
#          sequencing_region_cleaned + continent_group + (1 | bioproject_id)
#
# Uses NATIVE column names from metadata_oral_qc_dedup.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(MuMIn)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/alpha_diversity_v2")

# ============================================================================
# 1. Load data
# ============================================================================
cat("========================================\n")
cat("Step 1: Load and prepare data\n")
cat("========================================\n")

metadata <- read.table("../metadata_oral_qc_dedup.tsv", header = TRUE, sep = "\t",
                       fill = TRUE, comment.char = "", stringsAsFactors = FALSE, quote = "")

# Apply site mapping
site_map <- c("saliva" = "Saliva", "tongue" = "Tongue", "tooth" = "Tooth",
              "oral swab" = "Oral_swab", "oral mucosa" = "Oral_Mucosa",
              "gingiva" = "Gingiva", "other" = "Others")
metadata$site_cleaned <- ifelse(metadata$site_cleaned %in% names(site_map),
                                site_map[metadata$site_cleaned], metadata$site_cleaned)
# Fix age group encoding
metadata$age_cleaned_group <- gsub("＞", ">", metadata$age_cleaned_group)

cat(sprintf("Total samples: %d\n", nrow(metadata)))

# Select relevant columns (native names)
needed_cols <- c("Shannon", "sex_cleaned", "age_cleaned_group", "site_cleaned",
                 "sequencing_region_cleaned", "continent_group", "bioproject_id")
meta_sub <- metadata[, needed_cols]

# Remove missing
meta_sub <- meta_sub[complete.cases(meta_sub), ]
cat(sprintf("Complete cases: %d\n", nrow(meta_sub)))

# Filter: valid Sex (M/F), remove "Others" site
meta_sub <- meta_sub[meta_sub$sex_cleaned %in% c("M", "F"), ]
meta_sub <- meta_sub[!meta_sub$site_cleaned %in% c("Others", "Tonsils"), ]
meta_sub <- meta_sub[!is.na(meta_sub$Shannon), ]
cat(sprintf("After filtering: %d samples\n", nrow(meta_sub)))

# Set reference levels
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

# ============================================================================
# 2. Fit mixed effects model
# ============================================================================
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

# ============================================================================
# 3. Fixed effects
# ============================================================================
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

# Classify factor type
plot_effects$Factor <- NA
plot_effects$Factor[grepl("^sex_cleaned", plot_effects$Variable)] <- "Sex"
plot_effects$Factor[grepl("^age_cleaned_group", plot_effects$Variable)] <- "Age Group"
plot_effects$Factor[grepl("^site_cleaned", plot_effects$Variable)] <- "Site"
plot_effects$Factor[grepl("^sequencing_region_cleaned", plot_effects$Variable)] <- "16S Region"
plot_effects$Factor[grepl("^continent_group", plot_effects$Variable)] <- "Continent"

plot_effects$Label <- gsub("^(sex_cleaned|age_cleaned_group|site_cleaned|sequencing_region_cleaned|continent_group)",
                            "", plot_effects$Variable)

plot_effects$Sig <- ifelse(plot_effects$FDR < 0.001, "***",
                    ifelse(plot_effects$FDR < 0.01, "**",
                    ifelse(plot_effects$FDR < 0.05, "*", "NS")))

print(plot_effects[, c("Variable", "Estimate", "p_value", "FDR", "Sig")])
write.csv(plot_effects, "fixed_effects_results.csv", row.names = FALSE)

# ============================================================================
# 4. Model fit
# ============================================================================
cat("\n========================================\n")
cat("Step 4: Model fit statistics\n")
cat("========================================\n")
r2 <- r.squaredGLMM(model)
cat(sprintf("Marginal R2: %.4f\n", r2[1]))
cat(sprintf("Conditional R2: %.4f\n", r2[2]))
cat(sprintf("AIC: %.1f\n", AIC(model)))

# ============================================================================
# 5. Forest plot
# ============================================================================
cat("\n========================================\n")
cat("Step 5: Forest plot\n")
cat("========================================\n")

main_factors <- c("Sex", "Age Group", "Site", "16S Region", "Continent")

get_label <- function(lab, fac) {
  if (fac == "Sex" && lab == "M") return("Male")
  if (fac == "Age Group") {
    mapping <- c("0-6 months" = "Infancy (0-6 mo)",
                 "6-36 months" = "Toddler (6-36 mo)",
                 "3-12 years" = "Childhood (3-12 yr)",
                 "12-18 years" = "Adolescence (12-18 yr)",
                 ">60 years" = "Elderly (>60 yr)")
    return(mapping[lab])
  }
  if (fac == "Site" && lab == "Oral_swab") return("Oral Swab")
  return(lab)
}

plot_effects$ILab <- mapply(get_label, plot_effects$Label, plot_effects$Factor)

pf <- plot_effects[plot_effects$Factor %in% main_factors, ]
pf <- pf[!(pf$Factor == "16S Region" & pf$Label == ""), ]

colors <- c("***" = "#E41A1C", "**" = "#FF7F00", "*" = "#4DAF4A", "NS" = "gray50")

p <- ggplot(pf, aes(x = Estimate, y = reorder(ILab, Estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Sig), height = 0.2) +
  geom_point(aes(color = Sig, size = Sig)) +
  scale_color_manual(values = colors) +
  scale_size_manual(values = c("***" = 3, "**" = 2.5, "*" = 2.2, "NS" = 1.8)) +
  facet_grid(rows = vars(Factor), scales = "free_y", space = "free_y") +
  labs(x = "Effect Size (Coefficient)", y = NULL,
       title = "Shannon Diversity — Mixed Effects Model",
       subtitle = "Ref: Female / 18-60 yr / Saliva / V4 / Asia") +
  theme_bw(base_size = 11) +
  theme(strip.text.y = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "gray90"),
        legend.position = "bottom")

ggsave("shannon_forest_plot.pdf", p, width = 6, height = 12)
ggsave("shannon_forest_plot.png", p, width = 6, height = 12, dpi = 300)
cat("Saved: shannon_forest_plot.pdf/png\n")

# ============================================================================
# 6. Overview boxplots
# ============================================================================
cat("\n========================================\n")
cat("Step 6: Overview boxplots\n")
cat("========================================\n")

p_age <- ggplot(meta_sub, aes(x = age_cleaned_group, y = Shannon, fill = age_cleaned_group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon Diversity by Age Group") +
  theme_bw() + theme(legend.position = "none",
                      axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("shannon_age.pdf", p_age, width = 8, height = 5)
cat("Saved: shannon_age.pdf\n")

p_site <- ggplot(meta_sub, aes(x = site_cleaned, y = Shannon, fill = site_cleaned)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon Diversity by Oral Site") +
  theme_bw() + theme(legend.position = "none",
                      axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("shannon_site.pdf", p_site, width = 8, height = 5)
cat("Saved: shannon_site.pdf\n")

p_geo <- ggplot(meta_sub, aes(x = continent_group, y = Shannon, fill = continent_group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon Diversity by Continent") +
  theme_bw() + theme(legend.position = "none",
                      axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("shannon_continent.pdf", p_geo, width = 8, height = 5)
cat("Saved: shannon_continent.pdf\n")

p_region <- ggplot(meta_sub, aes(x = sequencing_region_cleaned, y = Shannon,
                                  fill = sequencing_region_cleaned)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon Diversity by 16S Region") +
  theme_bw() + theme(legend.position = "none",
                      axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("shannon_amplicon.pdf", p_region, width = 10, height = 5)
cat("Saved: shannon_amplicon.pdf\n")

cat("\n=== Alpha Diversity Analysis Complete ===\n")
