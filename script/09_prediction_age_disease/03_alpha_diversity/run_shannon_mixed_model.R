#!/usr/bin/env Rscript
# ============================================================================
# Alpha Diversity — Mixed Effects Model (Adapted for new data)
# Formula: Shannon ~ Sex + Age_Group + Site_Cleaned + Region_custom + Geography + (1 | BioProject)
#
# Key changes from old version:
#   - Shannon already in metadata (no separate alpha-diversity.tsv)
#   - Geography already in metadata (= continent_group)
#   - Data paths point to new_analysis/data_prep/
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

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/alpha_diversity")
source("../maaslin3/input_helpers.R")

cat("========================================\n")
cat("Step 1: Load new data\n")
cat("========================================\n")

# Alpha diversity uses the non-deduplicated oral QC metadata.
metadata <- load_standard_metadata(dedup = FALSE)

cat(sprintf("Total samples: %d\n", nrow(metadata)))

# Select relevant columns
meta_sub <- metadata[, c("Shannon", "Sex", "Age_Group", "Site_Cleaned",
                          "Region_custom", "Geography", "BioProject")]

# Remove samples with missing values
meta_sub <- meta_sub[complete.cases(meta_sub), ]
cat(sprintf("Complete cases: %d\n", nrow(meta_sub)))

# Filter: valid Sex (M/F), remove "Others" site
meta_sub <- meta_sub[meta_sub$Sex %in% c("M", "F"), ]
meta_sub <- meta_sub[!meta_sub$Site_Cleaned %in% c("Others", "Tonsils"), ]
meta_sub <- meta_sub[!is.na(meta_sub$Shannon), ]
cat(sprintf("After filtering: %d samples\n", nrow(meta_sub)))

# Set reference levels
meta_sub$Age_Group <- factor(meta_sub$Age_Group)
meta_sub$Age_Group <- relevel(meta_sub$Age_Group, ref = "18-60 years")

meta_sub$Site_Cleaned <- factor(meta_sub$Site_Cleaned)
meta_sub$Site_Cleaned <- relevel(meta_sub$Site_Cleaned, ref = "Saliva")

meta_sub$Region_custom <- factor(meta_sub$Region_custom)
meta_sub$Region_custom <- relevel(meta_sub$Region_custom, ref = "V4")

meta_sub$Geography <- factor(meta_sub$Geography)
meta_sub$Geography <- relevel(meta_sub$Geography, ref = "Asia")

meta_sub$Sex <- factor(meta_sub$Sex)
meta_sub$Sex <- relevel(meta_sub$Sex, ref = "F")

meta_sub$BioProject <- factor(meta_sub$BioProject)

cat("\nSample distribution:\n")
cat(sprintf("  Sex: F=%d M=%d\n", sum(meta_sub$Sex=="F"), sum(meta_sub$Sex=="M")))
cat(sprintf("  BioProjects: %d\n", length(unique(meta_sub$BioProject))))
cat(sprintf("  Age: ")); print(table(meta_sub$Age_Group))
cat(sprintf("  Site: ")); print(table(meta_sub$Site_Cleaned))
cat(sprintf("  Region: ")); print(table(meta_sub$Region_custom))
cat(sprintf("  Geography: ")); print(table(meta_sub$Geography))

# Save input data
write.table(meta_sub, "mixed_model_input_data.tsv", sep = "\t", quote = FALSE)

# ============================================================================
# Step 2: Fit mixed effects model
# ============================================================================
cat("\n========================================\n")
cat("Step 2: Fit model\n")
cat("========================================\n")

model <- lmer(Shannon ~ Sex + Age_Group + Site_Cleaned + Region_custom + Geography + (1 | BioProject),
              data = meta_sub, REML = TRUE)

saveRDS(model, "mixed_model_results.rds")

cat("\nModel Summary:\n")
print(summary(model))

# ANOVA
anova_results <- anova(model)
print(anova_results)
write.csv(anova_results, "model_anova.csv")

# ============================================================================
# Step 3: Extract and save fixed effects
# ============================================================================
cat("\n========================================\n")
cat("Step 3: Extract fixed effects\n")
cat("========================================\n")

fe <- as.data.frame(coef(summary(model)))
fe$Variable <- rownames(fe)
names(fe)[1:5] <- c("Estimate", "Std.Error", "df", "t_value", "p_value")

# 95% CI
ci <- confint(model, parm = "beta_", level = 0.95)
fe$CI_lower <- ci[, 1]
fe$CI_upper <- ci[, 2]
fe$FDR <- p.adjust(fe$p_value, method = "BH")

# Remove intercept for plotting
plot_effects <- fe[fe$Variable != "(Intercept)", ]

# Classify factor type
plot_effects$Factor <- NA
plot_effects$Factor[grepl("^Sex", plot_effects$Variable)] <- "Sex"
plot_effects$Factor[grepl("^Age_Group", plot_effects$Variable)] <- "Age_Group"
plot_effects$Factor[grepl("^Site_Cleaned", plot_effects$Variable)] <- "Site_Cleaned"
plot_effects$Factor[grepl("^Region_custom", plot_effects$Variable)] <- "Region_custom"
plot_effects$Factor[grepl("^Geography", plot_effects$Variable)] <- "Geography"

# Clean labels
plot_effects$Label <- gsub("^(Sex|Age_Group|Site_Cleaned|Region_custom|Geography)", "",
                            plot_effects$Variable)

# Significance
plot_effects$Sig <- ifelse(plot_effects$FDR < 0.001, "***",
                    ifelse(plot_effects$FDR < 0.01, "**",
                    ifelse(plot_effects$FDR < 0.05, "*", "NS")))

print(plot_effects[, c("Variable", "Estimate", "p_value", "FDR", "Sig")])
write.csv(plot_effects, "fixed_effects_results.csv", row.names = FALSE)

# ============================================================================
# Step 4: Model fit statistics
# ============================================================================
cat("\n========================================\n")
cat("Step 4: Model fit\n")
cat("========================================\n")
r2 <- r.squaredGLMM(model)
cat(sprintf("Marginal R²: %.4f\n", r2[1]))
cat(sprintf("Conditional R²: %.4f\n", r2[2]))
cat(sprintf("AIC: %.1f\n", AIC(model)))

# ============================================================================
# Step 5: Forest plot
# ============================================================================
cat("\n========================================\n")
cat("Step 5: Generate forest plot\n")
cat("========================================\n")

main_factors <- c("Sex", "Age_Group", "Site_Cleaned", "Region_custom", "Geography")

# Improved labels
get_label <- function(lab, fac) {
  if (fac == "Sex" && lab == "M") return("Male")
  if (fac == "Age_Group") {
    mapping <- c("0-6 months" = "Infancy (0-6 mo)",
                 "6-36 months" = "Toddler (6-36 mo)",
                 "3-12 years" = "Childhood (3-12 yr)",
                 "12-18 years" = "Adolescence (12-18 yr)",
                 ">60 years" = "Elderly (>60 yr)")
    return(mapping[lab])
  }
  if (fac == "Site_Cleaned" && lab == "Oral_swab") return("Oral Swab")
  return(lab)
}

plot_effects$ILab <- mapply(get_label, plot_effects$Label, plot_effects$Factor)

pf <- plot_effects[plot_effects$Factor %in% main_factors, ]
pf <- pf[!(pf$Factor == "Region_custom" & pf$Label == ""), ]

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
  theme(
    strip.text.y = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "bottom"
  )

ggsave("shannon_forest_plot.pdf", p, width = 6, height = 12)
ggsave("shannon_forest_plot.png", p, width = 6, height = 12, dpi = 300)
cat("Saved: shannon_forest_plot.pdf/png\n")

# ============================================================================
# Step 6: Generate overview plots
# ============================================================================
cat("\n========================================\n")
cat("Step 6: Overview plots by factor\n")
cat("========================================\n")

# Shannon by Age
p_age <- ggplot(meta_sub, aes(x = Age_Group, y = Shannon, fill = Age_Group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon Diversity by Age Group") +
  theme_bw() + theme(legend.position = "none",
                      axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("shannon_age.pdf", p_age, width = 8, height = 5)
cat("Saved: shannon_age.pdf\n")

# Shannon by Site
p_site <- ggplot(meta_sub, aes(x = Site_Cleaned, y = Shannon, fill = Site_Cleaned)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon Diversity by Oral Site") +
  theme_bw() + theme(legend.position = "none",
                      axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("shannon_site.pdf", p_site, width = 8, height = 5)
cat("Saved: shannon_site.pdf\n")

# Shannon by Geography
p_geo <- ggplot(meta_sub, aes(x = Geography, y = Shannon, fill = Geography)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon Diversity by Continent") +
  theme_bw() + theme(legend.position = "none",
                      axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("shannon_geography.pdf", p_geo, width = 8, height = 5)
cat("Saved: shannon_geography.pdf\n")

# Shannon by Region (amplicon)
p_region <- ggplot(meta_sub, aes(x = Region_custom, y = Shannon, fill = Region_custom)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  labs(title = "Shannon Diversity by 16S Region") +
  theme_bw() + theme(legend.position = "none",
                      axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("shannon_amplicon.pdf", p_region, width = 10, height = 5)
cat("Saved: shannon_amplicon.pdf\n")

cat("\n=== Alpha Diversity Analysis Complete ===\n")
