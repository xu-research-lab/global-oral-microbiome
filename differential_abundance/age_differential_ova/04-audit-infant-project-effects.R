# Descriptive audit of the displayed 6-36-month HC saliva genera.
# Run from the project root. Does not modify metadata or refit the DA models.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})
source("differential_abundance/age_differential_ova/da_helpers.R")

base <- "differential_abundance/result/age_hc_ova/age_da_consensus"
fig_dir <- file.path(base, "figures_saliva_target_top10_globalpeak_maaslin2coef")
out_dir <- file.path(base, "diagnostics_infant_age_project")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(out_dir, name), sep = "\t")

meta <- load_age_metadata(hc_only = TRUE)
count_ids <- strsplit(readLines(COUNTS_TSV, n = 1), "\t", fixed = TRUE)[[1]][-1]
meta <- meta[meta$sample_id %in% count_ids & meta$age_cleaned_group %in% AGE_ORDER, ]
meta$age_cleaned_group <- factor(meta$age_cleaned_group, levels = AGE_ORDER)
stopifnot(nrow(meta) == 17360L, !anyDuplicated(meta$sample_id))
write_table(meta %>% group_by(age_cleaned_group) %>%
              summarise(n_samples = n(), n_projects = n_distinct(bioproject_id),
                        .groups = "drop"), "sample_counts_by_age.tsv")
project_counts <- meta %>% count(bioproject_id, age_cleaned_group, name = "n_samples")
write_table(project_counts, "sample_counts_by_project_age.tsv")
infant <- meta %>% filter(age_cleaned_group %in% AGE_ORDER[1:2])
write_table(as.data.frame.matrix(xtabs(~ bioproject_id + age_cleaned_group,
                                      droplevels(infant))) %>%
              tibble::rownames_to_column("bioproject_id"), "infant_project_counts.tsv")
write_table(infant %>% count(bioproject_id, age_cleaned_group, age_cleaned,
                             name = "n_samples"), "infant_recorded_age_counts.tsv")

selected <- fread(file.path(fig_dir, "figureB_selected_2method.tsv"), data.table = FALSE)
selected <- selected[selected$target == "6-36 months", ]
focus <- selected$feature
stopifnot(length(focus) == 10L)
patterns <- tempfile()
rows_file <- tempfile()
writeLines(paste0(focus, "\t"), patterns)
status <- system2("rg", c("-F", "-f", shQuote(patterns), shQuote(RELABUND_TSV)),
                  stdout = rows_file)
stopifnot(status == 0L)
rel <- fread(rows_file, header = FALSE, data.table = FALSE)
rel_header <- strsplit(readLines(RELABUND_TSV, n = 1), "\t", fixed = TRUE)[[1]]
names(rel) <- c("feature", rel_header[-1])
stopifnot(setequal(rel$feature, focus), !anyDuplicated(rel$feature),
          all(meta$sample_id %in% names(rel)))
x <- t(as.matrix(rel[, meta$sample_id, drop = FALSE]))
colnames(x) <- rel$feature
x <- x[, focus, drop = FALSE]
stopifnot(identical(rownames(x), meta$sample_id))
unlink(c(patterns, rows_file))
genus <- function(f) sub(".*g__", "", f)

age_stats <- function(keep, scenario) {
  younger <- keep & meta$age_cleaned_group == "0-6 months"
  older <- keep & meta$age_cleaned_group == "6-36 months"
  stopifnot(sum(younger) > 0, sum(older) > 0)
  a <- colMeans(x[younger, , drop = FALSE])
  b <- colMeans(x[older, , drop = FALSE])
  data.frame(scenario = scenario, feature = focus, genus = genus(focus),
             n_0_6 = sum(younger), n_6_36 = sum(older),
             mean_0_6_pct = 100 * a, mean_6_36_pct = 100 * b,
             prevalence_0_6_pct = 100 * colMeans(x[younger, , drop = FALSE] > 0),
             prevalence_6_36_pct = 100 * colMeans(x[older, , drop = FALSE] > 0),
             fold_mean_6_36_vs_0_6 = b / a,
             higher_mean_in_6_36 = b > a, row.names = NULL)
}
shared <- names(which(table(infant$bioproject_id, infant$age_cleaned_group)[, 1] > 0 &
                        table(infant$bioproject_id, infant$age_cleaned_group)[, 2] > 0))
within <- bind_rows(lapply(shared, function(p) age_stats(meta$bioproject_id == p, p)))
write_table(within, "within_project_age_comparison.tsv")
scenarios <- list(
  all_samples = rep(TRUE, nrow(meta)),
  excluding_PRJEB35824 = meta$bioproject_id != "PRJEB35824",
  excluding_PRJNA997934 = meta$bioproject_id != "PRJNA997934",
  excluding_both_major_projects = !meta$bioproject_id %in% c("PRJEB35824", "PRJNA997934"),
  excluding_PRJNA997934_age_zero = !(meta$bioproject_id == "PRJNA997934" & meta$age_cleaned == 0)
)
sensitivity <- bind_rows(lapply(names(scenarios), function(nm) age_stats(scenarios[[nm]], nm)))
write_table(sensitivity, "descriptive_exclusion_sensitivity.tsv")

# Verify that the audit recovers the raw values underlying the displayed heatmap.
heat <- fread(file.path(fig_dir, "figureB_heatmap_data_2method.tsv"), data.table = FALSE)
for (a in AGE_ORDER) {
  hh <- heat[heat$feature %in% focus & as.character(heat$age_lab) == a, ]
  if (nrow(hh) == 0L) hh <- heat[heat$feature %in% focus & grepl(a, heat$group_lab, fixed = TRUE), ]
  stopifnot(nrow(hh) == length(focus))
  got <- colMeans(x[meta$age_cleaned_group == a, , drop = FALSE])
  stopifnot(max(abs(got[hh$feature] - hh$mean_abundance)) < 1e-10)
}

full <- fread(file.path(base, "full_table.tsv"), data.table = FALSE)
coefs <- full[full$feature %in% focus & full$comparison %in%
                c("saliva_enrich_0-6_months", "saliva_enrich_6-36_months"), ]
write_table(coefs, "original_da_results_selected_genera.tsv")

# Original deposited age categories expose inadequate precision for half-year bins.
sra <- fread("data/projects/PRJNA997934/SraRunTable.csv", data.table = FALSE)
sra <- sra[sra$Run %in% meta$sample_id, ]
write_table(sra %>% count(age_cat, AGE, name = "n_samples"), "PRJNA997934_deposited_age_counts.tsv")

# A descriptive comparison, without inferential p-values: subject identities are
# unavailable for the main longitudinal study, so samples are not independent people.
plot_data <- within[within$scenario %in% c("PRJEB35824", "PRJNA997934"), ]
plot_data <- rbind(plot_data, sensitivity[sensitivity$scenario == "excluding_both_major_projects", ])
plot_data$panel <- factor(plot_data$scenario,
  levels = c("PRJEB35824", "PRJNA997934", "excluding_both_major_projects"),
  labels = c("PRJEB35824\n360 vs 272 samples",
             "PRJNA997934*\n313 vs 194 samples",
             "Excluding both projects\n113 vs 77 samples"))
plot_data$genus <- factor(plot_data$genus, levels = rev(genus(focus)))
p <- ggplot(plot_data, aes(y = genus)) +
  geom_segment(aes(x = mean_0_6_pct, xend = mean_6_36_pct, yend = genus),
               colour = "grey70", linewidth = 0.4) +
  geom_point(aes(x = mean_0_6_pct, colour = "0-6 months"), size = 2) +
  geom_point(aes(x = mean_6_36_pct, colour = "6-36 months"), size = 2) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.001, base = 10),
                     breaks = c(0, 0.01, 0.1, 1, 10)) +
  scale_colour_manual(values = c("0-6 months" = "#440154", "6-36 months" = "#21918C")) +
  facet_wrap(~panel, nrow = 1) +
  labs(x = "Mean relative abundance (%)", y = NULL, colour = "Recorded age group",
       caption = paste("*Deposited ages are integer years; the 0Y category cannot establish 0-6-month membership.",
                       "Descriptive means include zeros. PRJEB35824 includes repeated visits with unavailable subject IDs.", sep = "\n")) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic", colour = "black"),
        panel.grid.minor = element_blank(), legend.position = "top",
        strip.background = element_rect(fill = "grey95"),
        plot.caption = element_text(size = 7, hjust = 0),
        plot.margin = margin(5, 8, 5, 5))
ggsave(file.path(out_dir, "project_mean_abundance_audit.pdf"), p,
       device = cairo_pdf, width = 200, height = 115, units = "mm")
ggsave(file.path(out_dir, "project_mean_abundance_audit.png"), p,
       device = ragg::agg_png, width = 200, height = 115, units = "mm", dpi = 300, bg = "white")

cat("Overall:", nrow(meta), "samples;", n_distinct(meta$bioproject_id), "BioProjects\n")
cat("Infant groups:", nrow(infant), "samples;", n_distinct(infant$bioproject_id),
    "BioProjects;", length(shared), "projects span both recorded groups\n")
print(within %>% group_by(scenario) %>% summarise(n_positive = sum(higher_mean_in_6_36), .groups = "drop"))
print(sensitivity %>% group_by(scenario) %>% summarise(n_positive = sum(higher_mean_in_6_36), .groups = "drop"))
cat("Audit values agree with saved heatmap data; PNG exported at 300 dpi.\n")
