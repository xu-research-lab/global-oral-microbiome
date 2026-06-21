# ============================================================
# plot-phylum-relative-site-stacked-260610.R
# Stacked bar of phylum-level relative abundance
#   X: Age Group   Y: Relative Abundance
#   Facet: Saliva / Tooth / Tongue / Oral Mucosa
#   Fill : Top 9 phyla + Others, dominant phylum at bar bottom
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# Paths
abundance_file <- "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/phylum-relative/results/phylum-relative-table/phylum-relative-table.tsv"
metadata_file  <- "/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/metadata_oral_qc.tsv"
out_dir        <- "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/phylum-relative/results/phylum-site-stacked"

# Fixed orderings
# site_cleaned in metadata is lowercase; match lowercase, display title case
target_sites_lower   <- c("saliva", "tooth", "tongue", "oral mucosa")
target_sites_display <- c("Saliva", "Tooth", "Tongue", "Oral Mucosa")
site_label_map       <- setNames(target_sites_display, target_sites_lower)

age_order    <- c("0-6 months", "6-36 months", "3-12 years",
                  "12-18 years", "18-60 years", ">60 years")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("==========================================\n")
cat("Start:        ", format(Sys.time()), "\n")
cat("Abundance:    ", abundance_file, "\n")
cat("Metadata:     ", metadata_file, "\n")
cat("Output dir:   ", out_dir, "\n")
cat("Facet order:  ", paste(target_sites_display, collapse = " | "), "\n")
cat("Age order:    ", paste(age_order, collapse = " | "), "\n")
cat("==========================================\n\n")

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------
cat("[1/6] Reading relative abundance table (skip=1 for biom comment line)...\n")
rel <- fread(abundance_file, header = TRUE, sep = "\t", skip = 1)
setnames(rel, 1, "Phylum")
cat("  dims: ", nrow(rel), " phyla x ", ncol(rel) - 1, " samples\n", sep = "")

cat("\n[2/6] Reading metadata...\n")
meta <- fread(metadata_file, header = TRUE, sep = "\t")
cat("  dims: ", nrow(meta), " samples x ", ncol(meta), " columns\n", sep = "")

required_cols <- c("sample_id", "site_cleaned", "age_cleaned_group")
missing_cols  <- setdiff(required_cols, names(meta))
if (length(missing_cols) > 0) {
  stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ------------------------------------------------------------
# 2. Long format + strip taxonomy prefix + merge with metadata
# ------------------------------------------------------------
cat("\n[3/6] Melting + merging metadata...\n")

phylum_cols <- setdiff(names(rel), "Phylum")
rel_long <- data.table::melt(rel,
                             id.vars       = "Phylum",
                             measure.vars  = phylum_cols,
                             variable.name = "sample_id",
                             value.name    = "Abundance")
rel_long[, sample_id := as.character(sample_id)]

rel_long[, Phylum_clean := gsub("^d__Bacteria;\\s*p__", "", Phylum)]
rel_long[, Phylum_clean := gsub("^p__", "", Phylum_clean)]
rel_long <- rel_long[!is.na(Phylum_clean) & Phylum_clean != ""]

meta_sub <- meta[, .(sample_id, site_cleaned, age_cleaned_group)]
meta_sub[, sample_id := as.character(sample_id)]

merged <- merge(rel_long, meta_sub, by = "sample_id", all.x = FALSE)
cat("  merged samples: ", uniqueN(merged$sample_id), "\n", sep = "")

# ------------------------------------------------------------
# 3. Filter target sites and age groups
# ------------------------------------------------------------
cat("\n[4/6] Filtering target sites and age groups...\n")

merged <- merged[!is.na(site_cleaned) & site_cleaned != ""]
merged <- merged[!is.na(age_cleaned_group) & age_cleaned_group != ""]

# Normalize site_cleaned: trim + lowercase (metadata stores lowercase)
merged[, site_cleaned := tolower(trimws(as.character(site_cleaned)))]
# Normalize fullwidth ">" to ASCII ">" in age labels
merged[, age_cleaned_group := gsub("＞", ">", age_cleaned_group)]

# --- Diagnostics: list unique values present in the data ---
site_uniq <- sort(unique(merged$site_cleaned))
age_uniq  <- sort(unique(as.character(merged$age_cleaned_group)))

cat("  Unique site_cleaned values (n=", length(site_uniq), "):\n", sep = "")
for (s in site_uniq) cat("    [", s, "]\n", sep = "")
cat("  Unique age_cleaned_group values (n=", length(age_uniq), "):\n", sep = "")
for (a in age_uniq) cat("    [", a, "]\n", sep = "")

cat("  Targets requested:\n")
cat("    sites: ", paste(target_sites_lower, collapse = " | "), "\n", sep = "")
cat("    ages : ", paste(age_order,          collapse = " | "), "\n", sep = "")
cat("  Site matches:  ", sum(target_sites_lower %in% site_uniq), "/", length(target_sites_lower), "\n", sep = "")
cat("  Age  matches:  ", sum(age_order          %in% age_uniq),  "/", length(age_order),          "\n", sep = "")

merged <- merged[site_cleaned %in% target_sites_lower]
merged <- merged[age_cleaned_group %in% age_order]

# Map lowercase site to display label, then factor with desired display order
merged[, site_display := site_label_map[site_cleaned]]
merged[, site_display := factor(site_display, levels = target_sites_display)]
merged[, age_cleaned_group := factor(age_cleaned_group, levels = age_order)]

if (nrow(merged) == 0) {
  stop("No data left after filtering. ",
       "Check the diagnostic lists above and align target sites / age_order ",
       "with the actual values in the metadata columns.")
}

cat("  samples per site x age:\n")
print(merged[, .(n_samples = uniqueN(sample_id)),
             by = .(site_display, age_cleaned_group)][
               order(site_display, age_cleaned_group)])

# ------------------------------------------------------------
# 4. Top 9 phyla + Others, mean abundance per (site, age)
# ------------------------------------------------------------
cat("\n[5/6] Computing Top 9 phyla + Others...\n")

global_top9 <- merged[, .(Total = sum(Abundance, na.rm = TRUE)),
                      by = Phylum_clean][order(-Total)][
                        seq_len(min(9, .N))]$Phylum_clean

cat("  Top 9 phyla (by global total abundance, descending):\n")
for (i in seq_along(global_top9)) {
  cat("    ", i, ". ", global_top9[i], "\n", sep = "")
}
cat("    10. Others (all remaining phyla merged)\n")

merged[, Phylum_group := ifelse(Phylum_clean %in% global_top9,
                                Phylum_clean, "Others")]

plot_data <- merged[, .(Mean_Abundance = mean(Abundance, na.rm = TRUE)),
                    by = .(site_display, age_cleaned_group, Phylum_group)]

# Renormalize each (site, age) bar to 100%
plot_data[, Mean_Abundance := Mean_Abundance / sum(Mean_Abundance),
          by = .(site_display, age_cleaned_group)]

# ------------------------------------------------------------
# 5. Color palette + stacking/legend order
# ------------------------------------------------------------
color_palette <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999"
)
unified_colors <- setNames(color_palette[seq_along(global_top9)], global_top9)
unified_colors["Others"] <- "#CCCCCC"

# In ggplot2 stack, the LAST factor level is drawn at the BOTTOM of the bar.
# Putting the dominant phylum at the end places it at the bar bottom
# and at the bottom of the legend; Others appears at the top of the legend.
dominant_phylum <- global_top9[1]
other_phyla     <- setdiff(global_top9, dominant_phylum)
legend_order    <- c("Others", rev(other_phyla), dominant_phylum)

plot_data[, Phylum_group := factor(Phylum_group, levels = legend_order)]

cat("\n  Stacking + legend order (top -> bottom):\n")
for (i in seq_along(legend_order)) {
  pos <- if (i == 1) "top" else if (i == length(legend_order)) "bottom (dominant)" else "middle"
  cat("    ", i, ". ", legend_order[i], "  [", pos, "]\n", sep = "")
}

# ------------------------------------------------------------
# 6. Plot + save
# ------------------------------------------------------------
cat("\n[6/6] Plotting and saving...\n")

p <- ggplot(plot_data,
            aes(x = age_cleaned_group, y = Mean_Abundance, fill = Phylum_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.95) +
  facet_grid(. ~ site_display, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = unified_colors, name = "Phylum",
                    breaks = legend_order, labels = legend_order) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Age Group", y = "Relative Abundance") +
  theme_minimal() +
  theme(
    axis.text.x   = element_text(size = 26, face = "bold", colour = "black",
                                 angle = 45, hjust = 1, vjust = 1,
                                 margin = margin(t = 5)),
    axis.text.y   = element_text(size = 24, face = "bold", colour = "black"),
    axis.title.x  = element_text(size = 32, face = "bold", colour = "black",
                                 margin = margin(t = 14)),
    axis.title.y  = element_text(size = 32, face = "bold", colour = "black",
                                 margin = margin(r = 12)),
    legend.position   = "right",
    legend.title      = element_text(size = 30, face = "bold", colour = "black"),
    legend.text       = element_text(size = 28, face = "bold", colour = "black"),
    legend.key.size   = unit(1.6, "lines"),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    strip.text        = element_text(size = 32, face = "bold", colour = "black",
                                     margin = margin(6, 6, 6, 6)),
    strip.background  = element_rect(fill = "lightblue", color = "black",
                                     linewidth = 0.6),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.spacing     = unit(0.3, "lines"),
    axis.line    = element_line(color = "black", linewidth = 0.6),
    axis.ticks   = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

png_path <- file.path(out_dir, "phylum_site_age_stacked-0610.png")
pdf_path <- file.path(out_dir, "phylum_site_age_stacked-0610.pdf")
csv_path <- file.path(out_dir, "phylum_site_age_mean_abundance-0610.csv")

ggsave(png_path, plot = p, width = 22, height = 9, dpi = 300, limitsize = FALSE)
ggsave(pdf_path, plot = p, width = 22, height = 9, limitsize = FALSE)

# Summary table: rows = phylum, columns = site|age
plot_data[, site_age := paste(site_display, age_cleaned_group, sep = "|")]
summary_wide <- data.table::dcast(plot_data, Phylum_group ~ site_age,
                                  value.var = "Mean_Abundance")
fwrite(summary_wide, csv_path)

cat("\n==========================================\n")
cat("DONE\n")
cat("Output files:\n")
cat("  ", png_path, "\n")
cat("  ", pdf_path, "\n")
cat("  ", csv_path, "\n")
cat("End:          ", format(Sys.time()), "\n")
cat("==========================================\n")
