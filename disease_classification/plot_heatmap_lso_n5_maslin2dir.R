# RF LOSO heatmap: row set = the genera listed in n_5.csv ONLY (no top-10
# union). Cell direction always comes from the MaAsLin2 Disease_Status term
# in result/maaslin2_all_results.tsv (coef > 0 = higher in disease vs HC),
# regardless of significance -- every genus with a coefficient is filled by
# its direction. The asterisk marks the (disease, genus) pairs listed in
# n_5.csv, i.e. the DA-significant entries the shortlist records; genera not
# significant for a given disease in n_5 stay unstarred even when filled.
# Rows are NOT clustered: ordered by the count of n_5 disease-higher entries
# (descending), then n_5 HC-higher entries (ascending), so the bottom runs to
# the all-negative genera; ties break on total mean RF importance.
# Run: Rscript disease_classification/plot_heatmap_lso_n5_maslin2dir.R
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtext)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE))
} else if (file.exists("result/RF_feature_importances_lsoP.csv")) {
  "."
} else {
  "disease_classification"
}
RESULT_DIR <- file.path(SCRIPT_DIR, "result")
# TOP_N is the rank-label cutoff: digits are shown for full-set ranks <= TOP_N.
# Usage: Rscript plot_heatmap_lso_n5_maslin2dir.R [TOP_N]   (default 10)
args <- commandArgs(trailingOnly = TRUE)
TOP_N <- if (length(args)) as.integer(args[1]) else 10L
stopifnot(!is.na(TOP_N), TOP_N >= 1)
OUT_BASE <- file.path(RESULT_DIR, paste0("fi_heatmap_lso_genus_R_n5_maslin2dir_top", TOP_N))

normalize_disease <- function(x) {
  gsub("[^a-z0-9]", "", tolower(trimws(x)))
}
normalize_genus <- function(x) sub("^g__", "", trimws(x))

# ---------------------------------------------------------------- FI matrix
# Mean FI per (genus, comparison) straight from the lsoP fold importances.
fi <- read.csv(file.path(RESULT_DIR, "RF_feature_importances_lsoP.csv"),
               stringsAsFactors = FALSE)
fi$col <- sub("_(qiita_)?[^_]+$", "", fi$fold)
w <- fi %>%
  group_by(feature, col) %>%
  summarise(importance = mean(importance), .groups = "drop") %>%
  pivot_wider(names_from = col, values_from = importance)
imp_mat <- as.matrix(w[, -1])
rownames(imp_mat) <- w$feature
imp_mat <- imp_mat[, sort(colnames(imp_mat)), drop = FALSE]
rank_mat <- apply(imp_mat, 2, function(x) rank(-x, ties.method = "min", na.last = "keep"))
dimnames(rank_mat) <- dimnames(imp_mat)
# NOTE: no clustering anywhere -- rows get an explicit significance-count order
# further down, columns are ordered by n_5 significant count.

# Rank labels use the full-comparison rank (over all prevalence-standard
# genera); the top-10 rule applies only to labels and the row subset.
rank_top <- rank_mat

# ---------------------------------------------------------------- row set
# only the genera listed in n_5.csv (curated DA shortlist)
n5 <- as.data.frame(fread(file.path(RESULT_DIR, "n_5.csv"), stringsAsFactors = FALSE))
n5_map <- n5 %>%
  transmute(disease_key = normalize_disease(Disease),
            genus_key = normalize_genus(Genus),
            n5_direction = trimws(Direction)) %>%
  filter(n5_direction %in% c("+", "-", "=")) %>%
  distinct()
if (anyDuplicated(paste(n5_map$disease_key, n5_map$genus_key, sep = "|"))) {
  stop("Multiple n_5 rows for a disease + genus; resolve before plotting.")
}
n5_genera <- unique(n5_map$genus_key)

keep_rows <- normalize_genus(rownames(rank_top)) %in% n5_genera
rank_top <- rank_top[keep_rows, , drop = FALSE]
imp_mat <- imp_mat[keep_rows, , drop = FALSE]
missing_n5 <- setdiff(n5_genera, normalize_genus(rownames(rank_top)))
if (length(missing_n5)) stop("n_5 genera missing from the lsoP FI set: ",
                             paste(missing_n5, collapse = ", "))
cat(sprintf("rows kept: %d n_5 genera\n", nrow(rank_top)))

# Disease abbreviations supplied by the user; keys ignore spaces/punctuation.
DISEASE_ABBR <- c(
  acutemyeloidleukemia = "AML", alcoholusedisorder = "AUD",
  alzheimersdisease = "AD", aphthous = "AS", asthma = "Asthma",
  atopicdermatitis = "AtD", autoimmunehepatitis = "AIH",
  axialspondyloarthritis = "axSpA", breastcancer = "BC",
  celiacdisease = "Celiac Disease", chagasdisease = "CD",
  chronicobstructivepulmonarydisease = "COPD", copd = "COPD",
  clostridiumdifficileinfection = "CDI", colorectalcancer = "CRC",
  coronavirusdisease2019 = "COVID-19", covid19 = "COVID-19",
  dentalcaries = "DC", depression = "Depression", diabetes = "Diabetes",
  digestivetractcancer = "DTC", esophagealcancer = "EC",
  gastroesophagealrefluxdisease = "GERD", handfootandmouthdisease = "HFMD",
  headandneckcancer = "HNC", helicobacterpyloriinfection = "HPI",
  hypertension = "HTN", humanimmunodeficiencyvirus = "HIV", hiv = "HIV",
  inflammatoryboweldisease = "IBD", ibd = "IBD",
  intensivecareunit = "ICU", icu = "ICU", lungcancer = "LC",
  noncariouscervicallesions = "NCCL", nonmalignantbreastdisease = "NBD",
  obesity = "Obesity", obstructivesleepapnea = "OSA", orallichenplanus = "OLP",
  oralmucositis = "OM", oralsquamouscellcarcinoma = "OSCC", oscc = "OSCC",
  oralsubmucousfibrosis = "OSMF", oralverrucoushyperplasia = "OVH",
  pancreaticcancer = "Pancreatic Cancer", parkinsonsdisease = "Parkinson's Disease",
  periodontaldisease = "PD", posttraumaticstressdisorder = "PTSD",
  preterm = "PTB", prostatecancer = "PC", rheumatoidarthritis = "RA", ra = "RA",
  sjgrenssyndrome = "SS", sjogrenssyndrome = "SS", thyroidcancer = "TC"
)

# Diseases ordered by the number of n_5-significant entries (in the displayed
# comparisons) for that disease, descending; saliva/tooth columns of one
# disease share the same significance pattern (direction and n_5 entries are
# disease-level), so they stay grouped. Ties fall back to the reference AUC
# rule: descending number of sites, then median AUC, then name.
auc_data <- bind_rows(lapply(c("saliva", "tooth"), function(site) {
  read.csv(file.path(RESULT_DIR,
                     paste0("rf_loso_auc_all_diseases_genus_", site, ".csv")),
           stringsAsFactors = FALSE) %>%
    mutate(site = site)
})) %>%
  filter(!is.na(auc))
column_disease <- sub("^[^_]+_", "", colnames(rank_top))
n5_genera_shown <- normalize_genus(rownames(rank_top))
disease_sig <- n5_map %>%
  filter(genus_key %in% n5_genera_shown) %>%
  group_by(disease_key) %>%
  summarise(n_sig = n(), .groups = "drop")
disease_order <- auc_data %>%
  group_by(disease) %>%
  summarise(n_sites = n_distinct(site), median_auc = median(auc),
            .groups = "drop") %>%
  mutate(disease_key = normalize_disease(disease)) %>%
  left_join(disease_sig, by = "disease_key") %>%
  mutate(n_sig = coalesce(n_sig, 0L)) %>%
  arrange(desc(n_sig), desc(n_sites), desc(median_auc), disease) %>%
  pull(disease)
column_disease <- sub("^[^_]+_", "", colnames(rank_top))
missing_auc <- setdiff(column_disease, disease_order)
if (length(missing_auc)) {
  stop("Disease missing from reference AUC tables: ",
       paste(missing_auc, collapse = ", "))
}
ordered_columns <- unlist(lapply(disease_order, function(disease) {
  cols <- colnames(rank_top)[column_disease == disease]
  sites <- sub("_.*$", "", cols)
  cols[order(match(sites, c("saliva", "tooth")), na.last = TRUE)]
}), use.names = FALSE)
rank_top <- rank_top[, ordered_columns, drop = FALSE]
imp_mat <- imp_mat[, ordered_columns, drop = FALSE]
column_disease <- sub("^[^_]+_", "", colnames(rank_top))
column_sites <- sub("_.*$", "", colnames(rank_top))
column_abbr <- unname(DISEASE_ABBR[normalize_disease(column_disease)])
if (anyNA(column_abbr)) {
  stop("Disease abbreviation missing: ",
       paste(column_disease[is.na(column_abbr)], collapse = ", "))
}
column_labels <- ifelse(column_sites == "tooth",
                        paste0(column_abbr, " (tooth)"), column_abbr)
names(column_labels) <- colnames(rank_top)

# --------------------------------------- directions from maaslin2_all_results.tsv
# The disease term lives in rows with metadata == "Disease_Status" (one row per
# disease level vs the HC reference); coef > 0 means higher in that disease.
# The direction is disease-level (the model carries site_cleaned as a
# covariate), so it applies to every site column of that disease. Any genus
# with a disease-term coefficient is filled by its direction, significant or
# not -- significance is NOT required for colouring (stars come from n_5).
# Cells with no entry keep the empty background tile.
maslin2 <- fread(file.path(RESULT_DIR, "maaslin2_all_results.tsv"),
                 stringsAsFactors = FALSE, sep = "\t", header = TRUE)
maslin2 <- maslin2 %>%
  filter(metadata == "Disease_Status") %>%
  transmute(disease_key = normalize_disease(Disease),
            genus_key = normalize_genus(feature),
            coef = as.numeric(coef),
            qval = as.numeric(qval))
if (anyDuplicated(paste(maslin2$disease_key, maslin2$genus_key, sep = "|"))) {
  stop("Multiple Disease_Status rows for a disease + genus; resolve before plotting.")
}
direction_map <- maslin2 %>%
  mutate(direction = ifelse(is.finite(coef) & coef > 0, "+",
                            ifelse(is.finite(coef) & coef < 0, "-", "=")))

# --------------------------------------------------------- explicit row order
# No clustering. Count the n_5 entries (restricted to the displayed
# comparisons) per genus: disease-higher entries descending, then HC-higher
# entries ascending, so the tail runs from few negatives down to the
# all-negative genera. Ties fall back to total mean importance (descending)
# for stable output.
row_dir <- bind_rows(lapply(structure(as.list(colnames(rank_top)), names = colnames(rank_top)),
                            function(cmp) {
  data.frame(feature = rownames(rank_top), comparison = cmp,
             importance = as.numeric(imp_mat[, cmp]),
             stringsAsFactors = FALSE)
})) %>%
  mutate(disease_key = normalize_disease(sub("^[^_]+_", "", comparison)),
         genus_key = normalize_genus(feature)) %>%
  left_join(n5_map, by = c("disease_key", "genus_key")) %>%
  mutate(has_importance = is.finite(importance),
         in_n5 = !is.na(n5_direction)) %>%
  group_by(feature) %>%
  summarise(n_pos = sum(in_n5 & n5_direction == "+"),
            n_neg = sum(in_n5 & n5_direction == "-"),
            imp_sum  = sum(importance[has_importance], na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(n_pos), n_neg, desc(imp_sum))
new_rows <- row_dir$feature
rank_top <- rank_top[new_rows, , drop = FALSE]
imp_mat  <- imp_mat[new_rows, , drop = FALSE]
cat(sprintf("row order by n_5 significant-count: max n_pos = %d, %d rows with no n_5 entry here\n",
            max(row_dir$n_pos), sum(row_dir$n_pos == 0 & row_dir$n_neg == 0)))

long <- as.data.frame(rank_top) %>%
  mutate(feature = rownames(rank_top)) %>%
  pivot_longer(-feature, names_to = "comparison", values_to = "rank") %>%
  left_join(as.data.frame(imp_mat) %>%
              mutate(feature = rownames(imp_mat)) %>%
              pivot_longer(-feature, names_to = "comparison", values_to = "importance"),
            by = c("comparison", "feature")) %>%
  mutate(site = sub("_.*$", "", comparison),
         disease = sub("^[^_]+_", "", comparison),
         disease_key = normalize_disease(disease),
         genus_key = normalize_genus(feature)) %>%
  left_join(direction_map, by = c("disease_key", "genus_key")) %>%
  left_join(n5_map, by = c("disease_key", "genus_key")) %>%
  mutate(is_top10 = is.finite(rank) & rank >= 1 & rank <= TOP_N,
         has_importance = is.finite(importance),
         has_direction = !is.na(direction) & direction %in% c("+", "-", "="),
         in_consensus = has_direction,          # fill any coefficient entry
         significant = !is.na(n5_direction),    # asterisk = n_5 entry only
         signed_importance = case_when(
           !has_importance ~ NA_real_,
           !has_direction ~ NA_real_,  # no MaAsLin2 disease term -> no fill
           direction == "+" ~ importance,
           direction == "-" ~ -importance,
           direction == "=" ~ 0,
           TRUE ~ NA_real_
         ),
         label = ifelse(is_top10, sprintf("%d", as.integer(rank)), ""),
         # single centred cell label: rank digit with a trailing "*" when the
         # pair is a n_5 entry, or "*" alone when only significant
         label_cell = case_when(
           is_top10 & significant ~ paste0(label, "*"),
           is_top10 ~ label,
           significant ~ "*",
           TRUE ~ ""),
         xr = match(comparison, colnames(rank_top)),
         yr = match(feature, rev(rownames(rank_top))))

plot_data <- long %>%
  filter(has_importance) %>%
  mutate(comparison_label = unname(column_labels[comparison])) %>%
  select(comparison, comparison_label, site, disease, feature, importance,
         signed_importance, rank, is_top10, direction, in_consensus,
         significant, has_direction, has_importance, coef, qval, label)
write.csv(plot_data, paste0(OUT_BASE, "_data.csv"),
          row.names = FALSE, na = "")

importance_limit <- max(plot_data$importance)
row_labels <- normalize_genus(rev(rownames(rank_top)))
# Socransky red-complex genera highlighted red on the y axis:
# Porphyromonas_A = P. gingivalis, Tannerella = T. forsythia,
# and every Treponema row (T. denticola sits in Treponema_B_986141 in GTDB;
# the user asked to mark all Treponema entries).
RED_COMPLEX <- c("Porphyromonas_A", "Tannerella")
row_labels_html <- ifelse(row_labels %in% RED_COMPLEX |
                            grepl("^Treponema", row_labels),
                          paste0("<span style='color:#D65A50'>", row_labels, "</span>"),
                          row_labels)
# Fonts: real Arial (installed to ~/.fonts) resolves through fontconfig, so
# no systemfonts registration is needed and rich text keeps its italic face.
# Genus names are drawn via geom_richtext; axis/legend text uses family Arial.
p_heat <- ggplot(long, aes(xr, yr)) +
  geom_tile(fill = "#ECECEC", color = "white", linewidth = 0.3) +
  geom_tile(data = filter(long, has_importance, in_consensus),
            aes(fill = signed_importance),
            color = "white", linewidth = 0.3) +
  # centred rank digit ("17" or ""), default font -- one size for digits only
  geom_text(data = filter(long, is_top10 & has_importance),
            aes(label = label), color = "black", size = 2.3,
            show.legend = FALSE) +
  # n_5 asterisk in its OWN layer so it can be LARGER than the digits without
  # touching them (a single geom_text can only carry one size). Default font
  # as well: paired stars sit just right of the digit -- the offset widens with
  # the digit count so two-digit ranks ("13*") still clear the glyph -- with a
  # small downward nudge; lone stars get a larger downward nudge so their ink
  # sits on the tile centre (Liberation Sans "*" ink centre is 339.5 font
  # units ABOVE the advance-box centre -> 0.56 mm at size 3.4 = 0.18 tiles,
  # measured with fontTools, not eyeballed).
  geom_text(data = filter(long, significant & has_importance),
            aes(xr + ifelse(is_top10, 0.08 + 0.13 * nchar(label), 0),
                yr + ifelse(is_top10, -0.05, -0.18),
                label = "*"),
            color = "black", size = 3.4,
            show.legend = FALSE) +
  scale_fill_gradient2(low = "#397FB1", mid = "#FFFFFF", high = "#D65A50",
                       midpoint = 0, limits = c(-importance_limit, importance_limit),
                       na.value = "#ECECEC",
                       name = "Directional RF importance<br>
(&gt;0 higher in disease,<br>
&lt;0 higher in HC)",
                       guide = guide_colorbar(barheight = grid::unit(3.2, "cm"),
                                              barwidth = grid::unit(0.4, "cm"))) +
  scale_x_continuous(breaks = seq_len(ncol(rank_top)),
                     labels = unname(column_labels),
                     expand = expansion(add = 0.5)) +
  scale_y_continuous(breaks = NULL, expand = expansion(add = 0.5)) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  # y-axis genus labels: drawn as a geom_richtext layer instead of axis text.
  # ggplot2 4.x + this ggtext version renders element_markdown axis text as
  # plain text (raw <span> markup, no colour); geom_richtext keeps the
  # red-complex highlight working on agg_png and pdf devices.
  geom_richtext(data = data.frame(yr = seq_len(nrow(rank_top)),
                                  label = row_labels_html),
                aes(y = yr, label = label), inherit.aes = FALSE,
                x = 0.45, hjust = 1, vjust = 0.5,
                family = "Arial", fontface = "italic", size = 7 / .pt,
                fill = NA, label.color = NA, label.padding = grid::unit(0, "lines"),
                show.legend = FALSE) +
  theme_minimal(base_size = 8, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               size = 7, colour = "black", family = "Arial"),
    axis.text.y = element_blank(),
    legend.title = ggtext::element_markdown(size = 8, lineheight = 0.9,
                                            hjust = 0, colour = "black",
                                            family = "Arial"),
    legend.spacing.x = grid::unit(3, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = 8, family = "Arial"),
    legend.key.size = grid::unit(0.35, "cm"),
    # left margin reserves room for the richtext genus labels, which extend
    # beyond the panel (clip = "off"); tune with the longest genus name
    plot.margin = margin(3, 3, 3, 70),
    panel.grid = element_blank()
  )

FIG_WIDTH <- 6.4 * 0.9
FIG_HEIGHT <- 0.105 * nrow(rank_top) + 0.95
ggsave(paste0(OUT_BASE, ".png"), p_heat, width = FIG_WIDTH, height = FIG_HEIGHT,
       device = ragg::agg_png, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(paste0(OUT_BASE, ".pdf"), p_heat, width = FIG_WIDTH, height = FIG_HEIGHT,
       device = grDevices::cairo_pdf, limitsize = FALSE)
summary_by_comparison <- plot_data %>%
  group_by(comparison) %>%
  summarise(n_available = n(), n_top10 = sum(is_top10),
            n_dir = sum(has_direction & has_importance, na.rm = TRUE),
            n_sig = sum(significant & has_importance, na.rm = TRUE),
            n_case = sum(has_direction & has_importance & direction == "+", na.rm = TRUE),
            n_hc = sum(has_direction & has_importance & direction == "-", na.rm = TRUE),
            .groups = "drop")
print(as.data.frame(summary_by_comparison), row.names = FALSE)
cat("Column order:", paste(column_labels, collapse = " | "), "\n")
cat("Direction-coloured cells (all available):",
    sum(long$has_direction & long$has_importance, na.rm = TRUE),
    ";  asterisks (n_5 entries):", sum(long$significant & long$has_importance), "\n")
cat(sprintf("[write] %s.{png,pdf} and _data.csv\n", OUT_BASE))
cat(sprintf("[size] %.2f x %.2f inches; %d genera x %d comparisons\n",
            FIG_WIDTH, FIG_HEIGHT,
            nrow(rank_top), ncol(rank_top)))