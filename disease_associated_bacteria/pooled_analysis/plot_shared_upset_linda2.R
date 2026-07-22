#!/usr/bin/env Rscript
# UpSet for each method's significant genera across diseases (da_v2 linda2).
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ComplexUpset)
})

WORK_DIR <- "/home/zhangcaiyi/project/global_oral_diseae_da/da_v2"
OUT_DIR <- file.path(WORK_DIR, "method_comparison_linda2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

disease_short <- c(
  "Dental Caries" = "Dental Caries",
  "Head and Neck Cancer" = "H&N Cancer",
  "HIV" = "HIV", "OSCC" = "OSCC",
  "Periodontal Disease" = "Periodontal Dz",
  "RA" = "RA", "Breast Cancer" = "Breast Cancer",
  "IBD" = "IBD", "Lung Cancer" = "Lung Cancer"
)

make_upset <- function(sig_path, summary_path, title, out_prefix) {
  if (!file.exists(sig_path)) {
    cat("Skip missing:", sig_path, "\n"); return(invisible(NULL))
  }
  sig <- read.delim(sig_path, stringsAsFactors = FALSE, check.names = FALSE)
  sig <- sig %>% filter(!is.na(Genus), Genus != "")
  if (nrow(sig) == 0) {
    cat("No significant genera for", title, "\n"); return(invisible(NULL))
  }
  sig$DiseaseLabel <- ifelse(sig$Disease %in% names(disease_short),
                             disease_short[sig$Disease], sig$Disease)
  if (file.exists(summary_path)) {
    summ <- read.csv(summary_path, stringsAsFactors = FALSE)
    bp_col <- if ("Matched_BioProjects" %in% names(summ)) {
      "Matched_BioProjects"
    } else if ("Usable_BioProjects" %in% names(summ)) {
      "Usable_BioProjects"
    } else {
      NA_character_
    }
    if (!is.na(bp_col)) {
      bp_map <- setNames(summ[[bp_col]], summ$Disease)
      sig$DiseaseLabel <- vapply(seq_len(nrow(sig)), function(i) {
        d <- sig$Disease[i]
        lab <- sig$DiseaseLabel[i]
        nbp <- bp_map[[d]]
        if (is.null(nbp) || is.na(nbp)) lab else sprintf("%s [%dBP]", lab, nbp)
      }, character(1))
    }
  }
  counts_df <- sig %>%
    distinct(Disease, DiseaseLabel, Genus) %>%
    count(Disease, DiseaseLabel, name = "SignificantGenera") %>%
    arrange(desc(SignificantGenera), Disease)
  write.csv(counts_df, paste0(out_prefix, "_sig_counts.csv"), row.names = FALSE)

  disease_sets <- split(sig$Genus, sig$DiseaseLabel)
  disease_sets <- lapply(disease_sets, unique)
  disease_sets <- disease_sets[sapply(disease_sets, length) > 0]
  disease_sets <- disease_sets[order(sapply(disease_sets, length), decreasing = TRUE)]
  all_genera <- sort(unique(unlist(disease_sets)))
  mat <- as.data.frame(lapply(disease_sets, function(g) all_genera %in% g))
  rownames(mat) <- all_genera
  colnames(mat) <- names(disease_sets)
  mat[] <- lapply(mat, as.integer)
  set_cols <- colnames(mat)
  p <- upset(
    mat, set_cols, name = "Disease", width_ratio = 0.25, min_size = 1,
    sort_intersections_by = "cardinality", sort_sets = FALSE,
    base_annotations = list(
      "Shared genera count" = intersection_size(text = list(size = 3))
    ),
    set_sizes = (upset_set_size() + ylab("Significant genera per disease")),
    themes = upset_default_themes(text = element_text(size = 10))
  ) + ggtitle(title)
  w <- max(10, length(set_cols) * 0.9 + 4)
  ggsave(paste0(out_prefix, "_upset.pdf"), p, width = w, height = 7.2, limitsize = FALSE)
  ggsave(paste0(out_prefix, "_upset.png"), p, width = w, height = 7.2, dpi = 300,
         limitsize = FALSE)
  cat("Saved UpSet:", out_prefix, "\n")
}

make_upset(
  file.path(WORK_DIR, "maaslin3", "maaslin3_significant_genera_q0.05.tsv"),
  file.path(WORK_DIR, "maaslin3", "maaslin3_run_summary.csv"),
  "MaAsLin3 (LMM): Shared Disease-associated Genera",
  file.path(OUT_DIR, "maaslin3")
)
make_upset(
  file.path(WORK_DIR, "ancombc2", "ancombc2_significant_genera_q0.05.tsv"),
  file.path(WORK_DIR, "ancombc2", "ancombc2_run_summary.csv"),
  "ANCOM-BC2 (LMM): Shared Disease-associated Genera",
  file.path(OUT_DIR, "ancombc2")
)
make_upset(
  file.path(WORK_DIR, "linda2", "linda_significant_genera_q0.05.tsv"),
  file.path(WORK_DIR, "linda2", "linda_run_summary.csv"),
  "LinDA (LMM, no winsor): Shared Disease-associated Genera",
  file.path(OUT_DIR, "linda2")
)
