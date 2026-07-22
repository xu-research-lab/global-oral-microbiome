#!/usr/bin/env Rscript
# Three-way comparison: Genus + Direction, q <= 0.05 (da_v2 linda2).
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggVennDiagram)
  library(patchwork)
})

WORK_DIR <- "/home/zhangcaiyi/project/global_oral_diseae_da/da_v2"
OUT_DIR <- file.path(WORK_DIR, "method_comparison_linda2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MAASLIN_SIG <- file.path(WORK_DIR, "maaslin3", "maaslin3_significant_genera_q0.05.tsv")
ANCOMBC_SIG <- file.path(WORK_DIR, "ancombc2", "ancombc2_significant_genera_q0.05.tsv")
LINDA_SIG <- file.path(WORK_DIR, "linda2", "linda_significant_genera_q0.05.tsv")
RUN_SUMMARY <- file.path(WORK_DIR, "maaslin3", "maaslin3_run_summary.csv")

for (p in c(MAASLIN_SIG, ANCOMBC_SIG, LINDA_SIG, RUN_SUMMARY)) {
  if (!file.exists(p)) stop("Missing: ", p)
}

bp_by_disease <- {
  summ <- read.csv(RUN_SUMMARY, stringsAsFactors = FALSE)
  setNames(summ$Matched_BioProjects, summ$Disease)
}

m3 <- read.delim(MAASLIN_SIG, stringsAsFactors = FALSE, check.names = FALSE)
ac <- read.delim(ANCOMBC_SIG, stringsAsFactors = FALSE, check.names = FALSE)
ld <- read.delim(LINDA_SIG, stringsAsFactors = FALSE, check.names = FALSE)

short_to_full <- c(
  "Dental Caries" = "Dental Caries",
  "H&N Cancer" = "Head and Neck Cancer",
  "HIV" = "HIV", "OSCC" = "OSCC",
  "Periodontal Dz" = "Periodontal Disease",
  "RA" = "RA", "Breast Cancer" = "Breast Cancer",
  "IBD" = "IBD", "Lung Cancer" = "Lung Cancer"
)
normalize_disease <- function(x) {
  x <- as.character(x)
  ifelse(x %in% names(short_to_full), short_to_full[x], x)
}

prep_sig <- function(df, method) {
  if (!("Direction" %in% names(df))) stop("Missing Direction in ", method)
  if (!("Genus" %in% names(df))) stop("Missing Genus in ", method)
  df %>%
    mutate(
      DiseaseFull = normalize_disease(Disease),
      Direction = as.character(Direction),
      Genus = as.character(Genus),
      key = paste(Genus, Direction, sep = "|")
    ) %>%
    filter(!is.na(Genus), Genus != "", Direction %in% c("+", "-")) %>%
    group_by(DiseaseFull, Genus) %>%
    slice(1) %>%
    ungroup()
}

m3 <- prep_sig(m3, "MaAsLin3")
ac <- prep_sig(ac, "ANCOM-BC2")
ld <- prep_sig(ld, "LinDA2")

make_sets <- function(df) {
  lapply(split(df$key, df$DiseaseFull), function(x) unique(x[!is.na(x) & x != ""]))
}
m3_sets <- make_sets(m3)
ac_sets <- make_sets(ac)
ld_sets <- make_sets(ld)
diseases <- sort(Reduce(union, list(names(m3_sets), names(ac_sets), names(ld_sets))))

key_to_label <- function(keys) {
  if (length(keys) == 0) return("")
  paste(sort(keys), collapse = ";")
}

overlap_rows <- list()
conflict_rows <- list()
plot_list <- list()

for (d in diseases) {
  a <- if (d %in% names(m3_sets)) m3_sets[[d]] else character()
  b <- if (d %in% names(ac_sets)) ac_sets[[d]] else character()
  cset <- if (d %in% names(ld_sets)) ld_sets[[d]] else character()

  all3 <- Reduce(intersect, list(a, b, cset))
  m3_ac <- setdiff(intersect(a, b), cset)
  m3_ld <- setdiff(intersect(a, cset), b)
  ac_ld <- setdiff(intersect(b, cset), a)
  only_m3 <- setdiff(a, union(b, cset))
  only_ac <- setdiff(b, union(a, cset))
  only_ld <- setdiff(cset, union(a, b))

  overlap_rows[[d]] <- data.frame(
    Disease = d,
    n_maaslin3 = length(a), n_ancombc2 = length(b), n_linda = length(cset),
    n_all3 = length(all3),
    n_maaslin3_ancombc2 = length(m3_ac),
    n_maaslin3_linda = length(m3_ld),
    n_ancombc2_linda = length(ac_ld),
    n_only_maaslin3 = length(only_m3),
    n_only_ancombc2 = length(only_ac),
    n_only_linda = length(only_ld),
    all3 = key_to_label(all3),
    maaslin3_ancombc2 = key_to_label(m3_ac),
    maaslin3_linda = key_to_label(m3_ld),
    ancombc2_linda = key_to_label(ac_ld),
    only_maaslin3 = key_to_label(only_m3),
    only_ancombc2 = key_to_label(only_ac),
    only_linda = key_to_label(only_ld),
    stringsAsFactors = FALSE
  )

  g_m3 <- unique((m3 %>% filter(DiseaseFull == d))$Genus)
  g_ac <- unique((ac %>% filter(DiseaseFull == d))$Genus)
  g_ld <- unique((ld %>% filter(DiseaseFull == d))$Genus)
  shared_genus <- Reduce(union, list(
    intersect(g_m3, g_ac), intersect(g_m3, g_ld), intersect(g_ac, g_ld)
  ))
  for (g in shared_genus) {
    dirs <- c(
      MaAsLin3 = {
        x <- m3$Direction[m3$DiseaseFull == d & m3$Genus == g]
        if (length(x)) x[[1]] else NA_character_
      },
      ANCOMBC2 = {
        x <- ac$Direction[ac$DiseaseFull == d & ac$Genus == g]
        if (length(x)) x[[1]] else NA_character_
      },
      LinDA = {
        x <- ld$Direction[ld$DiseaseFull == d & ld$Genus == g]
        if (length(x)) x[[1]] else NA_character_
      }
    )
    present <- dirs[!is.na(dirs)]
    if (length(unique(present)) > 1) {
      conflict_rows[[paste(d, g, sep = "::")]] <- data.frame(
        Disease = d, Genus = g,
        Direction_MaAsLin3 = ifelse(is.na(dirs[["MaAsLin3"]]), "", dirs[["MaAsLin3"]]),
        Direction_ANCOMBC2 = ifelse(is.na(dirs[["ANCOMBC2"]]), "", dirs[["ANCOMBC2"]]),
        Direction_LinDA = ifelse(is.na(dirs[["LinDA"]]), "", dirs[["LinDA"]]),
        stringsAsFactors = FALSE
      )
    }
  }

  vlist <- list(MaAsLin3 = a, `ANCOM-BC2` = b, `LinDA (no winsor)` = cset)
  n_bp <- bp_by_disease[[d]]
  if (is.null(n_bp) || is.na(n_bp)) n_bp <- "?"
  panel_title <- sprintf("%s (BP=%s)", d, n_bp)
  p <- ggVennDiagram(vlist, label_alpha = 0) +
    scale_fill_gradient(low = "#F7FBFF", high = "#2171B5") +
    labs(title = panel_title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          legend.position = "none")
  plot_list[[d]] <- p
}

overlap_df <- bind_rows(overlap_rows) %>% arrange(desc(n_all3), Disease)
write.csv(overlap_df, file.path(OUT_DIR, "overlap_summary_3way.csv"), row.names = FALSE)

conflict_df <- if (length(conflict_rows) > 0) {
  bind_rows(conflict_rows) %>% arrange(Disease, Genus)
} else {
  data.frame(Disease = character(), Genus = character(),
             Direction_MaAsLin3 = character(), Direction_ANCOMBC2 = character(),
             Direction_LinDA = character(), stringsAsFactors = FALSE)
}
write.csv(conflict_df, file.path(OUT_DIR, "direction_conflicts.csv"), row.names = FALSE)

print(overlap_df[, c("Disease", "n_maaslin3", "n_ancombc2", "n_linda", "n_all3")])
cat(sprintf("Direction conflicts: %d\n", nrow(conflict_df)))

n <- length(plot_list)
ncol <- 3
nrow <- ceiling(n / ncol)
combined <- wrap_plots(plot_list, ncol = ncol) +
  plot_annotation(
    title = "da_v2: MaAsLin3 vs ANCOM-BC2 vs LinDA (no winsor)",
    subtitle = "q <= 0.05; overlap requires same Genus AND Direction",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )
w <- max(14, ncol * 4.5)
h <- max(9, nrow * 4.0)
ggsave(file.path(OUT_DIR, "venn_3way_per_disease.pdf"), combined,
       width = w, height = h, limitsize = FALSE)
ggsave(file.path(OUT_DIR, "venn_3way_per_disease.png"), combined,
       width = w, height = h, dpi = 300, limitsize = FALSE)
cat("Saved to ", OUT_DIR, "\n", sep = "")
