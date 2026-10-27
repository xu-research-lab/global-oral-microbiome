# Called only for --stage=fit/all; rendering needs none of the model packages.
set.seed(DA_SEED)
methods_to_run <- strsplit(method_arg, ",", fixed = TRUE)[[1]]
if (identical(methods_to_run, "all")) methods_to_run <- c("linda", "maaslin2", "ancombc2")
bad <- setdiff(methods_to_run, c("linda", "maaslin2", "ancombc2"))
if (length(bad)) stop("Unknown method: ", paste(bad, collapse = ", "))
for (m in methods_to_run) {
  pkg <- switch(m,
    linda = "MicrobiomeStat",
    maaslin2 = "Maaslin2",
    ancombc2 = "ANCOMBC"
  )
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
if ("ancombc2" %in% methods_to_run) suppressPackageStartupMessages(library(phyloseq))
sys.source(file.path(script_dir, "model.R"), envir = environment())

genus_counts <- load_genus_table("counts")
metadata <- if (analysis == "sex") load_sex_metadata(hc_only = TRUE) else load_age_metadata(hc_only = (scope == "hc"))
aligned <- align_samples(metadata, genus_counts)
metadata <- aligned$meta
genus_counts <- aligned$mat
rm(aligned)
if (analysis == "sex") {
  metadata <- metadata[metadata[[COL_SEX]] %in% c(SEX_REF, SEX_TGT), , drop = FALSE]
  genus_counts <- genus_counts[, metadata[[COL_SAMPLE]], drop = FALSE]
  if (any(table(factor(metadata[[COL_SEX]], levels = c(SEX_REF, SEX_TGT))) < MIN_SAMPLES_PER_SEX)) {
    stop("Pooled saliva needs at least ", MIN_SAMPLES_PER_SEX, " samples of each sex")
  }
  pairs_list <- list(list(
    site = SITE_SALIVA, continent = "pooled", g1 = SEX_REF, g2 = SEX_TGT,
    comparison_name = "saliva_pooled_F_vs_M"
  ))
} else {
  pairs_list <- make_age_ova_comparisons(metadata, SITE_ORDER)
}
if (!length(pairs_list)) stop("No valid comparisons")
if (nzchar(only_arg)) {
  nm <- vapply(pairs_list, `[[`, "", "comparison_name")
  idx <- if (grepl("^[0-9]+$", only_arg)) as.integer(only_arg) else match(only_arg, nm)
  if (length(idx) != 1L || is.na(idx) || idx < 1L || idx > length(pairs_list)) stop("Invalid --only: ", only_arg)
  pairs_list <- pairs_list[idx]
}
RUNNERS <- list(linda = run_pair_linda, maaslin2 = run_pair_maaslin2, ancombc2 = run_pair_ancombc2)
failed <- character()
for (m in methods_to_run) {
  for (i in seq_along(pairs_list)) {
    p <- pairs_list[[i]]
    # A completed fit does not need preprocessing again.
    if (file.exists(file.path(method_dirs[[m]], p$comparison_name, "analysis.done"))) next
    message(m, ": ", p$comparison_name)
    pp <- if (analysis == "sex") prep_pair(p$continent, p$g1, p$g2) else prep_ova(p$site, p$g2, p$other_groups)
    status <- RUNNERS[[m]](i, pp)
    if (status == "fail") failed <- c(failed, paste(m, p$comparison_name))
  }
}
if (length(failed)) stop("Failed fits: ", paste(failed, collapse = "; "))
