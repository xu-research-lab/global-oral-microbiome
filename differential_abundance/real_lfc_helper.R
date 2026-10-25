# ==============================================================================
# real_lfc_helper.R — shared by the da_v2 consensus bubble plot scripts.
#
# Real relative-abundance log2FC for a two-group contrast:
#   log2( (mean relabund in TARGET group + ps) / (mean relabund in REF group + ps) )
# i.e. the log2 ratio of group arithmetic means. ps = 1e-6 keeps genera that
# are completely absent in one group on a large-but-finite fold change instead
# of +/-Inf (agreed convention for all da_v2 bubble plots, 2026-09-21).
#
# Usage (in a plot script):
#   source("differential_abundance/real_lfc_helper.R")
#   lfc <- real_lfc_matrix(sdf, sdf$sex, ref = "F", tgt = "M")
# ==============================================================================
REAL_LFC_PS <- 1e-6

# s   : data.frame/matrix, samples x features of relative abundance
# g   : group vector/factor, length nrow(s)
# ref : reference group level (denominator, log2FC < 0 = enriched in ref)
# tgt : target group level (numerator, log2FC > 0 = enriched in tgt)
# returns a named numeric vector over colnames(s)
real_lfc_matrix <- function(s, g, ref, tgt, ps = REAL_LFC_PS) {
  g <- as.character(g)
  ir <- which(g == ref)
  it <- which(g == tgt)
  if (length(ir) == 0 || length(it) == 0)
    stop("real_lfc_matrix: empty group (ref='", ref, "' n=", length(ir),
         ", tgt='", tgt, "' n=", length(it), ")")
  m_ref <- colMeans(s[ir, , drop = FALSE])
  m_tgt <- colMeans(s[it, , drop = FALSE])
  out <- log2((m_tgt + ps) / (m_ref + ps))
  names(out) <- colnames(s)
  out
}