#!/usr/bin/env bash
# install_netmoss2.sh — download, patch and install the official NetMoss2
# (xiaolw95/NetMoss2, pinned commit 79da905 → v0.2.0) with a toy test that
# verifies the patch is live.
#
#   Usage:
#     bash install_netmoss2.sh [install_root]
#     RSCRIPT=/path/to/Rscript bash install_netmoss2.sh [install_root]
#
#   install_root (default ../external relative to this script) receives the
#   patched package source as <install_root>/NetMoss2 plus LOCAL_PATCH.diff
#   (the exact upstream->patched diff, generated without git involvement).
#
#   No git required: the pinned commit tarball is fetched from
#   codeload.github.com (reachable even where git-over-https is flaky;
#   override with UPSTREAM_SHA's equivalent pin if upstream moves).
#
# ── Why a patch is needed (see the header comment of netmoss_get_network.R) ──
# In the released directory-mode `getNetwork` (2023-09 source), the pooling
# accumulators iterate over SCALARS:
#     for (i in n1) / for (i in n2) / for (k in n1) / for (k in n2)
# where n1/n2 = length(list.files()) are scalars, so `for (k in n1)` runs once
# (k = n1): only the LAST study is pooled, all others silently dropped —
# contradicting the paper's all-cohort inverse-variance pooling.
# Fix: n1/n2 -> 1:n1/1:n2 at those four loop headers.
# Also `DESCRIPTION` imports `rsparcc`, which does not exist anywhere (CRAN,
# GitHub) and is never referenced by any R/ source — a phantom dependency that
# blocks source installation; it is deleted.
# Neither edit changes the algorithm the paper states. Toy test (2 matrices,
# gw1=0.5/n=30, gw2=-0.2/n=50) must return 0.1017185323
# (= (w_A*0.5 + w_B*(-0.2)) / (w_A + w_B), w_k = (n_k-1)/(1-r_k^2)).
# ────────────────────────────────────────────────────────────────────────────
set -euo pipefail

RSCRIPT=${RSCRIPT:-Rscript}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INSTALL_ROOT=${1:-"$SCRIPT_DIR/../external"}
UPSTREAM="https://github.com/xiaolw95/NetMoss2"
UPSTREAM_SHA="79da905efccfc8bec9c08a4f447bf1ddb12223f6"   # == v0.2.0 HEAD (2024)
PKGD="$(realpath -m "$INSTALL_ROOT/NetMoss2")"

if [ -d "$PKGD" ]; then
  echo "install_netmoss2.sh: patched source already at $PKGD — reuse (rm -rf to redo)"
else
  mkdir -p "$INSTALL_ROOT"
  # fetch pinned upstream snapshot (tarball: no git needed, works behind
  # firewalls where git-over-https is flaky)
  curl -fL --retry 5 --retry-delay 5 \
    "$UPSTREAM/archive/${UPSTREAM_SHA}.tar.gz" \
    -o "$INSTALL_ROOT/NetMoss2.tar.gz"
  tar -xzf "$INSTALL_ROOT/NetMoss2.tar.gz" -C "$INSTALL_ROOT"
  mv "$INSTALL_ROOT"/NetMoss2-* "$PKGD"
  rm -f "$INSTALL_ROOT/NetMoss2.tar.gz"

  cp "$PKGD/R/getNetwork.R" "$PKGD/R/getNetwork.R.orig"
  cp "$PKGD/DESCRIPTION"    "$PKGD/DESCRIPTION.orig"

  # patch 1: the four scalar loop headers (see rationale above)
  sed -i \
    -e 's/for (i in n1) {/for (i in 1:n1) {/' \
    -e 's/for (i in n2) {/for (i in 1:n2) {/' \
    -e 's/for (k in n1) {/for (k in 1:n1) {/' \
    -e 's/for (k in n2) {/for (k in 1:n2) {/' \
    "$PKGD/R/getNetwork.R"

  # patch 2: drop the phantom `rsparcc` dependency (blocks install, unused).
  # Upstream has it mid-list as "    rsparcc," (yaml follows), so a plain
  # line-delete is safe; the grep guard re-checks it was actually removed.
  DESCP="$PKGD/DESCRIPTION"
  sed -i '/^[ ]*rsparcc,[ ]*$/d' "$DESCP"
  $RSCRIPT -e "invisible(read.dcf('$DESCP'))" \
           || { echo "DESCRIPTION unreadable after patch"; exit 1; }
  /bin/grep -q 'rsparcc' "$DESCP" \
    && { echo "rsparcc still present in DESCRIPTION"; exit 1; } || true

  { diff -u "$PKGD/R/getNetwork.R.orig" "$PKGD/R/getNetwork.R" || true;
    diff -u "$PKGD/DESCRIPTION.orig"    "$PKGD/DESCRIPTION"    || true; } \
    > "$PKGD/LOCAL_PATCH.diff"
  rm -f "$PKGD/R/getNetwork.R.orig" "$PKGD/DESCRIPTION.orig"
  /bin/grep -c "^[+-][^+-]" "$PKGD/LOCAL_PATCH.diff" 2>/dev/null \
    && echo "── LOCAL_PATCH.diff written to $PKGD ──"
fi

# ── dependencies (WGCNA pulls most; preprocessCore/impute come from Bioconductor) ──
$RSCRIPT -e '
options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("igraph", "WGCNA", "ggplot2", "reshape2", "doBy", "WeightedROC",
          "pROC", "plotROC", "randomForest", "psych", "Rmisc", "ggforce",
          "ggpubr", "coin", "yaml")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bio <- c("preprocessCore", "impute")[
  !vapply(c("preprocessCore", "impute"), requireNamespace, logical(1), quietly = TRUE)]
if (length(bio)) BiocManager::install(bio, update = FALSE, ask = FALSE)'

# ── install + self-verifying toy test ──
$RSCRIPT -e "install.packages('$PKGD', repos = NULL, type = 'source')"

TMPD_R="$(mktemp -d)"
$RSCRIPT -e '
suppressPackageStartupMessages(library(NetMoss2))
D <- commandArgs(TRUE)[1]           # toy: 2 studies, gw1-g2 = 0.5 (n=30) / -0.2 (n=50)
dir.create(file.path(D, "data"), recursive = TRUE)
dir.create(file.path(D, "net"), recursive = TRUE)
mat17 <- function(m) matrix(format(c(m), digits = 17, trim = TRUE),
                            nrow(m), ncol(m), dimnames = dimnames(m))
stub <- function(n, f) {
  cat(paste(c("X", paste0("s", seq_len(n))), collapse = "\t"), "\n", file = f)
  cat(paste(c("x", rep("0", n)), collapse = "\t"), "\n", file = f, append = TRUE)
}
g <- c("g1", "g2", "g3", "g4")
mA <- diag(4); mA[1, 2] <- mA[2, 1] <- 0.5
mB <- diag(4); mB[1, 2] <- mB[2, 1] <- -0.2
dimnames(mA) <- dimnames(mB) <- list(g, g)
write.table(mat17(mA), file.path(D, "net", "A.txt"), sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = NA)
write.table(mat17(mB), file.path(D, "net", "B.txt"), sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = NA)
stub(30, file.path(D, "data", "A.txt")); stub(50, file.path(D, "data", "B.txt"))
r <- NetMoss2::getNetwork(file.path(D, "data"), file.path(D, "data"),
                          file.path(D, "net"), file.path(D, "net"))
stopifnot(isTRUE(all.equal(r[[1]], r[[2]], tolerance = 1e-12)))
p <- r[[1]][1, 2]
cat("toy pooled[1,2] =", format(p, digits = 10), "\n")
stopifnot(isTRUE(all.equal(p, 0.1017185322805387, tolerance = 1e-6)))
cat("TOY TEST PASS — patched getNetwork pools all studies\n")' "$TMPD_R"
rm -rf "$TMPD_R"

echo "install_netmoss2.sh: done — patched NetMoss2 at $PKGD"