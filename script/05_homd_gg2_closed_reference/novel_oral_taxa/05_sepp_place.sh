#!/usr/bin/env bash
# Phase 4: place non-chimeric candidate ASVs onto the Greengenes2 backbone via SEPP.
#
# GG2 is SEPP-native and unifies amplicon regions on one tree (a V4 and a V1-V3 ASV of the
# same organism cannot be sequence-clustered but CAN be co-placed). HOMD taxa are overlaid
# on the GG2 tree downstream (06) as labels, so a candidate placing in a HOMD-free
# neighborhood is a candidate novel oral taxon.
#
# Requires QIIME2 (q2-fragment-insertion) and the Greengenes2 SEPP reference.
# Cluster step; not executed locally.
set -euo pipefail

CANDIDATES="${1:?usage: 05_sepp_place.sh <nonchimeric_candidates.fasta> <gg2_sepp_reference.qza> <outdir> [threads]}"
GG2_SEPP_REF="${2:?path to greengenes2 SEPP reference (.qza)}"
OUTDIR="${3:?output directory}"
THREADS="${4:-16}"

mkdir -p "$OUTDIR"

# 1. FASTA -> QIIME2 artifact
qiime tools import \
  --input-path "$CANDIDATES" \
  --output-path "$OUTDIR/candidates.qza" \
  --type 'FeatureData[Sequence]'

# 2. Fragment insertion (SEPP) onto the GG2 backbone
qiime fragment-insertion sepp \
  --i-representative-sequences "$OUTDIR/candidates.qza" \
  --i-reference-database "$GG2_SEPP_REF" \
  --p-threads "$THREADS" \
  --o-tree "$OUTDIR/insertion_tree.qza" \
  --o-placements "$OUTDIR/insertion_placements.qza"

# 3. Export tree + placements for 06_define_phylotypes.py
qiime tools export --input-path "$OUTDIR/insertion_tree.qza" --output-path "$OUTDIR/tree_export"
qiime tools export --input-path "$OUTDIR/insertion_placements.qza" --output-path "$OUTDIR/placements_export"

echo "SEPP placement done."
echo "  tree:       $OUTDIR/tree_export/tree.nwk"
echo "  placements: $OUTDIR/placements_export/placements.json"
echo "Next: derive per-query nearest-reference + phylotype clusters (e.g. gappa examine"
echo "graft/assign on placements.json), then run 06_define_phylotypes.py."
