#!/usr/bin/env python3
"""
Parse SEPP placements JSON, extract per-query best placement,
cross-reference with HOMD, and output placements TSV for step 06.

Creates: pipeline/05_sepp/sepp_placements_parsed.tsv
Columns: query_id, edge_num, pendant_length, like_weight_ratio, likelihood,
         nearest_ref (from VSEARCH), nearest_ref_is_homd, cluster (tree edge group)
"""

import json, os, sys
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

PIPE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline")
SEPP_DIR = PIPE / "05_sepp"

def log(s): print(f"[parse-sepp] {s}", flush=True)

# ============================================================
# 1. Load SEPP placements JSON
# ============================================================
log("Loading SEPP placements.json...")
with open(SEPP_DIR / "placements_export/placements.json") as f:
    data = json.load(f)

log(f"  {len(data['placements']):,} placement entries")
log(f"  Fields: {data['fields']}")

# ============================================================
# 2. Extract per-query best placement
# ============================================================
log("Extracting per-query best placements...")

rows = []
missing_best = 0

for entry in data['placements']:
    # Query name
    nm = entry.get('nm', [])
    if not nm or not nm[0]:
        continue
    query_id = nm[0][0]  # e.g., "DEREP_926332"

    # Placement positions
    positions = entry.get('p', [])
    if not positions:
        continue

    # Best placement: highest like_weight_ratio (field index 2)
    best_idx = np.argmax([p[2] for p in positions])  # index 2 = like_weight_ratio
    best = positions[best_idx]

    distal_length, edge_num, lwr, likelihood, pendant_length = best

    rows.append({
        'query_id': query_id,
        'edge_num': int(edge_num),
        'distal_length': float(distal_length),
        'like_weight_ratio': float(lwr),
        'likelihood': float(likelihood),
        'pendant_length': float(pendant_length),
    })

df = pd.DataFrame(rows)
log(f"  Parsed {len(df):,} placements")
log(f"  Unique queries: {df['query_id'].nunique():,}")
log(f"  Median pendant_length: {df['pendant_length'].median():.4f}")
log(f"  Median like_weight_ratio: {df['like_weight_ratio'].median():.4f}")

# ============================================================
# 3. Load VSEARCH-based nearest_ref for cross-reference
# ============================================================
log("Loading VSEARCH nearest-ref mapping...")
vsearch = pd.read_csv(SEPP_DIR / "vsearch_placements.tsv", sep='\t', low_memory=False)
# Drop duplicates - keep best identity per query
vsearch_dedup = vsearch.sort_values('identity', ascending=False).drop_duplicates(subset='query_id', keep='first')
vsearch_map = vsearch_dedup.set_index('query_id')[['nearest_ref', 'nearest_ref_is_homd', 'identity']].to_dict('index')
log(f"  VSEARCH mapping: {len(vsearch_map):,} queries")

# Add nearest_ref info from VSEARCH (for HOMD check and grouping option)
df['nearest_ref'] = df['query_id'].map(lambda x: vsearch_map.get(x, {}).get('nearest_ref', ''))
df['nearest_ref_is_homd'] = df['query_id'].map(lambda x: vsearch_map.get(x, {}).get('nearest_ref_is_homd', None))
df['vsearch_identity'] = df['query_id'].map(lambda x: vsearch_map.get(x, {}).get('identity', np.nan))

n_with_ref = df['nearest_ref'].notna().sum()
n_homd = df['nearest_ref_is_homd'].sum()
log(f"  With nearest_ref: {n_with_ref:,} / {len(df):,}")
log(f"  nearest_ref IS HOMD: {n_homd:,}")

# ============================================================
# 4. Create cluster IDs: group by tree position (edge_num)
#    This is the KEY advantage of SEPP: phylotypes based on
#    actual tree position, not sequence identity
# ============================================================
log("Creating tree-position-based clusters...")
# Cluster by edge_num (same edge = co-placed = likely same organism)
# For multi-region unification: V4 and V1-V3 ASVs of same organism
# will be placed on same edge
df['cluster'] = 'EDGE_' + df['edge_num'].astype(str)

# Count clusters
n_clusters = df['cluster'].nunique()
log(f"  Clusters (unique edges): {n_clusters:,}")
log(f"  Median ASVs per cluster: {df.groupby('cluster').size().median():.0f}")

# Multi-ASV clusters (potential multi-region unification)
multi_asv = (df.groupby('cluster').size() > 1).sum()
log(f"  Clusters with >1 ASV: {multi_asv:,} (multi-region unification)")

# ============================================================
# 5. Save
# ============================================================
out_path = SEPP_DIR / "sepp_placements_parsed.tsv"
df.to_csv(out_path, sep='\t', index=False)
log(f"\nSaved: {out_path} ({os.path.getsize(out_path):,} bytes)")
log(f"  {len(df):,} rows x {len(df.columns)} columns")

# Summary
log(f"\n=== SEPP PLACEMENT SUMMARY ===")
log(f"Total queries placed: {len(df):,}")
log(f"  HOMD-gap (nearest ref NOT HOMD): {(~df['nearest_ref_is_homd'].fillna(True)).sum():,}")
log(f"  HOMD-proximal: {df['nearest_ref_is_homd'].fillna(False).sum():,}")
log(f"  Unique tree-edge clusters: {n_clusters:,}")
log(f"  Multi-ASV clusters: {multi_asv:,}")

# Pendant length distribution
for p in [10, 25, 50, 75, 90, 95]:
    val = np.percentile(df['pendant_length'], p)
    log(f"  pendant_length p{p}: {val:.4f}")

log("Done!")
