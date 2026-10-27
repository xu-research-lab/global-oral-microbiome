#!/bin/bash
#SBATCH --job-name=genus_dcntm
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/genus_decontam_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/genus_decontam_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=12:00:00

# ============================================================
# Genus-level decontam with ALL control samples
# Strategies A (within-study) + B (pooled)
# Uses existing SILVA annotations, no re-annotation against full SILVA
# ============================================================

set -e
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv

# Make sure scipy is available
pip install scipy -q 2>/dev/null || true

BASE_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"
PIPE="${BASE_DIR}/pipeline"
TEMP_DIR="${BASE_DIR}/control_quality_viz/temp"
OUTDIR="${PIPE}/03_decontam_genus"
mkdir -p "${OUTDIR}" "${BASE_DIR}/logs"

echo "=== Genus-level Decontam (ALL controls) ==="
echo "Started: $(date)"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os, subprocess, json
from collections import defaultdict

PIPE = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline"
TEMP_DIR = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/control_quality_viz/temp"
OUTDIR = f"{PIPE}/03_decontam_genus"
os.makedirs(OUTDIR, exist_ok=True)

ORAL_DADA2 = "/beegfs/dongbiao/oral"
CONTROL_DADA2 = [
    "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/control_dada2",
    "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/control_dada2_extra",
]
ORAL_MD = "/beegfs/huangxiaochang/global_oral/metadata_v2_merge/metadata-v2-merge-oral.csv"
VSEARCH = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"
CONTROL_CLASSIFICATION = f"{TEMP_DIR}/../control_classification.csv"

# ============================================================
# 1. Load existing SILVA annotations
# ============================================================
print("=== Step 1: Loading existing SILVA annotations ===")
cols = ['query', 'target', 'id', 'alnlen', 'mism', 'gaps',
        'qlo', 'qhi', 'tlo', 'thi', 'ql', 'tl', 'qcov', 'tcov']
hits = pd.read_csv(f"{TEMP_DIR}/silva_hits.txt", sep="\t", names=cols)
hits['true_id'] = (hits['alnlen'] - hits['mism'] - hits['gaps']) / hits['ql'] * 100
hits_best = hits[hits['true_id'] >= 80].drop_duplicates(subset='query', keep='first')

# Parse SILVA -> genus
silva_genus = {}
with open("/beegfs/dongbiao/silva_138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta") as f:
    for line in f:
        if line.startswith(">"):
            parts = line.strip().lstrip(">").split(" ", 1)
            acc = parts[0]
            if len(parts) > 1:
                ranks = parts[1].split(";")
                genus = ranks[5].strip() if len(ranks) > 5 else "Unknown"
            else:
                genus = "Unknown"
            silva_genus[acc] = genus

asv_genus_map = {}
for _, row in hits_best.iterrows():
    genus = silva_genus.get(row['target'], "Unknown")
    if genus != "Unknown":
        asv_genus_map[row['query']] = genus
print(f"  Derep ASVs with genus: {len(asv_genus_map):,}")

# Load derep sequences
with open(f"{TEMP_DIR}/asv_mapping.json") as f:
    asv_mapping = json.load(f)
asv_seqs = asv_mapping['asv_id_to_seq']
print(f"  Derep sequences: {len(asv_seqs):,}")

# Build derep FASTA for vsearch
derep_fa = f"{OUTDIR}/derep_asv_db.fa"
print("  Creating derep DB FASTA...")
with open(derep_fa, "w") as f:
    for did, seq in asv_seqs.items():
        f.write(f">{did}\n{seq}\n")

# ============================================================
# 2. Sample metadata
# ============================================================
print("\n=== Step 2: Sample metadata ===")
ctrl_df = pd.read_csv(CONTROL_CLASSIFICATION)
ctrl_samples = set(ctrl_df["sample"])
ctrl_sample_bp = dict(zip(ctrl_df["sample"], ctrl_df["bioproject"]))
blank_bps = set(ctrl_sample_bp.values())
print(f"  ALL control samples: {len(ctrl_samples):,} across {len(blank_bps)} BPs")

oral_md = pd.read_csv(ORAL_MD, low_memory=False)
oral_md = oral_md[oral_md["bioproject"].isin(blank_bps)]

md_rows = []
oral_set = set()
for _, row in oral_md.iterrows():
    sid, bp = row["sample_id"], row["bioproject"]
    if bp not in blank_bps: continue
    oral_set.add(sid)
    md_rows.append({"sample_id": sid, "bioproject": bp, "is_blank": sid in ctrl_samples})
for _, row in ctrl_df.iterrows():
    sid, bp = row["sample"], row["bioproject"]
    if sid not in oral_set:
        md_rows.append({"sample_id": sid, "bioproject": bp, "is_blank": True})

md_df = pd.DataFrame(md_rows)
n_ctrl = md_df["is_blank"].sum()
n_oral = (~md_df["is_blank"]).sum()
print(f"  Samples: {len(md_df)} ({n_oral} oral, {n_ctrl} blank), {md_df['bioproject'].nunique()} BPs")

# ============================================================
# 3. Per-study: vsearch local ASVs -> derep DB (97%), aggregate genus
# ============================================================
print("\n=== Step 3: Per-study ASV annotation + genus aggregation ===")
genus_counts = defaultdict(lambda: defaultdict(int))

for idx, bp in enumerate(sorted(blank_bps)):
    # Load oral + control FASTA
    def load_fasta(fa_path):
        result = {}
        if not os.path.exists(fa_path): return result
        with open(fa_path) as f:
            h, s = None, []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if h: result[h] = "".join(s).upper()
                    h = line.lstrip(">"); s = []
                else: s.append(line)
            if h: result[h] = "".join(s).upper()
        return result

    oral_seqs = load_fasta(os.path.join(ORAL_DADA2, bp, "results", "ASVs.fa"))
    ctrl_seqs = {}
    for cdir in CONTROL_DADA2:
        cf = os.path.join(cdir, bp, "results", "ASVs.fa")
        if os.path.exists(cf): ctrl_seqs = load_fasta(cf); break

    local_all = {**oral_seqs, **ctrl_seqs}
    if not local_all: continue

    # vsearch at 97% against derep DB
    tmp_q = f"/tmp/gq_{bp}.fa"; tmp_h = f"/tmp/gh_{bp}.txt"
    with open(tmp_q, "w") as f:
        for lid, seq in local_all.items():
            f.write(f">{lid}\n{seq}\n")

    subprocess.run([VSEARCH, "--usearch_global", tmp_q, "--db", derep_fa,
                    "--id", "0.97", "--maxaccepts", "1", "--maxrejects", "0",
                    "--userout", tmp_h, "--userfields", "query+target",
                    "--threads", "16"], capture_output=True)

    local_to_genus = {}
    if os.path.exists(tmp_h):
        with open(tmp_h) as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 2:
                    local_to_genus[p[0]] = asv_genus_map.get(p[1], "Unknown")

    for lid in local_all:
        if lid not in local_to_genus:
            local_to_genus[lid] = "Unknown"

    # Aggregate genus counts
    oral_ct = os.path.join(ORAL_DADA2, bp, "results", "ASVs_counts.tsv")
    ctrl_ct = None
    for cdir in CONTROL_DADA2:
        cf = os.path.join(cdir, bp, "results", "ASVs_counts.tsv")
        if os.path.exists(cf): ctrl_ct = cf; break

    def agg(ct_path):
        if not ct_path or not os.path.exists(ct_path): return
        ct = pd.read_csv(ct_path, sep="\t", index_col=0)
        for lid in ct.index:
            genus = local_to_genus.get(lid, "Unknown")
            if genus == "Unknown": continue
            for col in ct.columns:
                val = ct.loc[lid, col]
                if isinstance(val, pd.Series): val = val.iloc[0]
                val = int(val)
                if val > 0: genus_counts[genus][col] += val

    agg(oral_ct); agg(ctrl_ct)

    for f in [tmp_q, tmp_h]:
        if os.path.exists(f): os.remove(f)

    if (idx+1) % 5 == 0 or idx == 0:
        n_cls = sum(1 for v in local_to_genus.values() if v != "Unknown")
        print(f"  [{idx+1:2d}] {bp}: {n_cls}/{len(local_all)} ASVs classified")

print(f"\n  Total genera with data: {len(genus_counts):,}")

# ============================================================
# 4. Genus × sample matrix
# ============================================================
print("\n=== Step 4: Building genus x sample matrix ===")
all_genera = sorted(genus_counts.keys())
md_set = set(md_df["sample_id"])
all_samps = sorted(set(s for sc in genus_counts.values() for s in sc if s in md_set))

matrix = pd.DataFrame(0, index=all_genera, columns=all_samps, dtype=int)
for genus, sc in genus_counts.items():
    for samp, cnt in sc.items():
        if samp in matrix.columns:
            matrix.loc[genus, samp] += cnt

print(f"  Matrix: {len(all_genera)} genera x {len(all_samps)} samples")

# Save intermediate
matrix.to_csv(f"{OUTDIR}/genus_counts_matrix.tsv", sep="\t")
md_df.to_csv(f"{OUTDIR}/genus_decontam_metadata.tsv", sep="\t", index=False)

# ============================================================
# 5. Decontam strategies
# ============================================================
print("\n=== Step 5: Decontam strategies A + B ===")
from scipy.stats import fisher_exact

pa = matrix.set_index("feature_id") > 0

# --- A: within-study ---
print("\n--- Strategy A: within-study Fisher ---")
a_rows = []
for bp in sorted(md_df["bioproject"].unique()):
    bp_md = md_df[md_df["bioproject"]==bp]
    ccols = [c for c in pa.columns if c in set(bp_md[bp_md["is_blank"]]["sample_id"])]
    ocols = [c for c in pa.columns if c in set(bp_md[~bp_md["is_blank"]]["sample_id"])]
    if len(ccols) < 1 or len(ocols) < 2: continue
    for g in pa.index:
        cp = (pa.loc[g, ccols] > 0).sum() if ccols else 0
        op = (pa.loc[g, ocols] > 0).sum() if ocols else 0
        if cp == 0 and op == 0: continue
        try:
            _, p = fisher_exact([[cp, len(ccols)-cp], [op, len(ocols)-op]], alternative='greater')
        except: p = 1.0
        a_rows.append({"feature_id": g, "bioproject": bp, "ctrl_present": cp, "oral_present": op,
                       "ctrl_total": len(ccols), "oral_total": len(ocols),
                       "p_value": p, "contaminant": p < 0.05})

a_df = pd.DataFrame(a_rows)
a_df.to_csv(f"{OUTDIR}/decontam_A_within_study.tsv", sep="\t", index=False)
agg_a = a_df.groupby("feature_id")["contaminant"].sum()
for k in [1, 2, 3, 5]:
    print(f"  Flagged in >= {k} BP: {(agg_a >= k).sum():,}")

catalog_a = agg_a[agg_a >= 2].sort_values(ascending=False)
if len(catalog_a) > 0:
    print(f"  Top cross-lab (>=2 BP): {', '.join(catalog_a.head(15).index)}")

# --- B: pooled ---
print("\n--- Strategy B: pooled Fisher ---")
ccols_all = [c for c in pa.columns if c in set(md_df[md_df["is_blank"]]["sample_id"])]
ocols_all = [c for c in pa.columns if c in set(md_df[~md_df["is_blank"]]["sample_id"])]

b_rows = []
for g in pa.index:
    cp = (pa.loc[g, ccols_all] > 0).sum()
    op = (pa.loc[g, ocols_all] > 0).sum()
    if cp == 0 and op == 0: continue
    try:
        _, p = fisher_exact([[cp, len(ccols_all)-cp], [op, len(ocols_all)-op]], alternative='greater')
    except: p = 1.0
    cprev = cp/len(ccols_all); oprev = op/len(ocols_all)
    ratio = oprev/cprev if cprev > 0 else float('inf')
    b_rows.append({"feature_id": g, "ctrl_prevalence": round(cprev,4),
                   "oral_prevalence": round(oprev,4), "ctrl_present": cp, "oral_present": op,
                   "ctrl_total": len(ccols_all), "oral_total": len(ocols_all),
                   "p_value": p, "contaminant_p05": p<0.05, "contaminant_p01": p<0.01,
                   "oral_over_ctrl_ratio": round(ratio,1) if ratio!=float('inf') else "inf"})

b_df = pd.DataFrame(b_rows).sort_values("p_value")
b_df.to_csv(f"{OUTDIR}/decontam_B_pooled.tsv", sep="\t", index=False)
print(f"  Genera tested: {len(b_rows):,}")
print(f"  p<0.05: {b_df['contaminant_p05'].sum():,}")
print(f"  p<0.01: {b_df['contaminant_p01'].sum():,}")

print(f"\n  Top 25 contaminants (pooled p<0.05):")
for _, r in b_df[b_df["contaminant_p05"]].head(25).iterrows():
    print(f"  {r['feature_id']:30s}  ctrl={r['ctrl_prevalence']:.4f}  oral={r['oral_prevalence']:.4f}  ratio={r['oral_over_ctrl_ratio']}  p={r['p_value']:.2e}")

# ============================================================
# 6. Cross-check with HOMD-uncaptured candidates
# ============================================================
print(f"\n{'='*60}")
print(f"Step 6: Cross-check contaminant genera vs HOMD-uncaptured candidates")
print(f"{'='*60}")

cand = pd.read_csv(f"{PIPE}/02_chimera/candidates_all_nonchimeric.tsv", sep="\t", low_memory=False)
cand_genus = cand["derep_id"].map(asv_genus_map)

contam_a2 = set(agg_a[agg_a >= 2].index)
contam_b_p05 = set(b_df[b_df["contaminant_p05"]]["feature_id"])
# Union of both strategies
all_contam_genera = contam_a2 | contam_b_p05

print(f"Contaminant genera from A (>=2 BP): {len(contam_a2)}")
print(f"Contaminant genera from B (p<0.05): {len(contam_b_p05)}")
print(f"Union: {len(all_contam_genera)}")

n_known = cand_genus.notna().sum()
n_contam = cand_genus.isin(all_contam_genera).sum()
n_clean_known = (cand_genus.notna() & ~cand_genus.isin(all_contam_genera)).sum()
n_unknown = cand_genus.isna().sum()

print(f"\nNon-chimeric candidates: {len(cand):,}")
print(f"  Genus annotated: {n_known:,}")
print(f"  Genus UNKNOWN: {n_unknown:,}")
print(f"  In contaminant genera: {n_contam:,}")
print(f"  Clean (annotated, not contaminant): {n_clean_known:,}")

# Clean + unknown = surviving
surviving = n_clean_known + n_unknown
print(f"\n  >>> Surviving candidates (clean + unknown): {surviving:,} / {len(cand):,} ({surviving/len(cand)*100:.1f}%)")

if n_contam > 0:
    print(f"\n  Candidates affected (by genus):")
    for g, n in cand_genus[cand_genus.isin(all_contam_genera)].value_counts().head(15).items():
        a_flag = "A" if g in contam_a2 else ""
        b_flag = "B" if g in contam_b_p05 else ""
        print(f"    {g:30s}: {n:5d} candidates [{a_flag}{b_flag}]")

# Also report: top genera among candidates that are NOT contaminants
print(f"\n  Top 15 genera among CLEAN candidates:")
for g, n in cand_genus[cand_genus.notna() & ~cand_genus.isin(all_contam_genera)].value_counts().head(15).items():
    print(f"    {g:30s}: {n:5d} candidates")

# ============================================================
# 7. Summary
# ============================================================
print(f"\n{'='*60}")
print(f"SUMMARY")
print(f"{'='*60}")
print(f"Strategy A (within-study, >=2 BP): {len(contam_a2)} contaminant genera")
print(f"Strategy B (pooled, p<0.05):      {len(contam_b_p05)} contaminant genera")
print(f"Combined contaminant genera:        {len(all_contam_genera)}")
print(f"")
print(f"HOMD-uncaptured non-chimeric candidates: {len(cand):,}")
print(f"  Surviving after genus-level decontam:  {surviving:,} ({surviving/len(cand)*100:.1f}%)")
print(f"  Removed as contaminants:                {n_contam:,} ({n_contam/len(cand)*100:.1f}%)")

print(f"\n=== Outputs: {OUTDIR}/ ===")
for f in sorted(os.listdir(OUTDIR)):
    sz = os.path.getsize(os.path.join(OUTDIR, f))
    print(f"  {f} ({sz:,} bytes)")
print(f"\nDone at $(date)")
PYEOF

echo ""
echo "=== Finished at $(date) ==="
ls -lh "${OUTDIR}/"
