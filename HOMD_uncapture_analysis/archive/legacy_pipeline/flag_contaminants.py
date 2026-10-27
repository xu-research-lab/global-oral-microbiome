#!/usr/bin/env python3
"""Integrate contamination evidence: decontam + curated reference + oral-vs-blank.

Follows pipeline plan Phase 2 Step 5-8:
- Cross-lab contaminant catalog from decontam results
- Curated literature blacklist (contaminant_reference_genera.tsv)
- Record per-ASV evidence (in-blank, blacklist, low-abundance-cosmopolitan)
- Output: labeled candidates with keep/flag/remove action
"""

import pandas as pd
import json, os
from collections import defaultdict

PIPE = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline"
TEMP_DIR = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/control_quality_viz/temp"
REF_TSV = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/script/05_homd_gg2_closed_reference/contaminant_reference_genera.tsv"
OUR_CATALOG = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/contaminant_catalog/contaminant_asv_catalog.tsv"
OUTDIR = f"{PIPE}/05_contamination_flag"
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# 1. Load reference table
# ============================================================
print("=== Loading curated reference table ===")
ref = pd.read_csv(REF_TSV, sep="\t", comment="#")
auto_remove = set(ref[ref['classification'] == 'auto_remove']['taxon'].str.strip().str.lower())
flag_review = set(ref[ref['classification'] == 'flag_review']['taxon'].str.strip().str.lower())
oral_core = set(ref[ref['classification'] == 'oral_core']['taxon'].str.strip().str.lower())

# Build reference lookup
ref_lookup = {}
for _, row in ref.iterrows():
    taxon = row['taxon'].strip().lower()
    ref_lookup[taxon] = {
        'classification': row['classification'],
        'category': row['category'],
        'rank': row['rank'],
        'note': row['note'] if 'note' in row else '',
        'source': row['source'] if 'source' in row else '',
    }
print(f"  auto_remove: {len(auto_remove)}, flag_review: {len(flag_review)}, oral_core: {len(oral_core)}")

# ============================================================
# 2. Load decontam results
# ============================================================
print("\n=== Loading decontam results ===")
a_df = pd.read_csv(f"{PIPE}/03_decontam_genus/decontam_A_within_study.tsv", sep="\t")
b_df = pd.read_csv(f"{PIPE}/03_decontam_genus/decontam_B_pooled.tsv", sep="\t")

agg_a = a_df.groupby("feature_id")["contaminant"].sum()
contam_a1 = set(agg_a[agg_a >= 1].index)  # flagged in ≥1 BP
contam_a2 = set(agg_a[agg_a >= 2].index)  # flagged in ≥2 BP

contam_b = set()
if "contaminant" in b_df.columns:
    contam_b = set(b_df[b_df["contaminant"] == True]["feature_id"])
else:
    # decontam R output uses TRUE/FALSE
    b_df["contaminant_bool"] = b_df["contaminant"].astype(str).str.upper() == "TRUE"
    contam_b = set(b_df[b_df["contaminant_bool"]]["feature_id"])

# Build genus-level evidence
decontam_evidence = defaultdict(dict)
for _, row in b_df.iterrows():
    g = row['feature_id']
    decontam_evidence[g.lower()]['b_pooled'] = bool(row.get('contaminant', False))
    decontam_evidence[g.lower()]['b_pval'] = row.get('p', 1.0)
    decontam_evidence[g.lower()]['b_ctrl_prev'] = row.get('ctrl_prevalence', 0)
    decontam_evidence[g.lower()]['b_oral_prev'] = row.get('oral_prevalence', 0)

for g, n_bp in agg_a.items():
    decontam_evidence[g.lower()]['a_n_bps'] = int(n_bp)
    decontam_evidence[g.lower()]['a_flagged'] = int(n_bp) >= 2

print(f"  Decontam A (≥1 BP): {len(contam_a1)}, A (≥2 BP): {len(contam_a2)}, B (pooled): {len(contam_b)}")

# ============================================================
# 3. Load our cross-lab contaminant catalog (14,685 ASVs)
# ============================================================
print("\n=== Loading cross-lab contaminant catalog ===")
catalog = pd.read_csv(OUR_CATALOG, sep="\t", low_memory=False)
# Build sequence -> catalog info
catalog_seq = set(catalog['sequence'].str.upper())
# Also get contaminant genus list from catalog
catalog_contam_genera = set(catalog[catalog['is_contaminant_genus']]['genus'].str.lower())
catalog_oral_genera = set(catalog[catalog['is_oral_genus']]['genus'].str.lower())
print(f"  Catalog ASVs: {len(catalog):,}")
print(f"  Catalog contaminant genera: {len(catalog_contam_genera)}")

# ============================================================
# 4. Annotate candidates with SILVA genus
# ============================================================
print("\n=== Annotating candidates ===")
with open(f"{TEMP_DIR}/asv_mapping.json") as f:
    asv_seqs = json.load(f)['asv_id_to_seq']

silva_genus = {}
with open("/beegfs/dongbiao/silva_138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta") as f:
    for line in f:
        if line.startswith(">"):
            parts = line.strip().lstrip(">").split(" ", 1)
            if len(parts) > 1:
                ranks = parts[1].split(";")
                genus = ranks[5].strip() if len(ranks) > 5 else ""
            else: genus = ""
            if genus and genus != "Unknown":
                silva_genus[parts[0]] = genus

cols = ['query','target','id','alnlen','mism','gaps','qlo','qhi','tlo','thi','ql','tl','qcov','tcov']
hits = pd.read_csv(f"{TEMP_DIR}/silva_hits.txt", sep="\t", names=cols)
hits['true_id'] = (hits['alnlen'] - hits['mism'] - hits['gaps']) / hits['ql'] * 100
hits_best = hits[hits['true_id'] >= 80].drop_duplicates(subset='query', keep='first')

seq_to_genus = {}
for _, row in hits_best.iterrows():
    seq = asv_seqs.get(row['query'])
    if seq and row['target'] in silva_genus:
        seq_to_genus[seq] = silva_genus[row['target']]

cand = pd.read_csv(f"{PIPE}/02_chimera/candidates_all_nonchimeric.tsv", sep="\t", low_memory=False)

# Annotate: genus from SILVA (for reference table matching)
cand['silva_genus'] = [seq_to_genus.get(seq, "Unknown") for seq in cand['sequence']]

# Also check if sequence is in our cross-lab catalog
cand['in_catalog'] = [seq.upper() in catalog_seq for seq in cand['sequence']]

print(f"  Annotated: {(cand['silva_genus']!='Unknown').sum():,}/{len(cand):,}")

# ============================================================
# 5. Screen each candidate: oral_core protect > auto_remove > flag_review
# ============================================================
print("\n=== Screening candidates ===")
actions = []; reasons = []; evidence_sources = []

for _, row in cand.iterrows():
    genus = str(row['silva_genus']).lower()
    in_cat = row['in_catalog']

    evidence = []
    action = "keep"
    reason = ""

    # ---- Layer 1: Oral-core protection (always keep) ----
    if genus in oral_core:
        action = "keep"
        reason = "oral_core_protected"
        evidence.append("oral_core")

    # ---- Layer 2: Check reference table ----
    elif genus in auto_remove:
        action = "remove"
        reason = f"reference_auto_remove:{genus}"
        info = ref_lookup.get(genus, {})
        evidence.append(f"ref_auto_remove({info.get('category','unknown')})")

    elif genus in flag_review:
        action = "flag"
        reason = f"reference_flag_review:{genus}"
        info = ref_lookup.get(genus, {})
        evidence.append(f"ref_flag_review({info.get('category','unknown')})")

    # ---- Layer 3: Check decontam statistics ----
    if genus in decontam_evidence:
        de = decontam_evidence[genus]
        if de.get('a_flagged'):
            evidence.append(f"decontam_A(>=2BP)")
            if action == "keep":
                action = "flag"
                reason = f"decontam_cross_lab:{genus}"
        elif de.get('b_pooled', False) and de.get('b_pval', 1.0) < 0.01:
            evidence.append(f"decontam_B(p<0.01)")
            if action == "keep":
                action = "flag"
                reason = f"decontam_pooled:{genus}"

    # ---- Layer 4: Cross-lab catalog match ----
    if in_cat:
        evidence.append("in_contaminant_catalog")
        if action == "keep":
            action = "flag"
            reason = "in_cross_lab_contaminant_catalog"

    # ---- Layer 5: Unknown genus ----
    if genus == "unknown" and action == "keep":
        reason = "unknown_genus_unmatched"
        evidence.append("unknown_genus")

    actions.append(action)
    reasons.append(reason)
    evidence_sources.append(";".join(evidence) if evidence else "none")

cand['screen_action'] = actions
cand['screen_reason'] = reasons
cand['evidence'] = evidence_sources

# ============================================================
# 6. Summary
# ============================================================
print(f"\n{'='*60}")
print(f"SCREENING RESULTS")
print(f"{'='*60}")

for a in ["keep", "flag", "remove"]:
    sub = cand[cand['screen_action'] == a]
    print(f"\n  {a}: {len(sub):,} ({len(sub)/len(cand)*100:.1f}%)")
    for r, n in sub['screen_reason'].value_counts().head(10).items():
        print(f"    - {r}: {n:,}")

# Export
cand.to_csv(f"{OUTDIR}/candidates_screened.tsv", sep="\t", index=False)

keep = cand[cand['screen_action'] == 'keep']
flag = cand[cand['screen_action'] == 'flag']
remove = cand[cand['screen_action'] == 'remove']

keep.to_csv(f"{OUTDIR}/candidates_kept.tsv", sep="\t", index=False)
flag.to_csv(f"{OUTDIR}/candidates_flagged.tsv", sep="\t", index=False)
remove.to_csv(f"{OUTDIR}/candidates_removed.tsv", sep="\t", index=False)

# FASTA for survivors (keep + flag go to Phase 4)
survivors = cand[cand['screen_action'].isin(['keep', 'flag'])]
with open(f"{OUTDIR}/candidates_for_phase4.fasta", "w") as fh:
    for _, row in survivors.iterrows():
        fh.write(f">{row['derep_id']} action={row['screen_action']}|{row['screen_reason']}\n")
        seq = row['sequence']
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")

# Evidence summary
print(f"\n{'='*60}")
print(f"PHASE 4 INPUT")
print(f"{'='*60}")
print(f"Total candidates:               {len(cand):,}")
print(f"  keep (oral_core / clean):     {len(keep):,}")
print(f"  flag (need review):           {len(flag):,}")
print(f"  remove (contaminant):         {len(remove):,}")
print(f"  >>> For Phase 4 (keep+flag):  {len(survivors):,}")
print(f"       -> {OUTDIR}/candidates_for_phase4.fasta")

print(f"\n=== Outputs: {OUTDIR}/ ===")
for f in sorted(os.listdir(OUTDIR)):
    sz = os.path.getsize(os.path.join(OUTDIR, f))
    print(f"  {f} ({sz:,} bytes)")
print("Done!")
