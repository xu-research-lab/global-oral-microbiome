#!/bin/bash
#SBATCH --job-name=contam_cat
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/contam_cat_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/contam_cat_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=8:00:00

# ============================================================
# Phase 2 Step 5: Build cross-lab contaminant ASV catalog
# from 1,226 genuine control samples
# ============================================================

set -e
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv

BASE_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/06_control_catalog"
TEMP_DIR="${BASE_DIR}/work/control_quality_viz/temp"
CATALOG_DIR="${BASE_DIR}/work/contaminant_catalog"
mkdir -p "${CATALOG_DIR}"

echo "=== Building Cross-Lab Contaminant ASV Catalog ==="
echo "Started: $(date)"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import json, os, gzip
from collections import defaultdict, Counter

BASE_DIR = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/06_control_catalog"
TEMP_DIR = f"{BASE_DIR}/work/control_quality_viz/temp"
CATALOG_DIR = f"{BASE_DIR}/work/contaminant_catalog"

# ============================================================
# Load kept control sample list
# ============================================================
keep_list = pd.read_csv(f"{BASE_DIR}/work/control_quality_viz/genuine_controls_keep.csv")
keep_samples = set(keep_list['sample'].values)
keep_bps = set(keep_list['bioproject'].values)
print(f"Kept control samples: {len(keep_samples)} across {len(keep_bps)} BioProjects")

# Build sample -> bioproject + category map
sample_info = {}
for _, row in keep_list.iterrows():
    sample_info[row['sample']] = {
        'bioproject': row['bioproject'],
        'category': row['category'],
    }

# ============================================================
# Load SILVA taxonomy for ASV annotation
# ============================================================
print("Loading SILVA taxonomy...")
cols = ['query', 'target', 'id', 'alnlen', 'mism', 'gaps',
        'qlo', 'qhi', 'tlo', 'thi', 'ql', 'tl', 'qcov', 'tcov']
hits = pd.read_csv(f"{TEMP_DIR}/silva_hits.txt", sep="\t", names=cols)
hits['true_id'] = (hits['alnlen'] - hits['mism'] - hits['gaps']) / hits['ql'] * 100

silva_tax = {}
with open("/beegfs/dongbiao/silva_138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta") as f:
    for line in f:
        if line.startswith(">"):
            full_line = line.strip().lstrip(">")
            space_idx = full_line.find(" ")
            if space_idx > 0:
                silva_tax[full_line[:space_idx]] = full_line[space_idx+1:]

def parse_taxonomy(target_id):
    tax_str = silva_tax.get(target_id, "Unknown")
    ranks = [r.strip() for r in tax_str.split(";")]
    result = {}
    for i, rn in enumerate(['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
        result[rn] = ranks[i] if i < len(ranks) else 'Unknown'
    return result

# Build ASV -> taxonomy
asv_taxonomy = {}
for _, row in hits.iterrows():
    if row['true_id'] >= 80 and row['query'] not in asv_taxonomy:
        asv_taxonomy[row['query']] = parse_taxonomy(row['target'])

# Load ASV sequences
with open(f"{TEMP_DIR}/asv_mapping.json") as f:
    asv_mapping = json.load(f)
    asv_seq_map = asv_mapping['asv_id_to_seq']

print(f"  ASVs with taxonomy: {len(asv_taxonomy)}")
print(f"  Total ASVs: {len(asv_seq_map)}")

# Pre-build reverse index: sequence -> global ASV ID
seq_to_global_id = {seq: gid for gid, seq in asv_seq_map.items()}
print(f"  Built reverse sequence index for fast lookup")

# ============================================================
# Define contaminant indicator genera
# ============================================================
# Known reagent/lab contaminants (Salter 2014, Eisenhofer 2019)
CONTAMINANT_GENERA = {
    'Pseudomonas', 'Acinetobacter', 'Stenotrophomonas', 'Ralstonia',
    'Bradyrhizobium', 'Sphingomonas', 'Methylobacterium',
    'Cutibacterium', 'Propionibacterium', 'Staphylococcus',
    'Delftia', 'Acidovorax', 'Novosphingobium', 'Brevundimonas',
    'Aquabacterium', 'Curvibacter', 'Herbaspirillum', 'Massilia',
    'Pelomonas', 'Caulobacter', 'Chryseobacterium', 'Flavobacterium',
    'Sediminibacterium', 'Undibacterium', 'Enhydrobacter',
    'Micrococcus', 'Kocuria', 'Janthinobacterium', 'Comamonas',
    'Hydrogenophaga', 'Limnobacter', 'Variovorax', 'Sphingobium',
    'Paenibacillus', 'Bacillus', 'Enterobacter', 'Escherichia-Shigella',
    'Eperythrozoon', 'Mycoplasma',
}

ORAL_GENERA = {
    'Streptococcus', 'Prevotella', 'Veillonella', 'Haemophilus',
    'Neisseria', 'Rothia', 'Fusobacterium', 'Leptotrichia',
    'Porphyromonas', 'Actinomyces', 'Oribacterium', 'Capnocytophaga',
    'Parvimonas', 'Alloprevotella', 'Gemella', 'Granulicatella',
    'Selenomonas', 'Campylobacter', 'Tannerella',
    'Lachnoanaerobaculum', 'Peptostreptococcus', 'Schaalia',
    'Kingella', 'Eikenella', 'Aggregatibacter', 'Filifactor',
    'Cardiobacterium', 'Abiotrophia', 'Treponema',
}

# ============================================================
# Scan DADA2 results for kept samples, collect ASV statistics
# ============================================================
print("\nScanning DADA2 results for kept control samples...")

# ASV-level statistics across controls
asv_prevalence_samples = Counter()    # ASV ID -> n control samples
asv_prevalence_bioprojects = defaultdict(set)  # ASV ID -> set of BioProjects
asv_total_reads = Counter()           # ASV ID -> total reads in controls
asv_max_abundance = {}                # ASV ID -> max relative abundance in any control
asv_category_distribution = defaultdict(lambda: defaultdict(int))  # ASV -> category -> count

result_dirs = [
    f"{BASE_DIR}/work/control_dada2",
    f"{BASE_DIR}/work/control_dada2_extra"
]

total_kept_scanned = 0

for result_dir in result_dirs:
    if not os.path.isdir(result_dir):
        continue
    for bioproject in sorted(os.listdir(result_dir)):
        proj_path = os.path.join(result_dir, bioproject)
        if not os.path.isdir(proj_path):
            continue

        counts_file = os.path.join(proj_path, "results", "ASVs_counts.tsv")
        fasta_file = os.path.join(proj_path, "results", "ASVs.fa")
        if not os.path.exists(counts_file) or not os.path.exists(fasta_file):
            continue

        # Parse FASTA (local ASV header -> sequence)
        bioproject_asv_map = {}
        with open(fasta_file) as f:
            header = None; seq_parts = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header is not None:
                        bioproject_asv_map[header.lstrip(">")] = "".join(seq_parts).upper()
                    header = line; seq_parts = []
                else:
                    seq_parts.append(line)
            if header is not None:
                bioproject_asv_map[header.lstrip(">")] = "".join(seq_parts).upper()

        # Map local ASV sequence -> global ASV ID via pre-built reverse index
        seq_to_global = {}
        for local_id, seq in bioproject_asv_map.items():
            gid = seq_to_global_id.get(seq)
            if gid:
                seq_to_global[seq] = gid

        # Parse counts table
        counts_df = pd.read_csv(counts_file, sep="\t", index_col=0)

        for sample_name in counts_df.columns:
            # Check if this sample is in the keep list
            if sample_name not in keep_samples:
                continue

            total_kept_scanned += 1
            cat = sample_info[sample_name]['category']
            read_counts = counts_df[sample_name]
            nonzero = read_counts[read_counts > 0]
            total_reads = nonzero.sum()

            if total_reads == 0:
                continue

            for local_asv, count in nonzero.items():
                seq = bioproject_asv_map.get(local_asv)
                if seq is None:
                    continue
                global_id = seq_to_global.get(seq, local_asv)

                asv_prevalence_samples[global_id] += 1
                asv_prevalence_bioprojects[global_id].add(bioproject)
                asv_total_reads[global_id] += int(count)
                asv_category_distribution[global_id][cat] += 1

                # Track relative abundance in this sample
                rel_abund = count / total_reads
                if global_id not in asv_max_abundance or rel_abund > asv_max_abundance[global_id]:
                    asv_max_abundance[global_id] = rel_abund

print(f"Scanned {total_kept_scanned} kept samples out of {len(keep_samples)} expected")

# ============================================================
# Build catalog
# ============================================================
print("\nBuilding contaminant catalog...")

catalog_rows = []
n_oral_in_catalog = 0
n_contaminant_in_catalog = 0

for asv_id, n_samples in asv_prevalence_samples.items():
    n_bps = len(asv_prevalence_bioprojects[asv_id])
    total_reads = asv_total_reads[asv_id]
    max_abund = asv_max_abundance.get(asv_id, 0)
    seq = asv_seq_map.get(asv_id, '')

    # Taxonomy
    tax = asv_taxonomy.get(asv_id, {})
    genus = tax.get('genus', 'Unclassified')
    phylum = tax.get('phylum', 'Unclassified')
    domain = tax.get('domain', 'Unclassified')
    family = tax.get('family', 'Unclassified')
    order = tax.get('order', 'Unclassified')

    # Category distribution
    cat_dist = asv_category_distribution[asv_id]
    cat_str = ','.join([f"{c}:{n}" for c, n in sorted(cat_dist.items())])

    # Contamination flags
    is_contaminant_genus = genus in CONTAMINANT_GENERA
    is_oral_genus = genus in ORAL_GENERA
    if is_contaminant_genus:
        n_contaminant_in_catalog += 1
    if is_oral_genus:
        n_oral_in_catalog += 1

    # Confidence score: ASVs appearing in multiple independent labs are high-confidence
    if n_bps >= 3 and is_contaminant_genus:
        confidence = 'high'
    elif n_bps >= 2:
        confidence = 'medium'
    else:
        confidence = 'low'

    catalog_rows.append({
        'asv_id': asv_id,
        'sequence': seq,
        'seq_length': len(seq),
        'domain': domain,
        'phylum': phylum,
        'class': tax.get('class', 'Unclassified'),
        'order': order,
        'family': family,
        'genus': genus,
        'species': tax.get('species', 'Unclassified'),
        'n_control_samples': n_samples,
        'n_bioprojects': n_bps,
        'total_reads_in_controls': total_reads,
        'max_rel_abundance': round(max_abund * 100, 4),
        'is_contaminant_genus': is_contaminant_genus,
        'is_oral_genus': is_oral_genus,
        'contamination_confidence': confidence,
        'category_distribution': cat_str,
        'bioprojects': ','.join(sorted(asv_prevalence_bioprojects[asv_id])),
    })

catalog_df = pd.DataFrame(catalog_rows)
catalog_df = catalog_df.sort_values(['n_bioprojects', 'n_control_samples', 'total_reads_in_controls'],
                                     ascending=[False, False, False])

# ============================================================
# Summary statistics
# ============================================================
print("\n" + "="*70)
print("CONTAMINANT CATALOG SUMMARY")
print("="*70)
print(f"Total ASVs in catalog: {len(catalog_df):,}")
print(f"  ASVs from contaminant genera: {n_contaminant_in_catalog:,} ({n_contaminant_in_catalog/len(catalog_df)*100:.1f}%)")
print(f"  ASVs from oral genera: {n_oral_in_catalog:,} ({n_oral_in_catalog/len(catalog_df)*100:.1f}%)")

print(f"\nBy confidence level:")
for conf in ['high', 'medium', 'low']:
    n = (catalog_df['contamination_confidence'] == conf).sum()
    print(f"  {conf}: {n:,} ({n/len(catalog_df)*100:.1f}%)")

print(f"\nPrevalence distribution:")
for threshold in [1, 2, 3, 5, 10, 20, 50]:
    n = (catalog_df['n_control_samples'] >= threshold).sum()
    print(f"  ≥{threshold} control samples: {n:,} ASVs ({n/len(catalog_df)*100:.1f}%)")

print(f"\nCross-BioProject distribution:")
for threshold in [1, 2, 3, 5, 10]:
    n = (catalog_df['n_bioprojects'] >= threshold).sum()
    print(f"  ≥{threshold} BioProjects: {n:,} ASVs ({n/len(catalog_df)*100:.1f}%)")

print(f"\nTop 20 contaminant genera in catalog:")
for genus in sorted(CONTAMINANT_GENERA):
    sub = catalog_df[catalog_df['genus'] == genus]
    if len(sub) > 0:
        total_reads = sub['total_reads_in_controls'].sum()
        print(f"  {genus:30s}: {len(sub):6,d} ASVs, {total_reads:>12,d} reads, "
              f"in {sub['n_bioprojects'].max()} BPs max")

print(f"\nTop 20 oral genera in catalog (potential cross-contamination):")
oral_in_cat = catalog_df[catalog_df['is_oral_genus']]
oral_genus_summary = oral_in_cat.groupby('genus').agg(
    n_asvs=('asv_id', 'count'),
    total_reads=('total_reads_in_controls', 'sum'),
    max_bps=('n_bioprojects', 'max'),
    median_samples=('n_control_samples', 'median'),
).sort_values('total_reads', ascending=False)
for genus, row in oral_genus_summary.head(20).iterrows():
    print(f"  {genus:30s}: {int(row['n_asvs']):5,d} ASVs, {int(row['total_reads']):>10,d} reads, "
          f"max {int(row['max_bps'])} BPs, median {int(row['median_samples'])} samples")

# ============================================================
# Export catalog
# ============================================================
print("\n[Export] Writing catalog...")

# Full catalog
catalog_df.to_csv(f"{CATALOG_DIR}/contaminant_asv_catalog.tsv", sep="\t", index=False)
print(f"  -> contaminant_asv_catalog.tsv ({len(catalog_df):,} ASVs)")

# High-confidence subset (cross-lab, likely contaminant genera)
high_conf = catalog_df[catalog_df['contamination_confidence'] == 'high']
high_conf.to_csv(f"{CATALOG_DIR}/contaminant_catalog_high_confidence.tsv", sep="\t", index=False)
print(f"  -> contaminant_catalog_high_confidence.tsv ({len(high_conf):,} ASVs)")

# Oral ASVs found in controls (for awareness)
oral_in_controls = catalog_df[catalog_df['is_oral_genus']].sort_values(
    ['n_bioprojects', 'total_reads_in_controls'], ascending=[False, False])
oral_in_controls.to_csv(f"{CATALOG_DIR}/oral_asvs_in_controls.tsv", sep="\t", index=False)
print(f"  -> oral_asvs_in_controls.tsv ({len(oral_in_controls):,} ASVs)")

# FASTA for matching against HOMD-uncaptured candidates
with open(f"{CATALOG_DIR}/contaminant_catalog.fasta", "w") as f:
    for _, row in catalog_df.iterrows():
        f.write(f">{row['asv_id']} genus={row['genus']}|bps={row['n_bioprojects']}|samples={row['n_control_samples']}|conf={row['contamination_confidence']}\n")
        f.write(f"{row['sequence']}\n")
print(f"  -> contaminant_catalog.fasta ({len(catalog_df):,} sequences)")

# Summary statistics
summary = {
    'total_asvs': len(catalog_df),
    'n_kept_control_samples_used': total_kept_scanned,
    'n_bioprojects': len(keep_bps),
    'n_contaminant_genus_asvs': n_contaminant_in_catalog,
    'n_oral_genus_asvs': n_oral_in_catalog,
    'high_confidence': int((catalog_df['contamination_confidence'] == 'high').sum()),
    'medium_confidence': int((catalog_df['contamination_confidence'] == 'medium').sum()),
    'low_confidence': int((catalog_df['contamination_confidence'] == 'low').sum()),
    'asvs_in_1_sample': int((catalog_df['n_control_samples'] == 1).sum()),
    'asvs_in_2plus_samples': int((catalog_df['n_control_samples'] >= 2).sum()),
    'asvs_in_3plus_bps': int((catalog_df['n_bioprojects'] >= 3).sum()),
}
with open(f"{CATALOG_DIR}/catalog_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"  -> catalog_summary.json")

print(f"\n=== Catalog built at $(date) ===")
print(f"Output: {CATALOG_DIR}/")
for f in sorted(os.listdir(CATALOG_DIR)):
    size = os.path.getsize(os.path.join(CATALOG_DIR, f))
    print(f"  {f} ({size:,} bytes)")
PYEOF

echo ""
echo "=== Done at $(date) ==="
ls -lh "${CATALOG_DIR}/"
