#!/bin/bash
#SBATCH --job-name=sra2fastq
#SBATCH --output=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/sra2fastq_%j.out
#SBATCH --error=/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/sra2fastq_%j.err
#SBATCH --nodelist=cu09
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=24:00:00

# ============================================================
# Convert control sample SRA files to FASTQ.gz
#
# Uses SRA Toolkit fastq-dump (--split-files --gzip)
# Excludes mock community samples (ERR5002560-ERR5002567)
# Output organized by BioProject for downstream DADA2 processing
# ============================================================

set -e

FASTQ_DUMP="/home/huangxc/software/sratoolkit.2.11.0-centos_linux64/bin/fastq-dump"
SRA_DIR="/beegfs/dongbiao/oral/control_samples/sra"
FASTQ_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/control_fastq"
CONTROL_IDS="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/control_sample_ids.tsv"
MAX_PARALLEL=16

MOCK_SAMPLES="ERR5002560 ERR5002561 ERR5002562 ERR5002563 ERR5002564 ERR5002565 ERR5002566 ERR5002567"

mkdir -p "$FASTQ_DIR" "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs"

echo "=== SRA to FASTQ conversion started at $(date) ==="
echo "FASTQ-DUMP: $FASTQ_DUMP"
echo "SRA dir: $SRA_DIR"
echo "FASTQ output: $FASTQ_DIR"
echo "Max parallel: $MAX_PARALLEL"
echo ""

# ---- Step 1: Build task list ----
TASK_LIST="$FASTQ_DIR/conversion_tasks.txt"

python3 -c "
import pandas as pd
import os
import json

sra_dir = '$SRA_DIR'
fastq_dir = '$FASTQ_DIR'
mock = set('$MOCK_SAMPLES'.split())
ctrl = pd.read_csv('$CONTROL_IDS', sep='\t')

# Filter: strong categories only, downloaded, not mock
strong = ctrl[ctrl['category'].isin(['extraction_blank', 'blank', 'ntc_no_template'])]
sra_available = set(os.listdir(sra_dir))
strong = strong[strong['sample_id'].isin(sra_available) & ~strong['sample_id'].isin(mock)]

tasks = []
for _, row in strong.iterrows():
    sample_id = row['sample_id']
    bioproject = row['bioproject']
    category = row['category']
    sra_path = os.path.join(sra_dir, sample_id, f'{sample_id}.sra')
    out_dir = os.path.join(fastq_dir, bioproject)

    if os.path.exists(sra_path):
        # Check if already converted
        r1_exists = os.path.exists(os.path.join(out_dir, f'{sample_id}_R1.fastq.gz'))
        r2_exists = os.path.exists(os.path.join(out_dir, f'{sample_id}_R2.fastq.gz'))
        se_exists = os.path.exists(os.path.join(out_dir, f'{sample_id}.fastq.gz'))
        already_done = r1_exists or se_exists

        tasks.append({
            'sample_id': sample_id,
            'bioproject': bioproject,
            'category': category,
            'sra_path': sra_path,
            'out_dir': out_dir,
            'already_done': already_done,
        })

pending = [t for t in tasks if not t['already_done']]
done = [t for t in tasks if t['already_done']]

print(f'Total tasks: {len(tasks)}')
print(f'Already done: {len(done)}')
print(f'Pending: {len(pending)}')

# Write task list for bash processing
with open('$TASK_LIST', 'w') as f:
    for t in pending:
        f.write(f\"{t['sample_id']}\t{t['bioproject']}\t{t['category']}\t{t['sra_path']}\t{t['out_dir']}\n\")

# Per-project summary
from collections import Counter
bp_counts = Counter(t['bioproject'] for t in pending)
print(f'\nPending by BioProject:')
for bp, count in bp_counts.most_common():
    print(f'  {bp}: {count} samples')

# Count with existing BIOM data
biom_bps = {'PRJNA634162': 11, 'PRJNA431300': 4, 'PRJEB42371': 2}
for bp in biom_bps:
    if bp in bp_counts:
        print(f'  ** {bp}: {bp_counts[bp]} samples (HAS existing BIOM data!)')
"

N_PENDING=$(wc -l < "$TASK_LIST")
echo ""
echo "Pending tasks: $N_PENDING"

if [ "$N_PENDING" -eq 0 ]; then
    echo "All done, nothing to convert!"
    exit 0
fi

# ---- Step 2: Run fastq-dump in parallel ----
echo ""
echo "[convert] Starting parallel conversion (max $MAX_PARALLEL concurrent)..."

process_one() {
    local sample_id="$1"
    local bioproject="$2"
    local category="$3"
    local sra_path="$4"
    local out_dir="$5"

    mkdir -p "$out_dir"

    # Check if already done
    if [ -f "$out_dir/${sample_id}_R1.fastq.gz" ] || [ -f "$out_dir/${sample_id}.fastq.gz" ]; then
        echo "[skip] $sample_id ($bioproject) - already converted"
        return 0
    fi

    echo "[start] $sample_id ($bioproject / $category)"

    # Run fastq-dump --split-files: handles both PE and SE
    "$FASTQ_DUMP" \
        --split-files \
        --gzip \
        --skip-technical \
        --outdir "$out_dir" \
        "$sra_path" 2>&1 | tail -1

    # Detect output: paired-end -> _1.fastq.gz + _2.fastq.gz; single-end -> .fastq.gz
    # --split-files may also produce _pass_1 / _pass_2 with some SRA versions
    if [ -f "$out_dir/${sample_id}_1.fastq.gz" ] && [ -f "$out_dir/${sample_id}_2.fastq.gz" ]; then
        mv "$out_dir/${sample_id}_1.fastq.gz" "$out_dir/${sample_id}_R1.fastq.gz"
        mv "$out_dir/${sample_id}_2.fastq.gz" "$out_dir/${sample_id}_R2.fastq.gz"
        echo "[done] $sample_id - paired-end"
    elif [ -f "$out_dir/${sample_id}_1.fastq.gz" ]; then
        mv "$out_dir/${sample_id}_1.fastq.gz" "$out_dir/${sample_id}_R1.fastq.gz"
        echo "[done] $sample_id - single-end (or unpaired read 1)"
    elif [ -f "$out_dir/${sample_id}.fastq.gz" ]; then
        echo "[done] $sample_id - single-end"
    else
        # Try glob match for any naming variation
        FOUND=$(ls "$out_dir/${sample_id}"*.fastq.gz 2>/dev/null | head -5)
        if [ -n "$FOUND" ]; then
            echo "[done] $sample_id - found: $(echo $FOUND | tr '\n' ' ')"
        else
            echo "[FAIL] $sample_id - no output files found"
            return 1
        fi
    fi
}

export -f process_one
export FASTQ_DUMP="$FASTQ_DUMP"

# Process in parallel using xargs
cat "$TASK_LIST" | xargs -P $MAX_PARALLEL -I {} bash -c '
    IFS=$'"'"'\t'"'"' read -r sample_id bioproject category sra_path out_dir <<< "$1"
    process_one "$sample_id" "$bioproject" "$category" "$sra_path" "$out_dir"
' _ {}

# ---- Step 3: Summary ----
echo ""
echo "[summary] Conversion completed at $(date)"
echo ""

python3 -c "
import os
import pandas as pd

fastq_dir = '$FASTQ_DIR'
sra_dir = '$SRA_DIR'
mock = set('$MOCK_SAMPLES'.split())
ctrl = pd.read_csv('$CONTROL_IDS', sep='\t')

strong = ctrl[ctrl['category'].isin(['extraction_blank', 'blank', 'ntc_no_template'])]
sra_available = set(os.listdir(sra_dir))
strong = strong[strong['sample_id'].isin(sra_available) & ~strong['sample_id'].isin(mock)]

success = 0
failed = 0
total_bp = 0
for bp in strong['bioproject'].unique():
    bp_dir = os.path.join(fastq_dir, bp)
    if not os.path.isdir(bp_dir):
        continue
    fastq_files = [f for f in os.listdir(bp_dir) if f.endswith('.fastq.gz')]
    sample_ids = set()
    for f in fastq_files:
        # Extract sample ID: remove _R1/_R2 and .fastq.gz
        sid = f.replace('_R1.fastq.gz','').replace('_R2.fastq.gz','').replace('.fastq.gz','')
        sample_ids.add(sid)
    success += len(sample_ids)
    total_bp += 1
    expected = len(strong[strong['bioproject'] == bp])
    missed = expected - len(sample_ids)
    print(f'{bp}: {len(sample_ids)}/{expected} samples converted (PE files: {sum(1 for f in fastq_files if \"_R\" in f)})')
    if missed > 0:
        print(f'  WARNING: {missed} samples missing!')

print(f'\nTotal: {success} samples converted across {total_bp} BioProjects')
print(f'FASTQ output: {fastq_dir}')
"

echo ""
echo "=== All done at $(date) ==="
