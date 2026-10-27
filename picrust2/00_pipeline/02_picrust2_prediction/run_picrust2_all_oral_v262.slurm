#!/bin/bash
#SBATCH -J picrust2
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --exclude=cu01,cu02
#SBATCH -o logs/picrust2/picrust2.%j.out
#SBATCH -e logs/picrust2/picrust2.%j.err

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/home/zhangcaiyi/project/global_oral_picrust/close_ref4_v2}"

FASTA="${FASTA:-$PROJECT_ROOT/01.filter_input/oral_qc__dedup.nonzero.rep_seqs.fna}"
BIOM="${BIOM:-$PROJECT_ROOT/01.filter_input/oral_qc__dedup.nonzero.feature_table.biom}"
OUTDIR="${OUTDIR:-$PROJECT_ROOT/02.picrust2_out_v262}"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
OVERWRITE="${OVERWRITE:-0}"

cd "$PROJECT_ROOT"
mkdir -p logs/picrust2

echo "Job started at: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Node: $(hostname)"
echo "Workdir: $(pwd)"
echo "PROJECT_ROOT: $PROJECT_ROOT"
echo "FASTA: $FASTA"
echo "BIOM: $BIOM"
echo "OUTDIR: $OUTDIR"
echo "THREADS: $THREADS"
echo "OVERWRITE: $OVERWRITE"

if [ ! -s "$FASTA" ]; then
  echo "ERROR: FASTA not found or empty: $FASTA"
  exit 1
fi

if [ ! -s "$BIOM" ]; then
  echo "ERROR: BIOM not found or empty: $BIOM"
  exit 1
fi

# Conda activation can conflict with `set -u`, so disable nounset temporarily.
set +u

module load miniconda/4.9.2 || true

PICRUST2_ENV="/home/software/anaconda3/envs/picrust2"

if [ -f /home/software/anaconda3/etc/profile.d/conda.sh ]; then
  source /home/software/anaconda3/etc/profile.d/conda.sh
  conda activate "$PICRUST2_ENV"
else
  source activate "$PICRUST2_ENV"
fi

# Ensure env bin is first on PATH (module load miniconda can leave base ahead).
export PATH="$PICRUST2_ENV/bin:$PATH"
hash -r 2>/dev/null || true

set -u

echo "CONDA_PREFIX: ${CONDA_PREFIX:-NA}"
echo "Python: $(which python)"
python --version
echo "place_seqs.py: $(which place_seqs.py)"

echo "PICRUSt2:"
which picrust2_pipeline.py
picrust2_pipeline.py -h | head -n 10

echo "Conda packages:"
conda list | grep -E "^picrust2|^q2-picrust2|^biom-format" || true

# Fail fast if wrong place_seqs (base has a 2.5.2 stub that breaks 2.6.2).
if ! which place_seqs.py | grep -q "$PICRUST2_ENV/bin"; then
  echo "ERROR: place_seqs.py is not from picrust2 env: $(which place_seqs.py)"
  exit 1
fi

echo "Input files:"
ls -lh "$BIOM" "$FASTA"

echo "Optional input check:"
python scripts/check_picrust2_input.py \
  --biom "$BIOM" \
  --fasta "$FASTA" \
  --outdir "$PROJECT_ROOT/01.filter_input/check_before_picrust2" \
  --no-fail || true

if [ -d "$OUTDIR" ]; then
  if [ "$OVERWRITE" = "1" ]; then
    echo "OVERWRITE=1, removing old output:"
    echo "  $OUTDIR"
    rm -rf "$OUTDIR"
  else
    if [ -s "$OUTDIR/pathways_out/path_abun_unstrat.tsv.gz" ]; then
      echo "PICRUSt2 output already exists; skip."
      echo "Existing output: $OUTDIR/pathways_out/path_abun_unstrat.tsv.gz"
      exit 0
    else
      echo "Existing incomplete output directory found. Removing:"
      echo "  $OUTDIR"
      rm -rf "$OUTDIR"
    fi
  fi
fi

echo "Running PICRUSt2..."

picrust2_pipeline.py \
  -s "$FASTA" \
  -i "$BIOM" \
  -o "$OUTDIR" \
  -p "$THREADS" \
  --verbose

echo "Checking output..."

if [ -s "$OUTDIR/pathways_out/path_abun_unstrat.tsv.gz" ]; then
  echo "SUCCESS: $OUTDIR/pathways_out/path_abun_unstrat.tsv.gz"
else
  echo "ERROR: PICRUSt2 finished but path_abun_unstrat.tsv.gz was not found."
  echo "Output files:"
  find "$OUTDIR" -maxdepth 3 -type f | sort || true
  exit 1
fi

echo "Output files:"
find "$OUTDIR" -maxdepth 3 -type f | sort | head -200

echo "Job finished at: $(date)"
