#!/bin/bash
#SBATCH --job-name=dna_relax_smoke
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=%x_%j.out
set -euo pipefail

REPO=/project2/rohs_102/shewchuk/TF-conformation
MIN=${REPO}/stage3_minimize/stage3_minimize.py
IN=${REPO}/output/stage2_docked/tbp/1tgh_state_001.pdb
OUT=/scratch1/shewchuk/dna_relax_smoke
mkdir -p "$OUT"

source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate bioemu
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

echo "=== A: default (frozen DNA, no flag) ==="
python "$MIN" --input-pdb "$IN" --output-pdb "$OUT/A_frozen.pdb" \
    --scratch-dir "$OUT/tmpA"

echo "=== B: --dna-restraint-k 10.0 (same as protein k) ==="
python "$MIN" --input-pdb "$IN" --output-pdb "$OUT/B_dnak_same.pdb" \
    --scratch-dir "$OUT/tmpB" --dna-restraint-k 10.0

echo "=== C: --dna-restraint-k 0 (fully relaxed DNA) ==="
python "$MIN" --input-pdb "$IN" --output-pdb "$OUT/C_dnak0.pdb" \
    --scratch-dir "$OUT/tmpC" --dna-restraint-k 0

echo "=== ANALYSIS ==="
python "${REPO}/stage3_minimize/validation/analyze_dna_relax.py" "$IN" \
    A_frozen "$OUT/A_frozen.pdb" \
    B_dnak_same "$OUT/B_dnak_same.pdb" \
    C_dnak0 "$OUT/C_dnak0.pdb" | tee "$OUT/analysis.txt"
echo "=== DONE ==="
