#!/bin/bash
#SBATCH --job-name=pycurves_batch
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/pycurves_%A_%a.out
# One array task = one chunk of the worklist. pyCurves' XLA thread pool needs a
# compute node (login-node RLIMIT_NPROC crashes it); CPU-only JAX is sufficient.
set -eo pipefail
BASE="/project2/rohs_102/shewchuk/TF-conformation"
WL="${BASE}/analysis/dna_relax/pycurves_worklist.tsv"
PYC_ENV="/project2/rohs_102/shewchuk/conda/envs/pycurves"
CHUNK="${CHUNK:-20}"                     # rows per array task
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate "${PYC_ENV}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
start=$(( TASK * CHUNK + 1 )); end=$(( start + CHUNK - 1 ))
sed -n "${start},${end}p" "${WL}" | while IFS=$'\t' read -r tf cond state inpdb outpref; do
    [ -z "$inpdb" ] && continue
    mkdir -p "$(dirname "$outpref")"
    if [ -f "${outpref}.txt" ]; then echo "SKIP ${tf}/${cond}_${state} (done)"; continue; fi
    echo "RUN ${tf}/${cond}_${state} <- ${inpdb}"
    pycurves "$inpdb" --format curves --output-file "${outpref}.txt" \
        && echo "  OK ${outpref}.txt" || echo "  FAIL ${tf}/${cond}_${state}"
done
echo "task ${TASK} done (rows ${start}-${end})"
