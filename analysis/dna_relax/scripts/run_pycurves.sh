#!/bin/bash
# Run pyCurves on a DNA-protein PDB. pyCurves' JAX/XLA thread pool exceeds the
# login-node per-user RLIMIT_NPROC, so this MUST run on a compute node (srun/sbatch).
# CPU-only JAX in the pycurves env is sufficient for static structure analysis.
#   Usage: bash run_pycurves.sh <input.pdb> <out_prefix>
# On a compute node directly, or wrap in: srun -p rohs -A rohs_102 -c 4 -t 00:05:00 bash run_pycurves.sh ...
set -eo pipefail
IN="$1"; OUT="$2"
PYC_ENV="/project2/rohs_102/shewchuk/conda/envs/pycurves"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate "${PYC_ENV}"
pycurves "${IN}" --format curves --output-file "${OUT}.txt"
pycurves "${IN}" --format json  --output-file "${OUT}.json"
echo "pyCurves done -> ${OUT}.txt / ${OUT}.json"
