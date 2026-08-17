#!/bin/bash
#SBATCH --job-name="make_S_pymol"
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=00:30:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/analysis/figures/pymol/make_S_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/analysis/figures/pymol/make_S_%j.err

# Render S1/S2 structure figures on a compute node (login node has no GL and
# cmd.png segfaults). PyMOL edu license auto-loads from ~/.pymol/license.lic.
#
# Usage:
#   sbatch render_S.sbatch                      # default: ets1 lef1, modes S1,S2
#   sbatch render_S.sbatch ets1 lef1 --mode S1  # explicit pilots/mode
#
# Args after the script name are forwarded verbatim to make_S.py (after `--`).

set -eo pipefail
RUN_DIR=/project2/rohs_102/shewchuk/TF-conformation
mkdir -p "$RUN_DIR/analysis/figures/pymol"

source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate pymol
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "$RUN_DIR"
echo "host=$(hostname) start=$(date)"
echo "args: $@"

pymol -cq analysis/figscripts/make_S.py -- "$@"

echo "done=$(date)"
ls -la analysis/figures/pymol/S*_*.png 2>/dev/null || echo "no S PNGs found"
