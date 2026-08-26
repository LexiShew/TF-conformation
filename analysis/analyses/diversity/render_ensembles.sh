#!/bin/bash
#SBATCH --job-name="pymol_ensembles"
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=00:30:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/analysis/figures/pymol/render_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/analysis/figures/pymol/render_%j.err

set -eo pipefail
RUN_DIR=/project2/rohs_102/shewchuk/TF-conformation
mkdir -p "$RUN_DIR/analysis/figures/pymol"

source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate pymol
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "$RUN_DIR"
echo "host=$(hostname) start=$(date)"
pymol -cq analysis/analyses/diversity/render_ensembles.py -- all
echo "done=$(date)"
ls -la analysis/figures/pymol/*.png
