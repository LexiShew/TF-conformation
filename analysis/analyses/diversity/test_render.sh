#!/bin/bash
#SBATCH --job-name=pmtest
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=00:10:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/analysis/figures/pymol/test_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/analysis/figures/pymol/test_%j.err
set -eo pipefail
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh; conda activate pymol
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /project2/rohs_102/shewchuk/TF-conformation
echo "start $(date +%s)"
pymol -cq analysis/scripts/render_ensembles.py -- ets1
echo "end $(date +%s)"
ls -la analysis/figures/pymol/ets1_stack.png
