#!/bin/bash
#SBATCH --job-name=prepost
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/analysis/figures/pymol/prepost_min/render_%j.out
set -eo pipefail
cd /project2/rohs_102/shewchuk/TF-conformation/analysis/figscripts
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate pymol
mkdir -p /project2/rohs_102/shewchuk/TF-conformation/analysis/figures/pymol/prepost_min
# pilots passed as args; default all
PILOTS="${@:-all}"
pymol -cq render_prepost_min.py -- $PILOTS
echo "RENDER DONE"
