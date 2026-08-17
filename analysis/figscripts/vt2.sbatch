#!/bin/bash
#SBATCH -p rohs
#SBATCH -A rohs_102
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH -o /project2/rohs_102/shewchuk/TF-conformation/analysis/figures/pymol/_viewtest/vt2_%j.out
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh; conda activate pymol
pymol -cq /project2/rohs_102/shewchuk/TF-conformation/analysis/figscripts/view_test2.py 2>&1
echo VT2_END
