#!/bin/bash
#SBATCH --job-name=figs_interpret
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/figs_interpret_%j.log

source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/deeppbs

cd /project2/rohs_102/shewchuk/TF-conformation
python analysis/make_importance_figures.py
echo "Figure generation complete"
