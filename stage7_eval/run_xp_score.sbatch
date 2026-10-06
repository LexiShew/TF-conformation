#!/bin/bash
#SBATCH --job-name=xp_score
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/xp_score_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/xp_score_%j.err
set -eo pipefail
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/deeppbs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /project2/rohs_102/shewchuk/TF-conformation
python analysis/analyses/conf_vs_crystal_pwm/compare_crosspilot.py
