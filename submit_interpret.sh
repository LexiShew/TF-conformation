#!/bin/bash
#SBATCH --job-name=interpret_ets1
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/interpret_ets1_%j.log

set -eo pipefail

# Activate deeppbs env
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/deeppbs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# Run interpret
cd /project2/rohs_102/shewchuk/TF-conformation
python analysis/interpret_tfconf.py \
  --pilot ets1 \
  --out output/interpret_results \
  --arms baseline_ets1_fold0 augmented_ets1_fold0 augmented_ets1_fold0_dnarelax_s1

echo "Attribution run complete"
