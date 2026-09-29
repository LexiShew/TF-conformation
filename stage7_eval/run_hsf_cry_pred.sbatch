#!/bin/bash
#SBATCH --job-name=hsf_cry_pred
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --time=00:40:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/hsf_cry_pred_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/hsf_cry_pred_%j.err
set -eo pipefail
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/deeppbs
R=/project2/rohs_102/shewchuk/TF-conformation
cd $R/stage7_eval
python predict_ab_pwms.py --tf hsf \
  --id-file $R/stage7_eval/confcmp_idfiles/id_hsfcrystal_confcmp.txt \
  --combined-dir $R/output/stage5_aug/combined_assembly_hsf \
  --outputs-dir $R/output/stage6_train --fold 0 --kinds baseline \
  --pred-subdir predictions_confcmp
echo DONE
