#!/bin/bash
#SBATCH --job-name=ab_pwms
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/ab_pwms_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/ab_pwms_%j.err
set -eo pipefail
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/deeppbs
export OUTPUTS_DIR="/project2/rohs_102/shewchuk/TF-conformation/output/stage6_train"
cd /project2/rohs_102/shewchuk/TF-conformation/stage7_eval
for tf in csl engrailed err egr1 foxa hsf irf nfat runx; do
  echo "===== PWM regen: $tf ====="
  python predict_ab_pwms.py \
    --tf $tf \
    --id-file /project2/rohs_102/shewchuk/TF-conformation/stage7_eval/ab_idfiles/id_${tf}_ab.txt \
    --combined-dir /project2/rohs_102/shewchuk/TF-conformation/output/stage5_aug/combined_assembly_${tf} \
    --outputs-dir $OUTPUTS_DIR \
    --fold 0 \
    --kinds baseline augmented 2>&1 || echo "FAILED: $tf"
done
echo "ALL DONE"
