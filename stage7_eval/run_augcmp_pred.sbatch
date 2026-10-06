#!/bin/bash
#SBATCH --job-name=augcmp_pred
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48GB
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/augcmp_pred_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/augcmp_pred_%j.err
set -eo pipefail
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/deeppbs
R=/project2/rohs_102/shewchuk/TF-conformation
export OUTPUTS_DIR="$R/output/stage6_train"
cd $R/stage7_eval
nvidia-smi --query-gpu=name --format=csv,noheader || true

# Arm 1: augmented / frozen-DNA. Native inputs = stage5_aug features.
for tf in ets1 tbp hsf egr1 engrailed foxa csl err runx nfat lef1 irf; do
  echo "===== AUG-FROZEN: $tf ====="
  python predict_ab_pwms.py --tf $tf \
    --id-file $R/stage7_eval/confcmp_idfiles/id_${tf}_confcmp.txt \
    --combined-dir $R/output/stage5_aug/combined_assembly_${tf} \
    --outputs-dir $OUTPUTS_DIR --fold 0 --kinds augmented \
    --pred-subdir predictions_confcmp 2>&1 || echo "FAILED aug-frozen: $tf"
done

# Arm 2: augmented / relaxed-DNA. Run dirs are ..._fold0_dnarelax_s*, which
# predict_ab_pwms.py reaches by treating the fold as "0_dnarelax". Native
# inputs = stage5_aug_dnarelax features (same entry names, relaxed geometry).
for tf in ets1 tbp hsf egr1 engrailed foxa csl err runx nfat lef1 irf; do
  echo "===== AUG-RELAX: $tf ====="
  python predict_ab_pwms.py --tf $tf \
    --id-file $R/stage7_eval/confcmp_idfiles/id_${tf}_confcmp_relax.txt \
    --combined-dir $R/output/stage5_aug_dnarelax/combined_assembly_${tf} \
    --outputs-dir $OUTPUTS_DIR --fold 0_dnarelax --kinds augmented \
    --pred-subdir predictions_confcmp 2>&1 || echo "FAILED aug-relax: $tf"
done
echo "ALL DONE"
