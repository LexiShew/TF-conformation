#!/bin/bash
#SBATCH --job-name=crosspilot
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --time=14:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48GB
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/crosspilot_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/crosspilot_%j.err
set -eo pipefail
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/deeppbs
R=/project2/rohs_102/shewchuk/TF-conformation
export OUTPUTS_DIR="$R/output/stage6_train"
cd $R/stage7_eval
nvidia-smi --query-gpu=name --format=csv,noheader || true

PILOTS="csl egr1 engrailed err ets1 foxa hsf irf lef1 nfat runx tbp"

# Cross-pilot test: score TARGET pilot X's conformations + crystal under DONOR
# pilot Y's checkpoints (Y != X). Y's augmented arm was trained on Y's states,
# never on X's, so X is unseen at the structure level -- not merely at the
# entry level. Comparing augmented_Y vs baseline_Y on the SAME X entries
# isolates what augmentation buys on genuinely novel complexes.
for X in $PILOTS; do
  for Y in $PILOTS; do
    [ "$X" = "$Y" ] && continue
    for K in baseline augmented; do
      echo "===== target=$X donor=$Y kind=$K ====="
      python predict_ab_pwms.py --tf $Y \
        --id-file $R/stage7_eval/confcmp_idfiles/id_${X}_confcmp.txt \
        --combined-dir $R/output/stage5_aug/combined_assembly_${X} \
        --outputs-dir $OUTPUTS_DIR --fold 0 --kinds $K \
        --pred-subdir predictions_crosspilot 2>&1 | tail -3 || echo "FAILED $X/$Y/$K"
    done
  done
done
echo "ALL DONE"
