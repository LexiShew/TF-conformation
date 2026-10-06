#!/bin/bash
#SBATCH --job-name=hsf_aug_cry
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --time=00:40:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/hsf_aug_cry_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/hsf_aug_cry_%j.err
set -eo pipefail
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/deeppbs
R=/project2/rohs_102/shewchuk/TF-conformation
cd $R/stage7_eval
# HSF1 crystal through the augmented/frozen arm (its id list predated the
# 5D5U featurization, so this entry was never predicted for that arm).
python predict_ab_pwms.py --tf hsf \
  --id-file $R/stage7_eval/confcmp_idfiles/id_hsfcrystal_confcmp.txt \
  --combined-dir $R/output/stage5_aug/combined_assembly_hsf \
  --outputs-dir $R/output/stage6_train --fold 0 --kinds augmented \
  --pred-subdir predictions_confcmp
# keep the id list consistent with the other pilots for future reruns
grep -qx "5d5u_duplex_MA0486.2.jaspar.npz" $R/stage7_eval/confcmp_idfiles/id_hsf_confcmp.txt \
  || echo "5d5u_duplex_MA0486.2.jaspar.npz" >> $R/stage7_eval/confcmp_idfiles/id_hsf_confcmp.txt
echo DONE
