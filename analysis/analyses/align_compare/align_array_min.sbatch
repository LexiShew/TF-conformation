#!/bin/bash
#SBATCH --job-name=alignmin
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/analysis/analyses/align_compare/logs/min_%A_%a.out
set -eo pipefail
cd /project2/rohs_102/shewchuk/TF-conformation
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh; conda activate bioemu
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

# The task list (one line per state to minimize) is prebuilt at:
#   analysis/analyses/align_compare/tasklist.txt   with columns:  input_pdb  output_pdb
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" analysis/analyses/align_compare/tasklist.txt)
IN=$(echo "$LINE" | awk '{print $1}')
OUT=$(echo "$LINE" | awk '{print $2}')
mkdir -p "$(dirname "$OUT")"
python stage3_minimize/stage3_minimize.py \
  --input-pdb "$IN" --output-pdb "$OUT" \
  --scratch-dir /scratch1/shewchuk/deeppbs_min_tmp
echo "DONE task $SLURM_ARRAY_TASK_ID -> $OUT"
