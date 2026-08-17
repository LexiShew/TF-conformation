#!/bin/bash
#SBATCH --job-name="af3_tfconf"
#SBATCH --array=1-6
#SBATCH --partition=qcbgpu
#SBATCH --account=rohs_102
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100GB
#SBATCH --time=2:00:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/af3/slurm_output/slurm-%A_%a.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/af3/slurm_output/slurm-%A_%a.err

# AlphaFold3 for TF-conformation pilots: predict each crystal DBD-DNA complex.
# 10 structures/pilot (modelSeeds=[1,2] x 5 diffusion samples each).
# Crystal DBD construct sequences + crystal DNA duplex; egr1 includes 3x ZN.
# array index -> pilot (1-based): 1=egr1, 2=tbp, 3=ets1, 4=foxa, 5=lef1, 6=engrailed

set -eo pipefail

module load apptainer
module load gcc/13.3.0
module load cuda/12.6.3
module load nvhpc/24.5

export AF3_RESOURCES_DIR=/project2/rohs_102/share/alphafold3
export AF3_IMAGE=${AF3_RESOURCES_DIR}/image/alphafold3.sif
export AF3_CODE_DIR=${AF3_RESOURCES_DIR}/code
export AF3_MODEL_PARAMETERS_DIR=${AF3_RESOURCES_DIR}/weights
export AF3_DATABASES_DIR=${AF3_RESOURCES_DIR}/databases

export RUN_DIR=/project2/rohs_102/shewchuk/TF-conformation/af3
export AF3_INPUT_DIR=${RUN_DIR}/input
export AF3_OUTPUT_DIR=${RUN_DIR}/output

# array index -> input JSON (order matches header comment)
INPUTS=(egr1_1aay.json tbp_1tgh.json ets1_1k79.json foxa_1vtn.json lef1_2lef.json engrailed_3hdd.json)
JSON_FILE=${INPUTS[$((SLURM_ARRAY_TASK_ID - 1))]}
echo "[af3] task ${SLURM_ARRAY_TASK_ID} -> ${JSON_FILE}"

apptainer exec \
     --nv \
     --bind $AF3_INPUT_DIR:/root/af_input \
     --bind $AF3_OUTPUT_DIR:/root/af_output \
     --bind $AF3_MODEL_PARAMETERS_DIR:/root/models \
     --bind $AF3_DATABASES_DIR:/root/public_databases \
     --bind $AF3_CODE_DIR:/root/code \
     $AF3_IMAGE \
     python /root/code/alphafold3/run_alphafold.py \
     --json_path=/root/af_input/$JSON_FILE \
     --model_dir=/root/models \
     --db_dir=/root/public_databases \
     --output_dir=/root/af_output

echo "[af3] task ${SLURM_ARRAY_TASK_ID} (${JSON_FILE}) done"
