#!/bin/bash
#SBATCH --job-name=omm_probe
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/analysis/align_compare/omm_probe_%j.out
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh; conda activate bioemu
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
nvidia-smi -L
python -c "import openmm; print('OpenMM', openmm.__version__); print('platforms:', [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]); import openmm as mm; 
try:
    p=mm.Platform.getPlatformByName('CUDA'); print('CUDA OK')
except Exception as e: print('CUDA FAIL:', str(e)[:200])"
