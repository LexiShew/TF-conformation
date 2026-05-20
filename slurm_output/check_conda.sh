#!/bin/bash
#SBATCH -n 1
#SBATCH -p rohs
#SBATCH --account=rohs_102
#SBATCH --time=00:05:00
#SBATCH --output=check2.out

source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate bioemu
export CONDA_ROOT=/home1/shewchuk/.conda
env | grep -i conda

# Replicate what BioEmu does
python -c "
from bioemu.utils import get_conda_prefix
print('get_conda_prefix() returned:', get_conda_prefix())
"
