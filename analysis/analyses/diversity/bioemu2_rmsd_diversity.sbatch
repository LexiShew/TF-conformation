#!/bin/bash
#SBATCH --job-name=be2div
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24GB
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/logs/bioemu2/be2div_%j.out
#
# CPU-only: loads BioEmu2 state PDBs and does Kabsch Calpha RMSDs. No GPU.
#   sbatch --partition=qcbr --account=rohs_102 \
#     analysis/analyses/diversity/bioemu2_rmsd_diversity.sbatch
set -eo pipefail
TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
cd "${TFCONF_DIR}"
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate bioemu          # mdtraj 1.11 lives here, per af3_rmsd_diversity.py
python analysis/analyses/diversity/bioemu2_rmsd_diversity.py "$@"
