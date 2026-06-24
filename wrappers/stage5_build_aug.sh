#!/bin/bash
#SBATCH --job-name=build_aug
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/build_aug_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/build_aug_%j.err

set -eo pipefail
# TF-conformation is the authoritative pipeline root. Prefer an inherited
# TFCONF_DIR (exported by the launcher), else the cluster default.
TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
# shellcheck source=../lib/common.sh
source "${TFCONF_DIR}/lib/common.sh"
: "${TF_NAME:?TF_NAME must be set; sbatch --export=ALL,TF_NAME=<tf> ...}"
load_pilot_config "${TF_NAME}"

# shellcheck source=../stage5_build_aug/stage5_build_aug.sh
source "${TFCONF_DIR}/stage5_build_aug/stage5_build_aug.sh"
