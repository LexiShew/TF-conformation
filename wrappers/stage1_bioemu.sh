#!/bin/bash
#SBATCH --job-name=s1_bioemu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/s1_bioemu_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/s1_bioemu_%j.err
#
# Per-pilot Stage 1: generate BioEmu + HPacker conformations for THIS pilot's
# structure (all of its protein chains). Supersedes the old stage1_hpacker
# wrapper. To generate the WHOLE library at once instead, use
# stage1_bioemu/submit_stage1_bioemu.sh.

set -eo pipefail
TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
# shellcheck source=../lib/common.sh
source "${TFCONF_DIR}/lib/common.sh"
: "${TF_NAME:?TF_NAME must be set; sbatch --export=ALL,TF_NAME=<tf> ...}"
load_pilot_config "${TF_NAME}"

# BIOEMU_RAW_ROOT and PDB_ID come from common.sh + the pilot config; N_FRAMES is
# the per-pilot conformation count. stage1_bioemu.sh honors TFCONF_DIR.
bash "${TFCONF_DIR}/stage1_bioemu/stage1_bioemu.sh" \
    "${BIOEMU_RAW_ROOT}/${PDB_ID}_chains" "${N_FRAMES:-100}"
