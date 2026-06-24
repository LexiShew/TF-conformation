#!/bin/bash
#SBATCH --job-name=train_legacy
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --account=rohs_102
#SBATCH --partition=rohs
#SBATCH --gres=gpu:1
#SBATCH --array=0-0
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/train_legacy_%A_%a.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/train_legacy_%A_%a.err

set -eo pipefail
# TF-conformation is the authoritative pipeline root. Prefer an inherited
# TFCONF_DIR (exported by the launcher), else the cluster default.
TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
# shellcheck source=../lib/common.sh
source "${TFCONF_DIR}/lib/common.sh"
: "${TF_NAME:?TF_NAME must be set; sbatch --export=ALL,TF_NAME=<tf> ...}"
load_pilot_config "${TF_NAME}"

# shellcheck source=../stage6_train/train_legacy_aug.sh
source "${TFCONF_DIR}/stage6_train/train_legacy_aug.sh"
