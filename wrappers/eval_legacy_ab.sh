#!/bin/bash
#SBATCH --job-name=eval_ab
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --account=rohs_102
#SBATCH --partition=rohs
#SBATCH --gres=gpu:1
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/eval_ab_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/eval_ab_%j.err

set -eo pipefail
# TF-conformation is the authoritative pipeline root. Prefer an inherited
# TFCONF_DIR (exported by the launcher), else the cluster default.
TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
# shellcheck source=../lib/common.sh
source "${TFCONF_DIR}/lib/common.sh"
: "${TF_NAME:?TF_NAME must be set; sbatch --export=ALL,TF_NAME=<tf> ...}"
load_pilot_config "${TF_NAME}"

# shellcheck source=../stage7_eval/eval_legacy_ab.sh
source "${TFCONF_DIR}/stage7_eval/eval_legacy_ab.sh"
