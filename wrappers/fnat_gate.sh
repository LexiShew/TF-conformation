#!/bin/bash
#SBATCH --job-name=fnat_gate
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/fnat_gate_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/fnat_gate_%j.err

set -eo pipefail
TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
# shellcheck source=../lib/common.sh
source "${TFCONF_DIR}/lib/common.sh"
: "${TF_NAME:?TF_NAME must be set; sbatch --export=ALL,TF_NAME=<tf> ...}"
load_pilot_config "${TF_NAME}"

# shellcheck source=../fnat_gate/fnat_gate.sh
source "${TFCONF_DIR}/fnat_gate/fnat_gate.sh"
