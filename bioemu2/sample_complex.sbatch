#!/bin/bash
#SBATCH --job-name=be2
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:1
#
# NOTE ON THE GPU REQUEST. This deliberately asks for a generic gpu:1 rather
# than the gpu:rtx5000:1 the stage wrappers use, because BioEmu2's model card
# states "an 80 GB A100 being the smallest tested device" and an RTX 5000 is
# 16 GB. Check what the cluster actually has and pin it if there is something
# bigger:
#     sinfo -o '%P %G %m %N' | sort -u
#     sbatch --gres=gpu:a100:1 bioemu2/sample_complex.sbatch
# If it OOMs on a 16 GB card that is expected, not a bug in this script.
#
# Logs: the launcher below uses lib/common.sh's log_args, so filenames are
# timestamped the same way the pipeline's are.
#
# Usage:
#   sbatch --export=ALL,TF_NAME=ets1 bioemu2/sample_complex.sbatch
#   sbatch --export=ALL,TF_NAME=egr1,NUM_SAMPLES=50,EXCLUDE_LIGANDS=ZN \
#          bioemu2/sample_complex.sbatch
#
# Why this lives here and not in wrappers/: wrappers/ holds the SBATCH headers
# for the numbered pipeline DAG that run_pilot.sh submits. This is an
# exploratory parallel track outside that DAG, so it is self-contained in its
# own dir, matching how analysis/analyses/*/ keep their own .sbatch files.

set -eo pipefail

TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
# shellcheck source=../lib/common.sh
source "${TFCONF_DIR}/lib/common.sh"

: "${TF_NAME:?TF_NAME must be set; sbatch --export=ALL,TF_NAME=<tf> ...}"

# BioEmu2 needs its own env: the repo's 'bioemu' env is Python 3.11 and BioEmu2
# requires 3.12+. Overridable, like the other env vars in common.sh.
export BIOEMU2_ENV="${BIOEMU2_ENV:-bioemu2}"
conda activate "${BIOEMU2_ENV}"

NUM_SAMPLES="${NUM_SAMPLES:-100}"
SEED="${SEED:-1}"

EXTRA=()
[ -n "${EXCLUDE_LIGANDS:-}" ] && EXTRA+=( --exclude-ligands ${EXCLUDE_LIGANDS} )
[ "${NO_LIGANDS:-0}" = "1" ]  && EXTRA+=( --no-ligands )
[ -n "${CHECKPOINT:-}" ]      && EXTRA+=( --checkpoint "${CHECKPOINT}" )

echo "[bioemu2/${TF_NAME}] env=${BIOEMU2_ENV} samples=${NUM_SAMPLES} seed=${SEED}"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

python "${TFCONF_DIR}/bioemu2/sample_complex.py" \
    --tf-name "${TF_NAME}" \
    --num-samples "${NUM_SAMPLES}" \
    --seed "${SEED}" \
    "${EXTRA[@]}"

echo "[bioemu2/${TF_NAME}] sampling done; comparing against the existing pipeline"
conda activate "${DEEPPBS_ENV:-deeppbs}"   # fnat_gate + scipy live here
python "${TFCONF_DIR}/bioemu2/compare_to_pipeline.py" \
    --tf-name "${TF_NAME}" \
    --arms frozen dnarelax \
    --out-csv "${TFCONF_DIR}/output/bioemu2/${TF_NAME}_comparison.csv" || \
    echo "[bioemu2/${TF_NAME}] comparison skipped (no pipeline states?)"

echo "[bioemu2/${TF_NAME}] DONE"
