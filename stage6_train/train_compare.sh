#!/bin/bash
# lib/train_compare.sh — train baseline + augmented for a TF, one SLURM array task each.
#
# Reads SLURM_ARRAY_TASK_ID:
#   0 -> baseline (uses orig train${FOLD}.txt + assembly2024/)
#   1 -> augmented (uses train${FOLD}_aug_${TF_NAME}.txt + combined_assembly_${TF_NAME}/)
#
# Optional env var SEED: if set, run/config names get a "_s${SEED}" suffix so
# this same lib drives both the single-seed default pair (scripts/pipeline/run_pilot.sh) and
# the multi-seed paired comparison (scripts/pipeline/run_multiseed_pilot.sh).

require_var TF_NAME
require_var FOLD

conda activate deeppbs

# Self-contained: run this stage's co-located scripts, not a shared SCRIPTS_DIR.
STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${RUN_DIR}"

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID must be set}"

SEED_SUFFIX=""
if [ -n "${SEED:-}" ]; then
    SEED_SUFFIX="_s${SEED}"
fi

if [ "${TASK_ID}" -eq 0 ]; then
    RUN_NAME="baseline_${TF_NAME}_fold${FOLD}${SEED_SUFFIX}"
    CONFIG="config_baseline_${TF_NAME}_fold${FOLD}${SEED_SUFFIX}.json"
    TRAIN_FILE="./folds/train${FOLD}.txt"
    VALID_FILE="./folds/valid${FOLD}.txt"
elif [ "${TASK_ID}" -eq 1 ]; then
    RUN_NAME="augmented_${TF_NAME}_fold${FOLD}${SEED_SUFFIX}"
    CONFIG="config_augmented_${TF_NAME}_fold${FOLD}${SEED_SUFFIX}.json"
    TRAIN_FILE="./folds_aug/train${FOLD}_aug_${TF_NAME}.txt"
    VALID_FILE="./folds_aug/valid${FOLD}_${TF_NAME}.txt"
    # Augmented runs use the same valid set as baseline (no augmentation in valid)
    if [ ! -f "${VALID_FILE}" ]; then
        # Symlink the original valid file rather than duplicating it
        mkdir -p ./folds_aug
        ln -sf "$(readlink -f ./folds/valid${FOLD}.txt)" "${VALID_FILE}"
    fi
else
    echo "ERROR: unknown SLURM_ARRAY_TASK_ID=${TASK_ID} (expected 0 or 1)" >&2
    exit 1
fi

# Sanity check inputs exist
for f in "${CONFIG}" "${TRAIN_FILE}" "${VALID_FILE}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: required file not found: ${f}" >&2
        exit 1
    fi
done

echo "=== ${RUN_NAME} ==="
nvidia-smi --query-gpu=name --format=csv,noheader || true
echo "  config:      ${CONFIG}"
echo "  train:       ${TRAIN_FILE} ($(wc -l < "${TRAIN_FILE}") lines)"
echo "  valid:       ${VALID_FILE} ($(wc -l < "${VALID_FILE}") lines)"

python -W ignore "${STAGE_DIR}/driver.py" \
    "${TRAIN_FILE}" "${VALID_FILE}" \
    -c "${CONFIG}" \
    --balance unmasked \
    --eval_every 1 \
    --single_gpu \
    --run_name "${RUN_NAME}"

echo "=== Done ${RUN_NAME} ==="
