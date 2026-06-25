#!/bin/bash
# lib/train_legacy_aug.sh — train ONLY the augmented_legacy runs for a TF.
# Used for A/B testing: baselines from the prior multi-seed run are kept as-is
# (since they don't depend on Stage 3 protocol), and we pair them with new
# augmented runs trained on the legacy (--ignore-metals) augmenting set.
#
# Requires: SEED env var (1..5), LEGACY=1 (set by orchestrator)
#
# Reads SLURM_ARRAY_TASK_ID — should be 0 (we only run the augmented arm).
# Will be submitted as --array=0-0 (single-task array for uniformity).

require_var TF_NAME
require_var FOLD
require_var FOLDS_AUG_DIR
require_var OUTPUTS_DIR

if [ "${LEGACY:-0}" != "1" ]; then
    echo "ERROR: train_legacy_aug.sh requires LEGACY=1 to be set" >&2
    exit 1
fi

conda activate deeppbs

# Self-contained: run this stage's co-located scripts, not a shared SCRIPTS_DIR.
STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${RUN_DIR}"

SEED="${SEED:?SEED env var must be set (1..5)}"
SEED_SUFFIX="_s${SEED}"

RUN_NAME="augmented_legacy_${TF_NAME}_fold${FOLD}${SEED_SUFFIX}"
CONFIG="${OUTPUTS_DIR}/config_augmented_legacy_${TF_NAME}_fold${FOLD}${SEED_SUFFIX}.json"
TRAIN_FILE="${FOLDS_AUG_DIR}/train${FOLD}_aug_legacy_${TF_NAME}.txt"
VALID_FILE="${FOLDS_AUG_DIR}/valid${FOLD}_${TF_NAME}.txt"

# Augmented runs use the original valid set (read-only input under RUN_DIR/folds)
if [ ! -f "${VALID_FILE}" ]; then
    mkdir -p "${FOLDS_AUG_DIR}"
    ln -sf "$(readlink -f ./folds/valid${FOLD}.txt)" "${VALID_FILE}"
fi

# Sanity check
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
