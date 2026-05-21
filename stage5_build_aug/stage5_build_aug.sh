#!/bin/bash
# lib/stage5_build_aug.sh — build augmented train fold, combined assembly dir,
# and (default-seed) training configs for a TF. Fast pure-Python; runs inside
# the small SLURM job submitted by wrappers/stage5_build_aug.sh.
#
# Reads (via load_pilot_config): TF_NAME, FOLD, PWM_LABEL,
#   ORIG_FOLDS_DIR, STAGE4_DIR, FOLDS_AUG_DIR, ORIG_ASSEMBLY_DIR,
#   COMBINED_ASSEMBLY_DIR, SCRIPTS_DIR, RUN_DIR.

require_var TF_NAME
require_var FOLD
require_var PWM_LABEL
require_var STAGE4_DIR
require_var FOLDS_AUG_DIR
require_var COMBINED_ASSEMBLY_DIR

conda activate deeppbs

# 5a: build augmented train fold
python "${SCRIPTS_DIR}/build_augmented_fold.py" \
    --orig-train "${ORIG_FOLDS_DIR}/train${FOLD}.txt" \
    --stage4-dir "${STAGE4_DIR}/output" \
    --pwm-filter "${PWM_LABEL}" \
    --out-train  "${AUG_TRAIN_FOLD:-${FOLDS_AUG_DIR}/train${FOLD}_aug_${TF_NAME}.txt}"

# 5b: build combined assembly directory
python "${SCRIPTS_DIR}/build_combined_assembly.py" \
    --orig-dir "${ORIG_ASSEMBLY_DIR}" \
    --stage4-dir "${STAGE4_DIR}/output" \
    --pwm-filter "${PWM_LABEL}" \
    --out-dir "${COMBINED_ASSEMBLY_DIR}"

# 5c: build default-seed training configs (multiseed runs generate their own)
python "${SCRIPTS_DIR}/build_training_configs.py" \
    --tf-name "${TF_NAME}" \
    --combined-dir "${COMBINED_ASSEMBLY_DIR}" \
    --fold "${FOLD}" \
    --output-dir "${RUN_DIR}"

echo "[stage5/${TF_NAME}] DONE"