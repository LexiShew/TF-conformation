#!/bin/bash
# lib/eval_legacy_ab.sh — eval for the A/B test.
# Scans for THREE groups of checkpoints:
#   baseline_<tf>_fold<F>_s*           (original baselines, unchanged)
#   augmented_<tf>_fold<F>_s*          (current metal-preserving augmenting)
#   augmented_legacy_<tf>_fold<F>_s*   (legacy --ignore-metals augmenting)
#
# Output is named id_benchmark_<tf>_legacy_ab.json so it doesn't clobber the
# main eval. The downstream paired stats will show three pairs per seed:
#   baseline_s1 vs augmented_s1            (current)
#   baseline_s1 vs augmented_legacy_s1     (legacy A/B)
#   (and implicit augmented_s1 vs augmented_legacy_s1 by difference)

require_var TF_NAME
require_var FOLD
require_var TEST_PWM_FILTER
require_var TEST_FILTER_NAME

conda activate deeppbs

# Self-contained: run this stage's co-located scripts, not a shared SCRIPTS_DIR.
STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${RUN_DIR}"

declare -a COND_ARGS

scan_run_subdirs() {
    local outer="$1"
    [ -d "${outer}" ] || return 0
    for sub in "${outer}"/*/; do
        [ -d "${sub}" ] || continue
        if [ -f "${sub}Model.best.tar" ]; then
            local name
            name=$(basename "${sub}")
            COND_ARGS+=( --condition "${name}=${sub%/}" )
            echo "[eval-ab/${TF_NAME}] Found checkpoint: ${name}"
        fi
    done
}

# Three groups of outer dirs (each with multi-seed _s* suffix)
for kind in baseline augmented augmented_legacy; do
    for outer in "${OUTPUTS_DIR}/${kind}_${TF_NAME}_fold${FOLD}"_s*/; do
        [ -d "${outer}" ] || continue
        scan_run_subdirs "${outer%/}"
    done
done

if [ "${#COND_ARGS[@]}" -eq 0 ]; then
    echo "ERROR: No checkpoints found" >&2
    exit 1
fi

echo "[eval-ab/${TF_NAME}] Total conditions found: $(( ${#COND_ARGS[@]} / 2 ))"

# IMPORTANT: combined-dir for eval is just where features for the test set
# entries live — since id.txt entries are part of the ORIGINAL data, we point
# at the combined_assembly dir which contains them. (Doesn't matter whether
# we use _legacy or non-legacy version — the id.txt entries are identical
# in both.)
COMBINED_FOR_EVAL="${DATA_DIR}/combined_assembly_${TF_NAME}"
OUTPUT_JSON="${OUTPUTS_DIR}/id_benchmark_${TF_NAME}_legacy_ab.json"

echo "[eval-ab/${TF_NAME}] Running evaluator"
python "${STAGE_DIR}/evaluate_id_benchmark.py" \
    "${COND_ARGS[@]}" \
    --combined-dir "${COMBINED_FOR_EVAL}" \
    --filter "${TEST_PWM_FILTER}" \
    --filter-name "${TEST_FILTER_NAME}" \
    --output-json "${OUTPUT_JSON}"

echo "[eval-ab/${TF_NAME}] DONE"
echo "[eval-ab/${TF_NAME}] Results: ${OUTPUT_JSON}"
