#!/bin/bash
# lib/eval_benchmark.sh — run id.txt benchmark eval for a TF.
# Auto-discovers all checkpoints in:
#   DeepPBS_outputs/{baseline,augmented}_${TF_NAME}_fold${FOLD}/         (default-seed)
#   DeepPBS_outputs/{baseline,augmented}_${TF_NAME}_fold${FOLD}_s*/       (multi-seed, _s1, _s2, ...)
# This supports both single-default-seed runs and paired multi-seed pipelines.

require_var TF_NAME
require_var FOLD
require_var TEST_PWM_FILTER
require_var TEST_FILTER_NAME
require_var COMBINED_ASSEMBLY_DIR

conda activate deeppbs

# Self-contained: run this stage's co-located scripts, not a shared SCRIPTS_DIR.
STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${OUTPUT_ROOT}"

# Build --condition args from whatever subdirs exist
declare -a COND_ARGS

# scan_run_subdirs <kind> <outer_dir>
# Looks inside <outer_dir> for any subdir containing Model.best.tar.
scan_run_subdirs() {
    local kind="$1"
    local outer="$2"
    [ -d "${outer}" ] || return 0
    for sub in "${outer}"/*/; do
        [ -d "${sub}" ] || continue
        if [ -f "${sub}Model.best.tar" ]; then
            local name
            name=$(basename "${sub}")
            COND_ARGS+=( --condition "${name}=${sub%/}" )
            echo "[eval/${TF_NAME}] Found checkpoint: ${name}"
        fi
    done
}

for kind in baseline augmented; do
    # Default-seed outer dir (legacy/default)
    scan_run_subdirs "${kind}" "${OUTPUTS_DIR}/${kind}_${TF_NAME}_fold${FOLD}"
    # Multi-seed outer dirs (_s1, _s2, ...)
    for outer in "${OUTPUTS_DIR}/${kind}_${TF_NAME}_fold${FOLD}"_s*/; do
        [ -d "${outer}" ] || continue
        scan_run_subdirs "${kind}" "${outer%/}"
    done
done

if [ "${#COND_ARGS[@]}" -eq 0 ]; then
    echo "ERROR: No checkpoints found in ${OUTPUTS_DIR}/{baseline,augmented}_${TF_NAME}_fold${FOLD}/ or _s*/" >&2
    exit 1
fi

echo "[eval/${TF_NAME}] Total conditions found: $(( ${#COND_ARGS[@]} / 2 ))"

OUTPUT_JSON="${EVAL_OUT_DIR}/id_benchmark_${TF_NAME}.json"

echo "[eval/${TF_NAME}] Running evaluator"
python "${STAGE_DIR}/evaluate_id_benchmark.py" \
    "${COND_ARGS[@]}" \
    --combined-dir "${COMBINED_ASSEMBLY_DIR}" \
    --filter "${TEST_PWM_FILTER}" \
    --filter-name "${TEST_FILTER_NAME}" \
    --output-json "${OUTPUT_JSON}"

echo "[eval/${TF_NAME}] DONE"
echo "[eval/${TF_NAME}] Results: ${OUTPUT_JSON}"
