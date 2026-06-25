#!/bin/bash
# lib/stage3_recover.sh — re-attempt failed Stage 3 minimizations with a
# gentler ramp protocol. Skips any state that already has output.

require_var TF_NAME
require_var PDB_ID
require_var STAGE2_DIR
require_var STAGE3_DIR
require_var RECOVERY_RAMP_STAGES
require_var RECOVERY_STEPS_PER_STAGE
require_var N_FRAMES

conda activate "${BIOEMU_ENV:-bioemu}"

# Self-contained: run this stage's co-located scripts, not a shared SCRIPTS_DIR.
STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

n_recovered=0
n_attempted=0
n_still_failed=0

for i in $(seq 1 "${N_FRAMES}"); do
    STATE=$(printf "%03d" "${i}")
    INPUT="${STAGE2_DIR}/${PDB_ID}_state_${STATE}.pdb"
    OUTPUT="${STAGE3_DIR}/${PDB_ID}_state_${STATE}.pdb"

    if [ ! -f "${INPUT}" ]; then
        continue   # no Stage 2 input (rare)
    fi
    if [ -f "${OUTPUT}" ]; then
        continue   # already succeeded
    fi

    n_attempted=$((n_attempted + 1))
    echo "[stage3-recover/${TF_NAME}] Retrying state ${STATE} with gentler ramp"
    extra_args=()
    if [ "${STAGE3_IGNORE_METALS:-0}" = "1" ]; then
        extra_args+=( --ignore-metals )
    fi
    if python "${STAGE_DIR}/stage3_minimize.py" \
        --input-pdb "${INPUT}" \
        --output-pdb "${OUTPUT}" \
        --ramp-stages "${RECOVERY_RAMP_STAGES}" \
        --steps-per-stage "${RECOVERY_STEPS_PER_STAGE}" \
        "${extra_args[@]}"; then
        n_recovered=$((n_recovered + 1))
    else
        n_still_failed=$((n_still_failed + 1))
        echo "[stage3-recover/${TF_NAME}] State ${STATE} still failed"
    fi
done

echo "[stage3-recover/${TF_NAME}] DONE"
echo "  Attempted recovery: ${n_attempted}"
echo "  Recovered:          ${n_recovered}"
echo "  Still failed:       ${n_still_failed}"
