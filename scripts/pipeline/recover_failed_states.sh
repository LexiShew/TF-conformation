#!/bin/bash
# recover_failed_states.sh — array-parallel recovery of Stage-3 states that
# failed or timed out, ONE SLURM array task per missing state.
#
# Why this exists: stage3_recover.sh walks all N_FRAMES states SERIALLY in a
# single job and can hit its wall-clock when several states each need the
# expensive gentle ramp (this stalled the dux4 dimer trial: 7 stragglers x
# ~30 min serial ran into the 4 h limit, so its fnat gate never fired). This
# submits ONLY the missing states as a throttled array, so they recover in
# parallel and one slow state cannot sink the batch.
#
# Usage:
#   ./scripts/pipeline/recover_failed_states.sh <tf_name> [throttle]
#     throttle = max concurrent array tasks (default 4)
#
# After it completes, run the gate onward:
#   ./scripts/pipeline/run_pilot.sh <tf_name> 4 7
#
# Safety: idempotent (skips states that already have Stage-3 output) and
# ADDITIVE - it does not modify or depend on stage3_recover.sh, so the current
# pipeline behaviour is unchanged.

set -eo pipefail
export TFCONF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${TFCONF_DIR}"
# shellcheck source=../../lib/common.sh
source "${TFCONF_DIR}/lib/common.sh"

TF_NAME="${1:?Usage: $0 <tf_name> [throttle]}"
THROTTLE="${2:-4}"

load_pilot_config "${TF_NAME}"
require_var PDB_ID
require_var STAGE2_DIR
require_var STAGE3_DIR
require_var N_FRAMES

# States that have a Stage-2 input but NO Stage-3 output = failed/timed-out.
missing=()
for i in $(seq 1 "${N_FRAMES}"); do
    s=$(printf "%03d" "${i}")
    [ -f "${STAGE2_DIR}/${PDB_ID}_state_${s}.pdb" ] || continue
    [ -f "${STAGE3_DIR}/${PDB_ID}_state_${s}.pdb" ] && continue
    missing+=( "${i}" )
done

if [ "${#missing[@]}" -eq 0 ]; then
    echo "[recover/${TF_NAME}] No missing Stage-3 states - nothing to recover."
    exit 0
fi

ARRAY_SPEC=$(IFS=,; echo "${missing[*]}")
echo "[recover/${TF_NAME}] ${#missing[@]} missing state(s): ${ARRAY_SPEC}"
echo "[recover/${TF_NAME}] Submitting array --array=${ARRAY_SPEC}%${THROTTLE}"

JOB=$(sbatch --parsable \
    --array="${ARRAY_SPEC}%${THROTTLE}" \
    --export=ALL,TF_NAME="${TF_NAME}" \
    "${TFCONF_DIR}/wrappers/stage3_recover_array.sh")
echo "[recover/${TF_NAME}] Submitted recovery array: jobid=${JOB}"
echo "[recover/${TF_NAME}] When it finishes, continue with:"
echo "    ./scripts/pipeline/run_pilot.sh ${TF_NAME} 4 7"
