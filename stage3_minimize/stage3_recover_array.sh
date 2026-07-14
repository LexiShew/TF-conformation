#!/bin/bash
# lib/stage3_recover_array.sh — array-parallel recovery of failed Stage-3 states.
#
# One SLURM array task per state, using the gentle RECOVERY ramp. This is the
# parallel counterpart to stage3_recover.sh, which walks all N_FRAMES states
# SERIALLY in a single job and can hit its wall-clock when several states each
# need the expensive gentle ramp (this is what stalled the dux4 dimer trial:
# 7 stragglers x ~30 min serial ran into the 4 h limit). Running them as an
# array lets them recover in parallel; one slow state can no longer sink the
# whole batch.
#
# Idempotent: skips any state that already has Stage-3 output (identical guard
# to stage3_array.sh / stage3_recover.sh), so it only ever works on genuinely
# missing states. Best-effort: a state that still cannot be minimized logs
# loudly and exits 0, matching the serial recover's swallow-and-continue
# semantics (so a downstream afterok gate is not blocked by one stuck state).
#
# Submit via scripts/pipeline/recover_failed_states.sh, which computes the
# missing-state list and sets --array accordingly.

require_var TF_NAME
require_var PDB_ID
require_var STAGE2_DIR
require_var STAGE3_DIR
require_var RECOVERY_RAMP_STAGES
require_var RECOVERY_STEPS_PER_STAGE

conda activate "${BIOEMU_ENV:-bioemu}"

# Self-contained: run this stage's co-located scripts, not a shared SCRIPTS_DIR.
STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Use SLURM_ARRAY_TASK_ID directly, or first positional arg if testing locally.
TASK_ID="${SLURM_ARRAY_TASK_ID:-${1:-}}"
if [ -z "${TASK_ID}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID not set and no task ID arg given." >&2
    echo "Usage (local test): TF_NAME=dux4 bash stage3_recover_array.sh <task_id>" >&2
    exit 1
fi

# Zero-pad to 3 digits to match Stage 2 output convention.
STATE=$(printf "%03d" "${TASK_ID}")
INPUT="${STAGE2_DIR}/${PDB_ID}_state_${STATE}.pdb"
OUTPUT="${STAGE3_DIR}/${PDB_ID}_state_${STATE}.pdb"

if [ ! -f "${INPUT}" ]; then
    echo "[stage3-recover/${TF_NAME}] SKIP - input not found: ${INPUT}"
    exit 0
fi
if [ -f "${OUTPUT}" ]; then
    echo "[stage3-recover/${TF_NAME}] SKIP - output already exists: ${OUTPUT}"
    exit 0
fi

echo "[stage3-recover/${TF_NAME}] Recovering state ${STATE} with gentle ramp"
extra_args=()
if [ "${STAGE3_IGNORE_METALS:-0}" = "1" ]; then
    extra_args+=( --ignore-metals )
fi
if [ -n "${STAGE3_DNA_RESTRAINT_K:-}" ]; then
    extra_args+=( --dna-restraint-k "${STAGE3_DNA_RESTRAINT_K}" )
fi

# Best-effort: do NOT let a single stuck state fail the array task.
if python "${STAGE_DIR}/stage3_minimize.py" \
    --input-pdb "${INPUT}" \
    --output-pdb "${OUTPUT}" \
    --ramp-stages "${RECOVERY_RAMP_STAGES}" \
    --steps-per-stage "${RECOVERY_STEPS_PER_STAGE}" \
    "${extra_args[@]}"; then
    echo "[stage3-recover/${TF_NAME}] state ${STATE} RECOVERED"
else
    echo "[stage3-recover/${TF_NAME}] state ${STATE} STILL FAILED - leaving for inspection"
    exit 0
fi
