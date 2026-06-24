#!/bin/bash
# lib/stage3_array.sh — invoked by per-PDB SLURM array wrappers.
# One array task per BioEmu frame.

require_var TF_NAME
require_var PDB_ID
require_var STAGE2_DIR
require_var STAGE3_DIR
require_var RAMP_STAGES
require_var STEPS_PER_STAGE

conda activate bioemu

# Self-contained: run this stage's co-located scripts, not a shared SCRIPTS_DIR.
STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Use SLURM_ARRAY_TASK_ID directly, or first positional arg if testing locally.
TASK_ID="${SLURM_ARRAY_TASK_ID:-${1:-}}"
if [ -z "${TASK_ID}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID not set and no task ID arg given." >&2
    echo "Usage (local test): TF_NAME=dux4 bash stage3_array.sh <task_id>" >&2
    exit 1
fi

# Zero-pad to 3 digits to match Stage 2 output convention
STATE=$(printf "%03d" "${TASK_ID}")
INPUT="${STAGE2_DIR}/${PDB_ID}_state_${STATE}.pdb"
OUTPUT="${STAGE3_DIR}/${PDB_ID}_state_${STATE}.pdb"

if [ ! -f "${INPUT}" ]; then
    echo "[stage3/${TF_NAME}] SKIP — input not found: ${INPUT}"
    exit 0
fi

if [ -f "${OUTPUT}" ]; then
    echo "[stage3/${TF_NAME}] SKIP — output already exists: ${OUTPUT}"
    exit 0
fi

echo "[stage3/${TF_NAME}] Minimizing state ${STATE}"
extra_args=()
if [ "${STAGE3_IGNORE_METALS:-0}" = "1" ]; then
    extra_args+=( --ignore-metals )
fi
python "${STAGE_DIR}/stage3_minimize.py" \
    --input-pdb "${INPUT}" \
    --output-pdb "${OUTPUT}" \
    --ramp-stages "${RAMP_STAGES}" \
    --steps-per-stage "${STEPS_PER_STAGE}" \
    "${extra_args[@]}"
