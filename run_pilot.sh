#!/bin/bash
# run_pilot.sh — submit a complete TF pilot end-to-end with SLURM dependencies.
#
# Usage:
#   ./run_pilot.sh <tf_name> [stage_start] [stage_end]
#   stages: 1=hpacker, 2=redock, 3=minimize, 3r=recover, 4=preprocess,
#           5=build_aug, 6=train, 7=eval
#   Defaults: stage_start=1, stage_end=7
#
# Examples:
#   ./run_pilot.sh dux4              # full pipeline
#   ./run_pilot.sh dux4 4 7          # skip Stages 1-3 (already done)
#   ./run_pilot.sh dux4 6 7          # just retrain + reeval
#
# Dependency policy:
#   afterany — Stage 3 array → 3r, 3r → 4. Partial failure on
#              minimization is expected and recovery is designed to handle it.
#              Without this, a single task failure halts the whole pipeline.
#   afterok  — everywhere else; partial failure is fatal.

set -eo pipefail
JOBS_DIR="/project2/rohs_102/shewchuk/DeepPBS/run/jobs"
# shellcheck source=lib/common.sh
source "${JOBS_DIR}/lib/common.sh"

TF_NAME="${1:?Usage: $0 <tf_name> [stage_start] [stage_end]}"
STAGE_START="${2:-1}"
STAGE_END="${3:-7}"

load_pilot_config "${TF_NAME}"
require_var TF_NAME
require_var PDB_ID
require_var N_FRAMES
require_var COMBINED_ASSEMBLY_DIR

# Stage IDs we'll submit (so dependencies chain correctly)
PREV_JOB=""

# submit_stage [stage_id] [script] [dep_mode] [extra_sbatch_args...]
# dep_mode: "ok" (default, afterok) or "any" (afterany)
submit_stage() {
    local stage="$1"; shift
    local script="$1"; shift
    local dep_mode="${1:-ok}"; shift  # first remaining arg is dep mode

    local depend=()
    if [ -n "${PREV_JOB}" ]; then
        depend=( --dependency="after${dep_mode}:${PREV_JOB}" )
    fi

    local jobid
    jobid=$(sbatch --parsable \
        --export=ALL,TF_NAME="${TF_NAME}" \
        "${depend[@]}" \
        "$@" \
        "${script}")
    echo "[run_pilot] Stage ${stage} submitted: jobid=${jobid} depend=after${dep_mode}:${PREV_JOB:-none}"
    PREV_JOB="${jobid}"
}

stage_in_range() {
    local s="$1"
    [ "${s}" -ge "${STAGE_START}" ] && [ "${s}" -le "${STAGE_END}" ]
}

# Stage 1: HPACKER
if stage_in_range 1; then
    submit_stage 1 "${JOBS_DIR}/wrappers/stage1_hpacker.sh" ok
fi

# Stage 2: redock
if stage_in_range 2; then
    submit_stage 2 "${JOBS_DIR}/wrappers/stage2_redock.sh" ok
fi

# Stage 3: array of N_FRAMES tasks
if stage_in_range 3; then
    submit_stage 3 "${JOBS_DIR}/wrappers/stage3_array.sh" ok \
        --array="1-${N_FRAMES}%8"
fi

# Stage 3 recovery — uses afterany since the array is expected to have some
# failures (those are exactly what recovery is for).
if stage_in_range 3 && [ "${STAGE_END}" -ge 3 ]; then
    submit_stage "3r" "${JOBS_DIR}/wrappers/stage3_recover.sh" any
fi

# Stage 4: preprocess — uses afterany so a partial recovery still runs.
# Inside Stage 4, the script counts npz outputs and warns if anything is off.
if stage_in_range 4; then
    submit_stage 4 "${JOBS_DIR}/wrappers/stage4_preprocess.sh" any
fi

# Stage 5: build augmented fold (must succeed before training)
if stage_in_range 5; then
    submit_stage 5 "${JOBS_DIR}/wrappers/stage5_build_aug.sh" ok
fi

# Stage 6: training (array of 2)
if stage_in_range 6; then
    submit_stage 6 "${JOBS_DIR}/wrappers/train_compare.sh" ok
fi

# Stage 7: eval
if stage_in_range 7; then
    submit_stage 7 "${JOBS_DIR}/wrappers/eval_benchmark.sh" ok
fi

echo ""
echo "[run_pilot] Pipeline submitted. Monitor with:"
echo "  squeue -u \$USER"
echo "[run_pilot] Final results will be at:"
echo "  ${OUTPUTS_DIR}/id_benchmark_${TF_NAME}.json"