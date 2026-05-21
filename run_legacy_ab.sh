#!/bin/bash
# run_legacy_ab.sh — end-to-end A/B test: re-process a pilot through Stage 3
# WITHOUT the metal coordination cage (--ignore-metals), then train 5 paired
# augmented_legacy seeds, then eval against existing baselines AND existing
# augmented.
#
# Reuses the existing Stage 1 (HPACKER) and Stage 2 (docked) outputs since
# only Stage 3 protocol differs. Skips legacy Stage 2 reprocessing entirely:
# the legacy Stage 2 PDBs contain metals (carried from reference) but Stage 3
# --ignore-metals will strip them before minimization. This is equivalent to
# the pre-patch behavior where metals weren't in Stage 2 outputs at all.
#
# Usage:
#   ./run_legacy_ab.sh <tf_name> [start_stage]
#
# Stages:
#   3  Stage 3 with --ignore-metals  (writes to stage3_min_legacy/)
#   3r Stage 3 recovery
#   4  Stage 4 preprocess           (writes to stage4_npz_legacy/)
#   5  build legacy augmenting fold + combined_assembly_legacy_<tf>/
#   6  train 5 paired augmented_legacy_<tf>_s{1..5} runs
#   7  eval comparing baseline / augmented / augmented_legacy
#
# Default: start_stage=3 (full A/B run).

set -eo pipefail
JOBS_DIR="/project2/rohs_102/shewchuk/DeepPBS/run/jobs"
source "${JOBS_DIR}/lib/common.sh"

TF_NAME="${1:?Usage: $0 <tf_name> [start_stage]}"
START_STAGE="${2:-3}"
N_SEEDS=5

# Set LEGACY=1 BEFORE loading pilot config so derived paths get suffix
export LEGACY=1
load_pilot_config "${TF_NAME}"

echo "[legacy-ab/${TF_NAME}] Starting A/B run"
echo "  STAGE3_DIR:           ${STAGE3_DIR}"
echo "  STAGE4_DIR:           ${STAGE4_DIR}"
echo "  COMBINED_ASSEMBLY:    ${COMBINED_ASSEMBLY_DIR}"
echo "  AUG_TRAIN_FOLD:       ${AUG_TRAIN_FOLD}"
echo "  STAGE3_IGNORE_METALS: ${STAGE3_IGNORE_METALS}"

PREV_JOB=""
submit_stage() {
    local stage="$1"; shift
    local script="$1"; shift
    local dep_mode="${1:-ok}"; shift
    local depend=()
    if [ -n "${PREV_JOB}" ]; then
        depend=( --dependency="after${dep_mode}:${PREV_JOB}" )
    fi
    local jobid
    jobid=$(sbatch --parsable \
        --export=ALL,TF_NAME="${TF_NAME}",LEGACY=1 \
        "${depend[@]}" \
        "$@" \
        "${script}")
    echo "[legacy-ab] Stage ${stage} jobid=${jobid} depend=after${dep_mode}:${PREV_JOB:-none}"
    PREV_JOB="${jobid}"
}

stage_in_range() { [ "${1}" -ge "${START_STAGE}" ]; }

# --- Stage 3: array of N_FRAMES tasks ---
if stage_in_range 3; then
    submit_stage 3 "${JOBS_DIR}/wrappers/stage3_array.sh" ok \
        --array="1-${N_FRAMES}%8"
fi

# --- Stage 3r: recovery ---
if stage_in_range 3; then
    submit_stage "3r" "${JOBS_DIR}/wrappers/stage3_recover.sh" any
fi

# --- Stage 4: preprocess ---
if stage_in_range 4; then
    submit_stage 4 "${JOBS_DIR}/wrappers/stage4_preprocess.sh" any
fi

# --- Stage 5: build legacy aug + combined assembly ---
# We use stage5_build_aug.sh which honors AUG_TRAIN_FOLD env var (set by
# common.sh when LEGACY=1).
if stage_in_range 5; then
    submit_stage 5 "${JOBS_DIR}/wrappers/stage5_build_aug.sh" ok
fi

# --- Stage 5.5: build legacy training configs (this runs on the build_aug
# host since it's pure Python; piggyback as part of Stage 5 in a small extra job).
if stage_in_range 5; then
    BUILD_CFG_JOB=$(sbatch --parsable \
        --export=ALL,TF_NAME="${TF_NAME}",LEGACY=1 \
        --dependency="afterok:${PREV_JOB}" \
        --job-name=build_cfgs_legacy \
        --time=00:05:00 \
        --ntasks=1 --cpus-per-task=1 --mem=4GB \
        --account=rohs_102 --partition=rohs \
        --output="${LOGS_DIR}/build_cfgs_legacy_%j.out" \
        --error="${LOGS_DIR}/build_cfgs_legacy_%j.err" \
        --wrap="set -eo pipefail; source ${JOBS_DIR}/lib/common.sh; load_pilot_config ${TF_NAME}; conda activate deeppbs; python ${SCRIPTS_DIR}/build_legacy_training_configs.py --tf-name ${TF_NAME} --combined-dir ${COMBINED_ASSEMBLY_DIR} --fold ${FOLD} --seeds 1 2 3 4 5 --output-dir ${RUN_DIR}")
    echo "[legacy-ab] build_cfgs jobid=${BUILD_CFG_JOB} depend=afterok:${PREV_JOB}"
    PREV_JOB="${BUILD_CFG_JOB}"
fi

# --- Stage 6: train 5 paired augmented_legacy seeds ---
# Each seed is a separate sbatch (different SEED env var). All depend on
# build_cfgs completing. Track all of them so eval can depend on all.
TRAIN_JOBS=()
if stage_in_range 6; then
    for SEED in 1 2 3 4 5; do
        TRAIN_JOB=$(sbatch --parsable \
            --export=ALL,TF_NAME="${TF_NAME}",LEGACY=1,SEED="${SEED}" \
            --dependency="afterok:${PREV_JOB}" \
            --job-name="trleg_${TF_NAME}_s${SEED}" \
            "${JOBS_DIR}/wrappers/train_legacy_aug.sh")
        echo "[legacy-ab] train_legacy ${TF_NAME} seed=${SEED} jobid=${TRAIN_JOB}"
        TRAIN_JOBS+=( "${TRAIN_JOB}" )
    done
fi

# --- Stage 7: eval ---
# Depends on ALL 5 training jobs. Build the colon-separated dep list.
if stage_in_range 7; then
    if [ "${#TRAIN_JOBS[@]}" -gt 0 ]; then
        DEP_LIST=$(IFS=:; echo "${TRAIN_JOBS[*]}")
        EVAL_JOB=$(sbatch --parsable \
            --export=ALL,TF_NAME="${TF_NAME}" \
            --dependency="afterok:${DEP_LIST}" \
            "${JOBS_DIR}/wrappers/eval_legacy_ab.sh")
    else
        # No training submitted in this invocation — just run eval now
        EVAL_JOB=$(sbatch --parsable \
            --export=ALL,TF_NAME="${TF_NAME}" \
            "${JOBS_DIR}/wrappers/eval_legacy_ab.sh")
    fi
    echo "[legacy-ab] eval jobid=${EVAL_JOB}"
fi

echo ""
echo "[legacy-ab] Pipeline submitted for ${TF_NAME}."
echo "  Monitor: squeue -u \$USER"
echo "  Final results: ${OUTPUTS_DIR}/id_benchmark_${TF_NAME}_legacy_ab.json"