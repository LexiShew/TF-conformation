#!/bin/bash
# run_multiseed_pilot.sh — multi-seed paired comparison pilot.
#
# For a pilot TF (assuming Stages 1-5 have already been run), launches N
# paired-seed training pairs:
#   - baseline_<tf>_fold0_s1, augmented_<tf>_fold0_s1 (paired with seed=1)
#   - baseline_<tf>_fold0_s2, augmented_<tf>_fold0_s2 (paired with seed=2)
#   - ... up to s<N>
# Then runs the eval which auto-discovers all 2N checkpoints and produces
# paired statistics (bootstrap CIs and t-test p-values).
#
# Usage:
#   ./scripts/pipeline/run_multiseed_pilot.sh <tf_name> [n_seeds]
#
# Examples:
#   ./scripts/pipeline/run_multiseed_pilot.sh egr1 5      # 5 paired seeds for EGR1
#   ./scripts/pipeline/run_multiseed_pilot.sh tbp 5       # 5 paired seeds for TBP

set -eo pipefail
# TF-conformation repo root (this orchestrator's own dir); export so wrappers
# and stage scripts resolve it. TF-conformation is the authoritative pipeline.
export TFCONF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${TFCONF_DIR}"
# shellcheck source=lib/common.sh
source "${TFCONF_DIR}/lib/common.sh"

TF_NAME="${1:?Usage: $0 <tf_name> [n_seeds]}"
N_SEEDS="${2:-5}"

load_pilot_config "${TF_NAME}"
require_var TF_NAME
require_var FOLD
require_var COMBINED_ASSEMBLY_DIR

# Sanity check: data prep should already be done
if [ ! -d "${COMBINED_ASSEMBLY_DIR}" ]; then
    echo "ERROR: combined assembly dir not found: ${COMBINED_ASSEMBLY_DIR}" >&2
    echo "Run stages 1-5 first via ./scripts/pipeline/run_pilot.sh ${TF_NAME} 1 5" >&2
    exit 1
fi
AUG_TRAIN="${FOLDS_AUG_DIR}/train${FOLD}_aug_${TF_NAME}.txt"
if [ ! -f "${AUG_TRAIN}" ]; then
    echo "ERROR: augmented train file not found: ${AUG_TRAIN}" >&2
    exit 1
fi

echo "[multiseed/${TF_NAME}] Setup looks good"
echo "  combined dir:  ${COMBINED_ASSEMBLY_DIR}  ($(ls "${COMBINED_ASSEMBLY_DIR}" | wc -l) entries)"
echo "  aug train:     ${AUG_TRAIN}  ($(wc -l < "${AUG_TRAIN}") entries)"
echo "  n_seeds:       ${N_SEEDS}"

# -------------------- Step 1: Build N paired-seed configs --------------------
conda activate deeppbs
echo "[multiseed/${TF_NAME}] Building ${N_SEEDS} paired-seed configs"
for s in $(seq 1 "${N_SEEDS}"); do
    python "${TFCONF_DIR}/stage5_build_aug/build_training_configs.py" \
        --tf-name "${TF_NAME}" \
        --combined-dir "${COMBINED_ASSEMBLY_DIR}" \
        --fold "${FOLD}" \
        --seed "${s}" \
        --seed-suffix "_s${s}" \
        --output-dir "${RUN_DIR}"
done

# -------------------- Step 2: Submit training jobs --------------------
# 2N tasks: for each seed, one baseline + one augmented
echo "[multiseed/${TF_NAME}] Submitting ${N_SEEDS} paired training jobs"

declare -a TRAIN_JOB_IDS=()
for s in $(seq 1 "${N_SEEDS}"); do
    JOB_ID=$(sbatch --parsable \
        --job-name="train_${TF_NAME}_s${s}" \
        --output="${LOGS_DIR}/train_${TF_NAME}_s${s}_%A_%a.out" \
        --error="${LOGS_DIR}/train_${TF_NAME}_s${s}_%A_%a.err" \
        --export=ALL,TF_NAME="${TF_NAME}",SEED="${s}" \
        "${WRAPPERS_DIR}/train_compare.sh")
    TRAIN_JOB_IDS+=( "${JOB_ID}" )
    echo "  seed ${s}: train job ${JOB_ID}"
done

# -------------------- Step 3: Submit eval job depending on all training --------------------
DEPEND_STR=$(IFS=':'; echo "${TRAIN_JOB_IDS[*]}")
EVAL_JOB=$(sbatch --parsable \
    --dependency="afterok:${DEPEND_STR}" \
    --export=ALL,TF_NAME="${TF_NAME}" \
    "${WRAPPERS_DIR}/eval_benchmark.sh")
echo "[multiseed/${TF_NAME}] eval job ${EVAL_JOB} depends on training jobs"

echo ""
echo "[multiseed/${TF_NAME}] Pipeline submitted:"
echo "  ${N_SEEDS} paired training jobs: ${TRAIN_JOB_IDS[*]}"
echo "  1 eval job (auto-discovers all 2*${N_SEEDS} checkpoints): ${EVAL_JOB}"
echo ""
echo "Final results: ${EVAL_OUT_DIR}/id_benchmark_${TF_NAME}.json"
echo "Monitor with: squeue -u \$USER"
