#!/bin/bash
# run_pilot.sh — submit a complete TF pilot end-to-end with SLURM dependencies.
#
# Usage:
#   ./scripts/pipeline/run_pilot.sh <tf_name> [stage_start] [stage_end]
#   stages: 1=hpacker, 2=redock, 3=minimize, 3r=recover, 4=preprocess,
#           5=build_aug, 6=train, 7=eval
#   Defaults: stage_start=1, stage_end=7
#
# Examples:
#   ./scripts/pipeline/run_pilot.sh ndt80             # full pipeline
#   ./scripts/pipeline/run_pilot.sh ndt80 4 7         # skip Stages 1-3 (already done)
#   ./scripts/pipeline/run_pilot.sh ndt80 6 7         # just retrain + reeval
#   ARM=dnarelax ./scripts/pipeline/run_pilot.sh ndt80   # DNA-relaxation arm
#   CONF_SRC=bioemu2 ./scripts/pipeline/run_pilot.sh ets1 4 7
#                                                     # BioEmu2 complex frames;
#                                                     # stage them first with
#                                                     # bioemu2/stage_as_stage3.sh
#
# Dependency policy:
#   afterany — Stage 3 array → 3r only. Partial failure on minimization is
#              expected and recovery is designed to handle it; without this a
#              single task failure would halt the whole pipeline.
#   afterok  — everywhere else, INCLUDING the fnat gate (3g) → Stage 4. The gate
#              exits non-zero when no conformation clears FNAT_FLOOR, and that
#              MUST halt the pipeline so Stage 4 never trains on an empty set.
#              Recovery (3r) exits 0, so the gate's afterok on it still fires.

set -eo pipefail
# TF-conformation repo root (this launcher's own dir); export so the wrappers
# and stage scripts resolve it. Submit from here so log paths land in this
# repo's slurm_output/.
export TFCONF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${TFCONF_DIR}"
# shellcheck source=lib/common.sh
source "${TFCONF_DIR}/lib/common.sh"

TF_NAME="${1:?Usage: $0 <tf_name> [stage_start] [stage_end]}"
STAGE_START="${2:-1}"
STAGE_END="${3:-7}"

# ---------------------------------------------------------------------------
# Pre-flight gates. Both are cheap (<2s) and both exist because of real errors:
#
#  1. Config drift. config/pilots/*.sh are GENERATED from config/pilots.tsv.
#     If someone hand-edited a generated config, or edited the registry without
#     regenerating, the file on disk no longer matches the declared intent.
#  2. Label / subset correctness. The TBP pilot ran for months against NDT80's
#     motif (docs/tbp_label_error_downstream.md). validate_pilot.py cross-checks
#     PWM_LABEL against config/deeppbs_pdb_pwms.csv, catches test-set leakage,
#     empty on-target subsets, and pseudoreplicated ones.
#
# Set PILOT_SKIP_VALIDATE=1 to bypass (debugging only — never for a real run).
# ---------------------------------------------------------------------------
if [ "${PILOT_SKIP_VALIDATE:-0}" != "1" ]; then
    echo "[run_pilot] pre-flight: checking generated configs match config/pilots.tsv"
    if ! python3 "${TFCONF_DIR}/scripts/pipeline/gen_pilot_configs.py" --check; then
        echo "ERROR: pilot configs are out of date with config/pilots.tsv." >&2
        echo "       Run: python scripts/pipeline/gen_pilot_configs.py" >&2
        exit 1
    fi
    echo "[run_pilot] pre-flight: validating pilot ${TF_NAME}"
    if ! python3 "${TFCONF_DIR}/scripts/pipeline/validate_pilot.py" "${TF_NAME}" --quiet; then
        echo "ERROR: ${TF_NAME} failed validation (see above). Refusing to submit." >&2
        echo "       Override with PILOT_SKIP_VALIDATE=1 only if you know why." >&2
        exit 1
    fi
else
    echo "[run_pilot] WARNING: PILOT_SKIP_VALIDATE=1 — pre-flight checks skipped" >&2
fi

load_pilot_config "${TF_NAME}"
require_var TF_NAME
require_var PDB_ID
require_var N_FRAMES
require_var COMBINED_ASSEMBLY_DIR

echo "[run_pilot] arm=${ARM} (ARM=dnarelax for the DNA-relaxation variant)"
echo "[run_pilot] conformation source=${CONF_SRC:-bioemu1}"

# One timestamp for this whole submission, so every stage of this pilot run
# shares a log prefix and the files sort together:
#   slurm_output/20261015_143022_redock_123456.out
# log_args() lives in lib/common.sh; see the note there about SLURM having no
# date token and about this being SUBMIT time, not start time.
# Set TFCONF_LOG_TIMESTAMP=0 to fall back to the wrappers' un-stamped headers.
export SUBMIT_STAMP="${SUBMIT_STAMP:-$(log_stamp)}"
echo "[run_pilot] log prefix: ${SUBMIT_STAMP}_<stage>_<jobid>  (in ${LOGS_DIR})"

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

    # Timestamped log paths. Tag comes from the wrapper's own name
    # (stage2_redock.sh -> stage2_redock); array jobs need %A_%a not %j.
    # Passing --output on the command line overrides the wrapper's #SBATCH
    # header, which cannot carry a timestamp. See lib/common.sh:log_args.
    local tag
    tag=$(basename "${script}" .sh)
    # An array job needs %A_%a; a plain job needs %j. If we gave an array job
    # %j, every task would write to the SAME file and clobber the others.
    # Two ways a job can be an array, and both must be detected:
    #   1. the caller passes --array=   (stage3_array.sh gets 1-N%8 this way)
    #   2. the WRAPPER declares it in its own header
    #      (train_compare.sh has "#SBATCH --array=0-1": task 0 = baseline,
    #       task 1 = augmented, so %j would collide the two arms into one log)
    local arrflag=""
    case " $* " in *" --array="*) arrflag="--array" ;; esac
    if [ -z "${arrflag}" ] && grep -qE '^#SBATCH[[:space:]]+--array' "${script}"; then
        arrflag="--array"
    fi
    local logargs=()
    mapfile -t logargs < <(log_args_lines "${tag}" ${arrflag})

    local jobid
    # ARM is exported explicitly alongside TF_NAME. --export=ALL would carry it
    # anyway, but naming it makes the arm selection visible at the submit site
    # and survives a future change to a narrower --export list.
    jobid=$(sbatch --parsable \
        --export=ALL,TF_NAME="${TF_NAME}",ARM="${ARM}",CONF_SRC="${CONF_SRC}" \
        "${logargs[@]}" \
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

# Stage 1: BioEmu + HPacker conformation generation
if stage_in_range 1; then
    submit_stage 1 "${WRAPPERS_DIR}/stage1_bioemu.sh" ok
fi

# Stage 2: redock
if stage_in_range 2; then
    submit_stage 2 "${WRAPPERS_DIR}/stage2_redock.sh" ok
fi

# Stage 3: array of N_FRAMES tasks
if stage_in_range 3; then
    submit_stage 3 "${WRAPPERS_DIR}/stage3_array.sh" ok \
        --array="1-${N_FRAMES}%8"
fi

# Stage 3 recovery — uses afterany since the array is expected to have some
# failures (those are exactly what recovery is for).
if stage_in_range 3 && [ "${STAGE_END}" -ge 3 ]; then
    submit_stage "3r" "${WRAPPERS_DIR}/stage3_recover.sh" any
fi

# fnat gate (B7): the single structural-quality filter, Stage 3 -> Stage 4.
# Depends afterok on Stage 3 + its recover step (PREV_JOB = 3r, which exits 0),
# scores post-minimization fnat, and builds STAGE3_PASS_DIR. Runs whenever Stage
# 4 will run, so the pass dir is always fresh before preprocessing. If nothing
# clears FNAT_FLOOR the gate exits non-zero and the afterok edge halts Stage 4.
if stage_in_range 4; then
    submit_stage "3g" "${WRAPPERS_DIR}/fnat_gate.sh" ok
fi

# Stage 4: preprocess — afterok on the gate (NOT on recover), so an empty
# pass-list halts the pipeline. Stage 4 reads only STAGE3_PASS_DIR.
if stage_in_range 4; then
    submit_stage 4 "${WRAPPERS_DIR}/stage4_preprocess.sh" ok
fi

# Stage 5: build augmented fold (must succeed before training)
if stage_in_range 5; then
    submit_stage 5 "${WRAPPERS_DIR}/stage5_build_aug.sh" ok
fi

# Stage 6: training (array of 2)
if stage_in_range 6; then
    submit_stage 6 "${WRAPPERS_DIR}/train_compare.sh" ok
fi

# Stage 7: eval
if stage_in_range 7; then
    submit_stage 7 "${WRAPPERS_DIR}/eval_benchmark.sh" ok
fi

echo ""
echo "[run_pilot] Pipeline submitted. Monitor with:"
echo "  squeue -u \$USER"
echo "[run_pilot] Final results will be at:"
echo "  ${EVAL_OUT_DIR}/id_benchmark_${TF_NAME}${CONDITION_NAME_SUFFIX}.json"