#!/bin/bash
# lib/common.sh — sourced by every job script.
# Sets up paths, conda init, and sanity-checks required env vars.
# Caller must have already exported any TF-specific vars before sourcing this.

set -eo pipefail

# -------------------- Paths --------------------
export PROJECT_ROOT="/project2/rohs_102/shewchuk"

# --- Code / config / logs: rooted in THIS repo (TF-conformation is now the
# authoritative pipeline; nothing is sourced from DeepPBS/run/jobs anymore).
# Self-locating: this file lives at <repo>/lib/common.sh.
export TFCONF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export LIB_DIR="${TFCONF_DIR}/lib"
export CONFIG_DIR="${TFCONF_DIR}/config"
export PILOTS_DIR="${CONFIG_DIR}/pilots"
export WRAPPERS_DIR="${TFCONF_DIR}/wrappers"
export LOGS_DIR="${TFCONF_DIR}/slurm_output"

# Source-chain library: all protein chains for every DeepPBS structure.
# (Was deeppbs_pdbs/monomer_chains — monomers only — before the source_chains
# refactor; stage1_bioemu now samples every chain here.)
export BIOEMU_RAW_ROOT="${TFCONF_DIR}/structures/source_chains"

# Stage 1 (BioEmu + HPacker) conformation library. stage1_bioemu writes one dir
# per (PDB, source-chain) here as <PDB>_chain<X>_conformations/, each holding
# topology.pdb, samples.xtc, samples_sidechain_rec.{pdb,xtc}. This is the single
# canonical Stage 1 output location that Stage 2 reads from (see STAGE1_DIR).
export STAGE1_OUT_ROOT="${TFCONF_DIR}/structures/stage1_bioemu_output"

# --- Data / outputs / training trees: still live in the DeepPBS data trees on
# the cluster (large; not part of this repo). These are inputs/outputs, not code.
export REPO_DIR="${PROJECT_ROOT}/DeepPBS"
export RUN_DIR="${REPO_DIR}/run"
export DATA_DIR="${PROJECT_ROOT}/DeepPBS_data"
export ORIG_ASSEMBLY_DIR="${DATA_DIR}/deeppbsmar24/data/assembly2024"
export ORIG_FOLDS_DIR="${REPO_DIR}/run/folds"
export CONFORMATIONS_DIR="${REPO_DIR}/data/conformations"
export FOLDS_AUG_DIR="${REPO_DIR}/run/folds_aug"

# Pipeline outputs now land in this repo's output/ tree (was DeepPBS_outputs).
#   stage6_train/  — training run dirs (Model.best.tar, etc.); built into the
#                    training configs by stage5's build_training_configs.py and
#                    scanned by stage7's eval. OUTPUTS_DIR keeps that name so the
#                    write side (training) and read side (eval) stay in sync.
#   stage7_eval/   — benchmark JSON + per-entry prediction npz.
export OUTPUT_ROOT="${TFCONF_DIR}/output"
export OUTPUTS_DIR="${OUTPUT_ROOT}/stage6_train"
export EVAL_OUT_DIR="${OUTPUT_ROOT}/stage7_eval"

mkdir -p "${LOGS_DIR}" "${CONFORMATIONS_DIR}" "${FOLDS_AUG_DIR}" \
         "${OUTPUTS_DIR}" "${EVAL_OUT_DIR}"

# -------------------- Conda init --------------------
# In non-interactive SLURM shells, conda init isn't sourced automatically.
# Source it explicitly so `conda activate <env>` works.
CONDA_PREFIX_PATH="/apps/conda/miniforge3/24.11.3"
if [ -f "${CONDA_PREFIX_PATH}/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1091
    source "${CONDA_PREFIX_PATH}/etc/profile.d/conda.sh"
else
    echo "WARNING: conda hook not found at ${CONDA_PREFIX_PATH}" >&2
fi

# Some BioEmu internals look for these even when conda is properly active
export CONDA_ROOT="/home1/${USER}/.conda"
export HPACKER_PYTHONBIN="${CONDA_ROOT}/envs/hpacker/bin/python"
export HPACKER_REPO_DIR="${CONDA_ROOT}/envs/hpacker"

# -------------------- TF-specific config loader --------------------
# Caller must export TF_NAME before sourcing this if they want pilot config loaded.
load_pilot_config() {
    local tf_name="${1:?usage: load_pilot_config <tf_name>}"
    local cfg="${PILOTS_DIR}/${tf_name}.sh"
    if [ ! -f "${cfg}" ]; then
        echo "ERROR: pilot config not found: ${cfg}" >&2
        echo "Available configs:" >&2
        ls "${PILOTS_DIR}/"*.sh 2>/dev/null | sed 's/^/  /' >&2
        exit 1
    fi
    # shellcheck disable=SC1090
    source "${cfg}"

    # Guard the required pilot vars up front so a malformed config fails with a
    # clear message instead of a confusing downstream path error (S2). Note:
    # BINDING_CHAIN is the source-chain filename letter (the Stage 1 dir tag,
    # e.g. "A"), distinct from PROTEIN_CHAIN which is the 0-based mdtraj chainid
    # into the .cif — see B2 / the Stage 2 sequence-match guard that asserts the
    # two indexings actually point at the same protein.
    require_var PDB_ID
    require_var FOLD
    require_var PWM_LABEL
    require_var BINDING_CHAIN

    # Set derived paths once config has been read.
    # LEGACY=1 toggles to a parallel pipeline that uses --ignore-metals in
    # Stage 3 and writes everything to *_legacy/ paths. Used for A/B-testing
    # whether the metal-cage restraints affect downstream training outcomes.
    local suffix=""
    if [ "${LEGACY:-0}" = "1" ]; then
        suffix="_legacy"
        export STAGE3_IGNORE_METALS=1
        echo "[common] LEGACY=1: using --ignore-metals in Stage 3 and writing to *${suffix}/ paths"
    else
        export STAGE3_IGNORE_METALS=0
    fi

    export WORK_DIR="${CONFORMATIONS_DIR}/${TF_NAME}"
    export STAGE0_DIR="${WORK_DIR}/stage0_raw"

    # Stage 1 output = the binding chain's conformation dir in the new per-chain
    # library (B1). stage1_bioemu writes samples_sidechain_rec.{pdb,xtc} here;
    # Stage 2 reads them. BIOEMU_DIR is the same dir (raw BioEmu output lives
    # alongside the HPacker reconstruction).
    export BIOEMU_DIR="${STAGE1_OUT_ROOT}/${PDB_ID}_chain${BINDING_CHAIN}_conformations"
    export STAGE1_DIR="${BIOEMU_DIR}"
    export STAGE1_RELAX_PREFIX="samples_sidechain_rec"

    # STAGE2_DIR intentionally has NO legacy suffix: docking is metal-independent,
    # so the legacy A/B (Stage 3 --ignore-metals) reuses the same docked frames (S5).
    export STAGE2_DIR="${WORK_DIR}/stage2_docked"
    # Stages 3+ diverge in legacy mode
    export STAGE3_DIR="${WORK_DIR}/stage3_min${suffix}"
    export STAGE4_DIR="${WORK_DIR}/stage4_npz${suffix}"
    export REF_CIF="${BIOEMU_RAW_ROOT}/${PDB_ID}_chains/${PDB_ID}.cif"
    export COMBINED_ASSEMBLY_DIR="${DATA_DIR}/combined_assembly${suffix}_${TF_NAME}"
    # Augmented fold also gets a legacy suffix
    export AUG_TRAIN_FOLD="${FOLDS_AUG_DIR}/train${FOLD}_aug${suffix}_${TF_NAME}.txt"
    # And the conditioning suffix used by training-config output dirs
    export CONDITION_NAME_SUFFIX="${suffix}"

    # Verify the reference .cif exists (the source_chains refactor should have
    # copied each <PDB>.cif into its <PDB>_chains/ dir) — S3.
    if [ ! -f "${REF_CIF}" ]; then
        echo "ERROR: reference CIF not found: ${REF_CIF}" >&2
        echo "Expected the source_chains refactor to place ${PDB_ID}.cif in its chains dir." >&2
        exit 1
    fi

    # NB: STAGE1_DIR is Stage 1's OUTPUT — do not pre-create it here, so a
    # missing Stage 1 run surfaces as a clear "samples_* not found" in Stage 2.
    mkdir -p "${WORK_DIR}" "${STAGE2_DIR}" "${STAGE3_DIR}" "${STAGE4_DIR}/output"

    echo "[common] Loaded pilot config: ${TF_NAME} (PDB ${PDB_ID}, PWM ${PWM_LABEL})${suffix:+ [LEGACY]}"
}

# Validate that a required variable is set; bail with a helpful message if not.
require_var() {
    local name="$1"
    if [ -z "${!name:-}" ]; then
        echo "ERROR: required variable ${name} is not set" >&2
        echo "Did you load a pilot config? source it from ${PILOTS_DIR}/" >&2
        exit 1
    fi
}

# Auto-load pilot config if TF_NAME is set in the environment
if [ -n "${TF_NAME:-}" ]; then
    load_pilot_config "${TF_NAME}"
fi