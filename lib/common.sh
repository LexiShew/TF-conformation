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
export STAGE1_OUT_ROOT="${TFCONF_DIR}/output/stage1_bioemu"

# --- Base INPUTS: now vendored into this repo's data/ tree (totally isolated;
# nothing is read from the DeepPBS trees anymore). Overridable via env.
#   data/folds/       — baseline train/valid splits + id.txt (small; in git)
#   data/assembly2024 — base DeepPBS training npz (large; gitignored, copied in
#                       on the cluster: see data/README.md)
export DATA_ROOT="${TFCONF_DIR}/data"
export ORIG_FOLDS_DIR="${ORIG_FOLDS_DIR:-${DATA_ROOT}/folds}"
export ORIG_ASSEMBLY_DIR="${ORIG_ASSEMBLY_DIR:-${DATA_ROOT}/assembly2024}"

# Make the vendored deeppbs package (lib/deeppbs) and the model definition
# (stage6_train/models) importable from the repo, so nothing is pulled from
# DeepPBS/run via sys.path. The conda env still supplies the dependencies.
export PYTHONPATH="${TFCONF_DIR}/lib:${TFCONF_DIR}/stage6_train:${PYTHONPATH:-}"

# ALL generated pipeline data lands in this repo's output/ tree, per stage — the
# repo is totally isolated (no writes into the DeepPBS trees):
#   stage2_docked/ stage3_min/ stage4_npz/  — per-TF intermediate work
#   stage5_aug/      — augmented folds (folds_aug/) + combined assembly
#   stage6_train/    — training run dirs + the generated training configs (the
#                      run-dir path is baked into each config and scanned by eval)
#   stage7_eval/     — benchmark JSON + per-entry prediction npz
export OUTPUT_ROOT="${TFCONF_DIR}/output"
export OUTPUTS_DIR="${OUTPUT_ROOT}/stage6_train"   # training run dirs + configs
export EVAL_OUT_DIR="${OUTPUT_ROOT}/stage7_eval"   # benchmark json + predict npz
export FOLDS_AUG_DIR="${OUTPUT_ROOT}/stage5_aug/folds_aug"

mkdir -p "${LOGS_DIR}" "${OUTPUT_ROOT}" "${OUTPUTS_DIR}" "${EVAL_OUT_DIR}" \
         "${FOLDS_AUG_DIR}"

# -------------------- Conda init --------------------
# In non-interactive SLURM shells, conda init isn't sourced automatically.
# Source it explicitly so `conda activate <env>` works. Both the conda
# installation prefix and the per-user envs root are overridable via env, so the
# pipeline ports to another machine without editing scripts.
export CONDA_PREFIX_PATH="${CONDA_PREFIX_PATH:-/apps/conda/miniforge3/24.11.3}"
if [ -f "${CONDA_PREFIX_PATH}/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1091
    source "${CONDA_PREFIX_PATH}/etc/profile.d/conda.sh"
else
    echo "WARNING: conda hook not found at ${CONDA_PREFIX_PATH}" >&2
fi

# Per-stage conda envs, defined once here and overridable via the environment
# (e.g. BIOEMU_ENV=my-bioemu). Stages `conda activate "${BIOEMU_ENV}"` etc.
#   BIOEMU_ENV  — stages 1-3: BioEmu sampling, HPacker relax, OpenMM minimize
#   DEEPPBS_ENV — stages 4-7 + fnat gate: DeepPBS, 3DNA, biopython
#   HPACKER_ENV — side-chain reconstruction, invoked by BioEmu (not activated)
export BIOEMU_ENV="${BIOEMU_ENV:-bioemu}"
export DEEPPBS_ENV="${DEEPPBS_ENV:-deeppbs}"
export HPACKER_ENV="${HPACKER_ENV:-hpacker}"

# Some BioEmu internals look for these even when conda is properly active
export CONDA_ROOT="${CONDA_ROOT:-/home1/${USER}/.conda}"
export HPACKER_PYTHONBIN="${CONDA_ROOT}/envs/${HPACKER_ENV}/bin/python"
export HPACKER_REPO_DIR="${CONDA_ROOT}/envs/${HPACKER_ENV}"

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

    # ---------------- Pilot defaults ----------------
    # These were byte-identical in all 12 hand-written pilot configs, so they
    # live here instead of being repeated per file. A pilot that genuinely needs
    # to differ sets the value in its config/pilots.tsv row, which the generator
    # emits into its config — and because these use ${VAR:-default}, a config
    # value always wins.
    export FOLD="${FOLD:-0}"
    export RAMP_STAGES="${RAMP_STAGES:-0.1,0.3,0.5,0.7,1.0}"
    export STEPS_PER_STAGE="${STEPS_PER_STAGE:-500}"
    export RECOVERY_RAMP_STAGES="${RECOVERY_RAMP_STAGES:-0.05,0.1,0.2,0.4,0.7,1.0}"
    export RECOVERY_STEPS_PER_STAGE="${RECOVERY_STEPS_PER_STAGE:-1000}"

    # ---------------- Arm selection ----------------
    # ARM replaces the 11 *_dnarelax.sh config files, every one of which carried
    # the same two-line delta. The downstream machinery is unchanged: the
    # _dnarelax output suffix is still keyed on STAGE3_DNA_RESTRAINT_K being set
    # (see below), so ARM just sets those two variables.
    #
    #   ARM=frozen    (default) DNA held at the crystal pose
    #   ARM=dnarelax  soft tether + late release (the headline relaxed mode)
    #
    # Usage: ARM=dnarelax ./scripts/pipeline/run_pilot.sh <tf>
    export ARM="${ARM:-frozen}"
    case "${ARM}" in
        frozen)
            : # nothing to set; STAGE3_DNA_RESTRAINT_K stays unset
            ;;
        dnarelax)
            # Soft tether (not fully free): a stiffness floor against GBSA-driven
            # B-DNA melting, and the on-ramp to the Tier-2 cgDNA+ stiffness prior.
            export STAGE3_DNA_RESTRAINT_K="${STAGE3_DNA_RESTRAINT_K:-1.5}"
            # Keep DNA pinned at the protein k through the 4 clash-resolving
            # vdW-ramp stages (largest forces); release only for the final ramp
            # stage plus the full-vdW final minimization.
            export STAGE3_DNA_RELEASE_STAGE="${STAGE3_DNA_RELEASE_STAGE:-5}"
            ;;
        *)
            echo "ERROR: ARM=${ARM} is not valid. Use 'frozen' or 'dnarelax'." >&2
            exit 1
            ;;
    esac

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

    # DNA relaxation (Tier 1) is a SEPARATE parallel pipeline, keyed on the pilot
    # config having opted in via STAGE3_DNA_RESTRAINT_K. It appends "_dnarelax" to
    # the suffix so relaxed Stage-3+ outputs never overwrite the frozen-DNA
    # baseline (which the 5-seed benchmark depends on). Like LEGACY, STAGE2_DIR is
    # NOT suffixed — docking is DNA-restraint-independent, so relaxed runs reuse
    # the same docked frames. Composes with LEGACY (e.g. _legacy_dnarelax).
    if [ -n "${STAGE3_DNA_RESTRAINT_K:-}" ]; then
        suffix="${suffix}_dnarelax"
        echo "[common] DNA relaxation ON (STAGE3_DNA_RESTRAINT_K=${STAGE3_DNA_RESTRAINT_K}${STAGE3_DNA_RELEASE_STAGE:+, release stage ${STAGE3_DNA_RELEASE_STAGE}}): writing Stage 3+ to *${suffix}/ paths"
    fi

    # Stage 1 output = the binding chain's conformation dir in the new per-chain
    # library (B1). stage1_bioemu writes samples_sidechain_rec.{pdb,xtc} here;
    # Stage 2 reads them. BIOEMU_DIR is the same dir (raw BioEmu output lives
    # alongside the HPacker reconstruction).
    export BIOEMU_DIR="${STAGE1_OUT_ROOT}/${PDB_ID}_chain${BINDING_CHAIN}_conformations"
    export STAGE1_DIR="${BIOEMU_DIR}"
    export STAGE1_RELAX_PREFIX="samples_sidechain_rec"

    # STAGE2_DIR intentionally has NO legacy suffix: docking is metal-independent,
    # so the legacy A/B (Stage 3 --ignore-metals) reuses the same docked frames (S5).
    export STAGE2_DIR="${OUTPUT_ROOT}/stage2_docked/${TF_NAME}"
    # Stages 3+ diverge in legacy mode
    export STAGE3_DIR="${OUTPUT_ROOT}/stage3_min${suffix}/${TF_NAME}"
    # fnat gate (B7): the single structural-quality filter, Stage 3 -> Stage 4.
    # The gate scores post-minimization fnat and builds STAGE3_PASS_DIR (a mirror
    # of STAGE3_DIR holding only states with fnat >= FNAT_FLOOR); Stage 4 reads
    # ONLY that pass dir. FNAT_FLOOR is overridable globally or per-pilot config.
    export FNAT_FLOOR="${FNAT_FLOOR:-0.5}"
    export STAGE3_PASS_DIR="${STAGE3_DIR}_pass"
    export STAGE4_DIR="${OUTPUT_ROOT}/stage4_npz${suffix}/${TF_NAME}"
    export REF_CIF="${BIOEMU_RAW_ROOT}/${PDB_ID}_chains/${PDB_ID}.cif"
    export COMBINED_ASSEMBLY_DIR="${OUTPUT_ROOT}/stage5_aug${suffix}/combined_assembly_${TF_NAME}"
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
    mkdir -p "${STAGE2_DIR}" "${STAGE3_DIR}" "${STAGE4_DIR}/output" \
             "${COMBINED_ASSEMBLY_DIR%/*}"

    echo "[common] Loaded pilot config: ${TF_NAME} (PDB ${PDB_ID}, PWM ${PWM_LABEL})${suffix:+ [mode:${suffix#_}]}"
}

# -------------------- Timestamped SLURM log paths --------------------
# SLURM has NO filename pattern for date/time. Its tokens are %j (job id),
# %A/%a (array job/task), %x (job name), %u (user), %N (node) -- and a
# "#SBATCH --output=" header is a static directive, so $(date) inside one is
# NOT expanded. The only way to get a timestamp into the filename is to compute
# it in the shell and pass --output/--error on the sbatch COMMAND LINE, which
# overrides the header in the script.
#
# log_args <tag> [--array]
#   Echoes the --output/--error pair for a job, prefixed with a sortable
#   timestamp:  20261015_143022_<tag>_<jobid>.out
#   Pass --array to use %A_%a instead of %j for array jobs.
#
# Usage:
#   sbatch $(log_args redock) --export=ALL,TF_NAME=csl wrappers/stage2_redock.sh
#   mapfile -t LOG < <(log_args_arr train --array)   # when word-splitting matters
#
# IMPORTANT -- the stamp is SUBMIT time, not start time. A job that queues for
# six hours still carries its submit timestamp. That is deliberate: every job in
# one dependency chain is submitted together, so they share a prefix and sort as
# a group. If you need start time instead, there is no sbatch-side way to get
# it; the job itself would have to rename its own log.
#
# Set TFCONF_LOG_TIMESTAMP=0 to disable and fall back to the wrappers' own
# un-stamped headers.
log_stamp() {
    date +%Y%m%d_%H%M%S
}

log_args() {
    local tag="${1:?usage: log_args <tag> [--array]}"
    local pat="%j"
    [ "${2:-}" = "--array" ] && pat="%A_%a"
    if [ "${TFCONF_LOG_TIMESTAMP:-1}" = "0" ]; then
        return 0          # emit nothing; wrapper headers apply
    fi
    # One stamp per *submission batch*: SUBMIT_STAMP is set once by a launcher so
    # every stage of one pilot run shares it. Falls back to now if unset.
    local stamp="${SUBMIT_STAMP:-$(log_stamp)}"
    printf -- '--output=%s/%s_%s_%s.out --error=%s/%s_%s_%s.err' \
        "${LOGS_DIR}" "${stamp}" "${tag}" "${pat}" \
        "${LOGS_DIR}" "${stamp}" "${tag}" "${pat}"
}

# Array-safe variant: prints one argument per line so callers can read it into
# an array without relying on word splitting.
log_args_lines() {
    local tag="${1:?usage: log_args_lines <tag> [--array]}"
    local pat="%j"
    [ "${2:-}" = "--array" ] && pat="%A_%a"
    [ "${TFCONF_LOG_TIMESTAMP:-1}" = "0" ] && return 0
    local stamp="${SUBMIT_STAMP:-$(log_stamp)}"
    printf -- '--output=%s/%s_%s_%s.out\n' "${LOGS_DIR}" "${stamp}" "${tag}" "${pat}"
    printf -- '--error=%s/%s_%s_%s.err\n'  "${LOGS_DIR}" "${stamp}" "${tag}" "${pat}"
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