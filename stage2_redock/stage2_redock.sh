#!/bin/bash
# lib/stage2_redock.sh — Cα-align Stage 1 frames + carry reference DNA + metals across.
#
# Optional config knobs (set in config/pilots/<tf>.sh if needed):
#   MISMATCH_ACTION  — what to do if BioEmu/reference residue sequences differ:
#                      'fail' (default, safe), 'warn' (proceed with risk),
#                      'trim' (slice both to longest common-sequence run).
#                      Use 'trim' when a reference crystal has disordered loops
#                      that BioEmu modeled across (e.g. DUX4 5z6z).
#   MAX_MISMATCHES   — tolerance threshold before mismatch-action triggers
#                      (default 0; raise to e.g. 5 to allow occasional HIS-variant
#                      label noise).

require_var TF_NAME
require_var PDB_ID
require_var STAGE1_DIR
require_var STAGE2_DIR
require_var REF_CIF
require_var PROTEIN_CHAIN
require_var DNA_CHAINS

# This script is meant to be `source`d by wrappers/stage2_redock.sh (it uses
# `return` for --inspect-only and inherits the pilot env). Guard against running
# it directly, where load_pilot_config never ran and `return` would error (S6).
if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
    echo "ERROR: source this from wrappers/stage2_redock.sh (after common.sh +" >&2
    echo "       load_pilot_config); do not run it directly." >&2
    exit 1
fi

conda activate "${BIOEMU_ENV:-bioemu}"

# Resolve this stage's own directory so we run the co-located stage2_redock.py
# (self-contained — no dependency on a shared SCRIPTS_DIR).
STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Stage 1 (BioEmu+HPacker) output: the HPacker-reconstructed full-atom ensemble
# for the binding chain (B1). Filenames are samples_sidechain_rec.{pdb,xtc};
# STAGE1_DIR already encodes ${PDB_ID}_chain${BINDING_CHAIN}_conformations.
TRAJ="${STAGE1_DIR}/${STAGE1_RELAX_PREFIX:-samples_sidechain_rec}.xtc"
TOP="${STAGE1_DIR}/${STAGE1_RELAX_PREFIX:-samples_sidechain_rec}.pdb"
if [ ! -f "${TRAJ}" ] || [ ! -f "${TOP}" ]; then
    echo "ERROR: Stage 1 outputs not found:" >&2
    echo "  ${TRAJ}" >&2
    echo "  ${TOP}" >&2
    echo "Did Stage 1 (stage1_bioemu) run for ${PDB_ID} chain ${BINDING_CHAIN}?" >&2
    exit 1
fi

if [ ! -f "${REF_CIF}" ]; then
    echo "ERROR: Reference CIF not found: ${REF_CIF}" >&2
    echo "Expected at ${BIOEMU_RAW_ROOT}/${PDB_ID}_chains/${PDB_ID}.cif" >&2
    exit 1
fi

# Optionally inspect chains first if STAGE2_INSPECT=1 is set
if [ "${STAGE2_INSPECT:-0}" = "1" ]; then
    python "${STAGE_DIR}/stage2_redock.py" \
        --pdb-id "${PDB_ID}" \
        --ref "${REF_CIF}" \
        --inspect-only
    return 0  # if sourced; caller decides whether to continue
fi

# Build mismatch-action args from optional env vars
mismatch_args=()
if [ -n "${MISMATCH_ACTION:-}" ]; then
    mismatch_args+=( --mismatch-action "${MISMATCH_ACTION}" )
fi
if [ -n "${MAX_MISMATCHES:-}" ]; then
    mismatch_args+=( --max-mismatches "${MAX_MISMATCHES}" )
fi

# Alignment-mode passthrough (B6). Default is the engine default (interface);
# set STAGE2_ALIGN_MODE=all|per_domain in the pilot config or environment to run
# a baseline / multidomain diagnostic.
align_args=()
if [ -n "${STAGE2_ALIGN_MODE:-}" ]; then
    align_args+=( --align-mode "${STAGE2_ALIGN_MODE}" )
fi

echo "[stage2/${TF_NAME}] Re-docking ${PDB_ID} frames (chain ${BINDING_CHAIN})"
if [ -n "${MISMATCH_ACTION:-}" ]; then
    echo "[stage2/${TF_NAME}] Mismatch policy: action=${MISMATCH_ACTION} max=${MAX_MISMATCHES:-0}"
fi
[ -n "${STAGE2_ALIGN_MODE:-}" ] && echo "[stage2/${TF_NAME}] Align mode: ${STAGE2_ALIGN_MODE}"

python "${STAGE_DIR}/stage2_redock.py" \
    --pdb-id "${PDB_ID}" \
    --ref "${REF_CIF}" \
    --traj "${TRAJ}" \
    --top  "${TOP}" \
    --out-dir "${STAGE2_DIR}" \
    --protein-chain "${PROTEIN_CHAIN}" \
    --dna-chains "${DNA_CHAINS}" \
    "${mismatch_args[@]}" \
    "${align_args[@]}"

n_pdbs=$(ls "${STAGE2_DIR}"/${PDB_ID}_state_*.pdb 2>/dev/null | wc -l)
echo "[stage2/${TF_NAME}] DONE — wrote ${n_pdbs} structures to ${STAGE2_DIR}"