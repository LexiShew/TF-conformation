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

conda activate bioemu

if [ ! -f "${REF_CIF}" ]; then
    echo "ERROR: Reference CIF not found: ${REF_CIF}" >&2
    echo "Expected at ${BIOEMU_RAW_ROOT}/${PDB_ID}_chains/${PDB_ID}.cif" >&2
    exit 1
fi

# Optionally inspect chains first if STAGE2_INSPECT=1 is set
if [ "${STAGE2_INSPECT:-0}" = "1" ]; then
    python "${SCRIPTS_DIR}/stage2_redock.py" \
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

echo "[stage2/${TF_NAME}] Re-docking ${PDB_ID} frames"
if [ -n "${MISMATCH_ACTION:-}" ]; then
    echo "[stage2/${TF_NAME}] Mismatch policy: action=${MISMATCH_ACTION} max=${MAX_MISMATCHES:-0}"
fi

python "${SCRIPTS_DIR}/stage2_redock.py" \
    --pdb-id "${PDB_ID}" \
    --ref "${REF_CIF}" \
    --traj "${STAGE1_DIR}/${PDB_ID}_sidechain_rec.xtc" \
    --top  "${STAGE1_DIR}/${PDB_ID}_sidechain_rec.pdb" \
    --out-dir "${STAGE2_DIR}" \
    --protein-chain "${PROTEIN_CHAIN}" \
    --dna-chains "${DNA_CHAINS}" \
    "${mismatch_args[@]}"

n_pdbs=$(ls "${STAGE2_DIR}"/${PDB_ID}_state_*.pdb 2>/dev/null | wc -l)
echo "[stage2/${TF_NAME}] DONE — wrote ${n_pdbs} structures to ${STAGE2_DIR}"