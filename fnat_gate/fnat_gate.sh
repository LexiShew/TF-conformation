#!/bin/bash
# fnat_gate/fnat_gate.sh — the pipeline's SINGLE structural-quality gate,
# Stage 3 (minimize) -> Stage 4 (preprocess) (B7).
#
# Scores every POST-MINIMIZATION Stage 3 state (model protein vs the DNA carried
# in the same minimized file, i.e. use_model_dna — the correct coordinate frame)
# and builds a pass-only mirror dir ${STAGE3_PASS_DIR} of symlinks to the states
# whose interface reproduces the native protein-DNA contacts (fnat >= FNAT_FLOOR).
# Stage 4 reads ONLY that pass dir, so sub-floor conformations never become
# DeepPBS training data.
#
# Why here and not at Stage 2: minimization moves per-state fnat in BOTH
# directions, so a docked pose just under the floor can be rescued (or a good
# one loosened). The gate must therefore score the minimized pose. Stage 2
# carries EVERY docked state forward; it never drops, moves, or quarantines.
# This is the only fnat filter in the DAG.
#
# Knob (pilot config or environment):
#   FNAT_FLOOR — minimum native-contact fraction to keep a state (default 0.5).
#
# Fails loud: if no state clears the floor the gate exits non-zero, so the DAG's
# afterok edge halts Stage 4 for that pilot.
#
# Sourced by wrappers/fnat_gate.sh after common.sh + load_pilot_config.

require_var TF_NAME
require_var PDB_ID
require_var STAGE3_DIR
require_var REF_CIF

# score_stage3.py / interface_rmsd.py need biopython (Bio.PDB) + numpy (deeppbs env).
conda activate "${DEEPPBS_ENV:-deeppbs}"

STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FNAT_FLOOR="${FNAT_FLOOR:-0.5}"
# STAGE3_PASS_DIR is exported by common.sh; fall back to the same convention.
PASS_DIR="${STAGE3_PASS_DIR:-${STAGE3_DIR}_pass}"
PASS_LIST="${STAGE3_DIR}/${PDB_ID}_pass.txt"
CSV="${STAGE3_DIR}/${PDB_ID}_fnat.csv"

n_in=$(ls "${STAGE3_DIR}"/${PDB_ID}_state_*.pdb 2>/dev/null | wc -l)
if [ "${n_in}" -eq 0 ]; then
    echo "ERROR: no Stage 3 states in ${STAGE3_DIR}. Did Stage 3 run?" >&2
    exit 1
fi

echo "[fnat_gate/${TF_NAME}] scoring ${n_in} post-min states in ${STAGE3_DIR} (floor=${FNAT_FLOOR})"
# score_stage3.py measures fnat vs each model's OWN DNA (use_model_dna=True,
# hard-wired there), writes <pdb>_fnat.csv, and with --floor writes the pass-list
# (one "<pdb>_state_NNN" per surviving state, no .pdb suffix).
python "${STAGE_DIR}/score_stage3.py" \
    --ref "${REF_CIF}" \
    --dir "${STAGE3_DIR}" \
    --pdb-id "${PDB_ID}" \
    --floor "${FNAT_FLOOR}" \
    --out "${CSV}" \
    --pass-out "${PASS_LIST}"

# Rebuild the pass-only mirror dir from scratch each run (idempotent). Guard the
# rm so it can only ever fire on a path ending in _pass — never on STAGE3_DIR.
case "${PASS_DIR}" in
    *_pass) : ;;
    *) echo "ERROR: refusing to rm pass dir not ending in _pass: ${PASS_DIR}" >&2; exit 1 ;;
esac
rm -rf "${PASS_DIR}"
mkdir -p "${PASS_DIR}"

# Symlink each passing state into the mirror dir under its identical filename, so
# Stage 4's glob/labels are unchanged — only the set of inputs shrinks.
n_pass=0
if [ -f "${PASS_LIST}" ]; then
    while IFS= read -r state; do
        [ -n "${state}" ] || continue
        src="${STAGE3_DIR}/${state}.pdb"
        if [ -f "${src}" ]; then
            ln -sf "${src}" "${PASS_DIR}/${state}.pdb"
            n_pass=$((n_pass + 1))
        else
            echo "[fnat_gate/${TF_NAME}] WARNING: pass-list names a missing state: ${src}" >&2
        fi
    done < "${PASS_LIST}"
fi

echo "[fnat_gate/${TF_NAME}] DONE — ${n_pass}/${n_in} states pass fnat>=${FNAT_FLOOR}"
echo "[fnat_gate/${TF_NAME}] pass dir:          ${PASS_DIR}"
echo "[fnat_gate/${TF_NAME}] pass list:         ${PASS_LIST}"
echo "[fnat_gate/${TF_NAME}] per-state metrics: ${CSV}"

# Fail loud: an empty pass-list must halt the DAG (Stage 4 depends afterok on
# this gate) so no pilot trains on zero qualifying conformations.
if [ "${n_pass}" -eq 0 ]; then
    echo "ERROR: no Stage 3 states cleared fnat floor ${FNAT_FLOOR}; halting before Stage 4." >&2
    exit 1
fi

export STAGE3_PASS_DIR="${PASS_DIR}"
