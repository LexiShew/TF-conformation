#!/bin/bash
# fnat_gate/fnat_gate.sh — fnat data-integrity gate between Stage 2 and Stage 3 (B7).
#
# Scores every Stage 2 docked state (model protein vs the DNA carried in the same
# docked file, i.e. --use_model_dna) and moves states with fnat below FNAT_FLOOR
# into ${STAGE2_DIR}/rejected_fnat/, so only well-docked conformations flow into
# Stage 3 -> 5 and become training pairs. Produces a per-state fnat CSV + drop log.
#
# Knob (pilot config or environment):
#   FNAT_FLOOR  — minimum native-contact fraction to keep a state (default 0.5).
#
# Sourced by wrappers/fnat_gate.sh after common.sh + load_pilot_config.

require_var TF_NAME
require_var PDB_ID
require_var STAGE2_DIR
require_var REF_CIF

# interface_rmsd.py needs biopython (Bio.PDB) + numpy; both live in the deeppbs env.
conda activate "${DEEPPBS_ENV:-deeppbs}"

STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FNAT_FLOOR="${FNAT_FLOOR:-0.5}"
CSV="${STAGE2_DIR}/interface_metrics.csv"
DROPLOG="${STAGE2_DIR}/fnat_drops.log"

n_in=$(ls "${STAGE2_DIR}"/${PDB_ID}_state_*.pdb 2>/dev/null | wc -l)
if [ "${n_in}" -eq 0 ]; then
    echo "ERROR: no Stage 2 states in ${STAGE2_DIR}. Did Stage 2 run?" >&2
    exit 1
fi

echo "[fnat_gate/${TF_NAME}] scoring ${n_in} states in ${STAGE2_DIR} (floor=${FNAT_FLOOR})"
python "${STAGE_DIR}/interface_rmsd.py" \
    --ref "${REF_CIF}" \
    --states-dir "${STAGE2_DIR}" \
    --pdb-id "${PDB_ID}" \
    --use_model_dna \
    --floor "${FNAT_FLOOR}" \
    --out "${CSV}" \
    --fail "${DROPLOG}"

n_kept=$(ls "${STAGE2_DIR}"/${PDB_ID}_state_*.pdb 2>/dev/null | wc -l)
n_rej=$(ls "${STAGE2_DIR}"/rejected_fnat/${PDB_ID}_state_*.pdb 2>/dev/null | wc -l)
echo "[fnat_gate/${TF_NAME}] DONE — kept ${n_kept}, rejected ${n_rej} (floor ${FNAT_FLOOR})"
echo "[fnat_gate/${TF_NAME}] per-state metrics: ${CSV}"
echo "[fnat_gate/${TF_NAME}] drop log:          ${DROPLOG}"
