#!/bin/bash
# add_crystal_to_training.sh — process a crystal PDB into an npz with a custom
# PWM label and add it to a training fold.
#
# Uses the chain-split crystal files (1tgh_chainA_protein.pdb + 1tgh_dna.pdb)
# from the BioEmu monomer_chains directory, NOT the BioEmu-derived
# _docked.pdb (which is backbone-only and breaks PDB2PQR).
#
# Usage:
#   ./add_crystal_to_training.sh <pdb_id> <chain> <pwm_label> [fold]
#
# Example:
#   ./add_crystal_to_training.sh 1tgh A MA0343.1.jaspar 0

set -eo pipefail

PDB_ID="${1:?Usage: $0 <pdb_id> <chain> <pwm_label> [fold]}"
CHAIN="${2:?Usage: $0 <pdb_id> <chain> <pwm_label> [fold]}"
PWM_LABEL="${3:?Usage: $0 <pdb_id> <chain> <pwm_label> [fold]}"
FOLD="${4:-0}"

# Paths
PROJECT_ROOT="/project2/rohs_102/shewchuk"
REPO_DIR="${PROJECT_ROOT}/DeepPBS"
RUN_DIR="${REPO_DIR}/run"
ASSEMBLY_DIR="${PROJECT_ROOT}/DeepPBS_data/deeppbsmar24/data/assembly2024"
FOLDS_DIR="${RUN_DIR}/folds"

CRYSTAL_DIR="${PROJECT_ROOT}/TF-conformation/deeppbs_pdbs/monomer_chains/${PDB_ID}_chains"
if [ ! -d "${CRYSTAL_DIR}" ]; then
    echo "ERROR: monomer_chains directory not found: ${CRYSTAL_DIR}" >&2
    exit 1
fi

PROT_PDB="${CRYSTAL_DIR}/${PDB_ID}_chain${CHAIN}_protein.pdb"
DNA_PDB="${CRYSTAL_DIR}/${PDB_ID}_dna.pdb"

if [ ! -f "${PROT_PDB}" ]; then
    echo "ERROR: protein chain PDB not found: ${PROT_PDB}" >&2
    echo "Available chain files:" >&2
    ls "${CRYSTAL_DIR}"/*chain*protein*.pdb 2>/dev/null >&2 || true
    exit 1
fi
if [ ! -f "${DNA_PDB}" ]; then
    echo "ERROR: DNA PDB not found: ${DNA_PDB}" >&2
    exit 1
fi

# Set up a tmp processing dir
TMP_DIR="${PROJECT_ROOT}/DeepPBS/data/add_crystal/${PDB_ID}_${PWM_LABEL}"
rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}/pdb_input" "${TMP_DIR}/output"

# Concatenate protein + DNA chains into one PDB
echo "[add_crystal] Combining protein + DNA from chain-split files"
echo "  protein: ${PROT_PDB} ($(grep -c "^ATOM" "${PROT_PDB}") atoms)"
echo "  DNA:     ${DNA_PDB} ($(grep -c "^ATOM" "${DNA_PDB}") atoms)"

# Keep only ATOM records (drop CRYST1, HEADER, etc.); separate protein and DNA
# with a TER record so process_co_crystal sees distinct chains.
(
    grep "^ATOM" "${PROT_PDB}"
    echo "TER"
    grep "^ATOM" "${DNA_PDB}"
    echo "END"
) > "${TMP_DIR}/pdb_input/${PDB_ID}.pdb"

n_atoms=$(grep -c "^ATOM" "${TMP_DIR}/pdb_input/${PDB_ID}.pdb")
echo "[add_crystal] Combined PDB has ${n_atoms} atoms"

# Build input file
echo "${PDB_ID}.pdb,${PWM_LABEL}" > "${TMP_DIR}/input.txt"

cat > "${TMP_DIR}/process_config.json" << EOF
{
    "PDB_FILES_PATH": "${TMP_DIR}/pdb_input",
    "FEATURE_DATA_PATH": "${TMP_DIR}/output"
}
EOF

# Run process_co_crystal — must be from run/process/ for 3DNA paths
cd "${RUN_DIR}/process"
# shellcheck source=/dev/null
source ./proc_source.sh
rm -f *.pdb *.par *.pqr *.r3d *.dat *.log dna-dssr.json hstacking.pdb stacking.pdb 2>/dev/null

echo "[add_crystal] Running process_co_crystal.py"
python ../process_co_crystal.py \
    "${TMP_DIR}/input.txt" \
    "${TMP_DIR}/process_config.json" \
    2>&1 | tee "${TMP_DIR}/process.log"

rm -f *.pdb *.par *.pqr *.r3d *.dat *.log dna-dssr.json hstacking.pdb stacking.pdb 2>/dev/null

# Locate the output npz
EXPECTED="${PDB_ID}_${CHAIN}_${PWM_LABEL}.npz"
SRC="${TMP_DIR}/output/${EXPECTED}"
if [ ! -f "${SRC}" ]; then
    echo "ERROR: expected output not found: ${SRC}" >&2
    echo "Files in output dir:" >&2
    ls "${TMP_DIR}/output/" >&2 || true
    exit 1
fi

# Copy to assembly dir
DST="${ASSEMBLY_DIR}/${EXPECTED}"
if [ -f "${DST}" ]; then
    echo "WARNING: ${DST} already exists; will not overwrite"
else
    cp "${SRC}" "${DST}"
    echo "[add_crystal] Copied to ${DST}"
fi

# Add to training fold (with trailing-newline-safe append)
TRAIN_FILE="${FOLDS_DIR}/train${FOLD}.txt"
if grep -q "^${EXPECTED}$" "${TRAIN_FILE}"; then
    echo "[add_crystal] ${EXPECTED} already in ${TRAIN_FILE}"
else
    if [ -n "$(tail -c 1 "${TRAIN_FILE}")" ]; then
        echo "" >> "${TRAIN_FILE}"
    fi
    echo "${EXPECTED}" >> "${TRAIN_FILE}"
    echo "[add_crystal] Added ${EXPECTED} to ${TRAIN_FILE}"
fi

echo "[add_crystal] DONE"
echo "  npz:        ${DST}"
echo "  train file: ${TRAIN_FILE} (now has $(wc -l < "${TRAIN_FILE}") entries)"