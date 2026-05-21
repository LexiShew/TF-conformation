#!/bin/bash
# lib/stage4_preprocess.sh — convert minimized Stage 3 PDBs to npz training
# files via DeepPBS's process_co_crystal.py. Labels with PWM_LABEL.

require_var TF_NAME
require_var PDB_ID
require_var STAGE3_DIR
require_var STAGE4_DIR
require_var PWM_LABEL

conda activate deeppbs

# Build input.txt and process_config.json
mkdir -p "${STAGE4_DIR}/pdb_input" "${STAGE4_DIR}/output"
cd "${STAGE4_DIR}/pdb_input"
ln -sf "${STAGE3_DIR}"/${PDB_ID}_state_*.pdb . 2>/dev/null || true
n_inputs=$(ls "${STAGE4_DIR}/pdb_input"/*.pdb 2>/dev/null | wc -l)
echo "[stage4/${TF_NAME}] Linked ${n_inputs} structures from Stage 3"

if [ "${n_inputs}" -eq 0 ]; then
    echo "ERROR: No Stage 3 outputs found. Did Stage 3 run successfully?" >&2
    exit 1
fi

cd "${STAGE4_DIR}"
> input.txt
for pdb in pdb_input/*.pdb; do
    base=$(basename "${pdb}")
    echo "${base},${PWM_LABEL}" >> input.txt
done
echo "[stage4/${TF_NAME}] Wrote input.txt with $(wc -l < input.txt) entries"

cat > process_config.json << EOF
{
    "PDB_FILES_PATH": "${STAGE4_DIR}/pdb_input",
    "FEATURE_DATA_PATH": "${STAGE4_DIR}/output"
}
EOF

# CRITICAL: 3DNA looks for parameter files via paths relative to the cwd.
# Running from anywhere other than run/process/ causes silent SIGFPE crashes
# in Curves. Always cd to run/process/ before running, even if input/output
# paths are absolute. (Discovered the hard way during EGR1 pilot.)
cd "${RUN_DIR}/process"
# shellcheck disable=SC1091
source ./proc_source.sh

# Clean any leftover 3DNA temp files from a previous run
rm -f *.pdb *.par *.pqr *.r3d *.dat *.log dna-dssr.json hstacking.pdb stacking.pdb 2>/dev/null

echo "[stage4/${TF_NAME}] Running process_co_crystal.py (this takes ~25-40 min serial)"
python ../process_co_crystal.py \
    "${STAGE4_DIR}/input.txt" \
    "${STAGE4_DIR}/process_config.json" \
    2>&1 | tee "${STAGE4_DIR}/process.log"

# Clean up 3DNA tmp files we just created
rm -f *.pdb *.par *.pqr *.r3d *.dat *.log dna-dssr.json hstacking.pdb stacking.pdb 2>/dev/null

n_outputs=$(ls "${STAGE4_DIR}/output"/*.npz 2>/dev/null | wc -l)
echo "[stage4/${TF_NAME}] DONE — ${n_outputs} npz files in ${STAGE4_DIR}/output"

# Sanity: how many CONTACT COUNT lines did the log show? (one per success)
n_contacts=$(grep -c "CONTACT COUNT" "${STAGE4_DIR}/process.log" 2>/dev/null || echo 0)
echo "[stage4/${TF_NAME}] CONTACT COUNT lines in log: ${n_contacts}"
if [ "${n_contacts}" -ne "${n_outputs}" ]; then
    echo "[stage4/${TF_NAME}] WARNING: contact count (${n_contacts}) != npz count (${n_outputs})"
fi
