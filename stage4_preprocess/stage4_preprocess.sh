#!/bin/bash
# lib/stage4_preprocess.sh — convert minimized Stage 3 PDBs to npz training
# files via DeepPBS's process_co_crystal.py. Labels with PWM_LABEL.

require_var TF_NAME
require_var PDB_ID
require_var STAGE3_DIR
require_var STAGE4_DIR
require_var PWM_LABEL

conda activate "${DEEPPBS_ENV:-deeppbs}"

# Self-contained: run this stage's co-located scripts (process_co_crystal.py,
# proc_source.sh) and the vendored 3DNA toolchain in ../lib, not DeepPBS/run.
STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

# 3DNA writes scratch files into the cwd, so run from a dedicated writable work
# dir. proc_source.sh sets X3DNA to an absolute path, so 3DNA finds its
# parameter files regardless of cwd — which is what previously forced running
# from run/process/ and caused silent SIGFPE crashes in Curves when run
# elsewhere (discovered the hard way during the EGR1 pilot).
WORK3DNA="${STAGE4_DIR}/3dna_work"
mkdir -p "${WORK3DNA}"
cd "${WORK3DNA}"
# shellcheck disable=SC1091
source "${STAGE_DIR}/proc_source.sh"

# Clean any leftover 3DNA temp files from a previous run
rm -f *.pdb *.par *.pqr *.r3d *.dat *.log dna-dssr.json hstacking.pdb stacking.pdb 2>/dev/null

echo "[stage4/${TF_NAME}] Running process_co_crystal.py (this takes ~25-40 min serial)"
python "${STAGE_DIR}/process_co_crystal.py" \
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
