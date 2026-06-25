#!/bin/bash
# stage1_bioemu/submit_stage1_bioemu.sh — fan stage1_bioemu.sh out over every
# structure in the source-chains library, one SLURM job per <PDB>_chains/ dir.
#
# Each job runs BioEmu + HPacker on all protein chains of that structure (see
# stage1_bioemu.sh). Output lands flat in output/stage1_bioemu/ as
# <PDB>_chain<X>_conformations/. SLURM logs go to <repo>/slurm_output/.
#
# Usage:
#   ./submit_stage1_bioemu.sh [CHAINS_ROOT] [NUM_CONFORMATIONS]
#     CHAINS_ROOT        -- dir containing <PDB>_chains/ subdirs
#                           (default: structures/source_chains)
#     NUM_CONFORMATIONS  -- conformations per chain (default: 100)

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "${SCRIPT_DIR}")"

CHAINS_ROOT=${1:-"${REPO_DIR}/structures/source_chains"}
NUM_CONFS=${2:-100}

if [ ! -d "${CHAINS_ROOT}" ]; then
    echo "ERROR: chains root not found: ${CHAINS_ROOT}" >&2
    exit 1
fi

LOG_DIR="${REPO_DIR}/slurm_output"
mkdir -p "${LOG_DIR}"

# Submit from the repo root so SLURM_SUBMIT_DIR resolves the engine path.
cd "${REPO_DIR}"

shopt -s nullglob
n=0
for chains_dir in "${CHAINS_ROOT}"/*_chains/; do
    chains_dir=${chains_dir%/}
    pdb_id=$(basename "${chains_dir}" | cut -d'_' -f1)
    echo "Submitting stage1_bioemu: ${pdb_id} (${chains_dir}) x${NUM_CONFS} confs"
    sbatch --job-name="s1_${pdb_id}" \
        --output="${LOG_DIR}/stage1_bioemu_%j_%x.out" \
        --error="${LOG_DIR}/stage1_bioemu_%j_%x.err" \
        "${SCRIPT_DIR}/stage1_bioemu.sh" "${chains_dir}" "${NUM_CONFS}"
    n=$((n + 1))
done
echo "Submitted ${n} stage1_bioemu jobs from ${CHAINS_ROOT}"
