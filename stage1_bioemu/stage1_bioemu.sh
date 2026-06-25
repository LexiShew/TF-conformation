#!/bin/bash
# stage1_bioemu/stage1_bioemu.sh — BioEmu sampling + HPacker side-chain
# reconstruction for ONE <PDB>_chains/ directory, across ALL its protein chains.
#
# Supersedes the old stage1_hpacker, which only ran HPacker on pre-existing
# BioEmu samples for the single monomer chain. This stage does both the
# sampling AND the HPacker reconstruction, for every *_chain*_protein.pdb in
# the structure (multi-chain crystal complexes get every chain sampled).
#
# Reads:
#   <CHAINS_DIR>/<PDB>_chain<X>_protein.pdb   (one or more protein chains)
# Writes (per chain):
#   <CHAINS_DIR>/<PDB>_chain<X>_conformations/
#       topology.pdb, samples.xtc                 (BioEmu backbone-only)
#       samples_sidechain_rec.pdb, .xtc           (HPacker full-atom)
#
# DNA docking is NOT done here — that's stage2_redock, per complex.
#
# Usage (single structure, direct):
#   sbatch stage1_bioemu/stage1_bioemu.sh <CHAINS_DIR> [NUM_CONFORMATIONS]
# Usage (whole dataset): see submit_stage1_bioemu.sh, which fans this out over
#   every structures/source_chains/<PDB>_chains/ dir.

#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_output/stage1_bioemu_%j_%x.out
#SBATCH --error=slurm_output/stage1_bioemu_%j_%x.err

set -eo pipefail

CHAINS_DIR=$1
NUM_CONFORMATIONS=${2:-100}

if [ -z "${CHAINS_DIR}" ]; then
    echo "Usage: $0 <CHAINS_DIR> [NUM_CONFORMATIONS]" >&2
    echo "  CHAINS_DIR  -- a structures/source_chains/<PDB>_chains/ directory" >&2
    exit 1
fi
if [ ! -d "${CHAINS_DIR}" ]; then
    echo "ERROR: not a directory: ${CHAINS_DIR}" >&2
    exit 1
fi

# Activate conda in the (non-interactive) SLURM environment. Paths inherit from
# common.sh when run via the wrapper; inline defaults keep standalone runs working.
source "${CONDA_PREFIX_PATH:-/apps/conda/miniforge3/24.11.3}/etc/profile.d/conda.sh"
conda activate "${BIOEMU_ENV:-bioemu}"
export CONDA_ROOT="${CONDA_ROOT:-/home1/${USER}/.conda}"

# BioEmu / HuggingFace caches on scratch (large, regenerable). Override by
# exporting these before submitting.
export BIOEMU_CACHE_DIR="${BIOEMU_CACHE_DIR:-/scratch1/shewchuk/.bioemu_embeds_cache}"
export HF_HOME="${HF_HOME:-/scratch1/shewchuk/.cache/huggingface}"
mkdir -p "${BIOEMU_CACHE_DIR}" "${HF_HOME}"

# Repo root: prefer TFCONF_DIR (set by common.sh / the stage1 wrapper), else
# SLURM_SUBMIT_DIR (submit from repo root), else this script's parent dir.
REPO_DIR="${TFCONF_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}}"
GEN_SCRIPT="${REPO_DIR}/stage1_bioemu/generate_monomer_confs.py"

# All chains from all structures land flat in one output dir; each chain's
# <PDB>_chain<X>_conformations/ name is globally unique, so there are no
# collisions. Override with STAGE1_OUTPUT_DIR if needed.
OUTPUT_DIR="${STAGE1_OUTPUT_DIR:-${REPO_DIR}/structures/stage1_bioemu_output}"
mkdir -p "${OUTPUT_DIR}"

pdb_id=$(basename "${CHAINS_DIR}" | cut -d'_' -f1)
echo "[stage1_bioemu] ${pdb_id}: BioEmu + HPacker on all chains in ${CHAINS_DIR}"
echo "[stage1_bioemu] num_conformations=${NUM_CONFORMATIONS} cache=${BIOEMU_CACHE_DIR}"
echo "[stage1_bioemu] output -> ${OUTPUT_DIR}"

python "${GEN_SCRIPT}" \
    --chains-dir "${CHAINS_DIR}" \
    --all-chains \
    --reconstruct-sidechains \
    -n "${NUM_CONFORMATIONS}" \
    --cache-dir "${BIOEMU_CACHE_DIR}" \
    --output-dir "${OUTPUT_DIR}"

echo "[stage1_bioemu] ${pdb_id} DONE"
ls -d "${OUTPUT_DIR}"/${pdb_id}_chain*_conformations 2>/dev/null
