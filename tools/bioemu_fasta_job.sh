#!/bin/bash
# tools/bioemu_fasta_job.sh — SLURM GPU job: BioEmu sampling + HPacker
# side-chain reconstruction for every *.fasta in a directory.
#
# Normally submitted by tools/bioemu_from_pdb.sh. Each <NAME>.fasta writes
#   <OUTPUT_DIR>/<NAME>_conformations/
#       topology.pdb, samples.xtc             (BioEmu backbone-only)
#       samples_sidechain_rec.pdb, .xtc       (HPacker full-atom)
# Chains whose samples_sidechain_rec.xtc already exists are skipped, so the
# job is resumable.
#
# Usage:
#   sbatch tools/bioemu_fasta_job.sh <FASTA_DIR> [NUM_CONFORMATIONS] [OUTPUT_DIR]

#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p rohs,qcbgpu
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_output/bioemu_from_pdb_%j_%x.out
#SBATCH --error=slurm_output/bioemu_from_pdb_%j_%x.err

set -eo pipefail

FASTA_DIR=$1
NUM_CONFORMATIONS=${2:-100}
OUTPUT_DIR=${3:-"$(dirname "${FASTA_DIR}")"}

if [ -z "${FASTA_DIR}" ] || [ ! -d "${FASTA_DIR}" ]; then
    echo "Usage: $0 <FASTA_DIR> [NUM_CONFORMATIONS] [OUTPUT_DIR]" >&2
    exit 1
fi

# Same environment setup as stage1_bioemu/stage1_bioemu.sh
source "${CONDA_PREFIX_PATH:-/apps/conda/miniforge3/24.11.3}/etc/profile.d/conda.sh"
conda activate "${BIOEMU_ENV:-bioemu}"
export CONDA_ROOT="${CONDA_ROOT:-/home1/${USER}/.conda}"

export BIOEMU_CACHE_DIR="${BIOEMU_CACHE_DIR:-/scratch1/shewchuk/.bioemu_embeds_cache}"
export HF_HOME="${HF_HOME:-/scratch1/shewchuk/.cache/huggingface}"
mkdir -p "${BIOEMU_CACHE_DIR}" "${HF_HOME}" "${OUTPUT_DIR}"

echo "[bioemu_fasta_job] fastas=${FASTA_DIR} num_conformations=${NUM_CONFORMATIONS}"
echo "[bioemu_fasta_job] output -> ${OUTPUT_DIR}"

python - "${FASTA_DIR}" "${NUM_CONFORMATIONS}" "${OUTPUT_DIR}" "${BIOEMU_CACHE_DIR}" <<'EOF'
import glob
import os
import sys

from bioemu.sample import main as sample
from bioemu.sidechain_relax import main as sidechain_relax

fasta_dir, num_confs, output_dir, cache_dir = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]

fastas = sorted(glob.glob(os.path.join(fasta_dir, "*.fasta")))
if not fastas:
    sys.exit(f"No *.fasta files in {fasta_dir}")

for fasta in fastas:
    name = os.path.splitext(os.path.basename(fasta))[0]
    out_dir = os.path.join(output_dir, f"{name}_conformations")
    if os.path.exists(os.path.join(out_dir, "samples_sidechain_rec.xtc")):
        print(f"SKIP {name}: samples_sidechain_rec.xtc exists")
        continue
    os.makedirs(out_dir, exist_ok=True)

    print(f"=== {name}: BioEmu x{num_confs} -> {out_dir}", flush=True)
    sample(
        sequence=fasta,
        num_samples=num_confs,
        output_dir=out_dir,
        cache_embeds_dir=cache_dir,
        filter_samples=False
    )
    for npz in glob.glob(os.path.join(out_dir, "*.npz")):
        os.remove(npz)

    print(f"=== {name}: HPacker side-chain reconstruction", flush=True)
    sidechain_relax(
        pdb_path=os.path.join(out_dir, "topology.pdb"),
        xtc_path=os.path.join(out_dir, "samples.xtc"),
        outpath=out_dir,
        md_equil=False,
    )
    if not os.path.exists(os.path.join(out_dir, "samples_sidechain_rec.xtc")):
        sys.exit(f"HPacker finished but samples_sidechain_rec.xtc missing in {out_dir}")
EOF

echo "[bioemu_fasta_job] DONE"
ls -d "${OUTPUT_DIR}"/*_conformations 2>/dev/null
