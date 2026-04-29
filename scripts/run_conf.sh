#!/bin/sh

# SLURM batch script for generating monomer conformations with BioEmu.
#
# Runs generate_monomer_confs.py for a single <PDB>_chains/ directory on one GPU.
# Output is written to <CHAINS_DIR>/<PDB>_conformations/, with DNA docking
# applied when a *_dna.pdb is present in CHAINS_DIR.
#
# Typically submitted by submit_all.sh. To run directly, invoke from the
# project root so SLURM_SUBMIT_DIR points there:
#   sbatch scripts/run_conf.sh <CHAINS_DIR> <NUM_CONFORMATIONS>

#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_output/conf_%j_%x.out
#SBATCH --error=slurm_output/conf_%j_%x.err

CHAINS_DIR=$1
NUM_CONFORMATIONS=${2:-100}

if [ -z "$CHAINS_DIR" ]; then
  echo "Usage: $0 <CHAINS_DIR> [NUM_CONFORMATIONS]"
  exit 1
fi

python "${SLURM_SUBMIT_DIR}/scripts/generate_monomer_confs.py" \
  --chains-dir "$CHAINS_DIR" \
  -n "$NUM_CONFORMATIONS" \
  --reconstruct-sidechains
