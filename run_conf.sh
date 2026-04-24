#!/bin/sh

# SLURM batch script for generating monomer conformations with BioEmu.
#
# Runs generate_monomer_confs.py for a single protein on one GPU node.
# Typically submitted by submit_all_monomers.sh, but can be run individually.
#
# Usage: sbatch --job-name="${PDB_ID}_${NUM_CONFS}" run_conf.sh <PDB_ID> <NUM_CONFORMATIONS>
#   PDB_ID             -- 4-letter PDB identifier (case-insensitive)
#   NUM_CONFORMATIONS  -- number of conformations to generate

#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --time=24:00:00
#SBATCH --output=conf_%j_%x.out
#SBATCH --error=conf_%j_%x.err

PDB_ID=$1
NUM_CONFORMATIONS=$2

# check if PDB ID is provided
if [ -z "$PDB_ID" ]; then
  echo "Usage: $0 <PDB_ID> <NUM_CONFORMATIONS>"
  exit 1
fi

python "${SLURM_SUBMIT_DIR}/generate_monomer_confs.py" $PDB_ID $NUM_CONFORMATIONS

