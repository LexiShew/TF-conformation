#!/bin/bash

# Submit SLURM jobs to generate monomer conformations for all proteins in monomers/.
#
# Iterates over each monomers/<PDB>_chains/ subdirectory, finds *_chain*_protein.pdb
# files, and submits a separate SLURM job (via run_conf.sh) for each protein chain.
# Each job runs generate_monomer_confs.py to produce BioEmu conformations.
#
# Usage: ./submit_all.sh [NUM_CONFORMATIONS]
#   NUM_CONFORMATIONS -- number of conformations per protein (default: 100)

NUM_CONFS=${1:-100}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

for pdb_file in "$SCRIPT_DIR"/monomers/*_chains/*_chain*_protein.pdb; do
  [ -f "$pdb_file" ] || continue
  pdb_id=$(basename "$pdb_file" | cut -d'_' -f1)
  echo "Submitting: $pdb_id with $NUM_CONFS conformations"
  sbatch --job-name="${pdb_id}_${NUM_CONFS}" "$SCRIPT_DIR/run_conf.sh" "$pdb_id" "$NUM_CONFS"
done
