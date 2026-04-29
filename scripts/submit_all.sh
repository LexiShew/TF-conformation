#!/bin/bash

# Submit SLURM jobs to generate monomer conformations for all proteins in a directory.
#
# Iterates over each <MONOMER_DIRECTORY>/<PDB>_chains/ subdirectory and submits
# a separate SLURM job (via run_conf.sh) for each. Each job runs
# generate_monomer_confs.py to produce BioEmu conformations and dock DNA when
# a *_dna.pdb is present.
#
# Output is written to <MONOMER_DIRECTORY>/<PDB>_chains/<PDB>_conformations/.
# SLURM logs go to <PROJECT_ROOT>/slurm_output/.
#
# Usage: ./submit_all.sh <MONOMER_DIRECTORY> [NUM_CONFORMATIONS]
#   MONOMER_DIRECTORY  -- directory containing <PDB>_chains/ subdirs
#   NUM_CONFORMATIONS  -- number of conformations per protein (default: 100)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

MONOMER_DIRECTORY=${1:?"Usage: $0 <MONOMER_DIRECTORY> [NUM_CONFORMATIONS]"}
NUM_CONFS=${2:-100}

LOG_DIR="$PROJECT_ROOT/slurm_output"
mkdir -p "$LOG_DIR"

# Point BioEmu and HuggingFace caches at /scratch (large, regeneratable).
# Override either by exporting the env var before invoking this script.
export BIOEMU_CACHE_DIR="${BIOEMU_CACHE_DIR:-/scratch1/shewchuk/.bioemu_embeds_cache}"
export HF_HOME="${HF_HOME:-/scratch1/shewchuk/.cache/huggingface}"
mkdir -p "$BIOEMU_CACHE_DIR" "$HF_HOME"

# Run sbatch from PROJECT_ROOT so SLURM_SUBMIT_DIR is consistent for run_conf.sh.
cd "$PROJECT_ROOT"

shopt -s nullglob
for chains_dir in "$MONOMER_DIRECTORY"/*_chains/; do
  chains_dir=${chains_dir%/}
  pdb_id=$(basename "$chains_dir" | cut -d'_' -f1)
  echo "Submitting: $pdb_id ($chains_dir) with $NUM_CONFS conformations"
  sbatch --job-name="${pdb_id}_${NUM_CONFS}" \
    --output="$LOG_DIR/conf_%j_%x.out" \
    --error="$LOG_DIR/conf_%j_%x.err" \
    "$SCRIPT_DIR/run_conf.sh" "$chains_dir" "$NUM_CONFS"
done
