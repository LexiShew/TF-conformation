#!/bin/bash
#SBATCH --job-name=be2msa
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/logs/bioemu2/be2msa_%j.out
#
# Serial ColabFold MSA prefetch. No GPU: this phase is pure network + cif parsing,
# and holding a GPU while sleeping out a 429 backoff is exactly what we are fixing.
# Submit with a partition/account on the command line, e.g.
#   sbatch --partition=qcb --account=qcb_640 bioemu2/warm_msa_cache.sbatch
set -eo pipefail
TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
cd "${TFCONF_DIR}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export BIOEMU2_MSA_CACHE="${BIOEMU2_MSA_CACHE:-${TFCONF_DIR}/output/bioemu2/_msa_cache}"
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate "${BIOEMU2_ENV:-bioemu2}"
python "${TFCONF_DIR}/bioemu2/warm_msa_cache.py" \
    --worklist "${WORKLIST:-${TFCONF_DIR}/bioemu2/worklist/source_chains_survey.tsv}" \
    --report   "${BIOEMU2_MSA_CACHE}/warm_report.tsv"
