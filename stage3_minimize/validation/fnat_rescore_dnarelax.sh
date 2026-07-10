#!/bin/bash
#SBATCH --job-name=fnat_rescore_dnarelax
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/fnat_rescore_dnarelax_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/fnat_rescore_dnarelax_%j.err
#
# DIAGNOSTIC (brief step 3): score the relaxed Stage-3 states through the fnat
# metric and compare, state-by-state, to the frozen-DNA baseline. This does NOT
# filter, threshold, or build a pass dir — it only writes the per-state fnat CSV
# and prints the paired Δfnat, so we can inspect the distribution before deciding
# whether the FNAT_FLOOR needs re-tuning when DNA is allowed to move.
set -eo pipefail
export TFCONF_DIR=/project2/rohs_102/shewchuk/TF-conformation
cd "$TFCONF_DIR"
source lib/common.sh
TF_NAME=tbp_dnarelax
load_pilot_config tbp_dnarelax   # sets STAGE3_DIR=stage3_min_dnarelax/tbp, REF_CIF, PDB_ID

FROZEN_DIR="${OUTPUT_ROOT}/stage3_min/tbp"    # frozen-DNA baseline (comparison)
OUT_CSV="${STAGE3_DIR}/${PDB_ID}_fnat_dnarelax.csv"

conda activate "${DEEPPBS_ENV:-deeppbs}"

echo "=== fnat re-score: relaxed vs frozen baseline ==="
echo "relaxed dir : ${STAGE3_DIR}"
echo "frozen dir  : ${FROZEN_DIR}"
echo "ref         : ${REF_CIF}"
echo "floor (diag only, NOT applied): ${FNAT_FLOOR}"

# No --floor / --pass-out: this is diagnostic-only. Omitting --floor means
# score_stage3.py writes the per-state fnat CSV and the paired Δfnat comparison
# but does NOT compute or write any threshold-derived pass-list. (Its printed
# summary still reports pct>=0.5 as a reference; that is display only.)
python fnat_gate/score_stage3.py \
    --ref "${REF_CIF}" \
    --dir "${STAGE3_DIR}" \
    --pdb-id "${PDB_ID}" \
    --out "${OUT_CSV}" \
    --compare-dir "${FROZEN_DIR}"

echo "=== CSV written: ${OUT_CSV} ==="
echo "=== DONE ==="
