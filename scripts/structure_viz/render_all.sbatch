#!/bin/bash
#SBATCH --job-name=viz_ensemble
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --output=slurm_output/viz_ensemble_%j.out
# Render the crystal-vs-ensemble overlay for every pilot that has Stage-3 output.
# Usage: sbatch scripts/structure_viz/render_all.sbatch [n_states] [seed]
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
NSTATES="${1:-10}"
SEED="${2:-42}"
for TF in tbp egr1 foxa lef1 engrailed ets1; do
    if ls output/stage3_min/${TF}/*_state_*.pdb >/dev/null 2>&1; then
        echo "=== rendering $TF ==="
        bash scripts/structure_viz/render_ensemble.sh "$TF" "$NSTATES" "$SEED" || echo "  $TF FAILED"
    else
        echo "=== skip $TF (no Stage-3 output) ==="
    fi
done
echo "ALL DONE"
