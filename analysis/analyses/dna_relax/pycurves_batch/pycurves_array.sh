#!/bin/bash
#SBATCH --job-name=pycurves_batch
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/pycurves_%A_%a.out
# One array task = one chunk of the worklist. pyCurves' XLA thread pool needs a
# compute node (login-node RLIMIT_NPROC crashes it); CPU-only JAX is sufficient.
set -eo pipefail
BASE="/project2/rohs_102/shewchuk/TF-conformation"
WL="${BASE}/analysis/analyses/dna_relax/pycurves_worklist.tsv"
PYC_ENV="/project2/rohs_102/shewchuk/conda/envs/pycurves"
CHUNK="${CHUNK:-20}"                     # rows per array task
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate "${PYC_ENV}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
start=$(( TASK * CHUNK + 1 )); end=$(( start + CHUNK - 1 ))
sed -n "${start},${end}p" "${WL}" | while IFS=$'\t' read -r tf cond state inpdb outpref; do
    [ -z "$inpdb" ] && continue
    mkdir -p "$(dirname "$outpref")"
    # Run BOTH axis conventions (flags identical across crystal/frozen/relaxed):
    #   legacy    -> global bend + curvature + minimization (headline induced-fit metric)
    #   curvesplus-> Curves+ local base-pair axis params (literature-comparable)
    # Grooves come from either. NOTE: --ends is NOT used -- it requires equal-length
    # strands and errors out on the many asymmetric/overhang duplexes here (dux4,
    # egr1, ...). Groove coverage is fine without it.
    # No --visualization here (aggregation batch); viz is a separate representatives run.
    for conv in legacy curvesplus; do
        out="${outpref}_${conv}.json"
        if [ -f "$out" ]; then echo "SKIP ${tf}/${cond}_${state} ${conv} (done)"; continue; fi
        echo "RUN ${tf}/${cond}_${state} ${conv} <- ${inpdb}"
        pycurves "$inpdb" --format json --output-file "$out" \
            --axis-convention "$conv" \
            && echo "  OK $out" || echo "  FAIL ${tf}/${cond}_${state} ${conv}"
    done
done
echo "task ${TASK} done (rows ${start}-${end})"
