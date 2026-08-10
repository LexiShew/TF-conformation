#!/bin/bash
#SBATCH --job-name=af3_pycurves
#SBATCH --partition=rohs,qcbgpu
#SBATCH --account=rohs_102
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/af3/slurm_output/af3pyc_%A_%a.out
set -eo pipefail
BASE=/project2/rohs_102/shewchuk/TF-conformation
PYC_ENV=/project2/rohs_102/shewchuk/conda/envs/pycurves
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate "$PYC_ENV"
cd "$BASE/af3/af3_dna"
# one array task per pilot dir
PILOTS=(egr1_1aay tbp_1tgh ets1_1k79 foxa_1vtn lef1_2lef engrailed_3hdd csl_3brg dux4_5z6z err_1lo1 nfat_1a66 runx_1hjc hsf_5d5u irf_1if1)
PIL=${PILOTS[$((SLURM_ARRAY_TASK_ID-1))]}
echo "[af3pyc] $PIL"
for pdb in "$PIL"/*_dna.pdb; do
    [ -f "$pdb" ] || continue
    out="${pdb%.pdb}_legacy.json"
    if [ -f "$out" ]; then echo "SKIP $pdb"; continue; fi
    pycurves "$pdb" --format json --output-file "$out" --axis-convention legacy \
      && echo "OK $out" || echo "FAIL $pdb"
done
echo "[af3pyc] $PIL done"
