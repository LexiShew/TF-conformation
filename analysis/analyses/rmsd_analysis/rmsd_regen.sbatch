#!/bin/bash
#SBATCH --job-name=rmsd_regen
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=00:45:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/rmsd_analysis/rmsd_regen_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/rmsd_analysis/rmsd_regen_%j.err
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate pycurves
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /project2/rohs_102/shewchuk/TF-conformation/rmsd_analysis
echo "host=$(hostname) start=$(date)"
python compute_rmsds.py \
    --tfs csl dux4 egr1 engrailed err ets1 foxa hsf irf lef1 nfat runx tbp \
    --interface-dir interface_residues \
    --output-dir . 2>&1
echo "exit_pipe=${PIPESTATUS[0]} end=$(date)"
echo "=== coverage after ==="
tail -n +2 per_state_rmsds.csv | cut -d, -f1 | sort | uniq -c
