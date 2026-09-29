#!/bin/bash
#SBATCH --job-name=hsf_crystal_feat
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/hsf_crystal_feat_%j.out
#SBATCH --error=/project2/rohs_102/shewchuk/TF-conformation/slurm_output/hsf_crystal_feat_%j.err
set -eo pipefail
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/deeppbs
R=/project2/rohs_102/shewchuk/TF-conformation
W=$R/output/stage4_npz/hsf_crystal
mkdir -p $W/pdb_input $W/output $W/3dna_work
cp $R/structures/source_chains/5d5u_chains/5d5u_duplex.cif $W/pdb_input/5d5u_duplex.cif
cat > $W/input.txt <<EOF
5d5u_duplex.cif,MA0486.2.jaspar
EOF
cat > $W/process_config.json <<EOF
{
    "PDB_FILES_PATH": "$W/pdb_input",
    "FEATURE_DATA_PATH": "$W/output"
}
EOF
cd $W/3dna_work
source $R/stage4_preprocess/proc_source.sh
rm -f *.pdb *.par *.pqr *.r3d *.dat *.log dna-dssr.json hstacking.pdb stacking.pdb 2>/dev/null || true
python $R/stage4_preprocess/process_co_crystal.py $W/input.txt $W/process_config.json 2>&1 | tee $W/process.log
rm -f *.pdb *.par *.pqr *.r3d *.dat *.log dna-dssr.json hstacking.pdb stacking.pdb 2>/dev/null || true
echo "=== outputs ==="
ls -la $W/output/
