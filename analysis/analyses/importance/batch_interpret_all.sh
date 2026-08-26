#!/bin/bash
R=/project2/rohs_102/shewchuk/TF-conformation
PILOTS="csl egr1 engrailed err ets1 foxa lef1 nfat runx tbp hsf irf"

echo "Submitting attribution jobs for all 12 pilots..."
for P in $PILOTS; do
  cat > /tmp/submit_${P}.sh << 'EOFJOB'
#!/bin/bash
#SBATCH --job-name=interp_PILOT
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --output=LOGDIR/interpret_PILOT_%j.log

source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/deeppbs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

cd PROJDIR
python analysis/analyses/importance/interpret_tfconf_all.py --pilot PILOT --out output/interpret_results_all
echo "Complete: PILOT"
EOFJOB
  
  sed -i "s|PILOT|${P}|g; s|PROJDIR|$R|g; s|LOGDIR|$R/slurm_output|g" /tmp/submit_${P}.sh
  JID=$(sbatch --parsable /tmp/submit_${P}.sh)
  echo "  $P: job $JID"
  sleep 1
done
echo ""
echo "All jobs submitted. Track: squeue -u shewchuk | grep interpret"
