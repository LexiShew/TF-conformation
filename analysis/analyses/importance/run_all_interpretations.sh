#!/bin/bash
set -e
R=/project2/rohs_102/shewchuk/TF-conformation

echo "=== Complete Attribution Workflow for All 12 Pilots ==="
echo "Start: $(date)"
echo ""
echo "This will:"
echo "  1. Submit GPU jobs for all 12 pilots (baseline vs augmented-frozen vs augmented-relaxed)"
echo "  2. Wait for all jobs to complete (12-24 hours)"
echo "  3. Compile cross-pilot comparison table"
echo ""

echo "Step 1: Submit attribution jobs..."
bash $R/analysis/analyses/importance/batch_interpret_all.sh

echo ""
echo "Step 2: Waiting for GPU jobs to complete..."
while [ $(squeue -u shewchuk 2>/dev/null | grep interpret | wc -l) -gt 0 ]; do
  N=$(squeue -u shewchuk | grep interpret | wc -l)
  echo "  $(date): $N jobs still running"
  sleep 120
done

echo ""
echo "Step 3: Compile results from all pilots"
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/deeppbs

cd $R
python analysis/analyses/importance/compare_importance_all.py

echo ""
echo "=== Complete ==="
echo "Master summary: $R/output/interpret_results_all/all_pilots_importance_summary.csv"
echo "Individual results: $R/output/interpret_results_all/<pilot>/"
echo "End: $(date)"
