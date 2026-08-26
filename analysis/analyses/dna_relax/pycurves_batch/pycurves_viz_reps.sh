#!/bin/bash
#SBATCH --job-name=pyc_viz
#SBATCH -p qcb
#SBATCH -A rohs_102
#SBATCH -c 4
#SBATCH -t 00:30:00
#SBATCH -o /project2/rohs_102/shewchuk/TF-conformation/slurm_output/pyc_viz_%j.out
# Generate interactive HTML pyCurves viewers for the representative structures
# (median-bend state per TF per condition). Regenerates the JSON WITH
# --visualization (the aggregation batch omitted it to save space), then runs
# pycurves-viewer to emit standalone HTML.
#
# Reads a TSV worklist: <tf> <cond> <state> <inpdb>  (one rep per line)
set -eo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /project2/rohs_102/shewchuk/conda/envs/pycurves

BASE=/project2/rohs_102/shewchuk/TF-conformation
WL="$BASE/analysis/analyses/dna_relax/pycurves_viz/viz_worklist.tsv"
OUT="$BASE/analysis/analyses/dna_relax/pycurves_viz"
mkdir -p "$OUT"
cd "$OUT"

n_ok=0; n_fail=0
while IFS=$'\t' read -r tf cond state inpdb; do
    [ -z "$tf" ] && continue
    tag="${tf}_${cond}_${state}"
    json="$OUT/${tag}_curvesplus_viz.json"
    html="$OUT/${tag}.viewer.html"
    if [ ! -f "$inpdb" ]; then echo "MISS input $tag <- $inpdb"; n_fail=$((n_fail+1)); continue; fi
    echo "=== $tag ==="
    # 1) JSON with visualization geometry (curvesplus convention, matches batch)
    pycurves "$inpdb" --format json --output-file "$json" \
        --axis-convention curvesplus --visualization \
        && echo "  json OK" || { echo "  json FAIL $tag"; n_fail=$((n_fail+1)); continue; }
    # 2) standalone HTML viewer (point it at the actual structure)
    pycurves-viewer "$json" -o "$html" --structure "$inpdb" \
        && { echo "  html OK $html"; n_ok=$((n_ok+1)); } || { echo "  html FAIL $tag"; n_fail=$((n_fail+1)); }
done < "$WL"
echo "DONE: $n_ok html ok, $n_fail failed"
ls -la "$OUT"/*.viewer.html 2>/dev/null | wc -l
