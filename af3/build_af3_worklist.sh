#!/bin/bash
# Build a pyCurves worklist for the AlphaFold3 DNA predictions, then run it with the
# shared toolkit (scripts/pycurves/). Replaces the old standalone af3_pycurves.sh.
#
# Columns (tab-sep): TF  COND  STATE  INPUT_PDB  OUT_PREFIX   (COND is always "af3").
# Outputs land next to each PDB (e.g. .../seed-1_sample-0_dna_legacy.json), unchanged.
#
# Usage:
#   bash af3/build_af3_worklist.sh
#   WL=af3/af3_pycurves_worklist.tsv
#   N=$(( ($(wc -l < "$WL") + 19) / 20 ))            # CHUNK=20 rows/array task
#   CONVENTIONS=legacy sbatch --array=0-$((N-1)) scripts/pycurves/pycurves_array.sh "$WL"
# (CONVENTIONS=legacy preserves af3's legacy-only output; drop it for legacy+curvesplus.)
set -eo pipefail
BASE=/project2/rohs_102/shewchuk/TF-conformation
WL="${BASE}/af3/af3_pycurves_worklist.tsv"
: > "${WL}"
for pdb in "${BASE}"/af3/af3_dna/*/*_dna.pdb; do
    [ -f "$pdb" ] || continue
    dir=$(basename "$(dirname "$pdb")")          # <tf>_<pdbid>, e.g. egr1_1aay
    tf=${dir%%_*}                                 # egr1
    state=$(basename "$pdb" _dna.pdb)             # seed-1_sample-0
    printf "%s\t%s\t%s\t%s\t%s\n" "$tf" af3 "$state" "$pdb" "${pdb%.pdb}" >> "${WL}"
done
echo "worklist: ${WL} ($(wc -l < "${WL}") rows)"
cut -f1 "${WL}" | sort | uniq -c
