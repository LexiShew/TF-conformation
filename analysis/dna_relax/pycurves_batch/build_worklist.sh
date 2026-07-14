#!/bin/bash
# Build a worklist (one line per structure) for the pyCurves batch array.
# Columns (tab-sep): TF  COND  STATE  INPUT_PDB  OUT_PREFIX
# COND in {crystal,frozen,relaxed}. Crystal is one row per TF (the source DNA).
set -eo pipefail
BASE="/project2/rohs_102/shewchuk/TF-conformation"
OUTROOT="${BASE}/analysis/dna_relax/pycurves"
WL="${BASE}/analysis/dna_relax/pycurves_worklist.tsv"
declare -A PID=( [tbp]=1tgh [dux4]=5z6z [egr1]=1aay [engrailed]=3hdd [ets1]=1k79 [foxa]=1vtn [lef1]=2lef )
: > "${WL}"
for tf in "${!PID[@]}"; do
    pid=${PID[$tf]}
    # crystal (source DNA)
    cx="${BASE}/structures/source_chains/${pid}_chains/${pid}_dna.pdb"
    [ -f "$cx" ] && printf "%s\t%s\t%s\t%s\t%s\n" "$tf" crystal crystal "$cx" "${OUTROOT}/${tf}/crystal" >> "${WL}"
    # frozen + relaxed ensembles
    for cond in frozen relaxed; do
        if [ "$cond" = frozen ]; then dir="${BASE}/output/stage3_min/${tf}"
        else dir="${BASE}/output/stage3_min_dnarelax/${tf}"; fi
        for f in "${dir}"/${pid}_state_*.pdb; do
            [ -f "$f" ] || continue
            st=$(basename "$f" | sed -E 's/.*_state_([0-9]+)\.pdb/\1/')
            printf "%s\t%s\t%s\t%s\t%s\n" "$tf" "$cond" "$st" "$f" "${OUTROOT}/${tf}/${cond}_state_${st}" >> "${WL}"
        done
    done
done
echo "worklist: ${WL} ($(wc -l < "${WL}") rows)"
cut -f1,2 "${WL}" | sort | uniq -c
