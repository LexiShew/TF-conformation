#!/bin/bash
# build_pilot_worklist.sh - worklist of just the pilots, for the BioEmu2 batch.
#
# Why pilots are the default scope: extract_one_copy.py (--tf-name, required) and
# compare_to_pipeline.py both read chain selection from config/pilots/<tf>.sh, so
# ONLY a pilot can be reduced to one protein + its cognate duplex and fnat-scored.
# Without that reduction interface_rmsd.score() silently falls back to the
# REFERENCE DNA (finding F10, docs/codebase_review_2026-10-15.md) and reports a
# number that is arithmetic rather than a measurement. A non-pilot ensemble is
# therefore unscorable with what is in the repo today.
#
#   bash bioemu2/build_pilot_worklist.sh                 # -> worklist/pilots.tsv
set -eo pipefail
TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
cd "${TFCONF_DIR}"
SURVEY="${SURVEY:-bioemu2/worklist/source_chains_survey.tsv}"
OUT="${OUT:-bioemu2/worklist/pilots.tsv}"

head -1 "${SURVEY}" > "${OUT}"
for cfg in config/pilots/*.sh; do
    pid=$(bash -c "set -a; source '${cfg}' >/dev/null 2>&1; echo \${PDB_ID}" | tr 'A-Z' 'a-z')
    [ -n "${pid}" ] || continue
    awk -F'\t' -v q="${pid}" '$1==q && $9=="ok"' "${SURVEY}" >> "${OUT}"
done
# Keep the survey's small-first order rather than config glob order.
{ head -1 "${OUT}"; tail -n +2 "${OUT}" | sort -t"$(printf '\t')" -k8,8n; } > "${OUT}.tmp"
mv "${OUT}.tmp" "${OUT}"
echo "pilots worklist -> ${OUT} ($(($(wc -l < "${OUT}") - 1)) structures)"
cut -f1,8 "${OUT}" | tail -n +2 | tr '\n' ' '; echo
