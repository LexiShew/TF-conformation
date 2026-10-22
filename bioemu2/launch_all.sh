#!/bin/bash
# launch_all.sh - build a BioEmu2 worklist and submit the sampling array.
#
# Basic use
#   bash bioemu2/launch_all.sh                      # every ok structure, qcbgpu L40S
#   LIMIT=20 bash bioemu2/launch_all.sh             # smoke test: first 20 only
#   DRY_RUN=1 bash bioemu2/launch_all.sh            # print the sbatch, submit nothing
#
# Size-routed split across two partitions (what is actually in flight, 2026-10-22).
# The worklist is sorted small-first, so a positional split would send the LARGEST
# systems to whichever partition got the back half. Route by size instead: small
# systems fit a 16 GB rtx5000, the large tail needs the L40S's 46 GB.
#
#   # small systems -> rohs (lab's own partition, rtx5000 16 GB)
#   WORKLIST=bioemu2/worklist/source_chains_small.tsv TOKEN_MAX=4000 \
#     PARTITION=rohs GRES=gpu:rtx5000:1 MAX_CONCURRENT=8 \
#     bash bioemu2/launch_all.sh
#
#   # large tail -> qcbgpu (L40S 46 GB), with a longer clock
#   WORKLIST=bioemu2/worklist/source_chains_large.tsv TOKEN_MIN=4001 \
#     PARTITION=qcbgpu GRES=gpu:l40s:1 MAX_CONCURRENT=6 TIME_LIMIT=12:00:00 \
#     bash bioemu2/launch_all.sh
#
# Knobs
#   TOKEN_MIN / TOKEN_MAX  inclusive bounds on the survey's n_tokens_proxy (col 8)
#   PARTITION / GRES       sbatch overrides; the .sbatch header defaults to qcbgpu/l40s
#   TIME_LIMIT             sbatch --time override (header default 04:00:00)
#   MAX_CONCURRENT         the array's %N throttle (default 8)
#   SKIP_DONE              1 (default) drops structures whose manifest already
#                          reports NUM_SAMPLES frames, so the array is tight rather
#                          than spending a GPU slot per no-op skip
#   NUM_SAMPLES / SEED     sampling parameters (default 100 / 1)
#   WORKLIST               output path; MUST differ per concurrent array, because
#                          each task re-reads this file at runtime and a rewrite
#                          would repoint an already-queued array
set -eo pipefail

TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
cd "${TFCONF_DIR}"

SURVEY="${SURVEY:-bioemu2/worklist/source_chains_survey.tsv}"
WORKLIST="${WORKLIST:-bioemu2/worklist/source_chains_ok.tsv}"
OUT_ROOT="${OUT_ROOT:-${TFCONF_DIR}/output/bioemu2}"
MAX_CONCURRENT="${MAX_CONCURRENT:-8}"
NUM_SAMPLES="${NUM_SAMPLES:-100}"
SEED="${SEED:-1}"
TOKEN_MIN="${TOKEN_MIN:-0}"
TOKEN_MAX="${TOKEN_MAX:-100000000}"
SKIP_DONE="${SKIP_DONE:-1}"

[ -f "${SURVEY}" ] || {
    echo "ERROR: no survey at ${SURVEY}"
    echo "  run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \\"
    echo "       /project2/rohs_102/shewchuk/conda/envs/pycurves/bin/python \\"
    echo "       bioemu2/survey_source_chains.py --out ${SURVEY}"
    exit 1
}

mkdir -p "$(dirname "${WORKLIST}")" logs/bioemu2

# status==ok (both protein and DNA) within the token window, small-first order kept.
head -1 "${SURVEY}" > "${WORKLIST}"
awk -F'\t' -v lo="${TOKEN_MIN}" -v hi="${TOKEN_MAX}" \
    'NR>1 && $9=="ok" && $8+0>=lo+0 && $8+0<=hi+0' "${SURVEY}" >> "${WORKLIST}"

# Drop anything already sampled at this frame count. Parsed with awk rather than a
# python call per structure: 698 interpreter startups is slower than the sampling.
if [ "${SKIP_DONE}" = "1" ]; then
    awk -F'\t' -v root="${OUT_ROOT}" -v ns="${NUM_SAMPLES}" '
      NR==1 { print; next }
      {
        man = root "/" $1 "_complex/sample_manifest.json"
        have = -1
        while ((getline line < man) > 0)
            if (line ~ /"n_frames_written"/) { gsub(/[^0-9]/, "", line); have = line + 0 }
        close(man)
        if (have != ns + 0) print
      }' "${WORKLIST}" > "${WORKLIST}.tmp" && mv "${WORKLIST}.tmp" "${WORKLIST}"
fi

# Drop explicitly excluded pdb_ids (e.g. structures a still-running array task
# already owns, so two arrays never write the same output dir).
if [ -n "${EXCLUDE_FILE:-}" ] && [ -f "${EXCLUDE_FILE}" ]; then
    awk -F'\t' -v ex="${EXCLUDE_FILE}" '
      BEGIN { while ((getline l < ex) > 0) if (l != "") skip[l] = 1 }
      NR==1 { print; next }
      !($1 in skip)' "${WORKLIST}" > "${WORKLIST}.tmp" && mv "${WORKLIST}.tmp" "${WORKLIST}"
fi

# OFFSET drops the first N rows AFTER filtering, so one partition can take the
# head of the (small-first) list and another the tail without overlap.
if [ -n "${OFFSET:-}" ] && [ "${OFFSET}" -gt 0 ]; then
    { head -1 "${WORKLIST}"; tail -n +"$((OFFSET + 2))" "${WORKLIST}"; } > "${WORKLIST}.tmp" \
        && mv "${WORKLIST}.tmp" "${WORKLIST}"
fi

[ -n "${LIMIT:-}" ] && { head -n "$((LIMIT + 1))" "${WORKLIST}" > "${WORKLIST}.tmp" && mv "${WORKLIST}.tmp" "${WORKLIST}"; }

N=$(($(wc -l < "${WORKLIST}") - 1))
[ "${N}" -gt 0 ] || { echo "nothing to submit: worklist empty after filtering"; exit 0; }

echo "worklist : ${WORKLIST}  (${N} structures)"
echo "tokens   : ${TOKEN_MIN} .. ${TOKEN_MAX}  (survey n_tokens_proxy)"
echo "samples  : ${NUM_SAMPLES} per structure, seed ${SEED}"
echo "array    : 1-${N}%${MAX_CONCURRENT}"
echo "target   : partition=${PARTITION:-<header default>} gres=${GRES:-<header default>} time=${TIME_LIMIT:-<header default>}"

CMD=(sbatch --array="1-${N}%${MAX_CONCURRENT}")
[ -n "${PARTITION:-}"  ] && CMD+=(--partition="${PARTITION}")
[ -n "${GRES:-}"       ] && CMD+=(--gres="${GRES}")
[ -n "${TIME_LIMIT:-}" ] && CMD+=(--time="${TIME_LIMIT}")
CMD+=(--export="ALL,TFCONF_DIR=${TFCONF_DIR},WORKLIST=${TFCONF_DIR}/${WORKLIST},NUM_SAMPLES=${NUM_SAMPLES},SEED=${SEED}"
      bioemu2/sample_all.sbatch)

if [ "${DRY_RUN:-0}" = "1" ]; then printf '%q ' "${CMD[@]}"; echo; exit 0; fi
"${CMD[@]}"
