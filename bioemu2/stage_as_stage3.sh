#!/bin/bash
# stage_as_stage3.sh - present BioEmu2 frames to the pipeline as the Stage-3
# conformation set, so the fnat gate and Stages 4-7 run COMPLETELY UNMODIFIED.
#
# Why this shape. BioEmu2 samples the protein-DNA complex directly, replacing
# Stages 1-2 (and, given the raw-frame gate result, Stage 3 as well). The pipeline
# expects the conformations entering the gate to sit in STAGE3_DIR as
# <PDB_ID>_state_NNN.pdb; BioEmu2's states_1copy/ already uses exactly that
# naming. Symlinking them into STAGE3_DIR therefore needs no code change anywhere
# downstream: the gate scores them and builds STAGE3_PASS_DIR itself, exactly as
# it does for minimized BioEmu1 frames, so the two sources are gated identically.
#
# Symlinks, not copies: the frames are ~2.6 GB in total and output/ is gitignored.
#
# Requires CONF_SRC=bioemu2 so STAGE3_DIR resolves to output/stage3_min_be2/<tf>
# and nothing touches the BioEmu1 tree. Refuses to run otherwise.
#
# Usage:
#   CONF_SRC=bioemu2 bash bioemu2/stage_as_stage3.sh            # all eligible pilots
#   CONF_SRC=bioemu2 TF_LIST="ets1" bash bioemu2/stage_as_stage3.sh
set -eo pipefail

TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
cd "${TFCONF_DIR}"
source "${TFCONF_DIR}/lib/common.sh"

if [ "${CONF_SRC:-bioemu1}" != "bioemu2" ]; then
    echo "ERROR: run with CONF_SRC=bioemu2 - refusing to stage BioEmu2 frames into the BioEmu1 tree." >&2
    exit 1
fi

if [ -z "${TF_LIST:-}" ]; then
    TF_LIST=$(for c in config/pilots/*.sh; do basename "$c" .sh; done | tr '\n' ' ')
fi

for tf in ${TF_LIST}; do
    ( load_pilot_config "${tf}" >/dev/null 2>&1
      src="${TFCONF_DIR}/output/bioemu2/${PDB_ID}_complex/states_1copy"
      n_src=$( { ls "${src}"/${PDB_ID}_state_*.pdb 2>/dev/null || true; } | wc -l )
      if [ "${n_src}" -eq 0 ]; then
          echo "[stage_as_stage3/${tf}] SKIP - no states_1copy frames"
          exit 0
      fi
      mkdir -p "${STAGE3_DIR}"
      ln -sf "${src}"/${PDB_ID}_state_*.pdb "${STAGE3_DIR}/"
      n_dst=$( { ls "${STAGE3_DIR}"/${PDB_ID}_state_*.pdb 2>/dev/null || true; } | wc -l )
      echo "[stage_as_stage3/${tf}] ${PDB_ID}: ${n_dst} frames -> ${STAGE3_DIR#${TFCONF_DIR}/}"
    ) || echo "[stage_as_stage3/${tf}] FAILED (continuing)"
done
