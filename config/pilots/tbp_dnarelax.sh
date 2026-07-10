#!/bin/bash
# Pilot config: TBP / 1tgh — DNA-relaxation variant (Tier 1).
#
# Identical to tbp.sh but opts in to the DNA-relaxation pipeline. TF_NAME stays
# "tbp" so Stage 1 (BioEmu library) and Stage 2 (docked frames) are REUSED from
# the frozen-DNA run; common.sh appends a "_dnarelax" suffix so Stage 3+ outputs
# land in a parallel tree (output/stage3_min_dnarelax/tbp, etc.) and never
# overwrite the frozen-DNA baseline.
#
# Headline mode (agreed 2026-07-10): soft tether + late release.
#   - STAGE3_DNA_RESTRAINT_K=1.5  : DNA on a soft tether (not fully free) — a
#     stiffness floor against GBSA-driven B-DNA melting, and the natural on-ramp
#     to the Tier-2 cgDNA+ stiffness prior (a non-zero restraint by construction).
#   - STAGE3_DNA_RELEASE_STAGE=5  : keep DNA pinned at the protein k through the
#     4 clash-resolving vdW-ramp stages (largest forces), release to k=1.5 only
#     for the final ramp stage + the full-vdW final minimization.

# Inherit the full TBP config (TF_NAME, PDB_ID, chains, ramp, N_FRAMES, fold...).
_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${_here}/tbp.sh"

# --- DNA-relaxation opt-in (this is the whole delta vs tbp.sh) ---
export STAGE3_DNA_RESTRAINT_K=1.5
export STAGE3_DNA_RELEASE_STAGE=5

# fnat gate: re-score the relaxed structures through the EXISTING floor first as
# a diagnostic (see the 2026-07-10 log). Do NOT re-threshold yet — inspect the
# fnat distribution before deciding. Inherit FNAT_FLOOR default (0.5) for now.
