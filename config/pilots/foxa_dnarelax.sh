#!/bin/bash
# Pilot config: FOXA — DNA-relaxation variant (Tier 1).
# forkhead (rigid control)
#
# Identical to foxa.sh but opts in to the DNA-relaxation pipeline. TF_NAME is
# inherited from foxa.sh so Stage 1 (BioEmu library) and Stage 2 (docked frames)
# are REUSED from the frozen-DNA run; common.sh appends a "_dnarelax" suffix so
# Stage 3+ outputs land in a parallel tree and never overwrite the frozen-DNA
# baseline.
#
# Headline mode (agreed 2026-07-10): soft tether + late release.
#   - STAGE3_DNA_RESTRAINT_K=1.5  : DNA on a soft tether (stiffness floor vs
#     GBSA B-DNA melting; the on-ramp to the Tier-2 cgDNA+ stiffness prior).
#   - STAGE3_DNA_RELEASE_STAGE=5  : DNA pinned at protein k through the 4
#     clash-resolving vdW-ramp stages, released to k=1.5 for the final stage +
#     full-vdW final minimization.

# Inherit the full foxa config (TF_NAME, PDB_ID, chains, ramp, N_FRAMES, fold...).
_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${_here}/foxa.sh"

# --- DNA-relaxation opt-in (this is the whole delta vs foxa.sh) ---
export STAGE3_DNA_RESTRAINT_K=1.5
export STAGE3_DNA_RELEASE_STAGE=5

# fnat gate: re-score relaxed structures through the EXISTING floor as a
# diagnostic first; do NOT re-threshold yet. Inherit FNAT_FLOOR default (0.5).
