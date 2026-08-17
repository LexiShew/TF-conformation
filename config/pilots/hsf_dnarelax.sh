#!/bin/bash
# Pilot config: HSF1 - DNA-relaxation variant.
# HSF1 winged-helix (rigid control candidate)
#
# Identical to hsf.sh but opts in to the DNA-relaxation pipeline. TF_NAME is
# inherited so Stage 1 (BioEmu) and Stage 2 (docked frames) are REUSED from the
# frozen run; common.sh appends "_dnarelax" so Stage 3+ land in a parallel tree
# and never overwrite the frozen baseline.
#
# Headline mode: soft tether k=1.5 + late release (stage 5).
_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${_here}/hsf.sh"

# --- DNA-relaxation opt-in ---
export STAGE3_DNA_RESTRAINT_K=1.5
export STAGE3_DNA_RELEASE_STAGE=5
