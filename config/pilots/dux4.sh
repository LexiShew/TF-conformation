#!/bin/bash
# Pilot config: DUX4 / 5z6z

export TF_NAME="dux4"
export PDB_ID="5z6z"

# 4 of 5 DUX4 test entries in id.txt use this PWM
export PWM_LABEL="DUX4_HUMAN.H11MO.0.A"

# Regex includes both DUX4_HUMAN PWM and the MA0468.1 PWM that 6dfy_C uses
# (gives us bonus matched-vs-unmatched within-pilot replication)
export TEST_PWM_FILTER="DUX4|MA0468"
export TEST_FILTER_NAME="DUX4"

# Stage 2 chain layout for 5z6z reference cif.
# IMPORTANT: this is a guess — verify with auto-detection (see below) before
# launching, or pass --auto-detect-chains to stage2_redock.py.
# For typical homeodomain-on-DNA: protein and two DNA strands.
export PROTEIN_CHAIN=0     # ⚠ verify with: python scripts/inspect_chains.py --ref <ref.cif>
export DNA_CHAINS="1,2"

# Stage 2 mismatch handling:
# 5z6z has a 9-residue disordered loop (resSeq 82-92) that BioEmu modeled
# across. The new sequence-identity check in stage2_redock.py will detect
# this and refuse to proceed unless told how to handle it. 'trim' uses the
# longest common-sequence run for Kabsch alignment, dropping the gap region.
# (Note: doesn't help with DUX4's underlying problem of BioEmu sampling far
# from the bound state — but at least doesn't fail Stage 2 outright.)
export MISMATCH_ACTION="trim"
export MAX_MISMATCHES=0

export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

export N_FRAMES=98
export FOLD=0