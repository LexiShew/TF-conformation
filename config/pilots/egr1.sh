#!/bin/bash
# Pilot config: EGR1 / 1aay (the original pilot, preserved as a reference).
# Set ALL variables that the rest of the pipeline reads.

export TF_NAME="egr1"
export PDB_ID="1aay"

# PWM label used for augmentation npzs.  This MUST match what the held-out
# test entries use, otherwise transfer won't happen (cf. EGR1 pilot's
# MA0162.1 vs MA0162.4 finding).
export PWM_LABEL="EGR1_MOUSE.H11MO.0.A"

# Regex matching test entries we focus on in eval reports.
# This is for the "EGR1 entries (n=5)" report subset.
export TEST_PWM_FILTER="EGR1|MA0162"
export TEST_FILTER_NAME="EGR1"

# Stage 2 chain layout — for reference structure (cif).
# 1aay: chain A=DNA, B=DNA, C=protein, D-F=ZN, G-I=waters
# 0-based chainids; we want the protein chain and the two DNA strands.
export PROTEIN_CHAIN=2
export DNA_CHAINS="0,1"

# Stage 3 minimization parameters (defaults match the EGR1 pilot)
export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

# Number of frames in BioEmu xtc (drives Stage 3 array size)
export N_FRAMES=98

# Training fold to augment (0 = use train0/valid0)
export FOLD=0
