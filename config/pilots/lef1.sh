#!/bin/bash
# Pilot config: LEF-1 HMG-box / 2LEF  (NMR)
#
# Clean monomer: one protein chain over one DNA duplex.
# Chain layout from 2lef.cif (0-based mdtraj order; verify with --inspect-only):
#   [0] chain B : DNA (15)   [1] chain C : DNA (15)   [2] chain A : protein (86)
#
# NMR structure (multiple models). Stage 1 (BioEmu) samples from sequence, so the
# NMR vs crystal distinction doesn't matter for sampling; Stage 2 docks onto the
# reference's first model (stage2_redock.py uses ref.xyz[0]).

export TF_NAME="lef1"
export PDB_ID="2lef"

# PWM (config/deeppbs_pdb_pwms.csv): LEF1_MOUSE.H11MO.0.B (HOCOMOCO only, no
# JASPAR). Not in the id.txt benchmark — wiring placeholder (unused until Stage 4).
export PWM_LABEL="LEF1_MOUSE.H11MO.0.B"
export TEST_PWM_FILTER="LEF1"
export TEST_FILTER_NAME="LEF1"

# B2 — Stage 1 ensemble selector: 2lef_chainA_protein.pdb -> "A".
export BINDING_CHAIN="A"

# Stage 2 chain layout (0-based cif chainids). PROTEIN_CHAIN=2 is chain A.
export PROTEIN_CHAIN=2
export DNA_CHAINS="0,1"

# Stage 3 minimization parameters.
export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

export N_FRAMES=100
export FOLD=0

# export FNAT_FLOOR=0.5
