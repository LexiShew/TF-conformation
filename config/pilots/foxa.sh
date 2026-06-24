#!/bin/bash
# Pilot config: HNF-3γ / FOXA forkhead / 1VTN
#
# Clean monomer: one protein chain over one DNA duplex.
# Chain layout from 1vtn.cif (0-based mdtraj order; verify with --inspect-only):
#   [0] chain A : DNA (13)   [1] chain B : DNA (13)   [2] chain C : protein (102)

export TF_NAME="foxa"
export PDB_ID="1vtn"

# PWM (config/deeppbs_pdb_pwms.csv): MA1683.1.jaspar / FOXA3_HUMAN.H11MO.0.B.
# Not in the id.txt benchmark — wiring placeholder (unused until Stage 4).
export PWM_LABEL="MA1683.1.jaspar"
export TEST_PWM_FILTER="MA1683|FOXA"
export TEST_FILTER_NAME="FOXA"

# B2 — Stage 1 ensemble selector: 1vtn_chainC_protein.pdb -> "C".
export BINDING_CHAIN="C"

# Stage 2 chain layout (0-based cif chainids). PROTEIN_CHAIN=2 is chain C.
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
