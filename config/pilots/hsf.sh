#!/bin/bash
# Pilot config: Heat shock factor 1 / HSF1 winged-helix / 5D5U
#
# Clean monomer: one protein chain over one DNA duplex (verified via
# stage2_redock.py --inspect-only on 5d5u.cif):
#   chainid 1 : protein (100 res)   chainid 0 : DNA (12)   chainid 2: HOH
# HSF1 DBD is a monomeric winged-helix; single protein on the duplex.

export TF_NAME="hsf"
export PDB_ID="5d5u"

# PWM (config/deeppbs_pdb_pwms.csv): MA0486.2.jaspar / HSF1_HUMAN.H11MO.0.A.
# Not in the id.txt benchmark - wiring placeholder (unused until Stage 4).
# Same-family benchmark entries (n=3): 5hdn/7dcj HSF1_HUMAN.H11MO.0.A, MA0486.2.
export PWM_LABEL="MA0486.2.jaspar"
export TEST_PWM_FILTER="MA0486|HSF"
export TEST_FILTER_NAME="HSF"

# B2 - Stage 1 ensemble selector: 5d5u_chainB_protein.pdb -> "B".
export BINDING_CHAIN="B"

# Stage 2 chain layout (0-based cif chainids). PROTEIN_CHAIN=1 is chain B;
# DNA_CHAINS=0 is the single duplex it binds.
export PROTEIN_CHAIN=1
export DNA_CHAINS="0"

# Stage 3 minimization parameters.
export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

export N_FRAMES=100
export FOLD=0

# export FNAT_FLOOR=0.5
