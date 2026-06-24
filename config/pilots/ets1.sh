#!/bin/bash
# Pilot config: Ets-1 ETS domain / 1K79
#
# Chain layout from 1k79.cif (0-based mdtraj order; verify with --inspect-only):
#   [0] B : DNA (15)  [1] C : DNA (15)  [2] E : DNA (15)  [3] F : DNA (15)
#   [4] A : protein (104)                [5] D : protein (104)
#
# 1K79 has TWO Ets-1 copies, each on its OWN duplex (B/C and E/F): two separate
# monomer complexes in the asymmetric unit. We take chain A and its cognate
# duplex only, so the monomer guard (B3) sees a single protein on the selected
# DNA. VERIFY which duplex chain A binds (0,1 = B/C assumed); if A actually binds
# E/F, set DNA_CHAINS="2,3".

export TF_NAME="ets1"
export PDB_ID="1k79"

# PWM (config/deeppbs_pdb_pwms.csv): MA0098.2.jaspar / ETS1_MOUSE.H11MO.0.A.
# Not in the id.txt benchmark — wiring placeholder (unused until Stage 4).
export PWM_LABEL="MA0098.2.jaspar"
export TEST_PWM_FILTER="MA0098|ETS1"
export TEST_FILTER_NAME="ETS1"

# B2 — Stage 1 ensemble selector: 1k79_chainA_protein.pdb -> "A".
export BINDING_CHAIN="A"

# Stage 2 chain layout (0-based cif chainids). PROTEIN_CHAIN=4 is chain A;
# DNA_CHAINS is the single duplex chain A binds (NOT all four DNA chains, or the
# monomer guard would see chain D contacting the other duplex).
export PROTEIN_CHAIN=4
export DNA_CHAINS="0,1"

# Stage 3 minimization parameters.
export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

export N_FRAMES=100
export FOLD=0

# export FNAT_FLOOR=0.5
