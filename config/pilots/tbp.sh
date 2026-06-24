#!/bin/bash
# Pilot config: TBP / 1tgh
 
export TF_NAME="tbp"
export PDB_ID="1tgh"
 
# All 5 TBP test entries in id.txt use this PWM
export PWM_LABEL="MA0343.1.jaspar"
 
# Regex matching test entries in id.txt
export TEST_PWM_FILTER="MA0343"
export TEST_FILTER_NAME="TBP"
 
# Stage 1 ensemble selector (B2): source-chain filename letter
# (1tgh_chainA_protein.pdb -> "A") -> selects 1tgh_chainA_conformations/.
# Distinct from PROTEIN_CHAIN (the 0-based cif chainid).
export BINDING_CHAIN="A"

# Stage 2 chain layout for 1tgh.cif (verified via stage2_redock.py --inspect-only)
# 1tgh: protein on chain 2 (180 residues), DNA on chains 0,1 (12bp each strand)
export PROTEIN_CHAIN=2
export DNA_CHAINS="0,1"
 
# Stage 3 minimization parameters (defaults match EGR1 pilot)
export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000
 
# Number of frames in BioEmu xtc (drives Stage 3 array size)
export N_FRAMES=99
 
# Training fold to augment
export FOLD=0