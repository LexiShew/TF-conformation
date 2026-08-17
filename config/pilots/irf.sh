#!/bin/bash
# Pilot config: Interferon regulatory factor 1 / IRF1 winged-helix / 1IF1
#
# NOTE - MULTI-CHAIN RISK (verified via stage2_redock.py --inspect-only):
#   chainid 2 : protein (105 res)   chainid 3 : protein (104 res)
#   chainid 0 : DNA (26)            chainid 1 : DNA (26)   [ONE duplex, 2 strands]
# 1IF1 has TWO IRF1 copies binding TANDEM on a SINGLE duplex (the ISRE).
# Unlike ets1/1k79 (two proteins on two separate duplexes), here both proteins
# contact the one selected duplex, so the monomer guard (B3) will see the second
# IRF (chain 3) contacting DNA while we dock a single-protein BioEmu frame.
# This is the dux4-like regime. We select ONE protein (chain A / chainid 2) and
# the single duplex (both strands 0,1) and let the fnat gate decide survivorship;
# if 0 frames survive (as with dux4), this pilot is not viable in the monomeric
# pipeline. Diagnostic checkpoint: inspect stage3 fnat pass-rate before training.

export TF_NAME="irf"
export PDB_ID="1if1"

# PWM (config/deeppbs_pdb_pwms.csv): IRF1_MOUSE.H11MO.0.A.
# Not in the id.txt benchmark - wiring placeholder (unused until Stage 4).
# Same-family benchmark entries (n=4): 7jm4_{A,B,G,H} MA1419.1 (IRF4).
export PWM_LABEL="IRF1_MOUSE.H11MO.0.A"
export TEST_PWM_FILTER="MA1419|IRF"
export TEST_FILTER_NAME="IRF"

# B2 - Stage 1 ensemble selector: 1if1_chainA_protein.pdb -> "A".
export BINDING_CHAIN="A"

# Stage 2 chain layout (0-based cif chainids). PROTEIN_CHAIN=2 is chain A (longest);
# DNA_CHAINS="0,1" is the single tandem-bound duplex (both strands).
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
