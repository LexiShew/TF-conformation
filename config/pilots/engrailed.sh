#!/bin/bash
# Pilot config: engrailed homeodomain / 1hdd  (NEW-PILOT ACCEPTANCE TEST)
#
# This is the end-to-end test for the B1–B3 wiring: a monomeric homeodomain
# added ONLY via this config, run Stage 1 -> Stage 2 -> fnat with no manual path
# edits. 1hdd has TWO protein copies in the asymmetric unit (chains C and D),
# so it also exercises the chain selector (B2) and the monomer guard (B3).
#
# ====================================================================
# BEFORE LAUNCHING, fill the TODOs below. Inspect the reference layout:
#   source lib/common.sh
#   python stage2_redock/stage2_redock.py \
#       --ref structures/source_chains/1hdd_chains/1hdd.cif --inspect-only
# That prints each chain's 0-based chainid + protein/DNA classification. Pick:
#   - the protein chainid that binds DNA  -> PROTEIN_CHAIN
#   - that protein's cognate DNA strands  -> DNA_CHAINS
# and set BINDING_CHAIN to the matching source-chain FILENAME letter
# (structures/source_chains/1hdd_chains/1hdd_chain<X>_protein.pdb -> "<X>").
# Stage 2's sequence-match guard will fail loudly if BINDING_CHAIN (the Stage 1
# ensemble) and PROTEIN_CHAIN (the .cif chain) disagree.
# ====================================================================

export TF_NAME="engrailed"
export PDB_ID="1hdd"

# TODO: set the DeepPBS PWM label this TF's test entries use (from id.txt).
export PWM_LABEL="TODO_ENGRAILED_PWM"          # e.g. an H11MO / MA id
export TEST_PWM_FILTER="TODO"                  # regex over id.txt entries
export TEST_FILTER_NAME="engrailed"

# B2 — Stage 1 ensemble selector. 1hdd source chains are chainC + chainD.
# Use the copy that binds DNA in the reference; "C" is the usual choice.
export BINDING_CHAIN="C"

# Stage 2 chain layout (0-based mdtraj chainids into 1hdd.cif).
# TODO: verify with --inspect-only above. DNA_CHAINS must be the SINGLE binding
# site's strands — if you include both copies' DNA, the monomer guard (B3) will
# correctly refuse, since both C and D would then contact "the DNA".
export PROTEIN_CHAIN=0      # TODO: chainid of the BINDING_CHAIN protein in the cif
export DNA_CHAINS="1,2"     # TODO: its two cognate DNA strands

# Stage 3 minimization parameters (defaults match the other pilots)
export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

# Number of conformations Stage 1 generates / Stage 3 array size.
export N_FRAMES=100

# Training fold to augment.
export FOLD=0

# Optional: fnat gate floor (B7). Defaults to 0.5 if unset.
# export FNAT_FLOOR=0.5
