#!/bin/bash
# Pilot config: Interferon regulatory factor 1 / IRF1 winged-helix / 1IF1
#
# 1IF1 biological assembly (from CIF _pdbx_struct_assembly): ONE assembly,
# "tetrameric" (count 4) = 2 IRF1 + 1 palindromic 26-bp duplex (2 strands).
# CIF chain letters (what stage2 reads via --ref): A,B = DNA strands (one duplex,
# self-complementary ISRE), C,D = protein. split_chains uses the OPPOSITE letters
# (protein A/B, dna C/D) and emits ONE 1if1_dna.pdb holding both strands, because
# the two proteins share the single duplex — each over one half-site.
#
# Monomer strategy (per user, 2026): keep protein CIF-chain C (= split-file
# protein A) as the sole biological binder over the FULL duplex, and declare the
# second IRF1 (CIF chain D) a non-biological partner via MONOMER_ASSEMBLY=1.
# NOTE: chain D is a genuine cooperative ISRE partner, not crystal packing; this
# asserts IRF1 as an isolated monomer on a duplex whose 2nd half-site is empty.
# Chain layout verified via stage2_redock.py --inspect-only:
#   chainid 2 (C) : protein (105)   chainid 3 (D): protein (104)
#   chainid 0 (A) : DNA (26)        chainid 1 (B): DNA (26)  [one duplex]

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

# Stage 2 chain layout (0-based cif chainids). PROTEIN_CHAIN=2 is CIF chain C
# (= split-file protein A); DNA_CHAINS="0,1" is the full duplex (both strands).
export PROTEIN_CHAIN=2
export DNA_CHAINS="0,1"

# Declare configured protein chain the sole biological binder; excludes the
# second IRF1 (chainid 3) from the monomer guard. See header caveat.
export MONOMER_ASSEMBLY=1

# Stage 3 minimization parameters.
export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

export N_FRAMES=100
export FOLD=0

# export FNAT_FLOOR=0.5
