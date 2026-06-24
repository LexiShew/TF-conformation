#!/bin/bash
# Pilot config: engrailed homeodomain / 3HDD  (new-pilot acceptance test)
#
# End-to-end test for the B1–B3 wiring: a homeodomain added only via this config,
# run Stage 1 -> Stage 2 -> fnat with no manual path edits.
#
# Chain layout from 3hdd.cif (0-based mdtraj order; verify with
#   python stage2_redock/stage2_redock.py --ref structures/source_chains/3hdd_chains/3hdd.cif --inspect-only):
#   [0] chain C : DNA (21)     [1] chain D : DNA (21)
#   [2] chain A : protein (55) [3] chain B : protein (56)
#
# NOTE: 3HDD has TWO engrailed copies (chains A, B) over the single C/D duplex.
# The monomer guard (B3) refuses if BOTH contact the selected DNA. If it fires,
# either confirm only one copy binds (keep as-is) or pass --allow-multimer
# (STAGE2_ALIGN_MODE aside) — engrailed binds DNA as a monomer biologically; the
# second copy is crystallographic. 1HDD (chains C/D protein over A/B duplex) is
# an equivalent alternative.

export TF_NAME="engrailed"
export PDB_ID="3hdd"

# PWM (from config/deeppbs_pdb_pwms.csv). 3hdd is not in the id.txt benchmark, so
# this is a wiring placeholder — PWM_LABEL is unused until Stage 4, and there are
# no held-out engrailed test entries for Stage 6/7 transfer eval.
export PWM_LABEL="MA0220.1.jaspar"
export TEST_PWM_FILTER="MA0220"
export TEST_FILTER_NAME="engrailed"

# B2 — Stage 1 ensemble selector: source-chain filename letter
# (3hdd_chainA_protein.pdb -> "A") -> selects 3hdd_chainA_conformations/.
export BINDING_CHAIN="A"

# Stage 2 chain layout (0-based cif chainids). PROTEIN_CHAIN=2 is chain A.
export PROTEIN_CHAIN=2
export DNA_CHAINS="0,1"

# Stage 3 minimization parameters (match the other pilots).
export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

export N_FRAMES=100
export FOLD=0

# Optional fnat gate floor (B7); defaults to 0.5 if unset.
# export FNAT_FLOOR=0.5
