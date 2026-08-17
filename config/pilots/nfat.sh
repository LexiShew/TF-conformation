#!/bin/bash
# Pilot config: NFATc1 Rel-homology domain / 1A66
#
# Rel-homology DNA-binding domain (Pfam PF00554, P53-like clan CL0073 / E-set
# CL0159). TRUE MONOMER: the PDB entry 1A66 is titled "SOLUTION NMR STRUCTURE
# OF THE CORE NFATC1/DNA COMPLEX" -- a binary complex of the human NFATc1
# DNA-binding domain with the ARRE2 site (single protein chain + one duplex).
# Verified via web search of the PDB entry (not PDB chain count). The related
# crystal 1OWR (NFATc2) is titled "...BOUND MONOMERICALLY TO DNA" and
# corroborates monomeric DNA binding for the family; the DBD is notably flexible
# -- a good conformational-augmentation target.
#
# Family transfer bonus: 7 Rel-homology entries in the id.txt benchmark
# (NF-kB p50/p65, NFAT) share the fold.
#
# Chain layout from 1a66.cif (0-based mdtraj order; verify with --inspect-only):
#   [0] DNA (12)  [1] DNA (12)  [2] protein (178)
# Clean single protein chain -> auto-selected PROTEIN_CHAIN=2, DNA_CHAINS=0,1.
# No structural metals.

export TF_NAME="nfat"
export PDB_ID="1a66"

# PWM (config/deeppbs_pdb_pwms.csv): NFAC1_HUMAN.H11MO.0.B.
# 1a66 is a TRAINING entry (not in id.txt) -> PWM_LABEL wiring placeholder.
export PWM_LABEL="NFAC1_HUMAN.H11MO.0.B"
# Stage 6/7 transfer eval: held-out Rel-homology subset of id.txt (NF-kB + NFAT).
export TEST_PWM_FILTER="NFAC|NFKB|REL_|TF65|MA0105|MA0107|MA0152|MA0778|MA0518"
export TEST_FILTER_NAME="Rel"

# B2 - Stage 1 ensemble selector: 1a66_chainA_protein.pdb -> "A".
export BINDING_CHAIN="A"

export PROTEIN_CHAIN=2
export DNA_CHAINS="0,1"

export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

export N_FRAMES=100
export FOLD=0

# export FNAT_FLOOR=0.5
