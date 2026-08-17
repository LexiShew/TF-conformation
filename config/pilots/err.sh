#!/bin/bash
# Pilot config: Estrogen-related receptor beta (ERRbeta) / 1LO1
#
# Nuclear receptor DNA-binding domain (Pfam PF00105, clan CL0839 Tbcl_zf).
# TRUE MONOMER: verified against the primary literature via web search (not PDB
# chain count). The source paper reports "Monomeric complex of human orphan
# estrogen-related receptor-2 with DNA: a pseudo-dimer interface mediates
# extended half-site" recognition (PubMed 12654265). New fold vs the 7 pilots.
#
# Family transfer bonus: 7 nuclear-receptor entries in the DeepPBS id.txt
# benchmark (RXR/GR half-site readers) share the C4-Zn fold.
#
# Chain layout from 1lo1.cif (0-based mdtraj order; verify with
#   python stage2_redock/stage2_redock.py \
#     --ref structures/source_chains/1lo1_chains/1lo1.cif --inspect-only):
#   [0] DNA (13)  [1] DNA (13)  [2] protein (90)  [3] ZN  [4] ZN
# Clean single protein chain -> auto-selected PROTEIN_CHAIN=2, DNA_CHAINS=0,1.
# 2 structural Zn ions: Stage 3 auto-detects and cages them (default, no flag).

export TF_NAME="err"
export PDB_ID="1lo1"

# PWM (config/deeppbs_pdb_pwms.csv): ERR2_HUMAN.H11MO.0.A / MA0141.3.jaspar.
# 1lo1 is a TRAINING entry (not in id.txt) -> PWM_LABEL is a wiring placeholder,
# unused until Stage 4.
export PWM_LABEL="MA0141.3.jaspar"
# Stage 6/7 transfer eval: score the held-out nuclear-receptor subset of id.txt
# (RXR / GR / NR5A half-site readers). Motif-based so only NR-fold entries match.
export TEST_PWM_FILTER="RXRA|GCR_|NR5A|ESRR|MA0065|MA0113|MA0115|MA0114|MA0512|MA0007"
export TEST_FILTER_NAME="NR"

# B2 - Stage 1 ensemble selector: 1lo1_chainA_protein.pdb -> "A".
export BINDING_CHAIN="A"

# Stage 2 chain layout (0-based cif chainids).
export PROTEIN_CHAIN=2
export DNA_CHAINS="0,1"

# Stage 3 minimization parameters (match the other pilots).
export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

export N_FRAMES=100
export FOLD=0

# export FNAT_FLOOR=0.5
