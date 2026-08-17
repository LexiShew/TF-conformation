#!/bin/bash
# Pilot config: RBPJ / CSL (Suppressor of Hairless) / 3BRG
#
# CSL DNA-binding module (Pfam PF09271 LAG1 + PF09270 beta-trefoil; clans
# CL0073 / CL0066). TRUE MONOMER: RBPJ/CSL binds DNA through a single chain; the
# Notch intracellular domain and MAML co-activators assemble OFF the DNA-reading
# module. 3BRG is a mouse Su(H)-DNA complex (single protein chain + duplex;
# verified via web search of the primary literature, not PDB chain count).
#
# SOURCE CHOICE: use 3BRG (a TRAINING entry), NOT 3iag -- 3iag is itself in the
# id.txt benchmark test set, so using it as an augmentation source would leak
# held-out data into training.
#
# Family transfer bonus: 4 CSL entries in the id.txt benchmark (3iag, 3v79,
# 6dks, 6wqu) share the SuH/RBPJ fold and motif.
#
# Chain layout from 3brg.cif (0-based mdtraj order; verify with --inspect-only):
#   [0] DNA (15)  [1] DNA (15)  [2] protein (414)  [3..] EDO/HOH heterogens
# Clean single protein chain -> auto-selected PROTEIN_CHAIN=2, DNA_CHAINS=0,1.
# No structural metals.

export TF_NAME="csl"
export PDB_ID="3brg"

# PWM (config/deeppbs_pdb_pwms.csv): SUH_MOUSE.H11MO.0.A.
# 3brg is a TRAINING entry (not in id.txt) -> PWM_LABEL wiring placeholder.
export PWM_LABEL="SUH_MOUSE.H11MO.0.A"
# Stage 6/7 transfer eval: held-out CSL/SuH subset of id.txt.
export TEST_PWM_FILTER="SUH_|MA1116"
export TEST_FILTER_NAME="CSL"

# B2 - Stage 1 ensemble selector: 3brg_chainC_protein.pdb -> "C".
export BINDING_CHAIN="C"

export PROTEIN_CHAIN=2
export DNA_CHAINS="0,1"

export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

export N_FRAMES=100
export FOLD=0

# export FNAT_FLOOR=0.5
