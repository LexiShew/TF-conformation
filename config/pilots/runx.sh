#!/bin/bash
# Pilot config: RUNX1 / AML1 Runt domain / 1HJC
#
# Runt domain (Pfam PF00853, P53-like clan CL0073). TRUE MONOMER on DNA: the
# Runt domain contacts DNA on its own. In the RUNX-CBFbeta-DNA regulatory
# complex, CBFbeta makes NO DNA contact -- it allosterically enhances Runt DNA
# affinity (established in the primary literature; verified via web search).
# BioEmu samples the isolated Runt domain, which is exactly the DNA-reading unit.
# For THIS pilot the point is moot: 1hjc contains no CBFbeta at all (see below).
#
# SOURCE CHOICE: 1HJC is a Runt-domain/DNA complex with NO CBFbeta present --
# the asymmetric unit holds TWO identical Runt copies, each on its OWN duplex
# (verified: both protein chains have identical sequence GELVRTDSPNFLCSVLP...).
# No CBFbeta-stripping needed. We take ONE Runt copy + its cognate duplex, as in
# the ets1/1k79 two-copies-on-separate-duplexes pilot.
# 1HJC is a TRAINING entry (not in id.txt) -> no test-set leakage.
#
# Family transfer bonus: 5 RUNX1-motif entries in the id.txt benchmark
# (1hjb, 3wts, 3wtu, 3wtw, 4l0z) -- the single largest uncovered test bucket.
# NOTE: the RUNX-ETS composite crystals (3wts/3wtu/6vg8/6vge) also carry
# ETS-motif benchmark entries; the motif-based filter below selects ONLY the
# RUNX1 half so the ETS entries are not miscounted as RUNX transfer targets.
#
# Chain layout from 1hjc.cif (0-based mdtraj order; verify with
#   python stage2_redock/stage2_redock.py \
#     --ref structures/source_chains/1hjc_chains/1hjc.cif --inspect-only):
#   [0] protein A (118)  [3] protein D (118)
#   [1] DNA (16) [2] DNA (16) [4] DNA (16) [5] DNA (16)
# Contact analysis: protein chain 0 (auth A) binds DNA chains 1,2;
#                   protein chain 3 (auth D) binds DNA chains 4,5.
# Take chain 0 + duplex (1,2). Setting DNA_CHAINS to only this duplex keeps the
# monomer guard (B3) from seeing chain 3 contacting the other duplex.

export TF_NAME="runx"
export PDB_ID="1hjc"

# PWM (config/deeppbs_pdb_pwms.csv): RUNX1_MOUSE.H11MO.0.A.
export PWM_LABEL="RUNX1_MOUSE.H11MO.0.A"
# Stage 6/7 transfer eval: held-out RUNX1-motif subset of id.txt (Runt half only).
export TEST_PWM_FILTER="RUNX"
export TEST_FILTER_NAME="RUNX"

# B2 - Stage 1 ensemble selector: 1hjc_chainA_protein.pdb -> "A".
export BINDING_CHAIN="A"

# Stage 2 chain layout (0-based cif chainids). PROTEIN_CHAIN=0 is chain A;
# DNA_CHAINS is the single duplex chain A binds (1,2), NOT all four DNA chains.
export PROTEIN_CHAIN=0
export DNA_CHAINS="1,2"

# 1hjc holds two identical Runt copies, each on its OWN duplex (chain 0 -> DNA
# 1,2 ; chain 3 -> DNA 4,5). Within 5 A of duplex (1,2) the monomer guard also
# sees a couple of grazing atoms from chain 3 (crystal packing; contact map:
# chain3-DNA1=1, chain3-DNA2=2, vs chain0=35,52). Declare chain 3 non-biological
# so the guard keeps chain 0 as the sole binder -- same mechanism as engrailed.
export MONOMER_ASSEMBLY=1

export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

export N_FRAMES=100
export FOLD=0

# export FNAT_FLOOR=0.5
