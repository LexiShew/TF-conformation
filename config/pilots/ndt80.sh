#!/bin/bash
# Pilot config: NDT80 / 2euv
#
# Rationale (2026-10-01). This pilot exists because the TBP pilot was found to
# carry NDT80's PWM: config/pilots/tbp.sh set PWM_LABEL=MA0343.1.jaspar, but
# MA0343.1 is NDT80 (S. cerevisiae, UniProt P38830), not TBP. The five id.txt
# entries matching MA0343 (2etw/2euz/2evf/2evi/2evj) are all NDT80 structures
# (Pfam PF05224, NDT80/PhoG-like). See docs/tbp_label_audit.md.
#
# Here the label is correct by construction: 2EUV *is* NDT80, so structure,
# Y_pwm target, and the five on-target test anchors all finally agree.
#
# Selection criteria (all satisfied):
#   - monomeric: 1 protein chain (290 res) + real 2-strand duplex (13+12 nt),
#     verified via stage2_redock.py --inspect-only. No IRF-style missing strand.
#   - not in test set: 2euv absent from data/folds/id.txt.
#   - new Pfam family: PF05224 (NDT80/PhoG-like), clan CL0073 P53-like.
#     Family is new across all 13 existing pilots.
#   - within-family test anchors: 5 (2etw_A, 2euz_A, 2evf_A, 2evi_A, 2evj_A).
#   - bonus: 2euv is in NEITHER train*.txt NOR valid*.txt NOR id.txt, so
#     augmenting cannot double-count an entry the baseline already trains on.
#     (Sibling NDT80 structures 2euw/2eux/2evg/2evh/1mnn ARE in the folds —
#     that is why 2euv specifically was chosen over them.)

export TF_NAME="ndt80"
export PDB_ID="2euv"

# MA0343.1 = NDT80. Correct for this structure.
export PWM_LABEL="MA0343.1.jaspar"

# Regex matching the NDT80 test entries in id.txt (n=5)
export TEST_PWM_FILTER="MA0343"
export TEST_FILTER_NAME="NDT80"

# Stage 1 ensemble selector: source-chain filename letter
# (2euv_chainA_protein.pdb -> "A")
export BINDING_CHAIN="A"

# Stage 2 chain layout for 2euv.cif (verified via stage2_redock.py --inspect-only):
#   chainid 2: protein (290 res) | chainid 0,1: DNA (13,12 nt)
export PROTEIN_CHAIN=2
export DNA_CHAINS="0,1"

# Stage 3 minimization parameters (defaults, match TBP/EGR1 pilots)
export RAMP_STAGES="0.1,0.3,0.5,0.7,1.0"
export STEPS_PER_STAGE=500
export RECOVERY_RAMP_STAGES="0.05,0.1,0.2,0.4,0.7,1.0"
export RECOVERY_STEPS_PER_STAGE=1000

# Number of frames in BioEmu xtc (drives Stage 3 array size)
export N_FRAMES=99

# Training fold to augment
export FOLD=0
