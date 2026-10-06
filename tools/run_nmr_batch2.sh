#!/bin/sh
# Batch 2 — BioEmu ensembles for the screened apo/holo NMR pairs (20 of 22 staged).
# Excluded by construct screen: c02 (autoinhibited Ets-1 dN301, 140v96 res, seq_id 0.81),
#                               c12 (73v112 res, seq_id 0.78).
# One sbatch job per PDB; each job generates conformations for every protein chain it finds.

sh tools/bioemu_from_pdb.sh -F 1a66 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 1ahd 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 1f4s 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 1ig6 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 1nfa 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 1p6r 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 1vf9 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 1vfc 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2ahq 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2hoa 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2jx1 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2jyd 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2k1n 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2k9n 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2kdz 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2kei 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2kej 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2kek 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2khl 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2ko0 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2l7f 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2l7z 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2ld5 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2lkx 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2ltd 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2ltt 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2mxe 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2mxf 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2nbj 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2o8k 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2o9l 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2oeh 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2p7c 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 2ro4 100 structures/nmr/
sh tools/bioemu_from_pdb.sh -F 3alc 100 structures/nmr/
