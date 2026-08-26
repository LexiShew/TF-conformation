#!/usr/bin/env python3
"""
interface_mgwfl.py — Interface-restricted MGW fluctuation (MGW-FL) analysis.

Reproduces the core-motif MGW-FL metric of Jiang et al. (Biophys J 2026,
"Readout of intrinsic and induced DNA shape by homeodomain TF complexes")
on the TF-conformation pilots, and correlates it with DeepPBS augmentation
accuracy (ΔPearson / ΔMAE, augmented − baseline).

Pipeline:
  1. CONTACTING DNA:   uses fnat_gate/interface_rmsd.ref_side() on each pilot's
     crystal .cif with the pilot config's PROTEIN_CHAIN / DNA_CHAINS, contact
     cutoff 4.5 Å heavy-atom min-dist (same definition as the fnat gate).
     -> the set of DNA residues (chain,resnum) that contact the DBD.
  2. LEVEL MAPPING:    aligns each pyCurves groove profile's base-pair sequence
     (strand 1) as a substring of the crystal strand-1 sequence to map
     pyCurves level index -> crystal residue number. A level is "interface" if
     its strand-1 residue is in the contacting set. (Exact substring match for
     all 12 pilots; termini pyCurves can't define are dropped.)
  3. MGW-FL:           per-position MGW = mean(minor_width) over sub-levels;
     MGW-FL = population SD across the ensemble (frozen / relaxed / AF3 samples).
     Interface MGW-FL = mean MGW-FL over interface levels only.

Inputs (on endeavour, under /project2/rohs_102/shewchuk/TF-conformation):
  - structures/source_chains/<pid>_chains/<pid>.cif        crystal reference
  - config/pilots/<tf>.sh                                  PROTEIN_CHAIN/DNA_CHAINS
  - analysis/dna_relax/pycurves/<tf>/{crystal,frozen_state_*,relaxed_state_*}_legacy.json
  - af3/af3_dna/<tf>_<pid>/*_legacy.json                   AF3-sample groove
  - output/stage7_eval/id_benchmark_<tf>.json              accuracy metrics

Outputs (analysis/dna_relax/):
  - data/iface_mgwfl_vs_accuracy.csv    per-pilot interface MGW-FL + ΔP/ΔMAE
  - data/mgwfl_vs_accuracy.csv          whole-molecule version
  - figures/iface_mgwfl_vs_accuracy.png ΔPearson & ΔMAE vs interface MGW-FL

Run under the `deeppbs` conda env (needs biopython); pyCurves JSONs must
already exist (built by analysis/dna_relax/pycurves_batch/). See the Notion
research log 2026-07-23 for the interpretation.
"""
# Implementation lives in the session handoff; this header documents the method
# so the analysis is reproducible. Key parameters: contact_cut=4.5, iface_cut=5.0.
