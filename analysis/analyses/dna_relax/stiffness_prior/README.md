# Sequence-dependent stiffness prior for DNA relaxation (hexABC)

Second-generation DNA-relaxation restraint: replace the uniform Stage-3
`k_dna = 1.5` with a **per-base-pair-step** stiffness derived from the hexABC
MD ensemble, so the restraint encodes known sequence-dependent deformability.

## Motivation
The current relaxation tethers every DNA backbone atom to one soft harmonic
restraint (`k_dna = 1.5`), identical everywhere regardless of sequence. But DNA
stiffness is strongly sequence-dependent — TA steps are hinges, AT/GC steps are
rigid. A uniform restraint therefore over-constrains floppy steps and
under-constrains stiff ones.

## Data source
`/project2/rohs_102/share/HexABC_data/seq*/analyses/average/*_stiffness.json`
— per-step elastic constants (twist/roll/tilt/shift/slide/rise + `sum`/`product`)
from 380 hexABC MD sequences (20-mer duplexes).

## Table (`stiffness_table.json`)
Built by `build_stiffness_table.py`: for each dinucleotide step label, mean±std
of `sum_stiffness` (summed diagonal elastic constant) across all 380 sequences.

Stiffness ordering (summed, flexible → stiff):
`TA/TA (10.0) < CA/TG ~ TG/CA (11.0) < CG/CG (12.0) < ... < GC/GC (16.8) < AT/AT (18.2)`
— a ~1.8× range; TA steps floppiest (matches TATA-box hinge biology), AT/AT stiffest.

## Mapping (`seq_to_kdna.py`)
`seq_to_kdna(seq, table, k_base=1.5)` → per-step k array (len = len(seq)−1),
mean-normalized so the average step gets `k_base`. Floppy steps drop below
`k_base` (more freedom), stiff steps rise above it.

Example (tbp TATA site `CGTATATATACG`): TA steps → k≈1.09, AT/GC flanks → k≈1.98.

## Integration into Stage-3 (planned)
1. In `stage3_minimize.py`, read the DNA sequence from the input PDB.
2. `k_step = seq_to_kdna(seq, table, k_base=DNA_RESTRAINT_K)`.
3. Broadcast each step's k onto the P,C1' atoms of its two flanking nucleotides
   (a per-atom k array on the existing DNA `CustomExternalForce`, replacing the
   scalar global parameter `k_dna`).
4. New config flag `STAGE3_DNA_STIFFNESS_PRIOR=hexABC` (default off → uniform k,
   byte-identical to current behavior).

## Status
Prototype: table + mapping validated on TATA vs GC test sequences. NOT yet wired
into Stage-3 — the in-flight benchmarks test the uniform k=1.5 restraint; the
prior-based restraint is a second-generation variant to be benchmarked against
that baseline once current runs land.
