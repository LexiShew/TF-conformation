# DNA-relaxation structure analysis

Structure analysis of the Stage-3 DNA-relaxation feature (`--dna-restraint-k` +
`--dna-release-stage`, soft tether k=1.5 released at ramp stage 5), comparing
the frozen-DNA baseline (`output/stage3_min/<tf>`) against the relaxed pipeline
(`output/stage3_min_dnarelax/<tf>`), both referenced to the docked Stage-2 input
(`output/stage2_docked/<tf>`).

## Contents

```
scripts/
  batch_dna_shape.py        numpy geometry over a full state ensemble (RMSD, bend,
                            P-P fraying, per-residue disp, per-bp C1'-C1' width)
  run_pycurves.sh           run pyCurves on one PDB (MUST run on a compute node)
  parse_pycurves_grooves.py parse groove widths + overall bend from pyCurves .txt
figures/
  tbp_dna_shape.png         4-panel TBP ensemble figure (frozen vs relaxed)
data/
  tbp_dna_perstate.csv      per-state RMSD/bend/P-P gap
  tbp_dna_perres.csv        per-residue mean P displacement
  tbp_dna_perbp.csv         per-base-pair C1'-C1' width
  pycurves_tbp_state002.json crystal/frozen/relaxed bend + minor-groove (state 2)
```

## How to run

### Ensemble geometry (deeppbs env, CPU, numpy only)
```bash
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh; conda activate deeppbs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python scripts/batch_dna_shape.py \
  ../../output/stage2_docked/tbp \
  ../../output/stage3_min/tbp \
  ../../output/stage3_min_dnarelax/tbp \
  data/tbp_dna
```
Note: the fixed base-pair register in `batch_dna_shape.py` (REG) is the 1tgh TBP
TATA duplex (chain B 101-112 <-> chain C 124-113). Edit REG for other constructs.

### pyCurves helical parameters (pycurves env, CPU JAX, MUST be a compute node)
The login node's per-user thread limit crashes XLA's thread pool; run via srun/sbatch:
```bash
srun -p rohs -A rohs_102 -c 4 -t 00:05:00 \
  bash scripts/run_pycurves.sh ../../output/stage3_min_dnarelax/tbp/1tgh_state_002.pdb pycurves/tbp/tbp_state_002_relax
python scripts/parse_pycurves_grooves.py pycurves/tbp
```

## Key findings (TBP, 1tgh, state ensemble)

- **DNA moves ~1.7x more under relaxation**: median DNA-backbone RMSD vs docked
  0.72 A (relaxed, n=96) vs 0.41 A (frozen, n=94).
- **No excess fraying vs the frozen baseline**: max adjacent P-P gap tail is
  matched (1 state >8 A in each condition); no state reaches the >=9 A
  helix-unwind regime. The k=1.5 stiffness floor holds the duplex batch-wide.
- **Displacement localizes to the central TATA bases** (the TBP kink region),
  not the ends.
- **Relaxation bends toward the bound conformation** (pyCurves, state 2):
  overall bend UU crystal 79.7 / frozen 81.6 / relaxed 86.7; PP crystal 72.2 /
  frozen 68.0 / relaxed 73.5 -- relaxed PP-bend lands closest to crystal while
  frozen underbends. All three keep TBP's widened minor groove (~12 A vs 5.7 A
  canonical B-DNA). Consistent with induced-fit DNA deformation.

Caveat: `batch_dna_shape.py` metrics are numpy geometric proxies; pyCurves
provides the reference-grade helical parameters. The ensemble pyCurves run
(bend/groove distributions over all states, not just state 2) is the natural
next step.
