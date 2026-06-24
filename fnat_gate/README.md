# fnat_gate

Data-integrity gate between **Stage 2 (redock)** and **Stage 3 (minimize)** (B7).

BioEmu samples backbones broadly; many docked conformations are in the wrong
register relative to the DNA. Without a filter, those wrong-register states
become training pairs in Stage 5 and teach DeepPBS the wrong structure→PWM map.
This gate scores every docked state and drops the bad ones before they propagate.

## What it does

`fnat_gate.sh` runs `interface_rmsd.py --states-dir` over `${STAGE2_DIR}`:

- Scores each `${PDB_ID}_state_*.pdb` against the reference `.cif`, measuring
  **fnat** with `--use_model_dna` (model protein vs the DNA carried in the same
  docked file — the canonical frame check).
- Moves any state with `fnat < FNAT_FLOOR` (or that fails to score) into
  `${STAGE2_DIR}/rejected_fnat/`, so Stage 3's array simply skips them.
- Writes a per-state metrics CSV (`interface_metrics.csv`) and a drop log
  (`fnat_drops.log`) in `${STAGE2_DIR}`.

## Knob

| Var | Default | Meaning |
|---|---|---|
| `FNAT_FLOOR` | `0.5` | Minimum native-contact fraction to keep a state. Set in a pilot config or the environment to tune per-TF. |

## Where it runs

`run_pilot.sh` submits it as step **2g**, after Stage 2 and before Stage 3
(only when the run continues past Stage 3). The scorer (`interface_rmsd.py`)
is vendored here so the pipeline doesn't depend on the separate
`TF_conf_init_outputs/` analysis tree; that tree keeps its own copy for
post-hoc analysis. Requires `biopython>=1.79` + `numpy` (the `deeppbs` env).
