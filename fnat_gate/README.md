# fnat_gate

The pipeline's **single** structural-quality gate, between **Stage 3 (minimize)**
and **Stage 4 (preprocess)** (B7).

BioEmu samples backbones broadly; many docked conformations are in the wrong
register relative to the DNA. Without a filter, those wrong-register states
become training pairs in Stage 4/5 and teach DeepPBS the wrong structure→PWM map.
This gate scores every **post-minimization** state and keeps only the good ones.

## Why post-minimization (Stage 3), not at dock time (Stage 2)

Minimization moves per-state fnat in **both** directions — a docked pose just
under the floor can be lifted over it, and a good pose can be loosened. Scoring
at Stage 2 would (a) measure the wrong (pre-min) pose and (b) silently shrink the
denominator by discarding near-misses before minimization could rescue them.

So the rule is: **Stage 2 carries every docked state forward** (it may log an
fnat readout for diagnostics, but never drops/moves/quarantines), and the **only**
fnat filter runs here, on the minimized pose.

## What it does

`fnat_gate.sh` runs `score_stage3.py` over `${STAGE3_DIR}`:

- Scores each `${PDB_ID}_state_*.pdb` against the reference `.cif`, measuring
  **fnat** vs the model's own DNA (`use_model_dna=True` — the correct frame for a
  minimized complex).
- Writes a per-state metrics CSV (`${PDB_ID}_fnat.csv`) and a pass-list
  (`${PDB_ID}_pass.txt`, one `<pdb>_state_NNN` per surviving state).
- Rebuilds `${STAGE3_PASS_DIR}` (= `${STAGE3_DIR}_pass`) from scratch as a mirror
  of **symlinks** to the passing states (identical filenames). Stage 4 reads
  ONLY this dir.
- **Fails loud:** if no state clears the floor it exits non-zero, so the DAG's
  `afterok` edge halts Stage 4 for that pilot.

The rebuild is idempotent: re-running the gate `rm -rf`s the pass dir and
recreates it. The guard only ever removes a path ending in `_pass`.

## Knob

| Var | Default | Meaning |
|---|---|---|
| `FNAT_FLOOR` | `0.5` | Minimum native-contact fraction to keep a state. |

`FNAT_FLOOR` is overridable globally (export it before launching) **or per-pilot**
by setting it in `config/pilots/<tf>.sh`, e.g.:

```bash
export FNAT_FLOOR=0.45   # in config/pilots/foxa.sh
```

The per-pilot value wins because `common.sh` sets `FNAT_FLOOR="${FNAT_FLOOR:-0.5}"`
*after* sourcing the pilot config, and the gate honors the inherited value.
Lowering it enlarges the pass-list; raising it shrinks it.

## Where it runs

`scripts/pipeline/run_pilot.sh` submits it as the **Stage 3 → Stage 4** node:
it depends `afterok` on Stage 3 (and the Stage 3 recover step, if present), and
Stage 4 depends `afterok` on the gate. The scorer (`interface_rmsd.py`,
`score_stage3.py`) is vendored here. Requires `biopython>=1.79` + `numpy`
(the `deeppbs` env).
