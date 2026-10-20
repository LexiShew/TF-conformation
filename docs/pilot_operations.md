# Pilot operations — registry, generation, and keeping the cluster busy

How pilots are defined and run as of 2026-10-15. This replaces hand-written
`config/pilots/<tf>.sh` files and manual `run_pilot.sh` babysitting.

---

## The short version

```bash
# 1. add one row to the registry
$EDITOR config/pilots.tsv

# 2. generate its config (derives PWM_LABEL; refuses to guess)
python scripts/pipeline/gen_pilot_configs.py

# 3. check it
python scripts/pipeline/validate_pilot.py <tf> 

# 4. run it — or let the queue pick it up
./scripts/pipeline/run_pilot.sh <tf>
ARM=dnarelax ./scripts/pipeline/run_pilot.sh <tf>     # the relaxed arm

# keep N chains in flight, forever (cron)
*/30 * * * * cd /path/to/TF-conformation && \
  MAX_INFLIGHT=6 python scripts/pipeline/pilot_queue.py --fill >> slurm_output/queue.log 2>&1
```

---

## What changed and why

### Before

A pilot was a hand-written ~50-line shell file, plus a second near-identical file for the
relaxed arm. 23 files for 12 pilots. Measured across the 12 base configs:

| parameter | distinct values across 12 pilots |
|---|---|
| `FOLD`, `RAMP_STAGES`, `STEPS_PER_STAGE`, `RECOVERY_RAMP_STAGES`, `RECOVERY_STEPS_PER_STAGE` | **1 each — byte-identical everywhere** |
| `BINDING_CHAIN` / `PROTEIN_CHAIN` / `DNA_CHAINS` / `N_FRAMES` | 2–3 each |
| `TF_NAME` / `PDB_ID` / `PWM_LABEL` / `TEST_PWM_FILTER` / `TEST_FILTER_NAME` | 12 each |

And all 11 `*_dnarelax.sh` files carried the *same* two-line delta
(`STAGE3_DNA_RESTRAINT_K=1.5`, `STAGE3_DNA_RELEASE_STAGE=5`).

So a third of every config was boilerplate, and half the files existed to set two variables.
Worse: `PWM_LABEL` was typed by hand, which is how the TBP pilot ran for months against
NDT80's motif.

### After

| artifact | role |
|---|---|
| `config/pilots.tsv` | **the single source of truth.** One row per pilot. |
| `scripts/pipeline/gen_pilot_configs.py` | generates `config/pilots/<tf>.sh`; `--check` for CI |
| `scripts/pipeline/validate_pilot.py` | 10 correctness checks, run as a `run_pilot.sh` pre-flight gate |
| `scripts/pipeline/pilot_queue.py` | infers state from disk, keeps N chains in flight |
| `lib/common.sh` | holds the 5 constant defaults and the `ARM` switch |

**23 config files → 12 generated files + 1 registry.**

---

## The registry

Columns are documented in the header of `config/pilots.tsv`. The important design decision is
what is *not* a column:

- **`PWM_LABEL` is derived**, looked up in `config/deeppbs_pdb_pwms.csv` keyed on `pdb_id`.
  A motif belonging to a different protein cannot be written down.
- **The 5 constant stage-3 parameters are not columns** — they are `lib/common.sh` defaults
  (`${VAR:-default}`, so a registry value still wins if a pilot needs to differ).
- **The relaxed arm is not a second row.** It is `ARM=dnarelax`.

### The one escape hatch, and why it is still safe

Four structures (`1lo1`, `1k79`, `1vtn`, `5d5u`) carry **both** a HOCOMOCO and a JASPAR label.
The generator refuses to pick between them and fails. Those four pilots pin their historical
choice in an optional `pwm_label` column — and **a pin is still verified against the master
table**, so pinning chooses among *that structure's own* motifs and cannot smuggle in another
protein's. The pins exist because the original runs used JASPAR, and silently switching to
HOCOMOCO would change `Y_pwm` and invalidate those pilots' results.

This came up during migration: the generator reported 4 mismatches rather than quietly
overwriting. That is the behaviour working.

---

## Validation gates

`validate_pilot.py` runs 10 checks. Current state: **12 pilots, 0 errors, 7 warnings.**

| check | what it catches |
|---|---|
| **C1** | `PWM_LABEL` not a motif of `PDB_ID` — **the TBP error**; would have caught it on day one |
| **C2** | filter that is an invalid regex, matches 0 entries, or is pseudoreplicated |
| **C3** | on-target subset whose family disagrees with the pilot's own |
| **C4** | `ON_TARGET_N_EXPECTED` drift; malformed `ON_TARGET_USABLE` |
| **C5** | source PDB present in `id.txt` (test-set leakage) |
| **C6** | `BINDING_CHAIN` naming a file that does not exist |
| **C7** | `PROTEIN_CHAIN` overlapping `DNA_CHAINS`; inline comments on the value |
| **C8** | missing required vars; implausible `N_FRAMES` |
| **C9** | annotation invariant — motif+Pfam agree on a family but `family` disagrees |
| **C10** | **seeding (F6)** — `driver.py` must read `random_seed`; every training config must carry one; `baseline_X_sN` and `augmented_X_sN` must share it. Negative-tested: removing the read, or desyncing a pair, fails the gate. |

The 7 remaining warnings are real scientific caveats, not bugs, and are deliberately not
suppressed:

```
ets1   5 entries / 4 PDBs      foxa  5 entries / 2 PDBs
hsf    3 entries / 2 PDBs      nfat 15 entries / 7 PDBs
tbp    n=1 — no paired statistic possible
```

`ON_TARGET_USABLE=no` is a hard declaration that a pilot has no usable on-target test.
Three pilots carry it: **engrailed** (own motif absent from the benchmark), **lef1** (only
HMG-box in the benchmark, no family fallback), **tbp** (n=1, and that entry is a histone chain).

---

## The continuous queue

`pilot_queue.py` turns 24 (pilot × arm) chains into one cron-able command.

```
$ python scripts/pipeline/pilot_queue.py --status
PILOT      ARM       PROGRESS              NEXT STAGE      STATE
csl        frozen    #####...  5/8         3g fnat gate    idle
ndt80      frozen    ##......  2/8         3 minimize      idle
...
24 unit(s): 0 complete, 0 in flight, 24 idle   (max-inflight=4)
```

**State comes from the filesystem, not a database.** A stage is done if its output exists and
is non-empty, so the queue is restartable and crash-proof: delete an output dir and it gets
redone. There is nothing to corrupt and no state to reconcile.

Scheduling order:

1. `ARM=frozen` before `ARM=dnarelax` — frozen is the baseline everything pairs against.
2. Pilots closest to completion first — finish work in progress before starting new.
3. Registry order, as a stable tiebreak.

Safety properties:

- Never submits for a (pilot, arm) that already has jobs in `squeue`.
- Runs `validate_pilot.py` per candidate and **skips failures loudly** rather than submitting.
- Submits at most one chain per unit per invocation, so a cron loop ramps smoothly instead of
  dumping 24 chains at once.
- Refuses to submit at all if `squeue` is unavailable (it will not submit blind).
- `--dry-run` prints the plan.

Tune concurrency with `MAX_INFLIGHT` or `--max-inflight`. With the compute you have, 6–8 is
reasonable; each chain is 9 jobs including an `N_FRAMES`-wide Stage-3 array throttled at `%8`.

### It also surfaces stage-4 attrition

The queue prints a warning when npz count < gated-input count:

```
Stage-4 attrition detected:
  irf        frozen    5 npz  !! 73 lost (94%) in stage 4
```

This is the failure that made the IRF pilot unusable and went unnoticed for months, because
Stage 4 compared contact-count to npz-count but never inputs to outputs. Now a routine status
check shows it. (IRF is retired — see `retired_pilots.md` — but the detector stays.)

---

## Adding a new pilot

1. **Get the structure in.** `structures/source_chains/<PDB>_chains/` needs `<PDB>.cif`,
   `<PDB>_chain<X>_protein.pdb`, `<PDB>_dna.pdb`.

2. **Inspect the chain layout** — do not guess:
   ```bash
   python stage2_redock/stage2_redock.py \
     --ref structures/source_chains/<PDB>_chains/<PDB>.cif --inspect-only
   ```

3. **Check the structure is eligible.** The generator enforces the first of these; the rest
   are judgement:
   - not in `data/folds/id.txt` (leakage — generation *fails* if it is)
   - a true monomer on DNA (dux4 was retired for being a homodimer)
   - a real two-strand duplex that 3DNA resolves as **one** helix — irf was retired because
     1if1 resolves as two, and Stage 4 silently dropped 94% of its states
   - ideally in neither `train*.txt` nor `valid*.txt`, so augmenting cannot double-count an
     entry the baseline already trains on (`ndt80`/2euv was chosen on exactly this basis)

4. **Add the row**, run the generator, run the validator, and either launch it or let the
   queue pick it up.

If `on_target_usable=yes` but the filter matches nothing, **generation fails** rather than
producing a pilot whose subset Stage 7 would silently skip. That is how engrailed, foxa and
lef1 went unnoticed.

---

## Still outstanding

These are review findings the operational work above does not fix:

- **F6 — seeding is FIXED in code, but existing results are not.** `driver.py` now seeds
  python/numpy/torch/cuda from `C["random_seed"]`, and check **C10** fails the pipeline if
  that read is ever removed or if a baseline/augmented pair disagrees on its seed. What this
  does *not* do is retroactively pair the 264 checkpoints already on disk — those were run
  unseeded. **Any claim that rests on a paired baseline/augmented comparison needs a
  retrain.** See "Retraining after the seeding fix" below.
- **F9 — Stage 7 p-values use a normal approximation**, understating p by 2.6–15× at n=5.
  `analysis/common/stats.py` now has the correct implementation; `evaluate_id_benchmark.py`
  still needs to be switched over to it.
- **F7 — Stage 4 still does not fail on attrition.** The queue *reports* it; the stage itself
  should refuse to proceed past a threshold.


---

## Retraining after the seeding fix (F6)

The seeding bug is fixed in code, but **every existing checkpoint was trained unseeded**, so
no existing `baseline_X_sN` / `augmented_X_sN` pair shares an initialisation. Fixing the code
does not retroactively pair them.

What the fix now guarantees, verified empirically:

| property | status |
|---|---|
| same seed → identical parameter init | **yes** (hash-identical across runs) |
| different seed → different init | **yes** |
| same seed → identical shuffle order, *same* dataset size | **yes** |
| same seed → identical shuffle order, *different* dataset size | **no** — and this is the baseline-vs-augmented case, since the augmented set is larger. The arms share init, not shuffle order. The pre-fix comment claimed both. |
| unseeded runs differ from each other | **yes** — what all 264 historical runs did |

### What needs rerunning, in priority order

1. **Any pilot whose conclusion rests on a per-seed paired Δ.** That is most of the
   `output/stage7_eval/id_benchmark_*.json` results and the mechanism analysis's own-family
   table. The per-*entry* pairing inside a run was always valid; it is the seed-level pairing
   that was fictional.
2. **tbp** — needs a rerun regardless, because its npz carry NDT80's `Y_pwm`
   (`docs/tbp_label_error_downstream.md`). Do the seeding fix and the label fix in one pass
   rather than retraining twice.
3. **ndt80** — the replacement pilot; has never been trained.

Stage 6 onward is the expensive part, and Stages 1–5 are unaffected by seeding, so:

```bash
# per pilot, both arms, seeds 1-5 — stages 6 and 7 only
for s in 1 2 3 4 5; do
  ./scripts/pipeline/run_pilot.sh <tf> 6 7           # frozen
  ARM=dnarelax ./scripts/pipeline/run_pilot.sh <tf> 6 7
done
# or let the queue work through everything:
MAX_INFLIGHT=8 python scripts/pipeline/pilot_queue.py --fill
```

Before trusting the rerun, confirm the seeds actually took:

```bash
python scripts/pipeline/validate_pilot.py --all        # C10 must pass
grep -h "Seeded all RNGs" slurm_output/train_*.out | sort -u
```

A seeded run prints `Seeded all RNGs (python/numpy/torch/cuda) with random_seed=N`. If a log
instead prints the `WARNING: ... UNSEEDED` line, the config lacked `random_seed` and that run
must be discarded.

### Optional: bit-identical reruns

Parameter-init pairing does not require cuDNN determinism, so it is off by default (it costs
speed). For bit-identical reruns set `"deterministic_cudnn": true` in the training config, or
pass `--deterministic_cudnn`.
