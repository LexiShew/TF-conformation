# wrappers/

SLURM submission wrappers for the per-pilot pipeline. One per stage. Each holds
the stage's `#SBATCH` resource request, then:

1. resolves `TFCONF_DIR` (inherited from the launcher, else the cluster default),
2. `source`s `lib/common.sh` (sets paths + loads the pilot config from `TF_NAME`),
3. `source`s the stage's logic from its `stage*/` directory.

This replaces the old deployment under `DeepPBS/run/jobs/{lib,wrappers,scripts}`:
**TF-conformation is now the authoritative pipeline** — nothing is sourced from
DeepPBS anymore (only shared *data* trees, `DeepPBS_data/` and `DeepPBS_outputs/`,
are still read/written; see `lib/common.sh`).

| Wrapper | Sources | Notes |
|---|---|---|
| `stage1_bioemu.sh` | `stage1_bioemu/stage1_bioemu.sh` | BioEmu+HPacker on the pilot's chains. |
| `stage2_redock.sh` | `stage2_redock/stage2_redock.sh` | |
| `stage3_array.sh` | `stage3_minimize/stage3_array.sh` | submitter sets `--array=1-${N_FRAMES}%8` |
| `stage3_recover.sh` | `stage3_minimize/stage3_recover.sh` | |
| `stage4_preprocess.sh` | `stage4_preprocess/stage4_preprocess.sh` | |
| `stage5_build_aug.sh` | `stage5_build_aug/stage5_build_aug.sh` | |
| `train_compare.sh` | `stage6_train/train_compare.sh` | array `0-1` (baseline+augmented) |
| `train_legacy_aug.sh` | `stage6_train/train_legacy_aug.sh` | array `0-0` |
| `eval_benchmark.sh` | `stage7_eval/eval_benchmark.sh` | |
| `eval_legacy_ab.sh` | `stage7_eval/eval_legacy_ab.sh` | |

Driven by the launchers in `scripts/pipeline/` (`run_pilot.sh` for the full DAG,
`run_multiseed_pilot.sh`, and `run_legacy_ab.sh`). Do not submit these wrappers
by hand unless you set `--export=ALL,TF_NAME=<tf>` yourself.
