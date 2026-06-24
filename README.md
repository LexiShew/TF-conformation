# TF-conformation

Augment [DeepPBS](https://github.com/timkartar/DeepPBS) (structure → binding-
specificity PWM) with **conformational ensembles** of monomeric transcription-
factor–DNA complexes. For each TF we sample a backbone ensemble with
[BioEmu](https://github.com/microsoft/bioemu), rebuild side chains with HPacker,
dock each conformation onto the crystal DNA, filter by interface fidelity, and
feed the survivors into DeepPBS training.

TF-conformation is the **authoritative pipeline repo**: it vendors the DeepPBS
package and the 3DNA toolchain under `lib/`, so nothing is sourced from a
separate deployment. Only the large training data/outputs still live in the
shared `DeepPBS_data/` and `DeepPBS_outputs/` trees on the cluster.

## Pipeline (per TF)

A pilot is one `config/pilots/<tf>.sh`. `scripts/pipeline/run_pilot.sh` submits
the whole SLURM DAG with dependencies:

| Stage | Dir | Does |
| --- | --- | --- |
| 1 | `stage1_bioemu/` | BioEmu sampling + HPacker side-chain reconstruction, for **every** protein chain of the structure → `structures/stage1_bioemu_output/<PDB>_chain<X>_conformations/` |
| 2 | `stage2_redock/` | Interface-aligned Kabsch dock of each frame onto the crystal DNA (carries DNA + structural metals across); monomer guard |
| 2g | `fnat_gate/` | Score each docked state's `fnat` and drop sub-floor states before they become training data (`FNAT_FLOOR`, default 0.5) |
| 3 | `stage3_minimize/` | OpenMM minimization with a metal-coordination cage |
| 4 | `stage4_preprocess/` | DeepPBS `process_co_crystal.py` → per-state `.npz` features (3DNA via vendored `lib/`) |
| 5 | `stage5_build_aug/` | Build the augmented training fold + combined assembly + train configs |
| 6 | `stage6_train/` | Train baseline vs augmented DeepPBS models (paired) |
| 7 | `stage7_eval/` | Benchmark eval; paired bootstrap/t-test stats |

Each stage dir is self-contained: its `.sh` runs the co-located `.py`. The
`wrappers/` dir holds the SBATCH wrapper per stage (resource header + `source
lib/common.sh` + `source` the stage logic). `lib/common.sh` resolves all paths
(self-locates `TFCONF_DIR`) and loads the pilot config.

## Repository layout

```
config/pilots/<tf>.sh          per-TF config (PDB_ID, BINDING_CHAIN,
                               PROTEIN_CHAIN, DNA_CHAINS, PWM_LABEL, FOLD, ...)
lib/                           common.sh + vendored deeppbs pkg + 3DNA toolchain
wrappers/                      one SBATCH wrapper per stage
stage1_bioemu/ ... stage7_eval/ per-stage logic + engine scripts
fnat_gate/                     interface-fidelity gate (vendored interface_rmsd.py)
scripts/pipeline/              run_pilot.sh, run_multiseed_pilot.sh, run_legacy_ab.sh
scripts/{analysis,pymol,classification,pdb_prep,maintenance,deprecated}/
                               standalone utilities (not part of the DAG)
structures/                    on-disk, gitignored: source_chains/ (per-chain
                               inputs) + stage1_bioemu_output/ (ensembles)
```

## Running

```bash
# Full pipeline for a pilot (stages 1–7):
./scripts/pipeline/run_pilot.sh <tf_name>

# A subrange (e.g. skip Stages 1–3, already done):
./scripts/pipeline/run_pilot.sh <tf_name> 4 7

# Multi-seed paired comparison (after stages 1–5):
./scripts/pipeline/run_multiseed_pilot.sh <tf_name> 5
```

### Adding a new TF

1. Make sure `structures/source_chains/<PDB>_chains/` exists (`<PDB>.cif`,
   `<PDB>_chain<X>_protein.pdb`, `<PDB>_dna.pdb`).
2. Inspect the reference chain layout:
   ```bash
   python stage2_redock/stage2_redock.py --ref structures/source_chains/<PDB>_chains/<PDB>.cif --inspect-only
   ```
3. Write `config/pilots/<tf>.sh` (copy an existing one). Set `BINDING_CHAIN`
   (the source-chain filename letter, picks the Stage 1 ensemble) and
   `PROTEIN_CHAIN`/`DNA_CHAINS` (0-based cif chainids). Stage 2's sequence-match
   guard fails loudly if `BINDING_CHAIN` and `PROTEIN_CHAIN` disagree.
4. `./scripts/pipeline/run_pilot.sh <tf_name>`.

See `config/pilots/engrailed.sh` for a worked new-pilot template (1hdd).

## Requirements

- Conda envs: `bioemu` (Stages 1–3), `deeppbs` (Stages 4–7 + the fnat gate,
  needs `biopython>=1.79`), `hpacker` (side-chain reconstruction).
- A SLURM cluster with a GPU partition. The wrappers request the `rohs`
  partition / `rohs_102` account — edit the `#SBATCH` headers in `wrappers/` and
  the cluster paths in `lib/common.sh` for other clusters.
- The `deeppbs` Python package (vendored at `lib/deeppbs/`, installed into the
  conda env) and the 3DNA toolchain (vendored at `lib/`).

See `PIPELINE_FIXES.md` for the current fix spec and design invariants.
