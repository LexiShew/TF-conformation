# bioemu2 — sampling protein–DNA complexes directly

An **exploratory parallel track**, not a pipeline stage. It tests whether
[BioEmu2](https://github.com/microsoft/bioemu2) can sample a protein–DNA
**complex** well enough to replace stages 1–3, and it is kept outside the
numbered DAG so nothing in `run_pilot.sh` depends on it.

Status as of 2026-10-21: **untested on this cluster.** See *Blockers* below.

---

## Why this is interesting, and not just a dependency bump

BioEmu2 is a different model from the BioEmu this repo uses, and two of its
changes bear directly on the project's weakest point.

**It models nucleic acids.** Stage 2 exists *because* BioEmu 1 cannot: the
pipeline samples a protein chain from its sequence, then Kabsch-docks every
frame onto the crystal DNA, then minimizes with the DNA backbone restrained.
That is exactly the limitation item #1 of `docs/concerns.txt` identifies:

> you're varying the protein conformation but pairing each frame with the same
> crystal-pose DNA. The protein–DNA interface — which is what shapes the PWM —
> is not actually being augmented … it's much closer to "protein ensemble around
> a fixed binding pose" than "ensemble of binding modes."

BioEmu2 can sample the complex, which attacks that criticism instead of working
around it. That is the one experiment here worth running.

**It is all-atom.** HPacker drops out — one fewer conda env, one fewer
reconstruction step, and no more "side-chain atoms in `samples.xtc` are
placeholders, not physically meaningful positions"
(`stage1_bioemu/generate_monomer_confs.py`).

**Ligands by CCD code.** Stage 3's metal-coordination cage is a workaround for
AMBER not knowing Zn²⁺. BioEmu2 takes `entity_type: ligand, ccd_code: ZN`
directly. `--exclude-ligands ZN` exists so you can A/B whether it matters —
relevant for `egr1` (C2H2 zinc finger) and `runx`.

---

## Blockers — check these before anything else

**1. Python.** `config/bioemu.yml` pins **3.11.14**; BioEmu2 needs **3.12+**.
A new env is required, not an upgrade. `bioemu2/environment.yml` has it.

**2. GPU memory.** The stage wrappers request `--gres=gpu:rtx5000:1` — a Quadro
RTX 5000, **16 GB**, compute capability 7.5. BioEmu2's model card states *"an
80 GB A100 being the smallest tested device."* That is 5× below their tested
floor. Your systems are small (180–414 residues plus a duplex, well under the
stated 10k-token ceiling) so it may fit, but you are off-road. Check for
something bigger:

```bash
sinfo -o '%P %G %m %N' | sort -u
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv
```

An OOM on a 16 GB card is the expected outcome, not a bug in these scripts.

---

## Install

```bash
conda env create -f bioemu2/environment.yml
conda activate bioemu2
git clone https://github.com/microsoft/bioemu2 ~/src/bioemu2
cd ~/src/bioemu2 && python -m pip install ".[cuda]"      # cc >= 7.5
```

Smoke-test the upstream install before touching this repo's data —
adenylate kinase is their canonical open/closed example:

```bash
bioemu2 sample-equilibrium --pdb-id 1AKE --out /tmp/adk
```

---

## Run

```bash
# interactive, one pilot
python bioemu2/sample_complex.py --tf-name ets1 --num-samples 100

# batch
sbatch --export=ALL,TF_NAME=ets1 bioemu2/sample_complex.sbatch

# does modelling the zinc change the ensemble?
sbatch --export=ALL,TF_NAME=egr1,EXCLUDE_LIGANDS=ZN bioemu2/sample_complex.sbatch

# compare against the existing pipeline, same metrics
python bioemu2/compare_to_pipeline.py --tf-name ets1 --arms frozen dnarelax
```

`ets1` (1k79) is the right first target: it is the strongest pilot and the one
whose own-family effect is closest to surviving correction.

### Outputs — `output/bioemu2/<pdb>_complex/`

| file | what |
|---|---|
| `topology.cif`, `samples.xtc` | BioEmu2 native, via its `save_frames` |
| `topology.pdb` | same topology as PDB, for the mdtraj/3DNA/PyMOL tools here |
| `states/<pdb>_state_NNN.pdb` | per frame, named so the **existing** fnat gate reads them |
| `physicality.json` | BioEmu2's own per-frame validity report |
| `sample_manifest.json` | checkpoint, seed, entity counts, timings |

This tree does **not** feed `stage2_redock`. That is the point: the track
replaces stages 1–2 rather than feeding them.

---

## What these three scripts deliberately do *not* do

Most of what a new sampling track appears to need already exists in this repo or
in BioEmu2. Reused, not rewritten:

| needed | reused from |
|---|---|
| system definition | BioEmu2 `System.from_cif()` on the crystal `.cif` this repo already keeps at `structures/source_chains/<PDB>_chains/<PDB>.cif` — protein chains, both DNA strands and structural metals come across in one call, so there is no YAML to build |
| cif → pdb | BioEmu2 returns an `mdtraj.Trajectory`, so `traj[0].save_pdb()` is a one-liner |
| structure QC | BioEmu2 `check_physicality()` |
| fnat / interface RMSD | `fnat_gate/interface_rmsd.py` and `score_stage3.py:score_dir` — called directly, so these numbers are produced the same way the pipeline's own gate produces them |
| effect size (CLES) | `analysis/analyses/nmr_rmsd/compute_pair_rmsd.py:cles` |
| BH-FDR | `analysis/common/stats.py` |
| cluster paths, conda envs, timestamped log names | `lib/common.sh` (`log_args`, `LOGS_DIR`, `*_ENV`) |
| chain-layout inspection | `stage2_redock/stage2_redock.py --inspect-only` |
| protein sequence extraction | not needed with `from_cif`; `stage1_bioemu` has `get_protein_sequence()` if it ever is |

The per-frame files are named `<pdb>_state_NNN.pdb` **specifically** so
`fnat_gate` globs them unchanged. That naming is load-bearing, not cosmetic.

---

## Reading the comparison honestly

Two traps, both printed by `compare_to_pipeline.py` at the end of every run:

**Lower BioEmu2 fnat is not automatically worse.** The frozen arm pairs every
frame with crystal-pose DNA, so its fnat is high partly by construction — the
DNA was never free to move. BioEmu2 samples the DNA too, so its fnat is an
actual measurement. What would be decisive is BioEmu2 reproducing the native
interface *while also* reproducing the crystal's DNA deformation — which stages
1–3 cannot test at all.

**Do not read state frequencies as populations.** BioEmu2's model card states
sampling frequency *"must not be interpreted as a calibrated Boltzmann
probability."* This repo's analysis reads frequencies off ensembles (best-state
counts, ensemble spread, bend IQR). That caveat applies to BioEmu 1 too but is
stated more sharply here, and it should be pinned down before any
ensemble-frequency claim leans on it.

**One live risk the comparison checks for.** `interface_rmsd.score()` silently
falls back to the *reference* DNA when a model's nucleotide count differs from
the reference (finding **F10** of `docs/codebase_review_2026-10-15.md`). The
existing pipeline carries crystal DNA so it never trips this; a BioEmu2 ensemble
samples the DNA and can. `compare_to_pipeline.py` counts and reports those
states rather than letting them pass as model-frame measurements.

---

## Files

```
bioemu2/
  README.md                 this
  environment.yml           Python 3.12 env (config/bioemu.yml is 3.11)
  sample_complex.py         System.from_cif -> sample -> project layout
  compare_to_pipeline.py    same-metric comparison vs stage3 states
  sample_complex.sbatch     SLURM wrapper; sources lib/common.sh
```
