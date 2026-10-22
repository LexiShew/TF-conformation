# bioemu2 — sampling protein–DNA complexes directly

An **exploratory parallel track**, not a pipeline stage. It tests whether
[BioEmu2](https://github.com/microsoft/bioemu2) can sample a protein–DNA
**complex** well enough to replace stages 1–3, and it is kept outside the
numbered DAG so nothing in `run_pilot.sh` depends on it.

Status as of 2026-10-22: **running, with a known size ceiling.** Both blockers below are resolved - the
`bioemu2` env exists at `/project2/rohs_102/shewchuk/conda/envs/bioemu2`
(Python 3.12), and an L40S (48 GB) on `qcbgpu` runs these systems fine
despite the model card's 80 GB A100 floor. A full sweep of all 698
`structures/source_chains/` entries is in flight - see *Batch: all
source_chains* below.

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


---

## Batch: all source_chains

The three scripts above are per-pilot. `survey_source_chains.py`,
`sample_all.sbatch` and `launch_all.sh` extend the track to **every** structure
in `structures/source_chains/` (698 of them) without touching the pilot path.

```bash
# 1. size every deposited assembly (~1 min, login node; needs gemmi => pycurves env)
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /project2/rohs_102/shewchuk/conda/envs/pycurves/bin/python \
  bioemu2/survey_source_chains.py --out bioemu2/worklist/source_chains_survey.tsv

# 2. build the worklist and submit the array
bash bioemu2/launch_all.sh                    # all ok structures, 8 concurrent
LIMIT=20 bash bioemu2/launch_all.sh           # smoke test: first 20 only
MAX_CONCURRENT=4 NUM_SAMPLES=50 bash bioemu2/launch_all.sh
DRY_RUN=1 bash bioemu2/launch_all.sh          # print the sbatch, submit nothing
```

**Why a survey before the array.** `System.from_cif()` imports the whole
deposited assembly and BioEmu2's model card states a ~10k-token ceiling, so
deposited size predicts whether a task takes minutes, takes hours, or OOMs.
One minute of login-node CPU turns an opaque 698-job queue into a worklist
ordered small-first, with the risky tail identified up front:

| tokens_proxy (protein res + nt) | structures |
|---|---|
| median | 276 |
| > 2000 | 70 |
| > 4000 | 62 |
| > 6000 | 35 |
| > 10000 (past the stated ceiling) | 11 |

`n_tokens_proxy` is a lower bound used for **ordering and flagging only** - it
counts polymer residues, not BioEmu2's real tokenization. Do not quote it as a
token count.

**Why `gpu:l40s:1` on `qcbgpu`.** Not the `gpu:rtx5000:1` the numbered stage
wrappers use, and not `sample_complex.sbatch`'s generic `gpu:1`. The one
verified run on this cluster (job 6537427, 1k79/ets1, 100 samples) completed in
9:44 with 5.7 GB peak RSS on an L40S, which settles the model card's "80 GB
A100 is the smallest tested device" question in practice for systems this size.
L40S is also the biggest card this account reaches in quantity (22 across
`qcbgpu`), so it is the throughput choice. `MAX_CONCURRENT` defaults to 8
rather than all 22 because `qcbgpu` is shared with the whole qcb community.

**Resumability.** `sample_all.sbatch` skips any structure whose
`sample_manifest.json` already reports the requested frame count, so a
part-finished array can be resubmitted wholesale. Per-structure outcomes land
in `output/bioemu2/_status/<pdb>.status` (`ok <secs>` / `failed rc=<n>` /
`skipped`), which is what to aggregate rather than re-reading 698 SLURM logs.

### Known caveat this batch inherits

Some deposited entries hold only **one strand** of a self-complementary duplex,
with crystallographic symmetry generating the partner (task 1, `2dgc`, imported
as 1 protein + a single 19 nt DNA entity). BioEmu2 then samples protein + ssDNA,
which is not the complex this track is about. The survey's `n_poly_chains`
column is what to filter on if that matters for a given analysis; the batch
does not filter on it, because a single-strand entry is still a valid BioEmu2
run, just not a duplex one.

### Post-processing is still pilot-only

`extract_one_copy.py` requires `--tf-name` and reads chain selection from
`config/pilots/<tf>.sh` via `interface_rmsd.ref_side`, and
`compare_to_pipeline.py` needs stage3 arms to compare against. Neither runs for
a bare source_chains entry - the batch generates ensembles only. Scoring the
non-pilot 685 needs a config-free way to pick the cognate protein+duplex copy,
which does not exist in the repo yet.

---

## Two hard limits found while running the batch (2026-10-22)

Both were discovered by running it, and both change how the batch must be driven.

### 1. The MSA step rate-limits before the GPU does any work

`equilibrium` mode calls `get_msa(system)` before diffusion, which queries the
**public** ColabFold server. Running 14 concurrent sampling tasks (8 on one
partition + 6 on another) produced 20 immediate failures:

```
requests.exceptions.HTTPError: 429 Client Error: Too Many Requests
    for url: https://api.colabfold.com/ticket/msa
```

Every one died ~31 s in, before touching the GPU. So the array's concurrency
ceiling is set by a public web API, not by available GPUs. bioemu2's own retry
(`stop_after_attempt` + `wait_fixed(5)`) is far too tight for a 429.

**Fix: fetch MSAs serially, once, outside the GPU jobs.**

```bash
sbatch --partition=qcbr --account=rohs_102 bioemu2/warm_msa_cache.sbatch
```

`warm_msa_cache.py` walks the worklist on one CPU with real exponential backoff
(30/60/120/300/600 s) and fills the cache. It is resumable — a cached system
returns in well under a second — so re-running after a partial pass costs
minutes. Dedup is free: the cache key is a sha256 over the system's *unique
protein sequence set*, so structures sharing a protein (`5iy6`-`5iyd`,
`6cnb`-`6cnd`, the `1a1*` zinc-finger series) collapse to one entry and one
request. Empirically ~0.6 fetches per structure.

The cache lives at `output/bioemu2/_msa_cache`, and
`~/.cache/bioemu2/colabfold_msa` is a **symlink** to it so jobs that did not set
`BIOEMU2_MSA_CACHE` still share one cache. `sample_all.sbatch` also staggers task
starts by `(task_id % 8) * 15 s`, which only matters on a cold cache.

### 2. A CUDA OOM ceiling near ~4000 tokens — not a partition problem

Every structure above 5000 tokens failed identically after ~35 min of real work:

```
File ".../bioemu2/biomoa/modules/conditioning.py", line 121, in forward
  pair_method = pair_method_per_token[:, :, None] + pair_method_per_token[:, None, :]
torch.OutOfMemoryError: CUDA out of memory. Tried to allocate 55.52 GiB.
GPU 0 has a total capacity of 44.39 GiB ... this process has 41.28 GiB in use
```

| pdb | tokens_proxy | outcome on L40S (46 GB) |
|---|---|---|
| 6v92 | 5121 | OOM |
| 6kw4 | 5193 | OOM |
| 5iyd | 5573 | OOM |
| 5iyc | 5737 | OOM |
| 6cnf | 5861 | OOM |
| 6cnc | 5869 | OOM |

This is the **N² pair representation in the trunk embedder**, computed once per
system. Three consequences worth being precise about:

- **Lowering `--num-samples` does not help.** The allocation is per-system, not
  per-frame. Only the token count moves it.
- **A bigger card barely helps.** The 5869-token system needs ~97 GB
  (41 in use + 55.5 requested), so an 80 GB A100 fails too. Only an H200-class
  card has a chance, and the `ce_hopper`/`isi` H200 partitions reject
  `rohs_102` — the one reachable A100 is on `qcbgpu` with a ~1 week queue.
- **The threshold used for the split (4000) is an estimate, not a measurement.**
  Scaling the failed allocation by N² puts the L40S ceiling near ~3800-4000
  tokens, so structures approaching 4000 are marginal. The batch runs
  small-first, so the true cutoff will be read off the failures at the tail of
  `output/bioemu2/_status/`.

**Practical scope: the ~636 structures at or below 4000 tokens.** The 62 above it
are out of reach on hardware this account can get, and were split into
`worklist/source_chains_large.tsv` rather than deleted so they are easy to retry
if H200 access appears.

### 3. rtx5000 works but costs 3-4x, because of a missing kernel cache

`rohs`'s Quadro RTX 5000 (16 GB, compute capability 7.5) *does* run the small
systems — 16 GB is not the binding constraint for them. But every task logs:

```
Error reading system-wide triton tuning cache file:
  .../cuequivariance_ops/triton/cache/fused_sigmoid_gated_dual_gemm_forward_kernel_wrapper.7.5.json
```

`cuequivariance` ships no tuned Triton kernels for cc 7.5 and never persists an
autotune (`~/.cache/cuequivariance-triton/` stays empty across tasks), so the
cost is paid every time, not once:

| GPU | model load | 100 samples, ~104-token system |
|---|---|---|
| L40S (cc 8.9, qcbgpu) | 3-28 s | 3.9-5.2 min |
| Quadro RTX 5000 (cc 7.5, rohs) | 309 s | 14.7-17.4 min |

`rohs` capacity compounds it: of its two GPU nodes one sits in `drain`, so the
practical concurrency is 3-4 tasks, against 8 on `qcbgpu`. Net throughput for
the small bulk is roughly 8-10x lower. The batch therefore keeps the bulk on
`qcbgpu` and leaves a 60-structure trickle on `rohs` so the lab's own cards are
not idle.

### Current split

| array | worklist | target | throttle |
|---|---|---|---|
| bulk, small | `source_chains_qcbgpu.tsv` (546) | `qcbgpu` `gpu:l40s:1` | `%8` |
| trickle | `source_chains_rohs.tsv` (60) | `rohs` `gpu:rtx5000:1` | `%3` |
| out of reach | `source_chains_large.tsv` (62) | needs >80 GB | not submitted |

Resource request was trimmed to `--cpus-per-task=4 --mem=32GB` after measuring
5.7 GB peak RSS; the job is GPU-bound and the original 8 CPU / 96 GB request was
blocking task packing on shared nodes.

**Finish with a sweep.** `launch_all.sh` skips completed structures, so after the
arrays drain, re-running it picks up the 429 casualties and any other stragglers:

```bash
WORKLIST=bioemu2/worklist/sweep.tsv TOKEN_MAX=4000 \
  PARTITION=qcbgpu GRES=gpu:l40s:1 MAX_CONCURRENT=8 bash bioemu2/launch_all.sh
```
