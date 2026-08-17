# Chapter 9 — Stage 3: Energy Minimization (the physics stage)

> **What this chapter is.** Stage 3 is the pipeline's one physics-based refinement step. Stage 2
> hands over docked poses that are geometrically plausible but energetically strained; Stage 3
> relaxes that strain into a structure DeepPBS can featurize without reading clashes as chemistry.
> Every choice below serves one contract: **relax the strain, preserve the pose.**
>
> **Engine:** `stage3_minimize/stage3_minimize.py` (OpenMM), 409 lines, the script's own header calls
> it *"Stage 3 v5: DNA-aware minimization that preserves metal-coordination geometry through pairwise
> sidechain restraints (no metal in the simulated system)."* Driven on the cluster by
> `stage3_array.sh` (one SLURM array task per BioEmu frame) with `stage3_recover.sh` /
> `stage3_recover_array.sh` for failed states.
>
> **Sourcing.** Everything below is grounded in the actual engine and wrapper source read from
> `endeavour:/project2/rohs_102/shewchuk/TF-conformation/stage3_minimize/` and the pilot configs in
> `config/pilots/`. Line-level details (defaults, thresholds, force construction) are verbatim from
> `stage3_minimize.py`. The engine source is saved alongside this chapter as
> `src_stage3_minimize.py` for reference.
>
> **Each item is structured:** *why it's necessary → how it's used → how to run it → caveats / known issues.*

---

## MD-concepts primer (read first)

- **Force field** — a formula for the potential energy of a set of atoms given their positions: harmonic springs for bonds and angles, periodic terms for torsions, and non-bonded electrostatics + van der Waals.
- **Energy minimization** — rolls the structure *downhill* on that energy surface to the nearest low-energy geometry. It is **not** dynamics: no heating, no time evolution, no sampling of new conformations.
  - The engine constructs a `LangevinIntegrator(300 K, 1/ps, 1 fs)` only because OpenMM requires an integrator to build a `Simulation` object. It calls `sim.minimizeEnergy(...)` and **never** calls `sim.step(...)`, so the temperature and friction never actually act on the system.
- **The goal, stated narrowly:** take a clashy docked pose and settle it into the nearest physically reasonable one **without letting it wander away from the bound geometry we trust.**

---

## 9.0 — The pipeline shape (per state)

The engine processes **one state (one BioEmu frame) per invocation**, in six steps:

1. **PDBFixer** standard repair — *and*, before any heterogens are stripped, record every structural metal and which side-chain atoms coordinate it.
2. **Add hydrogens** (`Modeller`, forcefield = ff14SB + GB-Neck2).
3. **Phase 0** — hydrogen-only minimization with heavy atoms frozen (mass = 0).
4. **Phase 1** — vdW ramp (σ = 0.1 → 1.0) with backbone restraints *and* the metal coordination-cage restraints active.
5. **Phase 2** — standard final minimization at full vdW.
6. **Strip hydrogens, write** the minimized PDB.

Each state's array task reads `${STAGE2_DIR}/${PDB_ID}_state_${STATE}.pdb` and writes
`${STAGE3_DIR}/${PDB_ID}_state_${STATE}.pdb` (state = 3-digit zero-padded task id). The task is
**idempotent** — it skips if the output already exists or the input is missing — so re-submitting
only fills gaps. `N_FRAMES` in the pilot config sets the array size (e.g. tbp = 99, egr1 = 98).

---

## 9.1 — Why minimize a docked pose at all?

**Why it's necessary**
- Stage 2 hands over a BioEmu conformer Kabsch-fit onto the crystal DNA. The *global* placement is right, but the newly formed interface carries atomic clashes, stretched bonds, and strained contacts — artifacts of dropping an independently sampled protein onto DNA it never saw during sampling.
- DeepPBS reads *interface geometry*. Featurize a raw docked pose and those clashes get encoded as features — the model would learn simulation artifacts, not recognition chemistry.

**How it's used**
- Minimization removes clashes and relieves local strain, producing a structure that is *physically self-consistent* at the interface.
- **The exact clash diagnostic (from the engine).** `count_clashes()` counts **heavy-atom protein–DNA pairs closer than 0.24 nm (2.4 Å)**, and reports the minimum protein–DNA heavy-atom distance. Protein = the 20 standard residues (+ His protonation variants HID/HIE/HIP); DNA = DA/DG/DC/DT; hydrogens excluded. It is printed at `Initial`, after every ramp stage, and at `Final`.
- **Output tree.** Frozen (default) → `output/stage3_min/<tf>/<pdb>_state_*.pdb`; DNA-relaxed → the parallel `output/stage3_min_dnarelax/<tf>/`. The per-pilot fnat CSV (`<pdb>_fnat.csv`, columns `pdb_id, state, n_iface_res, n_segments, seq_ident, iRMSD_global, iRMSD_seg_max, iRMSD_seg_mean, fnat`) is what the fnat gate (Ch. 8) and the F/I figure series read — scored *post*-minimization.

**How to run it**
- Third `afterok` edge of the per-TF DAG: `run_pilot.sh <tf>` (stage 3). Under the hood each state is a SLURM array task running `stage3_array.sh`, which `conda activate`s the **`bioemu`** env and calls `stage3_minimize.py`.
- Direct single-state invocation (for debugging):
  ```bash
  conda activate bioemu
  python stage3_minimize/stage3_minimize.py \
      --input-pdb  output/stage2_docked/<tf>/<pdb>_state_001.pdb \
      --output-pdb output/stage3_min/<tf>/<pdb>_state_001.pdb \
      --ramp-stages 0.1,0.3,0.5,0.7,1.0 --steps-per-stage 500
  ```

**Caveats / known issues**
- Minimization is not neutral w.r.t. the fnat gate — it moves per-state fnat both ways (can lift a near-miss over the 0.5 floor or loosen a good pose). That is exactly why the gate is scored *after* this stage. Downstream the pass-rate metric is literally "fraction of Stage-3 states with `fnat ≥ 0.5`."
- "Physically self-consistent" is *local*, not global — the strain is gone, not the pose made correct. Correctness is still Stage 2's interface alignment's responsibility.
- **PDBFixer is told `missingResidues = {}`** — it repairs missing *atoms* within present residues but does **not** rebuild whole missing residues (disordered loops stay absent). And `removeHeterogens(keepWater=False)` drops waters too.

---

## 9.2 — The force field and the solvent model

**Why it's necessary**
- Minimization needs an energy function — the force field defines what "downhill" means.
- A solvent model is needed because a protein–DNA interface is highly charged; minimizing in vacuum would let electrostatics collapse the structure unphysically.

**How it's used (verbatim from the engine)**
- `forcefield = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")` — **AMBER ff14SB** for protein and DNA, **GB-Neck2** implicit solvent (water as a continuous dielectric).
- `nonbondedMethod = app.NoCutoff` — implicit solvent has no periodic box, so all pairwise interactions are computed exactly.
- **Constraints differ by phase:** Phase 0's system is built with `constraints=None`; the Phase 1/2 system uses `constraints=app.HBonds` (H-bond lengths constrained). Both use `rigidWater=False` (there is no water after `removeHeterogens`).
- **Platform:** tries `CUDA` with `{"CudaPrecision": "mixed"}` first, falls back to `CPU` if CUDA is unavailable; the chosen platform is printed per run.

**Caveats / known issues**
- **The load-bearing scientific caveat — GBSA over-compaction.** GB-Neck2 has no explicit counter-ions, so it *over-stabilizes close ion pairs*. At a charged protein–DNA interface the relaxed geometry can come out **artificially compact**. This is the Part III caveat, and it originates *here* — a property of the solvent model, not a bug.
- ff14SB parameterizes standard residues/nucleotides only; anything non-standard (metals, modified bases, cofactors) is outside its parameter set — the direct cause of the metal problem in 9.4.

---

## 9.3 — Restraints: hold what you trust, relax what you don't

**Why it's necessary**
- Unrestrained, minimizing a strained complex under GBSA can slide the *whole protein off the DNA* to a lower-energy but meaningless pose. The energy minimum is not the bound pose.
- The docked backbone *is* the frame's identity — the specific BioEmu conformer being tested. Letting it float discards the very thing the augmentation experiment varies.

**How it's used (verbatim from the engine)**
- A **harmonic positional restraint** implemented as a `CustomExternalForce` with energy `0.5 * k * ((x-x0)² + (y-y0)² + (z-z0)²)` — a virtual spring pulling each restrained atom back toward reference `x0`; stiffness `k` (kcal/mol/Å²) sets how hard.
- **What gets pinned:** protein backbone atoms `("N","CA","C")` and — by default — DNA backbone atoms `("P","C1'")`, at **k = 10 kcal/mol/Å²** (`--restraint-k`, default 10.0).
- **The reference positions are the post-Phase-0 coordinates** (`relaxed_positions`). Because Phase 0 froze all heavy atoms, those heavy-atom coordinates equal the docked coordinates — so in practice the backbone is pinned to its docked position.
- **What stays free:** side chains and everything else — they carry recognition chemistry and need to repack against the docked DNA.
- **Evidence it works (analysis side).** The R-series figures measure whole-protein backbone RMSD (N, Cα, C, O — the Cα proxy, since `per_state_rmsds.csv` has no pure-Cα column) to crystal, Stage 2 vs Stage 3. R1: the two arms nearly coincide. R3: the per-state Δ (stage3 − stage2) has a **median slightly negative** — minimization nudges frames marginally *toward* the crystal, never away.

**Caveats / known issues**
- The restraint reference is the *docked* pose, not the crystal — restraints preserve whatever Stage 2 produced, error included. Garbage docking in → garbage held.
- `k = 10` is a fixed choice, not a swept parameter; no sensitivity study behind the exact value.

---

## 9.4 — The metal-coordination cage

**Why it's necessary**
- Many DBDs — C2H2 zinc fingers, GATA, homeodomain-Zn — depend on a **structural** Zn²⁺/Mg²⁺ that holds the fold together. (EGR1/1aay, for instance, has three ZN in the crystal — chains D–F.)
- Two problems collide: (1) `removeHeterogens` strips *all* non-standard residues, silently including the metal; and (2) even if kept, ff14SB can't parameterize a bare metal ion.
- Doing nothing → the coordination shell collapses into an **unphysical apo-Zn fold** that then gets featurized.

**How it's used — no metal in the simulated system (verbatim from the engine)**
- **Detect** structural metals *before* `removeHeterogens`: `STRUCTURAL_METALS = {ZN, MG, MN, FE, CA, CO, NI, CU}`.
- **Record the coordination shell:** ligand residues `{CYS, HIS, ASP, GLU, SER, THR, TYR}` (+ His protonation variants HID/HIE/HIP/HSD/HSE/HSP); ligand atoms `{SG, ND1, NE2, OD1, OD2, OE1, OE2, OG, OG1, OH}`; a side-chain heavy atom counts as coordinating if it is within **`--metal-coord-cutoff` = 3.0 Å** of the metal. Keep a cluster only if **≥ 2 ligands** (a lone contact is not a shell). Ligands are recorded by `(chain.id, residue.id, atom.name)` so they survive PDBFixer/Modeller re-indexing.
- **Cage it:** build a single `HarmonicBondForce` and, for **every pair** of coordinating atoms in a cluster, add a bond at that pair's **original separation** with **k = `--metal-cage-k` = 20 kcal/mol/Å²**. This holds the coordination geometry rigid without the ion ever being in the system.
- **Reusable pattern:** for any cofactor a force field can't describe, restrain the *geometry the cofactor imposes*, not the cofactor.

**How to run it**
- On by default. Cage stiffness `--metal-cage-k` (20), detection radius `--metal-coord-cutoff` (3.0 Å).
- **Off:** `--ignore-metals` skips detection entirely and lets `removeHeterogens` strip the metal like any heteroatom — the pre-patch / `legacy` behavior, and the arm used for the Ch. 19 metal-cage A/B. At the config level, set `STAGE3_IGNORE_METALS=1` (the wrapper appends `--ignore-metals`).
- **Diagnostics printed per run:** the number of metals and clusters found, each cluster's ligand list, the number of cage restraints, and after minimization the **cage drift** — `mean` and `max` of |final − equilibrium| pair distances (Å). Low drift = the cage held.

**Caveats / known issues**
- The cage preserves *geometry*, not *chemistry*: no ion, no charge, no true coordination energetics — a geometric splint, only as good as the coordination it recorded.
- The 3.0 Å + ≥2-ligand rule is a heuristic; an unusual/distorted shell could be mis-detected — check the printed cluster list and cage drift before trusting a metal-dependent frame.
- If a recorded ligand atom can't be re-located after H-addition it is dropped with a `WARNING`, and a cluster that falls below 2 ligands is discarded.
- Impact is uneven: the cage matters most where minimization is otherwise *least* constrained (TBP in Ch. 19); for well-packed folds the downstream effect is small.

---

## 9.5 — The three-phase ramp (why not just minimize once)

**Why it's necessary**
- Docked poses start *clashy*. Overlapping atoms feel enormous vdW repulsion, and a single hard minimization on such a start frequently **diverges** rather than settling. The fix: resolve clashes gently — soften the atoms, grow them back.

**How it's used — three phases (verbatim from the engine)**
- **Phase 0 — hydrogens only.** Freeze every heavy atom (`setParticleMass(i, 0.0)`), minimize just the added H's (`maxIterations=2000`). Fixes H placement before it can perturb the heavy-atom frame. *Because heavy atoms are frozen, the heavy-atom clash count is unchanged through Phase 0* — the clash-trajectory figure carries the initial value forward as the "H-min" point.
- **Phase 1 — vdW ramp.** Scale every atom's vdW radius σ through **`--ramp-stages` = 0.1, 0.3, 0.5, 0.7, 1.0**, `--steps-per-stage = 500` minimization iterations each. Implemented by rescaling each particle's σ on the `NonbondedForce` (charge and ε untouched) and `updateParametersInContext`. Backbone + cage restraints active throughout. PE, clash count, and min distance printed per stage.
- **Phase 2 — final minimization.** At full vdW (σ = 1.0), `--final-iterations = 10000`. Then strip H's and write with `keepIds=True`.
- The principle: softening-then-growing (soft-core-like) brings a badly clashing model to a stable state without blow-ups.

**How to read the log (for the Ch. 9.1/9.5 figure)**
- The per-state log lines the clash-trajectory figure parses (all tagged `[<basename>]`):
  - `Initial heavy-atom clashes: <N>, min_dist: <D> Å`
  - `Stage <i>/<n>: σ=<s>, PE=<…>, clashes=<N>, min_d=<D> Å`  (one per ramp stage)
  - `Final: PE=<…>, clashes=<N>, min_d=<D> Å`
  - and a closing `Trajectory: clashes <init> -> <final> | min_d <init> -> <final> Å`.
  - Rebuild: `python make_clash_trajectory.py --log <stage3_log.out> --out F.png` (any env with matplotlib; worked example FOXA 1vtn state 001).

**Caveats / known issues**
- The ramp is why Stage 3 is compute-heavy: 5 stages × 500 iters + a 10k final, per state, per pilot.
- Soft-core softening temporarily lets atoms interpenetrate — correct by design, but intermediate ramp-stage structures are *not* physically meaningful; only the Phase-2 output is.
- **A failing ramp stage `sys.exit(2)`; a failing final `sys.exit(3)`** — so a hard-failing state produces *no* output PDB and drops out of the denominator downstream (feeds the pass-rate spread and Ch. 13's per-condition denominator drift). This is what the recovery path (below) is for.

**Failure recovery — the gentler ramp**
- Failed states are re-run with a **gentler, longer 6-stage recovery ramp**: `RECOVERY_RAMP_STAGES = 0.05, 0.1, 0.2, 0.4, 0.7, 1.0`, `RECOVERY_STEPS_PER_STAGE = 1000` (both set in the pilot config). This is a finer low-σ approach (starts at 0.05, adds a 0.2 and 0.4 rung) with double the iterations per rung.
- Two recovery drivers: `stage3_recover.sh` walks all `N_FRAMES` states **serially** in one job; `stage3_recover_array.sh` runs them as a **SLURM array** (one task per missing state). The array version exists because the serial one can hit its 4 h wall-clock when several states each need the expensive gentle ramp — this is exactly what stalled the DUX4 dimer trial (7 stragglers × ~30 min serial). Both are idempotent and **best-effort**: a state that still fails after recovery logs loudly and exits 0, so it never blocks the downstream `afterok` gate.

---

## 9.6 — DNA treatment: frozen default vs relaxed variant

**Why it's necessary**
- The single most consequential knob in Stage 3. DNA restraint stiffness decides *what the augmentation experiment tests*:
  - stiff DNA → "an ensemble of protein poses **around a fixed binding site**";
  - soft DNA → "an ensemble of **genuine binding modes**" (protein *and* DNA co-adapting).
- It is the physical embodiment of the book's central question — hence its own subsection.

**How it's used (verbatim from the engine)**
- **Frozen DNA (default, `--dna-restraint-k` omitted).** Protein *and* DNA backbone share **one** `CustomExternalForce` at `k = --restraint-k = 10`, reproducing the original behavior byte-for-byte. Output: `output/stage3_min/<tf>/`.
- **Relaxation variant (`_dnarelax`).** When `--dna-restraint-k` is given, DNA backbone atoms are **split onto their own** `CustomExternalForce` (global parameter `k_dna` — two forces can't share a parameter name). Two sub-cases:
  - `k_dna > 0` → DNA on its own soft tether at that stiffness (**1.5 kcal/mol/Å²** for the headline mode);
  - `k_dna = 0` → DNA omitted from *all* restraint forces (fully free from stage 1, no DNA force built).
  Output: the parallel `output/stage3_min_dnarelax/<tf>/`.
- **Late release (`--dna-release-stage N`, requires `--dna-restraint-k`).** Hold DNA at the *protein* k for ramp stages `1..N-1` (protein settles first under the largest clash-resolving forces), then drop `k_dna` to the target at stage `N` and keep it there through the final minimization. Implemented by `setParameter("k_dna", …)` per stage on the live context. **In the `_dnarelax` config, N = 5 — which is the *last* of the 5 ramp stages**, so DNA stays pinned through all four clash-resolving stages and relaxes only for the final ramp rung + the full-vdW final minimization. (A guard re-applies the target `k_dna` before Phase 2 in case the release stage was set beyond the ramp.)

**How to run it**
- **Engine knobs:** `--restraint-k` (protein), `--dna-restraint-k` (DNA split-off), `--dna-release-stage`, `--metal-cage-k`, `--metal-coord-cutoff`, `--ramp-stages`, `--steps-per-stage`, `--final-iterations`, `--ignore-metals`, `--scratch-dir` (default `/scratch1/shewchuk/deeppbs_min_tmp`).
- **Config form — the whole delta from `tbp.sh` to `tbp_dnarelax.sh` is two lines:**
  ```bash
  source "${_here}/tbp.sh"          # inherits TF_NAME, PDB_ID, chains, ramp, N_FRAMES, FOLD…
  export STAGE3_DNA_RESTRAINT_K=1.5
  export STAGE3_DNA_RELEASE_STAGE=5
  ```
  `TF_NAME` stays the base name (`tbp`), so Stage 1 (BioEmu library) and Stage 2 (docked frames) are **reused**; `common.sh` appends a `_dnarelax` suffix so Stage 3+ never overwrites the frozen baseline. Launch with `run_pilot.sh <tf>_dnarelax 3 5` (start at stage 3 — stage 2 is shared).
- **Two independent seed-paired experiments.** The frozen and relaxed pipelines each retrain their *own* baseline, so absolute Pearson is **not** comparable across the two — the only cross-treatment quantity is the within-pipeline, seed-matched ΔPearson (grey = baseline, teal = augmented·frozen, green = augmented·relaxed).

**Caveats / known issues**
- **k_dna = 1.5 is a deliberate stiffness *floor*, not a principled value.** The config comment states its purpose: a floor against GBSA-driven B-DNA melting, and the natural on-ramp to a Tier-2 cgDNA+ sequence-dependent stiffness prior (non-zero by construction). The principled replacement is Ch. 20; the flat number is the interim.
- **GBSA fraying risk** is the reason the floor exists — releasing DNA toward `k_dna = 0` under implicit solvent (no counter-ions) risks the duplex fraying.
- **Coverage is partial.** The `_dnarelax` arm has data for 7 pilots, trained models for 2, benchmark eval for 1 (`tbp_dnarelax`, itself slightly worse than frozen TBP). The eval JSONs confirm the asymmetry — 10 standard `id_benchmark_<tf>.json` vs 9 `_dnarelax` (csl has no relaxed run). Treat all DNA-relaxation results as directional / preliminary (the Ch. 18 caveat). Supporting validation tooling lives in `stage3_minimize/validation/` (`analyze_dna_relax.py`, `fnat_rescore_dnarelax.sh`, `smoke_test_dna_relax.sh`, `smoke_results/`).
- **Forward links:** structural consequences (DNA bends toward the bound pose without fraying) → Ch. 18; principled stiffness priors → Ch. 20.

---

## Diagnostics — what a green Stage-3 result looks like

Confirm each before trusting a pilot's minimized ensemble:

- **Clashes resolved.** Per-state log shows heavy-atom protein–DNA clashes (< 2.4 Å pairs) falling from hundreds at `Initial` to near-zero at `Final`, min distance opening up. Rebuild with `make_clash_trajectory.py --log …`.
- **Pose preserved.** Whole-protein backbone RMSD to crystal barely changes Stage 2 → Stage 3 (R1); per-state Δ (R3) has a **median ≤ 0** (toward crystal). A clearly positive R3 median means minimization is walking frames *off* the bound pose — investigate before featurizing.
- **Cage held.** Small mean/max cage drift for metal folds; cage-ON render shows intact coordination geometry, not a collapsed apo-Zn shell.
- **Reachability sanity.** `compute_reachability.py --source stage3` reproduces the anchors to ~2 decimals: **ets1 d_min ≈ 0.87, tbp d_min ≈ 0.59**; `dux4 ≈ 10 Å` is the far-unreachable outlier (0 fnat survivors → never trains). The script prints a `VALIDATION` line; a mismatch means the wrong `--source` or a broken correspondence.
- **Non-empty pass list** into the fnat gate; an empty one exits non-zero and halts the DAG by design (Ch. 8).

---

## Runbook summary (Stage 3 at a glance)

| Engine flag | Default | What it controls |
|---|---|---|
| `--restraint-k` | 10.0 | protein backbone harmonic restraint (kcal/mol/Å²) |
| `--dna-restraint-k` | None (= use `--restraint-k`) | DNA-backbone restraint, split onto its own `k_dna` force; `1.5` = soft tether, `0` = fully free |
| `--dna-release-stage` | None | 1-based ramp stage at which `k_dna` drops from the protein k to the target (requires `--dna-restraint-k`) |
| `--metal-cage-k` | 20.0 | coordination-cage `HarmonicBondForce` stiffness |
| `--metal-coord-cutoff` | 3.0 | Å within which a side-chain heavy atom counts as coordinating |
| `--ramp-stages` | 0.1,0.3,0.5,0.7,1.0 | vdW σ soft-core ramp schedule |
| `--steps-per-stage` | 500 | minimization iters per ramp stage |
| `--final-iterations` | 10000 | Phase-2 final minimization iters |
| `--ignore-metals` | off | disable the cage (Ch. 19 A/B "cage-OFF" / `legacy` arm) |
| `--scratch-dir` | /scratch1/shewchuk/deeppbs_min_tmp | scratch workdir |

- **Config-level knobs (pilot `.sh`):** `RAMP_STAGES`, `STEPS_PER_STAGE`, `RECOVERY_RAMP_STAGES` (`0.05,0.1,0.2,0.4,0.7,1.0`), `RECOVERY_STEPS_PER_STAGE` (`1000`), `N_FRAMES` (array size), and for the relaxed arm `STAGE3_DNA_RESTRAINT_K` / `STAGE3_DNA_RELEASE_STAGE`; `STAGE3_IGNORE_METALS=1` for cage-off.
- **Inputs / outputs:** reads `output/stage2_docked/<tf>/<pdb>_state_NNN.pdb`; writes `output/stage3_min[_dnarelax]/<tf>/<pdb>_state_NNN.pdb` + the per-pilot `<pdb>_fnat.csv`.
- **Structural variants tracked downstream:** `metal_cage` (canonical) and `legacy` (cage-off) in `per_state_rmsds.csv`.
- **Idempotent + recoverable:** array tasks skip states that already have output; failed states re-run via the gentler recovery ramp (serial `stage3_recover.sh` or parallel `stage3_recover_array.sh`).

### Cluster / environment
- **Cluster:** endeavour, `rohs` partition / `rohs_102` account; wrappers hardcode `--gres=gpu:rtx5000:1`. Repo: `/project2/rohs_102/shewchuk/TF-conformation`. Conda base `/apps/conda/miniforge3/24.11.3`.
- **Stage-3 env:** `bioemu` (the array wrapper runs `conda activate ${BIOEMU_ENV:-bioemu}`); the OpenMM/PDBFixer stack lives there. Analysis diagnostics use `deeppbs` (matplotlib) and `pycurves` (mdtraj, for `compute_reachability.py`/`compute_rmsds.py`); PyMOL renders need `pymol`.
- **Login-node gotcha:** cap BLAS threads before any login-node python — `export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`. PyMOL `cmd.png()` segfaults on the login node (no GL) — render via `sbatch` on `rohs` with `ray=1`.

## Coverage (Stage-3 output discovered)
- **13 pilots:** `csl dux4 egr1 engrailed err ets1 foxa hsf irf lef1 nfat runx tbp`.
- Structural fidelity (R1/R2/R3, reachability): all 13. `dux4` has **0 fnat survivors** (the excluded dimer, shown as the 0% reference).
- Benchmark accuracy (baseline vs augmented): 12 pilots (all but `dux4`, no eval).
- DNA-relaxed Stage-3 output: 7 pilots have data, 2 trained, 1 (`tbp_dnarelax`) evaluated.

## Figures for this chapter
- *(create)* **Clash-count + potential-energy trajectory** — parse one representative state's engine log (`make_clash_trajectory.py`; FOXA 1vtn state 001). — 9.1 / 9.5.
- *(create)* **Zinc-finger coordination shell, cage-ON vs cage-OFF** (intact tetrahedral Zn vs collapsed apo-Zn), cage-drift annotated. — 9.4.
- *(create)* **Ramp schedule plot** — PE + clash count vs σ stage (0.1→1.0), gentle descent vs single-hard-minimization divergence. — 9.5.
- *(create)* **Release-schedule schematic** — DNA k held at the protein value through ramp stages 1–4, dropped to 1.5 at stage 5 + final min. — 9.6.
- **Existing (caption + place):** `R1_ca_rmsd_stages`, `R3_minimization_delta`, `rmsd_analysis/minimization_motion` (local relaxation, nudges toward crystal — shared with Ch. 15). Full DNA-shape consequences are Ch. 18.
