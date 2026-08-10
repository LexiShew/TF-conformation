# Chapters 10–11 — Stages 4–7: Featurization, Augmented Training Set, Paired Training, and Benchmark Evaluation

> **What these two chapters are.** Chapter 9 ended with a minimized, fnat-gated ensemble of
> physically self-consistent TF–DNA complexes. Chapters 10–11 turn that structural output into a
> trained model and a number. Chapter 10 (Stages 4–5) converts each passing structure into DeepPBS
> feature tensors and assembles the augmented training set. Chapter 11 (Stages 6–7) trains the
> matched baseline/augmented pair and scores both on the held-out benchmark. Every parameter below is
> verbatim from the engine sources on endeavour (`stage4_preprocess/`, `stage5_build_aug/`,
> `stage6_train/`, `stage7_eval/`), not paraphrased from the outline.
>
> **Sourcing note.** Grounded in the actual scripts read off the cluster:
> `process_co_crystal.py` (180 lines), `build_augmented_fold.py`, `build_combined_assembly.py`,
> `build_training_configs.py`, `driver.py` (418 lines), `evaluate_id_benchmark.py` (426 lines), the
> four SLURM wrappers, and a sample `id_benchmark_<tf>.json`. Where a value is read from one of these,
> it is marked **[code]**; where it is inferred from the outline or general DeepPBS knowledge, **[outline]**.
>
> **The through-line.** Stages 1–3 answered "can we build a trustworthy ensemble?" Stages 4–7 answer
> "does training on it change what DeepPBS predicts?" The honest headline, which Chapter 11 states and
> Part IV interprets: on the 130-entry cross-benchmark the net effect is small and mostly negative —
> and that is the setup, not the verdict.

---

# Chapter 10 — Stages 4–5: Featurization and the augmented training set

Stage 4 is the bridge from structure to tensor: it runs DeepPBS's own `process_co_crystal.py` on each
fnat-passing minimized PDB and emits a per-state `.npz` of protein-graph and DNA-shape features. Stage 5
assembles those per-state tensors into an augmented training fold and writes the matched baseline/augmented
training configs. Neither stage does any physics or learning — they are careful bookkeeping, and the care is
the point: this is where the shared-PWM-label confound (the one Part III and Part V both return to) is
physically created.

## 10.0 — Where Stages 4–5 sit

- **Input:** `${STAGE3_PASS_DIR}` — the fnat gate's symlink mirror holding only the Stage-3 states that
  cleared `FNAT_FLOOR`. **[code]** The Stage-4 wrapper links `${PDB_ID}_state_*.pdb` from that directory and
  refuses to glob the full Stage-3 output: *"Never glob the full STAGE3_DIR here — sub-floor states must not
  become training data (B7)."* **[code]**
- **Output of Stage 4:** one `.npz` per passing state under `${STAGE4_DIR}/output/`, each named
  `<pdb>_state_<i>_<PWM_ID>.npz`. **[code]**
- **Output of Stage 5:** a per-TF *combined assembly* directory (symlinks to the original training set + the
  new augmentation tensors) and a pair of training-config JSONs (`config_baseline_<tf>_fold0.json`,
  `config_augmented_<tf>_fold0.json`). **[code]**
- **Environment:** `deeppbs` conda env (`conda activate ${DEEPPBS_ENV:-deeppbs}`). **[code]** This is the
  second of the three envs — Stage 1 is `bioemu`/`hpacker`, Stage 3 is `bioemu`, Stages 4–7 are `deeppbs`.
- **Fail-loud:** the Stage-4 wrapper exits non-zero if zero gated inputs are present
  (`ERROR: No gated Stage 3 inputs`), so an empty gate halts the DAG rather than training on nothing. **[code]**

## 10.1 — Stage 4: what `process_co_crystal.py` actually computes

**Why it's necessary.** DeepPBS does not read PDB files at train time; it reads `.npz` feature tensors. Stage 4
is the only step that runs the DeepPBS featurizer, and it must run on the *minimized* structure so the features
reflect the relaxed interface, not the clashy docked pose.

**How it's used — the feature tensor.** For each line `<pdb_file>,<pwm_id>` in `input.txt`, the engine: **[code]**
1. Loads the structure and splits it into protein and DNA (`splitEntities`). **[code]**
2. Cleans each moiety — `cleanProtein(protein, add_charge_radius=True)` repairs missing atoms and attaches
   charge/radius; `cleanDNA(dna, fix_modified_nucleotide_hetflags=True)` normalizes modified nucleotides. **[code]**
3. Runs `processDNA` and asserts **exactly one helical segment** — the notorious guard
   `assert len(dna_helices) == 1, "number of helices present is not 1"`; a violation prints
   `ERROR: helix count problem <n>` and skips the state. **[code]** This is the self-complementary / symmetry-mate
   DNA failure the fix spec (Ch. 14) addresses upstream.
4. Builds the DNA coarse-grained bead cloud (`makeDNACG`) and a **14-column DNA shape matrix** `X_dna`: **[code]**
   - 6 intra-bp parameters: `buckle, shear, stretch, stagger, propeller, opening`
   - 6 inter-bp parameters: `shift, slide, rise, tilt, roll, twist` (mean-padded to length)
   - 2 backbone parameters: `major_groove_3dna, minor_groove_3dna` (midpoint-averaged)
5. Aligns the DNA sequence to the supplied PWM (`computeYAndMask`) to produce the label tensors `Y_pwm`,
   `pwm_mask`, `dna_mask`, and an alignment score. **[code]**
6. Builds the protein graph (`makeProteinGraph`) with per-atom features: **[code]**
   `charge, radius, SASA (getAtomSASA), Atchley factors (getAchtleyFactors), circular variance (getCV, radius 7.5)`.
7. Counts protein–DNA contacts (`countContacts`) and writes everything to
   `np.savez_compressed(<pdb>_<pwm_id>.npz, V_dna, X_dna, X_dna_point, V_prot, X_prot, E_prot, Y_pwm, pwm_mask,
   dna_mask, aln_score, contacts, ...)`. **[code]**

**How to run it.** **[code]**
```bash
# via the pipeline (SLURM wrapper resolves paths + config, deeppbs env):
sbatch --export=ALL,TF_NAME=<tf> wrappers/stage4_preprocess.sh
# what the wrapper does internally:
#   ln -sf ${STAGE3_PASS_DIR}/${PDB_ID}_state_*.pdb  ${STAGE4_DIR}/pdb_input/
#   build input.txt:  <state.pdb>,<PWM_LABEL>   (one line per passing state)
#   process_config.json = {"PDB_FILES_PATH":"./pdb","FEATURE_DATA_PATH":"./npz"}
#   python process_co_crystal.py input.txt process_config.json
```

**Caveats / known issues.**
- **The helix-count guard is a silent skip, not an error.** A state whose DNA parses to 0 or >1 helical segments
  is dropped with a printed `ERROR:` line and `continue` — it never becomes training data and quietly shrinks the
  per-TF denominator. **[code]** Watch the Stage-4 log line count vs the fnat-pass count.
- **The whole per-line body is wrapped in `try: … except: continue`.** **[code]** Any unhandled exception on one
  state is swallowed and the loop moves on. This is robust for batch processing but means a systematically broken
  input can produce a near-empty output with no non-zero exit. Diagnose by comparing `ls output/*.npz | wc -l`
  against `wc -l input.txt`.
- **The PWM label is attached here, per line, from `PWM_LABEL`.** **[code]** Every state of a complex is written
  with the *same* `pwm_id` — the physical origin of the shared-label confound (10.2, Ch. 12).
- **Features are computed on the minimized structure**, so the GBSA over-compaction caveat (Ch. 9.2) and the
  rigid-DNA restraint (Ch. 9.6) both propagate into `X_dna` — the DNA-shape channels encode whatever geometry
  Stage 3 produced.

## 10.2 — Stage 5: assembling the augmented fold (and creating the confound)

**Why it's necessary.** DeepPBS trains from a text file listing `.npz` entries and a data directory holding them.
The baseline model trains on the original crystal-only assembly; the augmented model must train on *that same set
plus* the new per-state tensors — and the comparison is only fair if the two sets are identical except for the
added frames. Stage 5 builds exactly that.

**How it's used — three scripts.** **[code]**
1. **`build_combined_assembly.py`** — makes a per-TF directory of symlinks: every original `assembly2024/*.npz`
   plus the new Stage-4 tensors filtered by PWM substring. *"Each pilot gets its own combined dir so they don't
   interfere."* Reports `n_orig`, `n_new (filtered)`, `n_total`. **[code]**
2. **`build_augmented_fold.py`** — copies the original train-fold list and appends the new `.npz` filenames
   (sorted, filtered by the `--pwm-filter` substring). It **enforces a trailing newline** on the original file
   *"to avoid the line-concatenation bug we hit during the EGR1 pilot,"* and warns if any new entry already
   appears in the original (leakage / double-counting guard). **[code]**
3. **`build_training_configs.py`** — writes the matched baseline and augmented config JSONs (10.3). **[code]**

**How the arms are kept comparable.** The baseline config's `data_dir` points at the original
`assembly2024`; the augmented config's `data_dir` points at the combined assembly. *"Both configs are identical
except for `data_dir` and `output_path`."* **[code]** Everything else — architecture, optimizer, epochs, and the
random seed — is shared.

**How to run it.** **[code]**
```bash
sbatch --export=ALL,TF_NAME=<tf> wrappers/stage5_build_aug.sh   # 10 min, 4 GB, no GPU
# internally: build_combined_assembly.py → build_augmented_fold.py → build_training_configs.py
```

**Caveats / known issues.**
- **The PWM-substring filter is load-bearing and silent on mismatch.** Both `build_combined_assembly.py` and
  `build_augmented_fold.py` select new entries by `if pwm_filter in filename`. **[code]** A wrong or empty filter
  matches zero entries; `build_augmented_fold.py` exits non-zero (`no .npz files matched filter`) but
  `build_combined_assembly.py` only prints a `WARNING` — so a bad filter can produce a "combined" dir that is
  just the original set, silently training an "augmented" model on no augmentation. **[code]** Verify `n_new > 0`
  in the Stage-5 log.
- **This is where the shared-PWM-label confound is physically created.** **[code]** Every appended frame carries the
  one crystal PWM label. The augmented fold is therefore *N frames → one label*; "does the ensemble help" is
  entangled with "does averaging N frames onto one label just regularize." Chapters 12 and 21 return to this; the
  decisive control (per-frame labels) is Ch. 21.
- **Symlinks, not copies.** The combined assembly is symlinks into `assembly2024/` and the Stage-4 output. **[code]**
  Moving or pruning either source silently breaks the training set. Keep the Stage-4 output in place until training
  completes.
- **Leakage guard is a warning, not a stop.** If a new entry name collides with an original, you get a printed
  `WARNING` and duplicates in the fold, not an abort. **[code]**

## 10.3 — The per-TF training config (baseline vs augmented)

**Why it matters.** The config is the contract that makes the A/B fair. Read from `build_training_configs.py`'s
`base_config()`, the two configs differ in exactly two fields and share the rest: **[code]**

| field | baseline | augmented | note |
|---|---|---|---|
| `data_dir` | `assembly2024` (orig) | `combined_assembly_<tf>` | the only substantive difference |
| `output_path` | `…/baseline_<tf>_fold0<suffix>` | `…/augmented_<tf>_fold0<suffix>` | where the run lands |
| `random_seed` | `--seed` (paired) | **same `--seed`** | *"paired seed so baseline and augmented share parameter init and shuffle order"* |
| `epochs` | 50 | 50 | **[code]** |
| `batch_size` | 1 | 1 | **[code]** |
| `loss` | `soft_ce` | `soft_ce` | **[code]** |
| `condition` | `prot_shape` | `prot_shape` | protein-shape readout **[code]** |
| `best_state_metric` | `mae` on `validation`, goal `min`, threshold `1.0` | same | checkpoint-selection rule **[code]** |
| `optimizer` | adam, lr 1e-3, weight_decay 1e-4 | same | **[code]** |
| `labels_key` | `Y_pwm` | `Y_pwm` | **[code]** |

**The paired-seed design.** The config comment is explicit: *"CRITICAL: paired seed so baseline and augmented
share parameter init and shuffle order. The augmentation effect is then isolated from init/shuffle noise."* **[code]**
Multi-seed runs pass `--seed 1..5` with `--seed-suffix _sN`, producing `baseline_<tf>_fold0_s3` /
`augmented_<tf>_fold0_s3` pairs. **[code]**

**The statistical caveat that lives here (Ch. 13).** The *default* seed is 42 and the generated configs set
`no_random: false`. **[outline+code]** In the single-seed regime, baseline and augmented within a seed share init
and shuffle (good), but *across* the small number of seeds the augmentation effect is not cleanly separable from
seed noise — which is why the honest analysis (Ch. 16) uses seed as the unit of replication and why the pooled
"power-recovered" p-values are pseudoreplication (Ch. 13).

## Figures to add — Chapter 10
- *(create)* **Feature-tensor schematic** — what one `.npz` contains: the protein graph (`V_prot`, `X_prot`,
  `E_prot` with charge/SASA/Atchley/CV channels) and the DNA point cloud (`X_dna_point`) + 14-column shape matrix
  (6 intra-bp, 6 inter-bp, 2 groove). Anchors 10.1.
- *(create)* **The augmented-fold diagram** — 1 crystal + N passing frames all pointing at one shared PWM label,
  drawn to foreshadow the confound. Reused in Ch. 12. Anchors 10.2.
- *(create)* **Stage-4 yield funnel** — fnat-pass states → helix-guard survivors → written `.npz` per pilot, to
  make the silent-skip attrition visible. Data: Stage-4 log line counts vs fnat-pass counts.

---

# Chapter 11 — Stages 6–7: paired training and benchmark evaluation

Stage 6 trains the matched baseline/augmented pair with DeepPBS's `driver.py`; Stage 7 scores both checkpoints on
the held-out 130-entry `id.txt` benchmark and computes the paired ΔPearson statistics that Part IV interprets.
This is where the project's headline number is produced — and where the statistics must be read carefully, because
the evaluator computes *two* kinds of p-value and only one of them is honest.

## 11.0 — Where Stages 6–7 sit

- **Input:** the two config JSONs and the combined assembly from Stage 5. **[code]**
- **Stage 6 output:** two run directories, each with `Model.best.tar`, `config.json`, `scaler.pkl`, `run.log`,
  and `predictions/`. **[code]**
- **Stage 7 output:** one `id_benchmark_<tf>.json` per pilot under `output/stage7_eval/`, holding
  `results`, `paired_stats`, `common_entries`, `subset_entries`. **[code]**
- **Environment:** `deeppbs`, on the `rohs` GPU partition (`--gres=gpu:1`, account `rohs_102`). **[code]**

## 11.1 — Stage 6: training with `driver.py`

**Why it's necessary.** This is the only stage that fits model weights. The baseline learns TF→PWM from crystal
structures alone; the augmented model learns from crystal + ensemble. The comparison of their held-out accuracy is
the entire experiment.

**How it's used.** `driver.py` takes a train file, a validation file, and a config; it: **[code]**
- Loads train/validation `.npz` sets (`loadDataset`, `balance="unmasked"`, feature scaling on). **[code]**
- Builds the DeepPBS model (`models.model_v2.Model(nF_prot, nF_dna, condition="prot_shape")`). **[code]**
- Trains with Adam (lr 1e-3, weight_decay 1e-4 via `addWeightDecay`, which excludes biases and 1-D params from
  decay), cross-entropy loss, for `epochs` (50). **[code]**
- **Selects the best checkpoint by `best_state_metric`** — MAE on the *validation* set, minimized, threshold 1.0.
  **[code]** This is the statistical caveat of Ch. 13: the validation set is a general-TF set, not the target TF,
  so the selected checkpoint is "best on general validation," fair across arms but not "best for TF X."
- Writes `Model.best.tar`, the fitted `scaler.pkl`, and per-entry `predictions/*.npz`. **[code]**

**How to run it.** **[code]**
```bash
sbatch --export=ALL,TF_NAME=<tf> wrappers/train_compare.sh
#   #SBATCH --array=0-1  → array task 0 = baseline, task 1 = augmented (the matched pair)
#   #SBATCH --gres=gpu:1 --mem=32GB --time=04:00:00 --partition=rohs --account=rohs_102
# multi-seed: run_multiseed_pilot.sh drives seeds s1..s5, each a fresh 0-1 array pair
```

**Caveats / known issues.**
- **Checkpoint selection on the wrong distribution.** `best_state_metric_dataset = "validation"` (general TFs),
  not the target TF. **[code]** Fair across arms, but the absolute numbers aren't "the best achievable for TF X."
- **Seed handling.** The paired seed lives in the config (`random_seed`), not on the CLI; `no_random: false` means
  the run is seeded but training still has GPU-nondeterministic ops. **[code]** The A/B is protected within a seed
  (shared init/shuffle) but not made bit-reproducible.
- **`try/except` is absent here** — a training failure exits non-zero and the `afterok` DAG edge halts, which is
  the desired fail-loud behavior for the expensive stage. **[code]**
- **The model is `model_v2`, hard-wired `condition="prot_shape"`.** **[code]** Feature ablations (shape-only,
  groove, dnaseqInfo) exist as capability but are never exercised (Ch. 14½, Ch. 22).

## 11.2 — Stage 7: benchmark evaluation and the two kinds of p-value

**Why it's necessary.** Training produces weights; Stage 7 produces the number. It runs both checkpoints over the
fixed 130-entry `id.txt` benchmark, computes per-entry metrics, and — crucially — the *paired* Δ-statistics that
answer "did augmentation help."

**How it's used — `evaluate_id_benchmark.py`.** **[code]**
- Loads every `--condition name=run_dir` (baseline/augmented, all seeds), each scored on the same `id.txt`
  entries with the training run's own `scaler.pkl`. **[code]**
- **Hard-checks denominators** — `check_denominators()` intersects the entry sets across conditions and restricts
  all aggregates to the common set, printing a warning if any condition evaluated a different set. **[code]** This
  is the fix for the per-condition denominator drift caveat (Ch. 13): a NaN/failed entry in one arm shrinks the
  common set for *all* arms, so the Δ compares means over identical entries.
- Reports six metrics: `auroc, mae, pearsonr, ic_weighted_pcc, spearmanr, brier_multi`. **[code]**
- **Computes per-entry paired Δ (augmented − baseline)** with a bootstrap 95% CI (10,000 resamples, seeded 42)
  and a paired one-sample t-test on the deltas. **[code]**
- Writes `id_benchmark_<tf>.json` with `results` (per-condition per-entry metrics), `paired_stats` (the Δ stats
  per pair), `common_entries`, and `subset_entries`. **[code]**

**How to run it.** **[code]**
```bash
sbatch --export=ALL,TF_NAME=<tf> wrappers/eval_benchmark.sh   # 30 min, GPU, deeppbs env
# internally, for each seed pair:
#   python evaluate_id_benchmark.py \
#     --condition baseline_<tf>_fold0_sN=<run_dir> \
#     --condition augmented_<tf>_fold0_sN=<run_dir> \
#     --combined-dir <combined_assembly> --filter '<family regex>' \
#     --filter-name <family> --output-json id_benchmark_<tf>.json
```

**The two p-values — read this before quoting significance.** The eval JSON's `paired_stats` is computed over the
`--filter` **own-family subset** (`subset_entries`), not the full 130-entry benchmark, and for each
`baseline_X / augmented_X` pair it holds an entry-level paired t-test. There are **two ways to read the resulting
numbers, and they disagree**: **[code + data]**
- **The pooled key** (e.g. `augmented:csl_fold0`) reports the paired t-test over the entry×seed rows and can look
  very significant — csl's own-family subset shows **ΔPearson +0.199, p = 0.001**. **[data]**
- **The honest per-seed reading** (`_s1 … _s5`) treats *seed* as the experimental unit (augmentation is applied
  once per seed). The five csl own-family seed deltas are `−0.024, −0.043, +0.159, +0.162, +0.164` — **two of five
  negative**. **[data]**

The pooled p-value is **pseudoreplication**: it counts the own-family entries × 5 seeds as independent when the
independent replicate is the seed. **[code comment + Ch. 13]** The evaluator *provides* the per-seed keys precisely
so the seed-level test can be done; the honest number is the seed-level one. This is the fig9 cautionary case of
Ch. 13, live in the data. *(The authoritative per-pilot own-family effects — recomputed with seed as the unit over
6 seeds — are the Chapter 16 / `FINDINGS.md` table; the five-seed slice shown here is only to make the pooled-vs-
seed contrast concrete, not a competing effect estimate.)*

**Caveats / known issues.**
- **Pair-matching is by regex.** Conditions are paired by `^(augmented_legacy|baseline|augmented)_(.+)$`; anything
  not matching `baseline|augmented` is silently dropped from `paired_stats` with a printed warning. **[code]**
  A misnamed run vanishes from the summary without failing the job.
- **`Model.best.tar` missing → condition skipped.** A condition whose run dir lacks the checkpoint prints
  `WARNING: Model.best.tar not found` and is skipped, not errored. **[code]** Check that every seed's checkpoint
  exists before trusting a cross-seed aggregate.
- **The subset filter defines "own family."** `--filter '<regex>'` selects the `subset_entries` the paired stats
  are computed over — this is the own-family subset that Part IV shows is where the signal actually lives. A wrong
  regex silently changes what "the effect" means. **[code]**
- **dux4 has no eval JSON** — 0 fnat survivors means no Stage-4 output, no training, no eval. **[data]** 12 of 13
  pilots are benchmarked; the monomer-only scope (Ch. 14) excludes the dux4 dimer.

## 11.3 — What the benchmark actually returns (the honest headline)

On the **130-entry cross-benchmark**, averaged over the whole benchmark, the net augmentation effect is small and
mostly flat-to-negative — for essentially every pilot the cross-benchmark ΔPearson sits near zero and is not
significant at the seed level. **[data]** This is the number the pipeline produces if you read the eval as a single
scalar, and taken alone it reads as "augmentation doesn't work."

Chapter 11 deliberately does **not** tabulate per-pilot effect estimates here — there is exactly one authoritative
set of those, computed with seed as the unit of replication over the own-family subset, and it lives in
**Chapter 16 (from `analysis/mechanism/FINDINGS.md`)**. Quoting a second, differently-scoped set of numbers in
this chapter is how the pooled-vs-seed and cross-benchmark-vs-own-family readings get conflated. The only figures
Chapter 11 should quote are: (a) the pooled-vs-seed contrast for a single pilot as the pseudoreplication example
(csl own-family, §11.2), and (b) the qualitative cross-benchmark headline (near-zero, non-significant).

The point Chapter 11 hands to Part IV: **a flat cross-benchmark mean is not "it doesn't work" — it is the signal
that a hidden variable is splitting the families.** The per-pilot resolution of that null (the own-family effect,
ETS1 +0.111 p = 0.013 and the rest of the table) is Chapter 16; the hidden variable itself is Chapter 17.

## Figures to add — Chapter 11
- **Existing (caption + place):** `figure_scripts/fig3_box_pearson`, `fig4_box_mae` — per-entry metric
  distributions over the 130-entry benchmark, three arms. The "what the eval produces" reference.
- *(create)* **The two-p-values contrast** — the csl own-family example drawn twice: the pooled entry×seed t-test
  (ΔP +0.199, p = 0.001) beside the five per-seed deltas (`−0.024, −0.043, +0.159, +0.162, +0.164`; 2 of 5
  negative). The single clearest picture of the pseudoreplication trap; shared with Ch. 13.
- *(create)* **The paired-training schematic** — one seed → shared init/shuffle → baseline vs augmented diverge
  only by `data_dir`; the `--array=0-1` pair and the `_sN` seed stack drawn as the experimental design.
- *(create)* **Stage-6/7 flow** — config pair → array train → checkpoint-on-validation-MAE → id.txt eval →
  denominator-intersection → paired Δ. The runbook as one diagram.

---

## Consolidated runbook (Stages 4–7)

| stage | wrapper | engine | env | resources | key output |
|---|---|---|---|---|---|
| 4 | `wrappers/stage4_preprocess.sh` | `process_co_crystal.py` | deeppbs | CPU | `<pdb>_state_i_<PWM>.npz` |
| 5 | `wrappers/stage5_build_aug.sh` | `build_{combined_assembly,augmented_fold,training_configs}.py` | deeppbs | 4 GB, 10 min | combined assembly + config pair |
| 6 | `wrappers/train_compare.sh` | `driver.py` | deeppbs | GPU, 32 GB, 4 h, `--array=0-1` | `Model.best.tar` ×2 |
| 7 | `wrappers/eval_benchmark.sh` | `evaluate_id_benchmark.py` | deeppbs | GPU, 8 GB, 30 min | `id_benchmark_<tf>.json` |

**Engine flags worth knowing** (all **[code]**):
- Stage 4: `process_co_crystal.py <input.txt> <process_config.json> [--no_pwm] [--no_cleanp]`;
  `process_config.json = {"PDB_FILES_PATH", "FEATURE_DATA_PATH"}`.
- Stage 5: `build_augmented_fold.py --orig-train --stage4-dir --pwm-filter --out-train [--no-overwrite]`;
  `build_combined_assembly.py --orig-dir --stage4-dir --pwm-filter --out-dir`;
  `build_training_configs.py --tf-name --combined-dir --seed --output-dir [--seed-suffix]`.
- Stage 6: `driver.py <train> <valid> -c <config> [--single_gpu] [--epochs] [--batch_size]`.
- Stage 7: `evaluate_id_benchmark.py --condition name=dir (repeatable) --combined-dir --filter --filter-name
  --output-json [--bootstrap-iters 10000]`.
