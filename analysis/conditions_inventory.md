# TF-conformation — Conditions & Trials Inventory

Every experimental knob that has been exercised in the pipeline, with **what was actually
run** (has output on disk) vs **capability only** (script/flag present, no output tree).
Read-only survey of `/project2/rohs_102/shewchuk/TF-conformation`.

## A. Experimental axes (the "conditions")

| Axis | Settings exercised | Status |
|---|---|---|
| **Pilot / TF** | 11 structures (see §B) | RUN — 10 benchmarked, dux4 excluded |
| **Model arm** | baseline (crystal only) vs augmented (crystal + ensemble), always paired | RUN — every pilot |
| **DNA treatment** (Stage 3) | **frozen DNA** (default, backbone k=10) vs **DNA-relax** (soft tether k=1.5, released @ ramp stage 5) | RUN — frozen: all; dnarelax: data 6 / trained 2 / evaluated 1 |
| **Metal coordination cage** (Stage 3) | **cage ON** (default, sidechain-pair restraints k=20) vs **cage OFF** (`--ignore-metals`, "legacy") | **RUN in predecessor tree** — A/B executed for dux4/egr1/tbp (5 seeds) in `old_results` (May 2026); not re-run in current tree |
| **Training seeds** | single-seed vs **5 paired seeds** (s1–s5) | RUN — 5-seed for 6 pilots; single for csl/err/nfat/runx |
| **Cross-val fold** | FOLD=0 only | RUN — fold 0 exclusively |
| **Stage-2 alignment** | **interface** (default) · all/global · per_domain (diagnostic) | RUN — interface only in production; all/per_domain = capability |
| **fnat gate floor** | 0.5 (all pilots) | RUN — not varied |
| **DeepPBS feature set** (eval) | **full prot_shape** (used) · shape-only · groove · dnaseqInfo | **CAPABILITY ONLY** — 4 pred_configs exist, benchmark used full model only |
| **AF3 vs BioEmu** | AlphaFold3, 2 seeds × 5 samples/TF | RUN — 6 pilots, diversity comparison (not fed to training) |
| **DNA-relax smoke test** | A=frozen · C=dna_k=0 · D=tether+late-release (one TBP state) — B=dna_k=10 defined in script but **no output kept** | RUN — A/C/D only, validation/parameter-selection |

## B. Trial matrix — which conditions were run per pilot

| Pilot | PDB | Family | Std eval | 5-seed | DNA-relax data | DNA-relax trained | DNA-relax eval | AF3 |
|---|---|---|---|---|---|---|---|---|
| **egr1** | 1aay | C2H2 zinc finger | ✓ | ✓ (s1-s5) | ✓ | — | — | ✓ (2×5) |
| **engrailed** | 3hdd | Homeodomain | ✓ | ✓ (s1-s5) | ✓ | — | — | ✓ (2×5) |
| **ets1** | 1k79 | ETS | ✓ | ✓ (s1-s5) | ✓ | — | — | ✓ (2×5) |
| **foxa** | 1vtn | Forkhead | ✓ | ✓ (s1-s5) | ✓ | — | — | ✓ (2×5) |
| **lef1** | 2lef | HMG-box | ✓ | ✓ (s1-s5) | ✓ | ✓ | — | ✓ (2×5) |
| **tbp** | 1tgh | TBP / β-saddle | ✓ | ✓ (s1-s5) | ✓ | ✓ | ✓ (s1-s5) | ✓ (2×5) |
| **csl** | 3brg | CSL / Rel-like | ✓ | single | — | — | — | — |
| **err** | 1lo1 | Nuclear receptor | ✓ | single | — | — | — | — |
| **nfat** | 1a66 | Rel / NFAT | ✓ | single | — | — | — | — |
| **runx** | 1hjc | Runt | ✓ | single | — | — | — | — |
| **dux4** | 5z6z | Homeodomain (dimer) | — (0 survivors) | — | — | — | — | — |

✓ = run (output present) · — = not run · dux4 excluded (dimer, monomer-guard, 0 fnat survivors).

## C. Trained models actually on disk (`output/stage6_train/`)
239 subdirs. Distinct trained conditions (excluding `config_*`):
- `baseline_<tf>_fold0` / `augmented_<tf>_fold0` — 11 pilots, single-seed base run.
- `{baseline,augmented}_<tf>_fold0_s{1..5}` — 6 pilots × 5 seeds × 2 arms (egr1, engrailed, ets1, foxa, lef1, tbp).
- `{baseline,augmented}_lef1_fold0_dnarelax` + `{baseline,augmented}_tbp_fold0_dnarelax_s{1..5}` — the DNA-relax trained set.
- one stale `baseline_tbp_fold0_dnarelax_s3.incomplete_*` (aborted run).

## D. Completed evaluations (`output/stage7_eval/id_benchmark_*.json`)
11 benchmark files: **csl, egr1, engrailed, err, ets1, foxa, lef1, nfat, runx, tbp** (standard,
frozen DNA) + **tbp_dnarelax**. Each scores all discovered checkpoints on the 130-entry
general benchmark (pearsonr, spearmanr, auroc, ic_weighted_pcc, mae, ...).

## E. Conditions that exist as capability but have NOT been run *in the current tree*
1. **Stage-2 alignment baselines** — `all` (global) and `per_domain` modes plumbed through
   `STAGE2_ALIGN_MODE`, but production docking used `interface` only.
2. **DeepPBS feature ablations** — shape-only / groove / dnaseqInfo pred_configs exist; the
   benchmark used the full `prot_shape` model only.
3. **DNA-relax eval for the other 5 pilots** — data/configs built for egr1/engrailed/ets1/foxa/lef1,
   trained only for lef1 (+tbp); only tbp_dnarelax reached eval.
4. **Folds other than 0**, and **fnat floors other than 0.5** — never varied.

## F. Predecessor output tree — `/project2/rohs_102/shewchuk/TF_conf_init_outputs/old_results/`
A separate, **earlier** results tree (benchmark JSONs dated **May 2026**), superseded by the
current `TF-conformation/output/` (July 2026). It is where several conditions I earlier flagged
as "capability only" were in fact executed. Contents:

- **Metal-cage A/B (the key extra condition here).** `augmented_legacy_<tf>` = Stage 3 run with
  `--ignore-metals` (cage OFF), 5 seeds, for **dux4, egr1, tbp**, evaluated against the
  cage-ON `augmented_<tf>` and `baseline_<tf>` in `id_benchmark_<tf>_legacy_ab.json`
  (3 conditions × 5 seeds × 3 pilots). Mean Pearson over the 130-entry benchmark:

  | TF | baseline | augmented (cage ON) | augmented_legacy (cage OFF) |
  |---|---|---|---|
  | dux4 | 0.636 | 0.634 | 0.629 |
  | egr1 | 0.620 | 0.608 | 0.609 |
  | tbp  | 0.621 | 0.627 | 0.585 |

  *Reading: cage-ON augmented ≈ cage-OFF for dux4/egr1; for TBP cage-ON is clearly better than
  cage-OFF (0.627 vs 0.585) — the metal cage matters most where minimization is least
  constrained. Both still ≈ baseline.*

- **Pilots present:** dux4, egr1, tbp (5-seed each, standard + legacy) and a **wt1** pilot
  (WT1, C2H2 zinc finger) as `{baseline,augmented}_fold0` + `..._seed1` → `id_benchmark_wt1.json`
  — **two seed pairs, disagreeing in sign:** seed fold0 baseline 0.631 vs augmented 0.556
  (augmentation hurts), seed1 baseline 0.608 vs augmented 0.635 (augmentation helps). Mean over
  both seeds: baseline 0.619 vs augmented 0.595 — net slightly negative but seed-unstable, not a
  clean result.
- **`_legacy/` subdir** — an even earlier dux4/tbp baseline+augmented set.
- **`id_benchmark_comparison.json`** — the original single-seed baseline-vs-augmented comparison
  (the wt1 fold0 seed only: baseline 0.631 vs augmented 0.556; `_seed1` variant is a separate file).
- **`id_benchmark_merged.csv`** — 9,750-row long-format table merging all per-entry metrics
  across TFs/models/seeds/legacy flags (has a `legacy_ab` column).
- **`smoketest/`** — one early DeepPBS training smoke run (Model.2.tar, predictions).
- **`plots/`, `plots_by_family/`, `plots_egr1/`, `plots_pilot_family/`** — the earlier analysis
  figure sets (per-family metric summaries, egr1-vs-all, per-pilot family breakdowns).

**Bottom line:** `old_results` is the May-2026 predecessor of the current pipeline. It uniquely
holds (a) the **metal-cage on/off A/B** for dux4/egr1/tbp and (b) the original **wt1** pilot
(two seeds, sign-unstable) — neither of which exists in the current `TF-conformation/output/` tree.
