# TF-conformation — Results Inventory

**Repo:** `/project2/rohs_102/shewchuk/TF-conformation` (endeavour HPC)
**Survey date:** by Claude, read-only walk of the pipeline output tree.

## What the project is
Augments **DeepPBS** (protein structure → binding-specificity PWM) with **conformational
ensembles** of monomeric TF–DNA complexes. Per TF: BioEmu backbone sampling → HPacker
side chains → interface-aligned re-dock onto crystal DNA → fnat interface gate → OpenMM
minimization → DeepPBS `.npz` featurization → paired **baseline vs augmented** DeepPBS
training → benchmark eval on the general-130 held-out set.

## Pipeline stages (per pilot)
| Stage | Dir | Output |
|---|---|---|
| 1 | `stage1_bioemu` | BioEmu + HPacker ensemble per protein chain |
| 2 | `stage2_redock` | Kabsch dock of each frame onto crystal DNA |
| 2g | `fnat_gate` | drop states below fnat floor (0.5) |
| 3 | `stage3_minimize` | OpenMM minimization (+ optional DNA relaxation) |
| 4 | `stage4_preprocess` | per-state `.npz` DeepPBS features |
| 5 | `stage5_build_aug` | augmented fold + combined assembly + train configs |
| 6 | `stage6_train` | baseline vs augmented DeepPBS (paired, multi-seed) |
| 7 | `stage7_eval` | benchmark eval; per-entry metrics JSON |

## Pipeline versions / conditions found
- **standard** — the main pipeline (11 pilots).
- **dnarelax** — Stage-3 with a soft DNA tether (`STAGE3_DNA_RESTRAINT_K=1.5`, released at
  ramp stage 5) so DNA co-relaxes with the protein. Configs exist for **7 pilots** (dux4, egr1, engrailed, ets1, foxa, lef1, tbp); **only
  tbp_dnarelax has a completed eval** (dux4_dnarelax still yields 0 fnat survivors).
- **multi-seed (s1–s5)** — 6 pilots (egr1, engrailed, ets1, foxa, lef1, tbp) re-run with
  5 additional training seeds per arm; the other 5 pilots (csl, dux4, err, nfat, runx) are
  single-seed. tbp_dnarelax has 4–5 seeds per arm.
- **AF3 comparison** — AlphaFold3 predictions (2 seeds × 5 samples) for 6 pilots
  (egr1, engrailed, ets1, foxa, lef1, tbp) under `af3/output/`, used for
  ensemble-diversity-vs-BioEmu analysis.

## Pilots and results

See `results_inventory.csv` for the full machine-readable table. Summary:

| Pilot | Cond. | PDB | Family | Frames | Docked | fnat pass | npz | Seeds | Base r | Aug r | ΔPearson | ΔAUROC | ΔMAE |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| tbp | dnarelax | 1tgh | TBP / β-saddle | 99 |  | 96 | 96 | 5 | 0.634 | 0.600 | -0.034 | -0.018 | 0.043 |
| csl | standard | 3brg | CSL / Rel-like | 100 | 92 | 57 | 57 | 1 | 0.618 | 0.645 | 0.027 | 0.007 | -0.032 |
| dux4 | standard | 5z6z | Homeodomain (double) | 98 | 94 | 0 | 0 |  |  |  |  |  |  |
| egr1 | standard | 1aay | C2H2 zinc finger | 98 | 107 | 89 | 90 | 6 | 0.641 | 0.614 | -0.027 | -0.008 | 0.044 |
| engrailed | standard | 3hdd | Homeodomain | 100 | 97 | 86 | 85 | 6 | 0.622 | 0.618 | -0.004 | -0.008 | 0.021 |
| err | standard | 1lo1 | Nuclear receptor | 100 | 92 | 46 | 41 | 1 | 0.661 | 0.635 | -0.026 | 0.002 | 0.025 |
| ets1 | standard | 1k79 | ETS | 100 | 100 | 100 | 99 | 6 | 0.631 | 0.632 | 0.001 | -0.004 | 0.011 |
| foxa | standard | 1vtn | Forkhead | 100 | 111 | 68 | 73 | 6 | 0.630 | 0.628 | -0.002 | -0.006 | 0.015 |
| lef1 | standard | 2lef | HMG-box | 100 | 89 | 16 | 13 | 6 | 0.627 | 0.619 | -0.009 | -0.007 | 0.001 |
| nfat | standard | 1a66 | Rel/NFAT | 100 | 88 | 26 | 19 | 1 | 0.647 | 0.662 | 0.015 | 0.003 | -0.012 |
| runx | standard | 1hjc | Runt | 100 | 97 | 36 | 35 | 1 | 0.628 | 0.628 | -0.000 | 0.005 | 0.018 |
| tbp | standard | 1tgh | TBP / β-saddle | 99 | 97 | 94 | 94 | 6 | 0.637 | 0.621 | -0.015 | -0.009 | 0.025 |

*Metrics pooled across seeds (mean of per-seed means over the 130-entry benchmark; n=130
entries per model). ΔPearson>0 = augmentation helps. dux4 has 0 fnat survivors → no eval.*

## Key findings (from `analysis/`)
- **fnat gate yield tracks recognition-module rigidity, not interface size.** ETS1 & TBP
  pass ~100%; LEF1 (mobile HMG-box) 18%; dux4 0%. TBP (40 iface residues) and LEF1 (39)
  are near-identical in size but opposite in fidelity.
- **Augmentation net effect is small and mostly negative** on the general-130 benchmark.
  Net-positive pilots: csl (+0.027), nfat (+0.015), ets1 (+0.001). Most negative: egr1
  (−0.027), err (−0.026), tbp (−0.015). tbp_dnarelax is worse than standard tbp (−0.034).
- **Mechanism hypothesis** (`analysis/mechanism_and_roadmap.md`): augmentation helps in the
  *conformational-selection* regime (rigid protein reaches bound pose, little DNA
  deformation: ETS1, EGR1) and hurts under *induced-fit-on-DNA* (TBP kink, LEF1 bend),
  because free-protein frames carry no info about the DNA distortion DeepPBS holds fixed.
- **Caveats flagged by the author:** single crystal PWM label shared across all frames of a
  complex (regularization confound), small n, some single-seed. Motif-level family
  assignment corrects an earlier per-PDB mislabeling of ETS1–RUNX1 co-crystals.

## Where results live on the cluster
- `output/stage7_eval/id_benchmark_<tf>.json` — per-entry eval metrics (the results).
- `output/stage6_train/{baseline,augmented}_<tf>_fold0[_sN][_dnarelax]/` — trained models
  (`Model.best.tar`, `Model_metrics.json`, `predictions/`). 239 subdirs total.
- `output/stage{1..5}[_dnarelax]/` — intermediate ensembles/docked/minimized/npz/aug sets.
- `analysis/` — 6-pilot analysis: figures (F/I/R/S/P series), data CSVs, mechanism memo.
- `analysis/dna_relax/` — pyCurves DNA-bend analysis for the dnarelax condition.
- `rmsd_analysis/` — stage-wise Cα/interface RMSD analysis + plots.
- `af3/output/<tf>_<pdb>/` — AlphaFold3 ensembles + ranking scores.
- `deck/TF_conformation_deck.pptx` — presentation deck.
- `config/pilots/<tf>.sh` — per-pilot definitions; `docs/PIPELINE_FIXES.md`, `README.md`.
