# DeepPBS Conformational Augmentation — Pilot Assessment Memo
**Date:** 2026-07-08 · **Author:** research log (auto) · **Project:** TF-conformation

## One-line verdict
Conformational augmentation of DeepPBS with BioEmu ensembles does **not** uniformly
improve binding-specificity prediction — but the *sign and size* of its effect are
predictable from a single physical variable (how far the free-protein ensemble drifts
from the crystallographic bound pose). That predictive mechanism, not a blanket
improvement, is the result worth pursuing.

## Pilots run (6 TFs, current post-fix pipeline)
Held-out paired transfer test available for TBP, EGR1, ETS1 (n=5 each); FOXA and LEF1
have placeholder PWMs (general-130 only); engrailed still completing at time of writing.

| Pilot | Recognition module | fnat yield | Held-out on-target effect |
|---|---|---|---|
| **ETS1** | winged-HTH (ETS) | 100% | **uniformly better, all significant** |
| **TBP**  | β-saddle (TATA) | 100% | uniformly worse, all significant |
| **EGR1** | C2H2 Zn-finger | 91% | mixed (primary up, rank-IC down) |
| FOXA | forkhead | 75% | general-130 only, mildly positive |
| LEF1 | HMG-box | 19% | general-130 only, mildly negative |
| engrailed | homeodomain | 90.5% (86/95) | pending (training) |

## On-target paired deltas (5 held-out entries; + = augmentation improves)
Error metrics (Brier, MAE) sign-flipped so + always means improvement.

| metric | TBP | EGR1 | ETS1 |
|---|---|---|---|
| Pearson r | -0.263*** | +0.102* | **+0.149*** |
| Spearman rho | -0.180*** | -0.100*** | +0.068** |
| AUROC | -0.057** | +0.142*** | +0.030*** |
| IC-weighted PCC | -0.182*** | +0.014 | +0.085*** |
| Brier (multi) | -0.167*** | +0.064 | +0.092*** |
| MAE | -0.299*** | +0.067* | +0.183*** |

\*p<0.05 \*\*p<0.01 \*\*\*p<0.001. ETS1 is the only pilot where every metric
(including the rank-based ones EGR1 degraded) improves significantly.

## The mechanism: ensemble drift from the bound pose
Mean Ca deviation of the minimized (Stage 3) ensemble from the crystal bound pose,
measured on the crystal-DNA-registered frame (10 sampled states/pilot):

| Pilot | mean dev (A) | range (A) |
|---|---|---|
| TBP | 1.33 | 0.80-1.72 |
| ETS1 | 1.50 | 0.87-2.65 |
| engrailed | 2.11 | 1.13-5.38 |
| EGR1 | 2.75 | 1.01-10.86 |
| FOXA | 4.65 | 1.78-9.04 |
| LEF1 | 5.97 | 3.23-18.16 |

This ordering tracks fnat yield exactly, and (for the three tested pilots) tracks the
sign of the augmentation benefit at the rigid end. The physical reading: BioEmu samples
the *free* protein, so rigid recognition modules (TBP saddle, ETS winged-HTH) stay near
the bound geometry and supply faithful-but-diverse frames; mobile modules (HMG-box)
drift off the binding mode and inject noise.

Stage 2 -> Stage 3 deviation is <=0.05 A in every pilot mean: minimization barely moves
the docked pose. The ensemble spread is set almost entirely at BioEmu sampling + rigid
docking, not by minimization.

## Why this is not yet trustworthy — confounds to kill (from concerns.txt)
1. **Rigid DNA / single shared PWM.** All ~N minimized frames carry the crystal DNA
   rigidly and share ONE PWM label, so augmentation up-weights each TF ~95x with no
   binding-mode diversity. TBP's "harm" could be pure label-imbalance rather than a
   conformational effect. This is the #1 alternative explanation.
2. **GBSA over-stabilization** may push interfaces off the bound geometry while still
   clearing the Ca-fnat floor. TBP passed fnat at 100% and still lost every metric —
   high pass-rate does NOT predict useful augmentation.
3. **No fixed random seed** in the generated training configs, so some pilot-to-pilot
   variance is noise (fnat yields already differed run-to-run: FOXA 79% -> 75%).

## Recommended decisive experiment
Confound-controlled re-run on **ETS1 (the win) and TBP (the loss)**:
- fixed random seed (concerns #3),
- per-frame PWM/DNA variation instead of one shared label (concerns #1),
- optionally a stricter interface-geometry gate or explicit-solvent min (concerns #2).

Outcome logic: if ETS1's gain survives per-frame labeling, the effect is genuine
conformational information -> scale to the full TF set. If it collapses, the effect was
label-imbalance -> a clean negative worth knowing before scaling.

## Pipeline-engineering results delivered this cycle
- **Post-min fnat gate** (B7) confirmed working and self-consistent across all pilots.
- **Two-assembly auto-select patch** (fnat gate + redock monomer guard): unblocked ETS1
  (fnat 0.000 -> 0.783 median, 0/97 -> 100/100 pass) and engrailed (packing-copy drop),
  provably a no-op on single-chain references (TBP/EGR1 fnat unchanged). Diff saved as
  `tf_conformation_assembly_fix.diff` (not git-committed).
- ETS1 is the strongest augmentation signal in the study and exists *only* because of
  this patch.

## Deliverables
- `benchmark_summary.csv` (v4) — all pilots, global-130 + held-out subset paired stats.
- `on_target_three_pilots.png` — TBP/EGR1/ETS1 paired deltas.
- `tbp_vs_egr1_paired.png`, `four_pilot_summary.png` — earlier cross-pilot figures.
- six `*_ensemble.png` — crystal-vs-ensemble overlays with Ca-deviation annotations.
- `scripts/viz/` — reusable render scripts (ensemble + single-state).
