# Baseline DeepPBS: bound conformations vs the crystal's PWM prediction

All fnat-gated, energy-minimized BioEmu bound conformations were run through the
**baseline (no-augmentation) DeepPBS arm** — the same `baseline_<tf>_fold0_s1..s5`
checkpoints the project's A/B comparison uses — and each predicted PWM was scored
against the prediction that the *same checkpoint* makes for the crystal complex.

## What was run

| | |
|---|---|
| Pilots | 12 (dux4 excluded — 0 states survive the fnat gate) |
| Conformations scored | **701** distinct minimized states |
| Checkpoints | baseline fold-0 seeds s1–s5 (4 seeds for LEF1; s2 absent) |
| Predictions | 3,611 PWM comparisons (states + crystals + ensemble means) |
| Crystal reference | per-pilot co-crystal entry, same PWM label, same checkpoint |
| HSF1 | 5D5U is absent from the upstream DeepPBS dataset, so the symmetry-completed duplex was featurized fresh (230 contacts) and predicted through the same checkpoints |
| Skipped | 30 of 3,641 (29 with <4 overlapping PWM columns after alignment; 1 LEF1 seed with no crystal prediction) |

Predictions are compared **in PWM coordinates**. DeepPBS returns a per-DNA-position
prediction plus a mask marking which positions aligned to the reference motif; the
masked prediction positions correspond in order to the masked columns of the
experimental PWM. Lifting every entry into that shared frame means conformations
whose alignment window differs in length from the crystal's are still comparable,
on the intersection of covered columns. Strand orientation is resolved against the
crystal before scoring.

## Headline result

**A single bound conformation usually does not reproduce the crystal's predicted
motif.** Across all 701 states, mean Pearson *r* to the crystal PWM is
**0.55** (median 0.60). Only **25%** of conformations reach *r* ≥ 0.8;
**14%** fall at or below *r* ≤ 0.2, and **7%** are anticorrelated with the
crystal prediction. The spread is wide and TF-dependent: Engrailed (0.84 ± 0.18)
and TBP (0.81 ± 0.23) are largely crystal-like, while HSF1 (0.27 ± 0.39),
NFAT (0.33), LEF1 (0.35) and ERR (0.38) are not.

**Averaging the ensemble recovers most of the crystal's motif.** The
ensemble-mean PWM reaches *r* = **0.73** against the crystal on average — a
**+0.20** gain over the median single conformation, winning in **58 of 59**
TF × seed pairs (Wilcoxon p = 2.5 × 10⁻¹¹). Pilots with more surviving states
average closer to the crystal (Spearman ρ = 0.60, p = 0.039), so the residual
gap is partly a sampling-depth effect.

**Conformational PWMs are flatter for 10 of the 12 pilots.** Mean information
content is **7.47 bits** per conformation vs **8.62 bits** for the crystal — a
1.14-bit deficit (Wilcoxon p = 1.5 × 10⁻⁶). The loss is largest for EGR1
(−2.69 bits), Ets-1 (−2.53) and NFAT (−2.37). Two pilots run the other way:
FoxA conformations are sharper than its crystal (+0.87 bits) and HSF1 is
essentially tied (+0.03), so both sit on or above the diagonal in panel (b). Flattening is what you expect if the
predicted specificity is diluted by interface geometry that has drifted from the
crystallographic contact set — the conformations are being read as less
discriminating, not as recognizing something different.

**The crystal still tracks the experimental motif best (11 of 12 pilots).** Against the
JASPAR/HOCOMOCO PWM, the crystal prediction averages *r* = 0.67, the
ensemble-mean 0.56, and a single conformation 0.40 (crystal vs single
conformation, Wilcoxon p = 3.4 × 10⁻¹⁰). Ensemble averaging closes most but not
all of that gap (mean deficit −0.10, p = 1.9 × 10⁻⁵); it is TF-specific — near
zero for Engrailed, LEF1, ERR, TBP and EGR1, but −0.20 to −0.31 for NFAT, Ets-1,
HSF1 and CSL.

**Conformation matters more than checkpoint seed.** For the same state, the SD
across the five seeds averages 0.16; for the same seed, the SD across
conformations averages 0.25. The variability is a property of the structures, not
of training noise.

Whether the pilot's crystal sits in the fold-0 training set does not explain the
ranking: in-training pilots average *r* = 0.52 vs 0.49 for held-out ones
(Mann-Whitney p = 0.76). Note however that the crystal's *agreement with the
experimental PWM* is optimistic for the 7 in-training pilots (CSL, EGR1,
Engrailed, FoxA, IRF1, LEF1, NFAT), so the crystal-vs-experiment column should be
read as an upper bound for those.

## Figures

![Conformation vs crystal PWM agreement]({{artifact:art_29f83f10-c89b-4bcf-8203-79a3faa9dc33}})

**(a)** Per-conformation *r* to the crystal PWM, every state × seed; bar = median
conformation, diamond = ensemble-averaged PWM. **(b)** Information content,
conformation mean vs crystal (points below the dashed diagonal lost information
relative to the crystal; FoxA and HSF1 did not); open markers are pilots whose crystal is in the
fold-0 training set. **(c)** Agreement with the experimental PWM for crystal,
ensemble average, and mean single conformation; the crystal ranks highest for every
pilot except TBP, where the ensemble average edges it by 0.005.

![Exemplar motif logos]({{artifact:art_a7281747-1df1-452d-9c43-351d566e0e5a}})

Experimental PWM, crystal prediction, and the best/worst single conformation
(seed 1) for four pilots. EGR1's worst state (*r* = 0.18) loses the poly-G
signature entirely; TBP's best state (*r* = 0.99) is indistinguishable from the
crystal's TATA readout.

## What this implies for the augmentation

Ensemble averaging behaving as a denoiser — +0.20 toward the crystal, 58/59 wins
— supports treating the conformer set as repeated noisy observations of one
binding mode rather than as distinct specificity states. But the per-state
scatter is not uniform noise: the 14% of states at *r* ≤ 0.2 and 7% anticorrelated
are contributing motifs that disagree with the crystal, and they enter training
unweighted. Two things worth testing next:

1. **Weight or filter states by predicted-PWM agreement** (or by fnat, if the two
   correlate) instead of the current flat fnat floor. The per-state CSV has the
   agreement scores needed to check whether the low-*r* tail is the same set of
   structures the fnat gate already ranks worst.
2. **Look at HSF1, NFAT, ERR and LEF1 specifically.** These four have the lowest
   crystal agreement. LEF1 (12 states) and NFAT (17) also have among the thinnest
   surviving ensembles, so low agreement there may be sampling depth rather than
   genuine conformational disagreement — but HSF1 (96 states) and ERR (41) are
   well-sampled, so their disagreement needs another explanation.

## Reproduce

On `endeavour`, in `/project2/rohs_102/shewchuk/TF-conformation`:

- `stage7_eval/confcmp_idfiles/id_<tf>_confcmp.txt` — entry lists (states + crystal)
- `stage7_eval/run_confcmp_pred.sbatch` — baseline inference, writes `predictions_confcmp/` in each run dir
- `stage7_eval/run_hsf_crystal_feat.sbatch`, `run_hsf_cry_pred.sbatch` — the 5D5U crystal
- `analysis/analyses/conf_vs_crystal_pwm/compare_pwms.py` — scoring; also writes per-TF/per-seed `pwms_<tf>_s<n>.npz` holding the aligned PWM stacks
