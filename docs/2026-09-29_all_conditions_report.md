# All three DeepPBS conditions: conformation PWMs vs the crystal's PWM

Extends the baseline-only comparison to both augmented arms. Every fnat-gated
conformation is now scored under three conditions, each **self-consistently** —
a conformation's predicted PWM is compared to the crystal PWM produced by the
*same checkpoint*, so a difference reflects the model, not a shifted reference.

| Condition | Run dirs | Features | States scored |
|---|---|---|---|
| no augmentation | `baseline_<tf>_fold0_s1..s5` | `data/assembly2024` (train), `stage5_aug` (inputs) | 3,493 |
| augmented, frozen DNA | `augmented_<tf>_fold0_s1..s5` | `stage5_aug` | 3,488 |
| augmented, relaxed DNA | `augmented_<tf>_fold0_dnarelax_s1..s5` | `stage5_aug_dnarelax` | 3,345 |

10,682 comparisons total across 12 pilots; 69 entries skipped (67 with <4
overlapping PWM columns, 2 missing a crystal prediction). The relaxed arm is
scored against its own relaxed-geometry features — same entry names, different
directory — since that is what it was trained on.

## The leakage caveat comes first

**The augmented models were trained on these exact conformations.** Every state
in the frozen arm's comparison appears in `train0_aug_<tf>.txt`, and likewise for
the relaxed arm — 100% in-sample, verified per entry rather than assumed
(`train_membership.json`). Their state-level agreement with the crystal is
therefore *not* a generalization measurement, and the large gap in panel (a) is
partly memorization by construction.

Crystal entries divide into three classes (`crystal_leakage.json`):

| Class | Pilots | Meaning |
|---|---|---|
| held out | ERR, Ets-1, HSF1, RUNX1, TBP | crystal entry in no training file |
| exact | CSL, EGR1, Engrailed, IRF1, LEF1, NFAT | crystal entry itself in training |
| structure only | FoxA | same complex under a *different* PWM label (`1vtn_C_FOXA3_HUMAN` in training vs `1vtn_C_MA1683.1` scored) |

This corrects the earlier report, which flagged FoxA as in-training on a
PDB-level match; at entry level its scored label was never trained on.

## Results

**Both augmented arms track their own crystal more closely for 11 of the 12
pilots.** Mean state-level *r* rises from **0.50** (no augmentation) to **0.75**
(frozen) and **0.73** (relaxed) — +0.25 and +0.23 paired per TF × seed (Wilcoxon
p = 4.7 × 10⁻¹⁰ and 9.7 × 10⁻¹⁰). Restricting to the 5 held-out pilots the gain
persists at +0.23 and +0.19 (p = 3.2 × 10⁻⁵, 1.2 × 10⁻⁴), so it is not purely an
artifact of the crystal being in training — though the states themselves still are.

**TBP is the one exception, and it is a held-out pilot.** Both augmented arms track
TBP's crystal *less* closely than the baseline does (0.76 frozen and 0.75 relaxed
vs 0.81 baseline, Δ = −0.051 and −0.066) — the only pilot where augmentation hurts,
and one of the 5 whose crystal is in no training set. TBP also has by far the
lowest-information crystal motif of the set (5.3 bits averaged over conditions, vs
7.8–12.8 for the others) and the weakest baseline crystal-vs-experiment agreement
of any pilot (*r* = 0.13, vs 0.32–0.92). Notably augmentation does help TBP's
*crystal* readout substantially (*r* = 0.13 → 0.64 frozen, 0.53 relaxed) even as it
hurts the state-level agreement — so the exception is specific to how augmented
models read TBP's conformations, not to TBP overall. Held-out state-level gain is
therefore 4 of 5 pilots, not 5 of 5.

**The crystal readout improves, and this is the cleaner evidence.** The crystal
entry is a *separate* structure from the conformations, so for the 5 held-out
pilots nothing about that comparison is in any training set. There, agreement
with the experimental PWM goes from **0.55** to **0.83** (frozen, +0.28,
p = 1.5 × 10⁻⁶) and **0.76** (relaxed, +0.21, p = 1.5 × 10⁻⁵). Augmentation
improves the model's reading of a structure it never saw.

**Augmented motifs are sharper, not flatter.** Mean conformational information
content rises from 7.50 to 8.67 bits (frozen) and 8.80 (relaxed). The baseline's
1.14-bit deficit against its crystal largely closes — the flattening reported
earlier is a property of the unaugmented model, not of the conformations.

**Frozen vs relaxed DNA is a coin flip.** Mean Δ*r* = +0.011 favouring frozen,
Wilcoxon p = 0.36 over 59 TF × seed pairs. Per pilot the two straddle the
diagonal: ERR favours frozen (+0.17), LEF1 and FoxA favour relaxed (−0.07,
−0.03). Whatever the DNA-relaxation protocol changes, it does not systematically
change how closely conformations reproduce the crystal motif.

**HSF1 stays the outlier.** Its states reach only *r* ≈ 0.61 under either
augmented arm despite 96 surviving conformations — the highest count of any
pilot with the lowest agreement. Combined with its unchanged information content
in the baseline, HSF1 looks like a genuine structural disagreement rather than
thin sampling.

![Three-condition comparison]({{artifact:art_3a01b120-c5d1-4d84-b1fe-ea95d8ce1eb5}})

**(a)** Per-pilot state-level agreement, mean ± 95% CI over 5 training seeds;
`*` marks pilots whose crystal entry was seen in training, and TBP — the only
pilot where baseline beats both augmented arms — is called out. **(b)** Paired
Δ*r* against no augmentation, one point per pilot, filled = crystal held out;
TBP is the sub-zero point in both columns.
**(c)** Crystal-prediction agreement with the experimental PWM by condition.
**(d)** Conformational information content. **(e)** Frozen vs relaxed, per pilot.

## What to take from this

The honest reading: augmentation clearly improves how DeepPBS reads a *crystal*
structure it never trained on (+0.28 on held-out pilots), and that is the result
worth quoting. The state-level agreement jump is inflated by the states being
in-sample and should not be quoted as generalization without a
leave-conformations-out retrain.

Two follow-ups that would sharpen this:

1. **Retrain with a held-out conformation split** — e.g. train on states 1–50,
   score 51–100. Without it, the state-level numbers here cannot separate
   memorization from improved structure reading. TBP is the natural control: if
   augmentation genuinely hurts there rather than reflecting its unusually flat
   crystal motif, a conformation-split retrain should reproduce the negative Δ.
2. **Check whether the frozen/relaxed tie survives per-family stratification.**
   The aggregate is flat, but ERR's +0.17 and LEF1's −0.07 suggest the protocol
   may matter for specific binding modes rather than uniformly.

## Reproduce

On `endeavour`, in `/project2/rohs_102/shewchuk/TF-conformation`:

- `stage7_eval/run_augcmp_pred.sbatch` — both augmented arms, 12 pilots each
- `stage7_eval/run_hsf_aug_cry.sbatch` — HSF1 crystal for the frozen arm (its
  entry list predated the 5D5U featurization)
- `stage7_eval/confcmp_idfiles/id_<tf>_confcmp_relax.txt` — relaxed-arm entry lists
- `analysis/analyses/conf_vs_crystal_pwm/compare_all_conditions.py` — scoring
- `analysis/analyses/conf_vs_crystal_pwm/check_leakage.py` — the leakage classification
