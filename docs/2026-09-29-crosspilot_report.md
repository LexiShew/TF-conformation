# Cross-pilot test: the +0.28 crystal gain does not survive

A reviewer argued that the held-out crystal result in `all_conditions_report.md`
is near-duplicate leakage rather than generalization: the augmented model trained
on ~90 perturbed copies of the same complex carrying the crystal's own PWM as
label, then was tested on that complex. **The critique is correct on its central
claim, and the decisive test now confirms it.**

## The test

For target pilot X and donor pilot Y (Y ≠ X), score X's conformations and X's
crystal under Y's baseline and Y's augmented checkpoints. Y's augmented arm
trained on Y's conformations only, so X's complex is unseen at the *structure*
level — not merely absent from an entry list. The paired quantity

    Δ(X, Y) = r_augmented_Y(X) − r_baseline_Y(X)

holds the scored structures fixed and varies only whether the model saw that
complex. All 12 × 11 ordered pairs × 2 arms × 5 seeds = 264 inference blocks,
1,416 scored pairs, no failures.

## Result: the crystal gain vanishes

| test structure's relation to training | Δ*r* conformations→crystal | Δ*r* crystal→experiment |
|---|---|---|
| own crystal (~90 near-copies in training) | **+0.246** | **+0.165** |
| unseen complex (cross-pilot donor) | **+0.025** | **−0.008** |

Mann-Whitney within vs cross: p = 2.2 × 10⁻¹⁷ (state) and 1.2 × 10⁻⁸ (crystal).
The crystal-readout gain is **indistinguishable from zero** on unseen complexes
(Wilcoxon vs 0, p = 0.24; 4 of 12 targets positive), and it collapses for 10 of
the 12 pilots individually. The remaining two, CSL and EGR1, had a slightly
*negative* within-pilot delta to begin with.

**So the +0.28 I reported as "clean" was not clean.** The reviewer's inversion of
my recommendation was right: the state-level number was inflated *and flagged*,
while the crystal number was inflated and read as trustworthy. That was the more
dangerous of the two, and it is now retracted.

A small residual survives on the state→crystal readout: +0.025 cross-pilot,
Wilcoxon vs 0 p = 3.4 × 10⁻⁵, positive for 10 of 12 targets. It is an order of
magnitude below the within-pilot +0.246 and is consistent in size with the
+0.03–0.08 own-family effects measured elsewhere in the project. That is the
defensible effect size.

## Where the critique's mechanism was wrong

The critique attributed the leakage to the state features being *near-identical*
to the crystal's — "near-identical by construction" in the DNA-shape channels.
Measured in the model's own input space, they are not. Comparing each state's
distance to its own crystal against the distance of unrelated training entries to
that same crystal:

- median ratio **0.81** — states sit only ~19% closer than an arbitrary
  unrelated training structure
- **ERR (1.07) and NFAT (1.34) are *further* from their own crystal** than
  unrelated entries are
- TBP is the tightest at 0.48, still nowhere near duplicate

So the effect is real but the proposed mechanism is not the geometry of the input
features. Whatever the augmented model memorizes, it is not because the inputs
are near-copies — the label is shared across ~90 training entries while the
inputs genuinely differ. Label-side memorization, not input-side duplication.

## What the critique got right beyond the main claim

**The variance result holds, and it is the cleaner version of the claim.** SD
across conformations drops from 0.235 (baseline) to 0.180 (frozen) and 0.172
(relaxed) — reduced in **11 of 12 pilots**, Wilcoxon p = 0.027 and 0.021. HSF1 is
the sole exception (0.386 → 0.519), which is why its mean barely moves.
Augmentation buys consistency across conformers even where it does not buy
accuracy.

**TBP's diagnosis is right in direction, weaker in evidence.** TBP baseline:
crystal→experiment *r* = 0.13 while state→crystal = 0.81 — the model predicts
much the same wrong motif everywhere and the states agree with it. Under
augmentation the crystal readout rises to 0.64 while state agreement falls to
0.76. State-vs-crystal agreement is only interpretable conditional on the crystal
prediction being good. But across the 12 pilots that conditional is not
statistically visible (Spearman ρ = 0.38, p = 0.22), and the three pilots with an
unreliable anchor average *r* = 0.51 against 0.50 for the rest — so it is a
correct caution about individual pilots, not a systematic confound.

**Follow-up #1 as I wrote it was wrong.** "Train on states 1–50, score 51–100"
holds out no information; the exchangeable unit is the complex, not the state.
The cross-pilot design replaces it and needed no retraining.

## What is left standing

| claim | status |
|---|---|
| crystal readout improves (+0.28 held-out) | **retracted** — −0.008 cross-pilot |
| state-level agreement improves | inflated; +0.025 is the clean estimate |
| conformational variance shrinks | **holds** — 11/12 pilots, p = 0.027 |
| motifs sharpen, not flatten | holds as a model property, unchanged |
| frozen vs relaxed is a coin flip | unchanged, but measured only within in-sample arms |

![Leakage gradient]({{artifact:art_b0921d2b-7338-4204-b478-fa30eaf029d8}})

**(a)** Augmentation Δ*r* against structural distance from the training frames,
both readouts. **(b)** Per-pilot collapse of the crystal-readout gain.
**(c)** Distance from each pilot's states to its own crystal, relative to
unrelated training entries — values near 1 mean the states are not near-copies.

## Reproduce

On `endeavour`, in `/project2/rohs_102/shewchuk/TF-conformation`:

- `stage7_eval/run_crosspilot.sbatch` — the 264-block cross-pilot sweep,
  writing `predictions_crosspilot/` in each donor run dir
- `analysis/analyses/conf_vs_crystal_pwm/compare_crosspilot.py` — paired scoring
- `analysis/analyses/conf_vs_crystal_pwm/check_nearduplicate.py` — the
  distance-ratio measurement
- `analysis/analyses/conf_vs_crystal_pwm/check_feature_similarity.py` — per-key
  feature deviation
