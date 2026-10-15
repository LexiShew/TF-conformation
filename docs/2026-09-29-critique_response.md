# Response to the critique — steps 1–3 run

**Run:** 2026-09-29, from data already on disk. **Three of the note's premises do not hold; one of its
recommendations should be dropped; and the running order should change.**

Headline: **step 1 passes cleanly, step 2 comes out in the note's own "diagnostic" branch, and step 3's
substitute endpoint does not measure the same thing.** Plus one new data-integrity finding (irf).

---

## Step 1 — mean vs IQR separability and partial correlations. **PASSES.**

The note's stated worry: *"mean-vs-IQR plausibly is at ρ > 0.8."* It is not.

| pair | Spearman | Pearson |
|---|---|---|
| relaxed median bend vs relaxed IQR | +0.345 (p = 0.30) | +0.267 |
| relaxed mean bend vs relaxed IQR | +0.482 (p = 0.13) | +0.298 |
| frozen median bend vs frozen IQR | +0.273 (p = 0.42) | +0.463 |
| crystal (static) bend vs relaxed IQR | +0.173 (p = 0.61) | — |

Excluding irf (see the integrity note below) the median-vs-IQR relation is +0.224. **Bend magnitude and
bend fluctuation are separable in this panel.** The collinearity that would have invalidated the axis
isn't there.

Partial Spearman of bend IQR against relaxed own-family ΔP, controlling one covariate at a time
(rank-residualized; unadjusted ρ = −0.782, p = 0.0045):

| covariate controlled | covariate ~ IQR | partial ρ | p | n |
|---|---|---|---|---|
| relaxed mean bend | +0.482 | -0.802 | 0.0053 | 11 |
| relaxed median bend | +0.345 | -0.790 | 0.0065 | 11 |
| crystal bend (static) | +0.173 | -0.783 | 0.0074 | 11 |
| ensemble size | -0.452 | -0.778 | 0.0080 | 11 |
| protein reachability | +0.273 | -0.778 | 0.0080 | 11 |
| free-state spread | +0.173 | -0.774 | 0.0086 | 11 |
| own-family subset size | +0.184 | -0.773 | 0.0087 | 11 |
| interface MGW fluctuation | +0.321 | -0.777 | 0.0137 | 10 |

**Every partial lands between −0.773 and −0.802 — the axis does not weaken under any single control, and
several *strengthen* it slightly.** The reverse test is the decisive one:

- relaxed **median** bend vs effect, alone: ρ = −0.136 (p = 0.69); controlling IQR: ρ = **+0.229** (p = 0.53)
- relaxed **mean** bend vs effect, alone: ρ = −0.191 (p = 0.57); controlling IQR: ρ = **+0.340** (p = 0.34)

Rank regression, effect ~ IQR + median bend (n = 11, R² = 0.63): IQR β = −0.834 (t = −3.65, **p = 0.0065**);
median bend β = +0.152 (t = 0.66, p = 0.53).

**Conclusion: the axis is fluctuation, not magnitude.** The note was right that this had to be checked
first and right that everything downstream depended on it. It just passes — so nothing is wasted.

One caveat the note would want stated: these are *single*-covariate partials at n = 11, which is exactly
the regime it flags as underpowered. They rule out each covariate acting alone; they cannot rule out a
combination. That limitation is real and survives step 1.

---

## Step 2 — within-pilot mixed model. **RUNS, and lands in the note's diagnostic branch.**

Built as specified: baseline checkpoints only (so no memorization contamination), per-state
|bend − crystal bend| as the fixed effect, pilot as random effect. Joined
`pycurves_all_perstructure.csv` to `per_state_all_conditions.csv` on (tf, arm, state) — **every state with
a PWM prediction has curves geometry, 701 unique states, 3,493 state×seed rows across 12 pilots.**

| model | β (per SD of bend deviation) | se | z | p |
|---|---|---|---|---|
| random intercept | -0.0108 | 0.0043 | -2.49 | 0.0128 |
| random intercept + slope | -0.0287 | 0.0152 | -1.89 | 0.0584 |
| random int (vs experiment) | -0.0076 | 0.0040 | -1.90 | 0.057 |
| random int+slope (vs experiment) | -0.0208 | 0.0125 | -1.66 | 0.0973 |

Random intercept alone gives β = −0.0108 (p = 0.013) — the predicted sign. **But adding the random slope
the note specifically asked for drops it to p = 0.058**, and that is the correct specification: with
random intercepts only, the 701 states are treated as more independent than they are.

The per-pilot slopes are the real answer:

| pilot | states | slope (Δr per °) | p | Spearman ρ |
|---|---|---|---|---|
| runx | 35 | -0.0344 | 0.3696 | -0.119 |
| engrailed | 85 | -0.0176 | 0.0002 | -0.193 |
| err | 41 | -0.0129 | 0.0035 | -0.354 |
| ets1 | 99 | -0.0100 | 0.1298 | -0.095 |
| hsf | 96 | -0.0085 | 0.3733 | -0.066 |
| tbp | 94 | -0.0085 | 0.0155 | -0.076 |
| egr1 | 90 | -0.0070 | 0.0048 | -0.166 |
| lef1 | 12 | +0.0026 | 0.4823 | -0.154 |
| nfat | 17 | +0.0034 | 0.0647 | +0.061 |
| foxa | 70 | +0.0036 | 0.5140 | +0.086 |
| csl | 57 | +0.0080 | 0.4759 | +0.058 |

**7 of 11 negative, sign test p = 0.55.** Effect sizes are tiny (≈ −0.01 r per degree). csl, foxa, nfat
and lef1 run the wrong way. Second level: per-pilot slope vs bend IQR gives ρ = +0.30 (frozen) / +0.59
(relaxed, p = 0.074) — the **wrong sign** for a mediation story, where a higher-IQR pilot should show a
more negative within-pilot slope.

The note pre-registered how to read this: *"If the within-pilot slope is flat while the across-pilot
correlation is −0.78, that's diagnostic on its own — it would mean bend fluctuation isn't acting through
individual states at all, and the across-pilot result is a pilot-level property (i.e. probably a
confound)."*

That is approximately what happened, and it deserves to be taken at face value. Two readings remain open:

1. **The pessimistic one the note names.** Bend fluctuation is a pilot-level property; the across-pilot
   −0.78 may be confounded by something not in step 1's covariate list.
2. **A structural alternative.** The endpoint is a *pilot-level* quantity by construction — own-family
   ΔPearson is one number per pilot per arm, obtained by training on the whole ensemble. There is no
   per-state ΔP, so step 2 is not the same question with more n; it substitutes per-state PWM *agreement*
   for per-state augmentation *effect*. A mechanism that operates through the aggregate label
   ("one crystal PWM mislabels a fluctuating ensemble") predicts exactly this: an ensemble-level effect
   with no per-state signature.

Reading 2 is not a rescue — it is a testable claim, and the honest statement is that step 2 **did not
confirm the mechanism and moved the burden onto the intervention.** Which is the note's own conclusion,
reached one step earlier than it expected.

---

## Step 3 — swap the endpoint to PWM dispersion. **Do not do this as specified.**

The note proposes per-state PWM dispersion across the ensemble as an inference-only, leakage-free
substitute, on the grounds it "measures the same underlying thing." It does not, on two counts.

**It does not reproduce the mechanism correlation.** Bend IQR vs PWM dispersion (baseline checkpoints,
per-pilot, dispersion = IQR of per-state r):

| axis | vs dispersion (SD) | vs dispersion (IQR) |
|---|---|---|
| frozen bend IQR | −0.018 (p = 0.96) | −0.273 (p = 0.42) |
| relaxed bend IQR | +0.082 (p = 0.81) | −0.100 (p = 0.77) |

Restricting to pilots with ≥30 states doesn't help (−0.03 to −0.38, all n.s.). And dispersion only weakly
tracks the endpoint it would replace (ρ = +0.28 to +0.42, n.s. — and the *wrong sign*: more dispersion
goes with more positive ΔP).

**It is substantially a restatement of mean agreement.** Dispersion vs mean per-state r: ρ = **−0.776**
(p = 0.003). A pilot whose states agree poorly with the crystal also disagrees among themselves — most of
what dispersion measures is already in `r_mean`, which is a quantity the project has.

So the proposed swap would replace a noisy-but-hypothesis-linked endpoint with a cheap one that is
uncorrelated with the axis and largely redundant with an existing metric. **Keep ΔPearson as the
endpoint.** The memorization problem is real but it is fixed by using baseline checkpoints (or cross-pilot
ones), not by changing what is measured.

The note's underlying instinct — get an inference-only endpoint so added pilots don't need training arms —
is right and worth pursuing. It just needs an endpoint validated against ΔP first, and this one isn't.

---

## New: a data-integrity finding on irf

`irf`'s relaxed-arm `bend_uu` is **156–179° across all 93 states** (median 169°) while its crystal bend is
**16.9°**. `bend_uu` is an unsigned axis angle, so values near 180° are a register/parse failure, not
extreme bending. `FINDINGS.md` already drops irf from the MGW-FL tables "for a groove-parse artifact" —
the same pilot has the same class of problem in the bend channel, and nobody has flagged it there.

Consequences:
- irf's bend **median/mean are unusable**. Its **IQR (5.5°) is unaffected** because quartiles are robust to
  the offset — which is why the headline axis is not corrupted.
- Any analysis using `cry_bend_uu`-vs-ensemble-bend deviation for irf is invalid. That includes irf's row
  in step 2 (excluded there anyway, 5 states).
- Several other pilots have single-state maxima near 145–170° (err 165, nfat 170, runx 146, foxa 144,
  dux4 168) — isolated frames with the same failure. The IQR-based axis is robust to them; a mean- or
  SD-based axis would not be. **This is an argument for keeping IQR rather than switching to SD.**

Recommend: add a range assertion to the pyCurves aggregation (flag any state whose `bend_uu` deviates
from the pilot median by more than ~90°) and re-inspect irf's register before using its bend anywhere.

---

## On "would more pilots help?"

The note's table asserts n ≈ 30 is where covariate control "becomes answerable." Checking that with a
Fisher-z power calculation for a partial correlation with one covariate:

| n | power at true ρ = 0.78 | power at true ρ = 0.45 |
|---|---|---|
| 11 | 0.79 | 0.25 |
| 20 | 0.99 | 0.49 |
| 30 | 1.00 | 0.70 |
| 40 | 1.00 | 0.83 |
| 50 | 1.00 | 0.91 |

**If the true effect is as large as observed, n = 11 already has 79% power and n = 20 has 99%** — so more
pilots are not needed to detect *this* effect. n ≈ 30 is the right target only if you expect the true
effect to be substantially more modest (ρ ≈ 0.45), which after step 2 is the prudent assumption. So the
note's number is defensible, but for the opposite reason to the one given: not to gain power on a large
effect, but to stay powered if the effect is smaller than the n = 11 estimate (which regression to the
mean says it probably is).

The note's cost argument is sound and independent of this: the candidate pool is the existing 130-entry
benchmark, and inference-only pilots avoid training arms. But per step 3, the inference-only endpoint
still needs to be found.

---

## Revised order of work

The note's steps 1–3 were correctly identified as the cheap, high-information ones. Having run them:

1. **Step 1 — done.** Axis is fluctuation, not magnitude, and survives every single-covariate control.
2. **Step 2 — done, negative-to-ambiguous.** Report it. It is the most honest thing in the analysis and
   it constrains the claim: the mechanism, if real, is ensemble-level, not per-state.
3. **Step 3 — drop as specified.** Keep ΔP; fix leakage by checkpoint choice. If an inference-only
   endpoint is wanted, first find one that correlates with ΔP at n = 11.
4. **irf register fix + range assertion.** Hours, and it is a correctness issue.
5. **The k sweep is now the priority, not step 4.** Step 2 removed the within-pilot support, which makes
   the intervention the only remaining way to move this from correlation to cause. The note's motivation
   for it is strengthened, not weakened: three independent measurements say k = 1.5 barely moves the axis
   (IQR ρ = +0.927 frozen-vs-relaxed, mean 5.46° → 5.79°), so k ∈ {0.5, 1.5, 5.0, ∞} is a genuine
   dose–response and the high-IQR pilots should move most.
6. **Then scale pilots**, with the target set by the ρ ≈ 0.45 row (n ≈ 40) rather than ρ ≈ 0.78.

The sequence-design variant remains the strongest test for the reason the note gives — protein fixed, DNA
flexibility varied — and it is the one experiment that would survive the confounding ceiling it correctly
identifies.

## Reproduce

Inputs (both already on disk, no new compute):
`analysis/analyses/dna_relax/data/pycurves_all_perstructure.csv`,
`analysis/analyses/conf_vs_crystal_pwm/per_state_all_conditions.csv`,
`analysis/analyses/mechanism/data/selffam_effects.csv`.
Join key (tf, arm, state) with `state` cast to int; `arm` maps condition baseline/aug_frozen → frozen,
aug_relax → relaxed. Mixed models: `statsmodels.formula.api.mixedlm`, REML, groups = tf,
`re_formula="~z_bend_dev"` for the random-slope fits.
