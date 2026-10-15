# v5 gating report — Part 1, item 2, item 4

**Run:** 2026-09-29, from data on disk. All three of the three-hour items complete.

**Bottom line.** Part 1 lands in v5's **third** outcome row, not the first — and the power gain is
0.20 → 0.24, which does not change any decision. Item 4's pooled-arm recommendation **weakens** the
claim rather than strengthening it, for a reason the reliability arithmetic makes explicit. Item 2's
positive control **fails**: no sequence-level flexibility proxy predicts measured bend IQR on the
existing pilots, which is a material problem for Part 3 as designed. And the item-2 work turned up a
**confound nobody has tested — duplex length at ρ = −0.83 against the axis** — which the headline
nevertheless survives.

---

## Part 1 — Entry-vs-seed variance decomposition. **Third outcome row; the gain is negligible.**

Crossed entry × seed random effects on paired entry-level ΔP (1,046 rows), per pilot, relaxed arm:

| pilot | entries | seed var | entry var | residual | current per-seed SD | floor (entries → ∞) | reducible |
|---|---|---|---|---|---|---|---|
| hsf | 3 | 1% | 5% | 94% | 0.0805 | 0.0148 | 82% |
| runx | 5 | 10% | 0% | 89% | 0.0572 | 0.0347 | 39% |
| egr1 | 24 | 4% | 15% | 81% | 0.0466 | 0.0338 | 27% |
| tbp | 6 | 33% | 0% | 67% | 0.1178 | 0.1017 | 14% |
| ets1 | 10 | 18% | 22% | 60% | 0.0631 | 0.0545 | 14% |
| engrailed | 11 | 18% | 17% | 64% | 0.0969 | 0.0843 | 13% |
| csl | 5 | 35% | 16% | 50% | 0.0581 | 0.0512 | 12% |
| foxa | 5 | 52% | 3% | 45% | 0.0811 | 0.0749 | 8% |
| nfat | 15 | 40% | 0% | 60% | 0.0908 | 0.0867 | 5% |
| err | 10 | 73% | 0% | 27% | 0.1420 | 0.1394 | 2% |

**A correction to the framing first.** v5 reasoned from the audit table that "egr1 has the most entries
(24) and the lowest SEM (0.018); hsf has 3 entries at 0.036" — but the relation is not monotone: err has
10 entries and the *highest* SD (0.134), tbp has 6 at 0.118, csl has 5 at 0.058. Entry count does not
order the noise.

**Why not.** The seed-to-seed SD of an entry-mean is √(var_seed + var_resid/n_entries). Only the second
term shrinks with more entries. The decomposition says **var_seed dominates in 7 of 10 pilots** — err
73%, irf 87%, foxa 52%, nfat 40%, tbp 33%, csl 35%. Entry expansion helps only where var_resid is large
relative to var_seed: **hsf (82% reducible), runx (39%), egr1 (27%)** and essentially nobody else.
(Model check: predicted per-seed SD correlates with the empirical one at ρ = +0.955.)

**Consequence for the power numbers.** Taking the median over ten pilots:

| scenario | median SEM | ratio to between-pilot SD | attenuation | interaction power (n = 12) |
|---|---|---|---|---|
| current entry sets | 0.0361 | 0.691 | 0.715 | **0.201** |
| all pilots expanded to 24 entries | 0.0301 | 0.575 | 0.776 | 0.232 |
| entries → ∞ (the floor) | 0.0289 | 0.553 | 0.788 | 0.239 |

**Even infinite eval entries buy 0.20 → 0.24.** v5 called this "the cheapest power gain in the project"
and "more important than any other hour available"; it is neither, because the noise is training
stochasticity, not eval small-sample. **The power ceiling stands and no power calculation below it
changes.** Part 2 can be timestamped without waiting on this.

Worth doing anyway for **hsf** specifically — 3 entries, 82% reducible, and hsf is the pilot with three
unexplained anomalies (§4.3). Sharpening its ΔP is cheap and might inform that. Not a project-level fix.

---

## Item 4 — Pooled-arm correlation. **Do not adopt it. It widens the CI through zero.**

| axis | effect | ρ | p | boot 95% CI | CI excludes 0 | ρ w/o tbp |
|---|---|---|---|---|---|---|
| relaxed axis | relaxed ΔP | -0.745 | 0.0133 | [-1.000, -0.156] | yes | -0.650 |
| relaxed axis | frozen ΔP | -0.527 | 0.1173 | [-0.925, +0.236] | **no** | -0.433 |
| relaxed axis | pooled ΔP | -0.648 | 0.0425 | [-1.000, +0.013] | **no** | -0.550 |
| frozen axis | relaxed ΔP | -0.721 | 0.0186 | [-0.976, -0.139] | yes | -0.633 |
| frozen axis | frozen ΔP | -0.564 | 0.0897 | [-0.925, +0.195] | **no** | -0.400 |
| frozen axis | pooled ΔP | -0.661 | 0.0376 | [-1.000, +0.038] | **no** | -0.533 |
| pooled axis | relaxed ΔP | -0.697 | 0.0251 | [-0.951, -0.065] | yes | -0.600 |
| pooled axis | frozen ΔP | -0.588 | 0.0739 | [-0.961, +0.170] | **no** | -0.433 |
| pooled axis | pooled ΔP | -0.636 | 0.0479 | [-0.987, +0.059] | **no** | -0.500 |

**Pooled axis × pooled ΔP: ρ = −0.636, p = 0.048, CI [−0.987, +0.059] — the CI now contains zero**,
where the relaxed/relaxed cell's does not. v5 estimated "≈ −0.71"; the actual value is −0.636, and
dropping tbp takes it to −0.500.

**Why pooling doesn't help, precisely.** v5's reasoning was: *"If the arms measure the same thing with
equal precision, averaging their ΔP per pilot uses both estimates instead of the one that happened to
correlate better."* The premise is right and the conclusion requires one more condition — that the two
arms' **measurement errors be independent**. They are not:

- observed arm–arm ΔP correlation: **+0.770**
- reliability implied by the measured SEMs: 0.707 (frozen), 0.696 (relaxed)
- arm–arm correlation expected under independent errors: √(0.707 × 0.696) = **0.701**

**The observed +0.770 exceeds the independent-error ceiling of 0.701**, so the errors are positively
correlated — unsurprisingly, since the two arms share own-family entries, share seed indices, and are
built from the same BioEmu frames with only the DNA restraint differing. Pooling correlated errors gives
much less than the √2 reduction, and the between-pilot SD confirms it: 0.0518 frozen, 0.0523 relaxed,
**0.0509 pooled** — a 2% change, i.e. no meaningful averaging benefit.

So pooling forfeits the relaxed cell's tighter CI without buying precision. **Keep relaxed/relaxed as
primary and report frozen/frozen beside it**, exactly as the cleaned-headline table already does.

**The standing anomaly stays open, and is now sharper.** The arms are physically near-identical, equally
precise, and their errors are correlated — and relaxed still correlates at −0.745 against frozen's
−0.564. With differential noise eliminated as an explanation and error-sharing established, there is no
measurement-side story left. The remaining possibilities are a genuine arm difference too small to
detect in the Δr = +0.011 readout test, or chance at n = 10. Both should be stated; neither can be
distinguished here.

---

## Item 2 — Shape-tool flexibility vs measured bend IQR. **The positive control fails.**

**First, a tooling finding.** There is no DNA-shape flexibility predictor installed: DNAshapeR is absent
from the R library, no Python shape package resolves in the `deeppbs` env, and `lib/cg_coefficients/`
contains a single Wigner-3j matrix file (`w3j_matrices-lmax=5-version=0.6.0.gz`) used by the cgDNA-based
`*_dnarelax.sh` restraint setup — not a sequence→flexibility predictor. **Part 3's manipulation has no
validated tool behind it in this stack.**

So I ran the test with established sequence-level proxies instead, extracted from the docked reference
structures (the DNA carried through from each pilot's crystal):

| pilot | PDB | length (bp) | GC frac | flexible YR/CG steps | longest A/T tract | relaxed bend IQR |
|---|---|---|---|---|---|---|
| egr1 | 1aay | 11 | 0.73 | 0.30 | 0 | 6.73 |
| nfat | 1a66 | 12 | 0.42 | 0.18 | 6 | 9.59 |
| hsf | 5d5u | 12 | 0.50 | 0.09 | 2 | 4.83 |
| tbp | 1tgh | 12 | 0.33 | 0.55 | 8 | 8.52 |
| foxa | 1vtn | 13 | 0.46 | 0.17 | 3 | 6.60 |
| err | 1lo1 | 13 | 0.62 | 0.25 | 2 | 5.76 |
| ets1 | 1k79 | 15 | 0.47 | 0.29 | 4 | 3.14 |
| csl | 3brg | 15 | 0.40 | 0.14 | 3 | 4.27 |
| runx | 1hjc | 16 | 0.56 | 0.27 | 2 | 2.36 |
| engrailed | 3hdd | 21 | 0.29 | 0.30 | 6 | 3.91 |

| predictor | vs relaxed IQR | vs frozen IQR |
|---|---|---|
| GC fraction | −0.067 (p = 0.85) | +0.018 (p = 0.96) |
| flexible YR/CG step fraction | +0.061 (p = 0.87) | +0.237 (p = 0.51) |
| longest A/T tract | +0.204 (p = 0.57) | +0.253 (p = 0.48) |
| AA+TT step fraction | −0.109 (p = 0.76) | −0.304 (p = 0.39) |
| **duplex length** | **−0.827 (p = 0.0031)** | **−0.772 (p = 0.0089)** |

**No sequence-composition proxy predicts measured ensemble bend IQR.** v5 called this "a free
retrospective validation of the manipulation" that "de-risks the whole experiment" — it did its job, and
the answer is negative. Three readings, and they are not equally likely:

1. **The proxies are too crude.** Plausible; a real cgDNA/DNAshapeR estimate is not GC fraction. But
   the negative result means Part 3 **cannot be scoped until a real predictor is installed and validated
   the same way.** Question 4 must be re-run with the actual tool before sequences are designed.
2. **n = 10 with 11–21 bp duplexes is too narrow a range** to show a sequence effect.
3. **Measured ensemble bend IQR is not intrinsic DNA flexibility.** It is the fluctuation of a duplex
   that is backbone-restrained toward a crystal and whose protein is varying — closer to a property of
   the restraint geometry than of the sequence. If so, Part 3's manipulation and Part 3's readout are
   measuring different things, which is the serious version of this result.

**Recommendation: keep Part 3 in scope but move question 4 to a blocking prerequisite**, and answer it
with an installed predictor. If a real tool also fails to predict measured IQR on these ten pilots,
reading 3 is the live one and the sequence-design experiment needs rethinking before it is costed.

### The new confound

Duplex length predicts the axis at **ρ = −0.83**, and this has a mechanical explanation: `bend_uu` is a
whole-duplex end-to-end axis angle, so a shorter duplex converts the same absolute end displacement into
a larger angle. Short duplexes look more flexible by construction.

**The headline survives it:**

- unadjusted: ρ = −0.745 (p = 0.013)
- controlling duplex length: **ρ = −0.743 (p = 0.022)**
- length vs relaxed ΔP: +0.451 (p = 0.19); controlling IQR: −0.443 (p = 0.23)
- rank regression ΔP ~ IQR + length: IQR β = −1.181 (**p = 0.022**), length β = −0.536 (p = 0.23), R² = 0.64

So length drives the *axis* strongly but not the *endpoint*, which is the pattern that leaves a partial
intact. Two consequences:

- **Add duplex length to the covariate list and to the test family** (now 20 tests; bend IQR's q moves
  to roughly 0.18–0.19). Claim 4's "q ≈ 0.15–0.17" becomes **q ≈ 0.15–0.19** depending on counting.
- **§4.2's interface restriction becomes more important, not less.** A footprint-restricted IQR is
  computed over a fixed number of base pairs, which removes the length artifact by construction. It is
  now the natural fix for two problems at once.

---

## Consolidated effect on v5

| v5 item | result | change |
|---|---|---|
| Part 1 variance decomposition | third outcome row; seed var dominates 7/10 | **does not gate anything**; power 0.20 → 0.24 at best; do hsf only |
| Item 2 shape-tool validation | **fails**; no predictor installed; no proxy works | Part 3 question 4 becomes blocking; install a real predictor first |
| Item 4 pooled-arm | **CI crosses zero** (−0.636, [−0.987, +0.059]) | do not adopt; keep relaxed/relaxed primary |
| — new: duplex length | ρ = −0.83 vs axis; headline survives at −0.743 | add to covariates and family; strengthens the case for §4.2 |
| Part 2 timestamping | unblocked | Part 1 was the stated gate and it does not bind |
| §4.2 interface restriction | promoted | fixes the length artifact by construction |

**Suggested reorder of Part 5:** §4.2 (interface restriction) moves from item 7 to item 2 — it now
addresses the length confound, is a prerequisite for k = ∞, and can still retract Claim 4. Part 1 drops
to a low-priority hsf-only task. Part 3's question 4 becomes item 3 and blocking.

---

## Reproduce

- Part 1: `perseed_perentry.csv`, `selffam == True`; pair augmented−baseline per (tf, dna, seed, entry);
  `smf.mixedlm("dP ~ 1", g, groups=ones, vc_formula={"entry":"0+e","seed":"0+s"})` per pilot; per-seed SD
  of the entry-mean = √(var_seed + var_resid/n_entries).
- Item 4: mean the two arms' ΔP and the two arms' IQR per pilot; 10,000 bootstrap resamples, seed
  20260929. Independent-error ceiling = √(rel_frozen × rel_relaxed), rel = 1/(1 + (SEM/between-SD)²).
- Item 2: DNA sequences from `output/stage2_docked/<tf>/<pdb>_state_001.pdb`, residues DA/DT/DG/DC,
  longest chain. Length = bp of that chain.
