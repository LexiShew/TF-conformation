# Analysis plan v5 — 2026-09-29

**Supersedes** v4. Written after the v4 audit report, which closed §1.1–§1.3 and established a
structural limit that reshapes the rest of the project.

**The finding that organizes this document.** At pilot level, with the measured effect size
(ρ ≈ −0.745), the measured noise (seed SEM ≈ 0.66 × between-pilot SD), and n = 12 as a hard ceiling,
**no affordable design is confirmatory.** The k sweep has ~21% power for its discriminating prediction
and ~34% at eight times the seed budget. A held-out sign test on 15 pilots has ~26%. Eighty percent
power would need ~60 pilots. This is not a scoping problem to be optimized around — it is a property
of the experimental system, and it belongs in the pre-registration and in the talk.

Once stated plainly the decision changes shape. Don't choose the least-underpowered pilot-level design.
Either accept Claim 4 as permanently hypothesis-generating and say so, or **move the unit of analysis
off the pilot.** This plan does both: it re-scopes the sweep as descriptive and honest, and it scopes
the experiment that actually escapes n = pilots.

---

## Part 0 — Ledger

### Settled

| finding | evidence |
|---|---|
| A crystal is a point estimate of a distribution | conformation SD 0.235 vs seed SD 0.16, n = 701, baseline checkpoints |
| Shared-label augmentation teaches conformational invariance | across-conformer SD 0.235 → 0.180 / 0.172, 11 of 12 pilots, p = 0.027 |
| The gain decays with distance from training frames | +0.246 own complex → +0.03–0.11 own family → +0.025 unseen → ≈0 cross-benchmark |
| Memorization is label-side, not input-side | median feature-distance ratio 0.81; ERR 1.07, NFAT 1.34 |
| The bend axis is fluctuation, not magnitude | reverse partial flips sign (−0.136 → +0.229); rank regression IQR β = −0.834, p = 0.0065 |
| The two arms are near-equivalent | axis ρ = +0.927; mean IQR 5.46° → 5.79°; readout Δr = +0.011, p = 0.36 |
| Bend fluctuation does not act per-state | within-pilot slopes 7/11, sign test p = 0.55 |
| irf's bend channel is broken and unrecoverable | frozen-vs-crystal deviation 153.9° vs a 3.8° median; `bend_pp` fails identically across all 184 states |
| The own-family middle point is clean | 0 of 11 pilots' reference PDBs appear in their own-family subsets |
| The cross-pilot gain is not family-mediated | broad-class proxy flat — though all 132 pairs are cross-family by construction |
| No inference-only endpoint exists | dispersion null and redundant with `r_mean`; `r_mean` itself frozen-arm-only at p = 0.098 |
| **Seed noise is ⅔ of the between-pilot spread** | SEM/between-SD = 0.644 frozen, 0.661 relaxed; no arm asymmetry (Wilcoxon p = 0.63) |
| **No pilot-level design is confirmatory** | ~21% power for the interaction at 5 seeds; ~60 pilots for 80% |

### Retracted — do not quote

The +0.28 held-out crystal gain (−0.008 cross-pilot). The information-content deficit as a *biological*
claim. The SSO7D calibration. "The mechanism axis strengthens in the relaxed arm." v1's HSF leakage
control. The frame-count bias argument for reachability. "The only correlation whose bootstrap CI
excludes zero" as a significance claim.

### Corrections to v4

- **`r_mean` is not an endpoint.** v4 inverted +0.552 from the dispersion table without noting it was
  the *frozen* arm. Direct computation: −0.552 (frozen, p = 0.098), −0.321 (relaxed, p = 0.366). A
  frozen-arm-only, non-significant relation is not a validated substitute for ΔP.
- **v4's state-count confound does not exist.** State count vs ΔP is +0.006 (frozen, p = 0.99) and
  +0.152 (relaxed, p = 0.68). It does not predict the endpoint, so it cannot mediate anything. The
  partial is slightly *larger* than the raw correlation.
- **The test family grows to 18, not shrinks.** `r_mean` vs bend IQR is ρ = −0.030, so the axes are
  independent and `r_mean` × 2 arms are two genuinely new tests. **Claim 4's q becomes ≈ 0.15–0.17.**
- **Drop the disattenuated −0.89.** Disattenuation corrects error in y while assuming x is error-free,
  but bend IQR is itself estimated from ~90 states. A one-sided correction at n = 10 is model-dependent
  and will be read as a claim. Report the attenuation ratio (0.66) as a limitation and stop.
- **Methods note:** ets1 and runx share four PDBs (3wts, 3wtu, 3wtw, 4l0z — the RUNX–ETS composites),
  scoring different chains with different labels. Correct handling, but 4 of runx's 5 and 4 of ets1's 8
  own-family entries sit on shared structures, so those two pilots' effects are not fully independent.

### The standing anomaly is now worse, and should be said out loud

The arms are physically near-identical (axis ρ = +0.927, readout Δr = +0.011, p = 0.36) **and equally
precise** (median SEM 0.0334 vs 0.0346, Wilcoxon p = 0.63). So there is no remaining mechanism for why
relaxed-effect correlates at −0.745 and frozen-effect at −0.564. Differential noise was the last
candidate explanation and it is dead.

**Consequence: report the pooled-arm correlation as primary.** If the arms measure the same thing with
equal precision, averaging their ΔP per pilot uses both estimates instead of the one that happened to
correlate better. Compute it properly with bootstrap and LOO; the rough value is ≈ −0.71. It is harder
to attack than either cell and it removes a question nobody can answer.

---

## Part 1 — The free power fix (do this first, 1 hour)

**Nobody has decomposed the variance in ΔP.** Per-seed ΔP is a **mean over own-family entries**, and
those subsets range from 3 (hsf) to 24 (egr1). If entry-to-entry variance dominates seed-to-seed, the
per-seed SD is mostly small-sample noise **in the eval set**, not training stochasticity.

The audit's own table points this way: egr1 has the most entries (24) and the lowest SEM (0.018); hsf
has 3 entries at 0.036; err has 10 entries and 4 seeds at 0.067.

**Why it matters more than any other hour available.** Expanding own-family entry sets is **eval-side
only — no retraining.** There are 130 benchmark entries and most own-family subsets are in single
digits. If the family definition can be loosened (PFAM clan, or the structural classes used in the
§3.5 proxy), you buy precision with inference cost instead of training cost.

**What.** Mixed model on `perseed_perentry.csv` with entry and seed as crossed random effects; extract
variance components per pilot.

| result | reading |
|---|---|
| Entry variance dominates | **Cheapest power gain in the project.** Expand own-family subsets by loosening the family definition; SEM falls without a single retrain; every pilot-level analysis improves at once. Recompute all power numbers afterward. |
| Seed variance dominates | The noise is training stochasticity. More seeds is the only lever, and the audit's table already shows what that buys (21% → 31% at 4× the budget). The power ceiling stands. |
| Mixed, pilot-dependent | Expand entries for the small-subset pilots (hsf 3, irf 4, csl 5, runx 5, foxa 5) where the gain is largest; leave egr1 alone. |

**This gates every power calculation below.** Run it before timestamping Part 2.

---

## Part 2 — The k sweep, re-scoped as descriptive (pre-registration)

**Timestamp this section and do not edit it after the calibration starts.**

### 2.0 — Stated power limitation

> The interaction test has approximately **21% power** at 5 seeds and n = 12 pilots for a true ρ of 0.5,
> and approximately 31% at 20 seeds. A null result therefore **does not distinguish "no mechanism" from
> "underpowered."** The pre-registered reading of a null is **"uninformative on the interaction,"** not
> "mechanism wrong." The interaction will be reported as an **effect size with a confidence interval**,
> and no hypothesis test on it will be performed.

Pre-registering an underpowered test is defensible. Pre-registering it without recording that it is
underpowered is not.

### 2.1 — Calibration gate (2 pilots, 1 day, run first — and it is the well-powered part)

Generate levels 3 and 4 for one high-IQR pilot (tbp, 11.78) and one low (runx, 1.29). Measure cleaned
bend IQR at each level.

This is a **within-pilot geometry measurement on ~90 states**, not a 12-point correlation, so it is
adequately powered and worth running regardless of what happens to the rest of the sweep.

| result | action |
|---|---|
| IQR separates levels, high-IQR pilot moves more | Proceed. |
| IQR separates levels, both move similarly | Proceed, noting the interaction prediction is weakened before the endpoint is measured. |
| IQR does not separate levels | **Stop.** k has not manipulated the axis. Go straight to Part 3. One day spent, weeks saved. |
| Structures degrade at k = ∞ (fnat pass collapses) | Drop level 4; run {0, 1.5, 5.0}; state why. |

### 2.2 — Design

| level | k |
|---|---|
| 1 | 0.0 (frozen arm) |
| 2 | 1.5 (existing relaxed arm) |
| 3 | 5.0 (the manipulation) |
| 4 | ∞ (ceiling; likely end-fraying) |

**Pilot-level endpoint:** ΔPearson at {0, 5.0}, two arms per pilot. Paired per pilot, so each pilot is
its own control — this removes every between-pilot confound, which is the sweep's genuine advantage and
should be stated as such. It does not escape n = 12.

**Seed budget** contingent on Part 1.

### 2.3 — Stated predictions (reported as effect sizes, not tested)

1. **Direction.** Softer restraints → more conformational disagreement with a shared label → the
   augmentation effect degrades, most in pilots that already fluctuate.
2. **Interaction.** Δ(ΔP) between k = 0 and k = 5.0 correlates negatively with intrinsic bend IQR.
3. **A main effect alone is not confirmation** — it is consistent with structures simply getting worse.
4. **Monotonicity**, if four levels run.
5. **An interaction in the wrong direction is to be taken seriously**, not explained away. It would mean
   fluctuating DNA benefits from softer restraints, contradicting the shared-label reading.

### 2.4 — Mandatory quality controls at every level

fnat pass rate, iRMSD distribution, surviving-state count, crystal-referenced range assertion on
`bend_uu` (|state − crystal| > 90°; the median-referenced version cannot catch an irf-type failure).
**If fnat pass rate falls with k, any readout change is confounded with structure quality.** Report
beside the endpoint, not in an appendix.

### 2.5 — The state-level precursor (the only well-powered part of the sweep)

Per-state bend change × per-state readout change across k, pilot as random effect, **n ≈ 700 states with
x manipulated.** The same BioEmu frame is minimized under each k, so states correspond across levels and
the pairing is exact.

**Label it correctly.** The endpoint here cannot be ΔP — ΔP is irreducibly pilot-level (one shared
label, one ensemble, one number). At state level the endpoint is r-to-crystal, which measures *"does DNA
restraint strength change what the model reads,"* **not** *"does DNA fluctuation drive the augmentation
effect."* That is the same construct substitution that made step 2 not-the-same-question — with a
manipulated x instead of an observed one.

So: a **precursor test**, not the mechanism test. And the prior from step 2 (7/11 slopes, n.s.) is not
encouraging.

| result | reading |
|---|---|
| Restraint strength changes per-state readout **more in high-IQR pilots** | A real mechanistic link, manipulated x, n ≈ 700. The strongest result the sweep can produce, and it should be the sweep's headline rather than the 12-point interaction. |
| Uniform change across pilots | Restraints affect the readout generically. Check QCs; likely structure quality. |
| No change | Consistent with step 2. The readout is insensitive to DNA restraint strength at the state level, which constrains the mechanism to operate only through the aggregate label. |

---

## Part 3 — The sequence-design experiment (the only confirmatory route)

Previously filed as a months-away aspiration. It is now **the only design on the table that escapes
n = pilots**, and it should be scoped in parallel with the sweep rather than after it — the sweep's
result will not change whether this is needed.

### 3.1 — Design

**One TF. 20–40 designed DNA variants spanning intrinsic flexibility.**

- **Protein fixed** → every between-pilot confound vanishes at once. No family, no DNA length, no crystal
  quality, no ensemble size, no PWM source.
- **Flexibility manipulated by sequence** — A-tract versus GC-rich flanks — and predictable in advance
  with the DNA shape-prediction tools already in the stack.
- **n = sequences, not pilots.** A manipulated predictor with no between-system confounding is a
  completely different power regime from anything in this project.
- **ΔP is not needed.** The readout is a per-sequence prediction, which sidesteps the shared-label
  problem, the training arms, and the seed noise simultaneously.

### 3.2 — Pilot selection

**ets1** is the obvious candidate: highest own-family effect (+0.111), 100% apo pass, clean geometry
(frozen-vs-crystal deviation 0.1°, the best in the set), 99 surviving states, and good experimental
coverage in the family.

Second choice if ets1 has a problem: **engrailed** (85 states, r = 0.84 to crystal, clean geometry) —
but the homeodomain's induced-fit literature makes it a harder system to interpret.

### 3.3 — What to measure

For each designed sequence: predicted intrinsic flexibility (shape tools), BioEmu ensemble → dock →
minimize → ensemble bend IQR (measured), and the DeepPBS readout per state and per ensemble.

**The prediction:** measured bend IQR tracks predicted intrinsic flexibility (a positive control on the
manipulation), and the readout degrades with flexibility.

### 3.4 — Why this settles what nothing else can

Every result in the project so far is a correlation between a structural property and a model-derived
effect, across systems that differ in many correlated ways. The audit showed the confounding worry is
not even the binding one — *nothing measured predicts ΔP*, which makes the controls uninformative rather
than reassuring. What binds is selection and n.

This design removes both: the predictor is manipulated, the system is fixed, and n is set by how many
sequences you design.

### 3.5 — Scoping questions to answer this week

1. Can 20–40 sequences be designed that span a meaningful intrinsic-flexibility range while keeping the
   ets1 core motif fixed? (Flanks only, or does the core have to vary?)
2. What does one sequence cost end-to-end — BioEmu + Stages 1–4 + inference?
3. Is there experimental specificity data for any of the designed variants, or is this purely a
   model-behaviour experiment? (Either is publishable; they are different papers.)
4. Does the shape-prediction tool's flexibility estimate correlate with measured ensemble bend IQR on
   the existing 12 pilots? **This is a free retrospective validation of the manipulation** and should be
   run before designing anything.

Question 4 is an hour and it de-risks the whole experiment.

---

## Part 4 — Remaining items

### §4.1 — Size-matched donor control (1 training arm; still highest stakes)

The last live alternative explanation for +0.025 — the only clean positive generalization effect in the
project. Donor Y's augmented checkpoint saw BioEmu frames **and** ~90 extra training examples.

**Control:** Y's baseline plus ~90 additional crystal entries from the benchmark, no BioEmu. Score
targets X ≠ Y as in the cross-pilot design.

| result | reading |
|---|---|
| ≈ 0 | +0.025 is conformational information; Claim 2's accuracy component survives. |
| ≈ +0.025 | **Augmentation is just more data.** Claim 2 becomes purely about invariance, with no accuracy benefit at all. |
| > +0.025 | Conformations are *worse* than crystal entries at matched volume. Sharpens the shared-label argument. |

### §4.2 — Interface restriction (0.5 d)

Per-bp bend variance along the duplex; IQR over footprint base pairs only; re-run the correlation.
**Prerequisite for k = ∞**, which will worsen end-fraying. Collapse ⇒ the axis was end-fraying and
Claim 4 comes out.

### §4.3 — Open and unexplained

- **HSF's anomalies** — worst agreement (0.27 → 0.61), sole pilot where across-conformer SD *rises*
  under augmentation (0.386 → 0.519). Cleared of parse failure (deviation 11.0°). Needs a real
  explanation; the variance exception is worth a slide.
- **The +0.59** — within-pilot slope vs bend IQR, backwards under every reading. Unaffected by cleaning.
- **Ch. 19½ provenance** — `analysis/analyses/nmr_rmsd/` results still not located. Open since the start
  of the session, still the cheapest unresolved item. Until resolved, Ch. 19½ stays out of the talk.
- **Chapter-map reconciliation** — two outlines with incompatible numbering. Include the
  `froz_bend_iqr` → `bend_iqr_frozen` rename; the arm-ambiguity that caused C2 is likely elsewhere too.

### §4.4 — Dropped

**The held-out sign test on 10–15 pilots.** My earlier recommendation, and the arithmetic doesn't
support it: at true ρ ≈ 0.5, sign accuracy is ≈ 0.67, so n = 15 needs 12 correct for p < 0.05 and
P(≥12) ≈ **0.26** — the same power as the sweep, before attenuation. A confirmatory version needs 30–40
pilots and costs more than Part 3. Reconsider only if Part 1 materially cuts the noise.

---

## Part 5 — Order

| # | item | cost | why here |
|---|---|---|---|
| 1 | Entry-vs-seed variance decomposition | 1 h | may cut SEM for free; gates every power number |
| 2 | Shape-tool flexibility vs measured bend IQR on 12 pilots | 1 h | free validation of Part 3's manipulation |
| 3 | Expand own-family entry sets, if item 1 says to | eval-only | no retraining |
| 4 | Pooled-arm correlation with bootstrap + LOO | 1 h | removes a question nobody can answer |
| 5 | **Timestamp Part 2**, then §2.1 calibration | 1 d | well-powered; can cancel the sweep |
| 6 | §4.1 size-matched donor | 1 arm | highest stakes |
| 7 | §4.2 interface restriction | 0.5 d | prerequisite for k = ∞ |
| 8 | §2.5 state-level precursor | with the sweep | the only well-powered part |
| 9 | Part 3 scoping questions | this week | de-risks the confirmatory experiment |
| 10 | Full sweep, if the gate passes | weeks | descriptive, not confirmatory |
| 11 | Sequence-design experiment | months | the one that settles it |

Items 1, 2 and 4 total three hours and change how everything after them is run.

---

## Part 6 — The talk

**Claims 1–3 unchanged.** None depends on the pilot-level regression.

1. **A crystal is a point estimate of a distribution.** Conformation SD 0.235 vs seed SD 0.16, n = 701.
2. **Shared-label augmentation teaches conformational invariance.** SD 0.235 → 0.180, 11 of 12 pilots,
   p = 0.027, while accuracy on unseen complexes rises only +0.025. The model maps ~90 structurally
   distinct conformations — not near-copies, median feature distance 0.81 — onto one motif. The training
   objective working as specified. *Subject to §4.1.*
3. **The gain decays with distance from the training frames.** +0.246 → +0.03–0.11 → +0.025 → ≈0. Three
   designs, one monotone gradient, all three points now audited.

**Claim 4, final form:**

> Ensemble DNA-bend fluctuation is the largest of ~18 axis×arm correlations tested on ten pilots
> (ρ = −0.745 arm-matched and cleaned; pooled-arm ≈ −0.71; uncorrected p = 0.013; BH q ≈ 0.15–0.17). It
> survives every single-covariate control — and notably, **nothing** we measured predicts the
> augmentation effect, which makes those controls uninformative rather than reassuring. It does not
> operate through individual states (7/11 slopes, n.s.). Seed noise is two-thirds of the between-pilot
> spread, so the estimate is attenuated. At n = 12 pilots with this noise, roughly **60 pilots** would be
> needed for 80% power on the discriminating prediction: **no affordable pilot-level design can confirm
> this.** It is a hypothesis, and testing it requires moving the unit of analysis off the pilot.

**Volunteer, do not defend:** step 2's within-pilot null; the failed multiple-testing correction; the
in-sample structure of the augmented arms; the power ceiling. Each is what an audience would go looking
for, and finding them first is the most credibility-buying move available.

**Say out loud:** the binding sample size is the pilot count; the pilot set is one-per-family across
twelve families, so no pilot-pair design here can test family transfer; and the next experiment is
designed around escaping that, not around optimizing within it.

---

## References

- Mitra, R. *et al.* **Geometric deep learning of protein–DNA binding specificity.** *Nature Methods*
  **21**, 1674–1683 (2024). doi:10.1038/s41592-024-02372-w
- Jiang, Y., Shewchuk, A. M., Chiu, T.-P., Li, J., Kribelbauer-Swietek, J. F., Gompel, N. & Rohs, R.
  **Readout of intrinsic and induced DNA shape by homeodomain transcription factor complexes.**
  *Biophysical Journal* (2026). doi:10.1016/j.bpj.2026.03.036 — *"minor groove width reflects intrinsic
  DNA structure, whereas minor groove width fluctuations capture protein-induced stabilization and
  reshaping of the DNA."* The per-position profile format for §4.2, and the sequence-dependence result
  that underwrites Part 3's manipulation.
- Lewis, S. *et al.* **Scalable emulation of protein equilibrium ensembles with generative deep
  learning.** *Science* (2025). doi:10.1126/science.adv9817 — relevant only if the Ch. 19½ provenance
  question resolves.
