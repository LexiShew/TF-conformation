# Analysis plan v4 — 2026-09-29

**Supersedes** v3. Written after the v3 execution report, which closed four items and demoted one claim.

**What this document is.** Everything still to run, in order, each with what it decides and how to read
every outcome before you see it. Part 2 is a **pre-registration** for the k sweep and is meant to be
timestamped and left unedited once the calibration starts.

**The state in one paragraph.** Claims 1–3 (point-estimate, invariance, leakage gradient) rest on
n = 701 states or on three independently-measured points and are not in dispute. Claim 4 (the mechanism)
rests on ten pilots, is the best of ~16 axis×arm tests, survives every single-covariate control, does
not survive correction, and does not operate through individual states. It is a hypothesis. The k sweep
is its test, and three half-day audits determine whether that test is worth running and how to size it.

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
| Bend fluctuation does not act per-state | within-pilot slopes 7/11, sign test p = 0.55, ≈ −0.01 r per degree |
| irf's bend channel is broken and unrecoverable | frozen-vs-crystal deviation 153.9° against a 3.8° median; `bend_pp` fails identically across all 184 states |
| The cross-pilot gain is not family-mediated | broad-class proxy flat (+0.0017 vs +0.0183, p = 0.75) — but see §1.1 |
| Crystal-free PWM dispersion is not a usable endpoint | null against bend IQR (−0.079) and ΔP (+0.309); redundant with `r_mean` (−0.982) |

### Retracted — do not quote

The +0.28 held-out crystal gain (−0.008 cross-pilot). The information-content deficit as a *biological*
claim (reverses under augmentation, 7.50 → 8.67/8.80 bits). The SSO7D nonspecific calibration, which
supported that reading. "The mechanism axis strengthens in the relaxed arm" (axis arm is irrelevant).
v1's HSF-crystal leakage control (HSF is 5th worst). The frame-count bias argument for reachability
(~0.003 Å against a 0.103 Å gap). **"The only correlation whose bootstrap CI excludes zero," as a
significance claim** (q ≈ 0.07–0.15 on the honest family).

### Corrections to v3 itself

- **The range assertion must be crystal-referenced, not median-referenced.** A ">90° from the pilot
  median" rule cannot catch irf, whose entire distribution is displaced. Use
  |state `bend_uu` − crystal `bend_uu`| > 90°: flags all 184 irf states plus 34 isolated frames in
  dux4/err/foxa/nfat/runx.
- **§3.5 was not computable and its conclusion was overstated.** All twelve pilots are in twelve
  different families, so all 132 ordered pairs are cross-family *by construction*. The flat result is
  not evidence against an own-family effect; the design cannot test one. See §1.1 for the test that can.
- **§3.3 (HSF) is cleared of irf's failure mode** — frozen-vs-crystal deviation 11.0°, an order of
  magnitude below irf. Its anomalies need a real explanation, not a parse fix.
- **§1.3's contingency has fired.** No inference-only endpoint exists from dispersion, so the sweep runs
  with ΔP — **unless §1.2 below rescues it.**

### The cleaned headline

irf dropped, 34 outlier states removed:

| axis arm | effect arm | matched | ρ | p | boot 95% CI | LOO range |
|---|---|---|---|---|---|---|
| frozen | frozen | yes | −0.564 | 0.0897 | [−0.925, +0.205] | [−0.783, −0.400] |
| relaxed | relaxed | **yes** | **−0.745** | **0.0133** | [−1.000, −0.170] | [−0.933, −0.650] |

n = 10.

---

## Part 1 — Three audits that size the sweep (half a day total)

All three change how Part 2 is run. Do them first.

### §1.1 — Leakage-free own-family recompute (1 h)

**Why.** Claim 3's decay gradient has three points. The outer two are audited: +0.246 own-complex
(in-sample, verified per entry) and +0.025 cross-pilot (clean design). **The middle point — own-family
ΔP, +0.03 to +0.11 — has never been audited**, and it is the one doing the most interpretive work,
since it is the effect size the talk quotes as defensible.

**What.** Ch. 16's own-family ΔP is `augmented_Y` vs `baseline_Y` on the 130-entry benchmark restricted
to family(Y). Those benchmark entries are not Y's conformations, so the design is nearly clean already.
The one possible leak: **Y's own crystal entry sitting inside its own-family benchmark subset.** Check
per pilot; if present, drop it and recompute.

| result | reading |
|---|---|
| No pilot's crystal is in its own-family subset | The middle point is clean as measured. Claim 3's gradient stands fully audited. |
| Some pilots leak; recomputed ΔP is materially smaller | The gradient is steeper than reported and the "defensible effect size" needs revising downward toward +0.025. |
| Recomputed ΔP is unchanged | Same as the first row, with the check on record. |

**Secondary finding to surface regardless:** the pilot set is one-per-family across twelve families.
That was presumably deliberate for coverage breadth, and the cost is that **no pilot-pair design can
ever test family transfer.** Worth one line in the methods and one in the talk.

### §1.2 — `r_mean` as an endpoint, with the state-count partial (2 h)

**Why.** Buried in §2.2's table: dispersion vs own-family ΔP (frozen) = +0.552, p = 0.098, n = 10.
Since dispersion is −0.982 redundant with `r_mean`, that inverts to **`r_mean` vs ΔP at ρ ≈ −0.55** —
pilots whose states agree poorly with the crystal benefit more from augmentation. This is v2's §9,
which has never been run with its mandatory control.

**The mandatory control.** Surviving-state count correlates with crystal agreement at ρ = 0.60, and
state count **is** the fnat pass rate **is** the Ch. 15 rigidity ordering. Without partialling state
count you will rediscover rigidity and call it a mechanism.

**Two things ride on this.** First, `r_mean` is inference-only and already computed for every pilot — if
it survives, the four-level k sweep is affordable again. Second, it is the second-largest effect anyone
has found, so **it belongs in the multiple-testing family**, which means the current 16-test correction
is incomplete.

| result | reading |
|---|---|
| Survives the partial | A cheap inference-only endpoint exists. Four-level sweep is back on; new pilots no longer need training arms; add it to the test family and recompute q for everything. |
| Disappears after partialling | The augmentation gain is a rigidity story, not a conformational one. Worth knowing before spending 24 training arms — and it would reframe Claim 3's middle point. |
| Survives but is redundant with bend IQR | The two axes are one axis. Simplifies the mechanism claim and reduces the effective test family, which raises the surviving q. |

### §1.3 — Per-pilot ΔP standard error across seeds (afternoon)

**Why.** Relaxed ΔP rests on 4–6 seeds and nobody has computed its per-pilot SEM. Measurement error in
y attenuates a correlation, so if seed noise is a meaningful fraction of the between-pilot spread, the
sweep is underpowered before it starts.

**What.** A mixed model over entry × seed with the appropriate random effects, giving per-pilot CIs on
ΔP rather than a 5-point mean. Compare the median per-pilot SEM against the between-pilot SD of ΔP.

| result | reading |
|---|---|
| SEM small relative to between-pilot spread | ΔP is a reasonably precise y. Sweep at two levels × current seed count is adequately sized. |
| SEM comparable to between-pilot spread | **The observed ρ is attenuated and the sweep is underpowered.** Spend on *more seeds at two levels* rather than more levels — the sweep's power comes from precision in y, not from the number of k values. |
| SEM differs systematically between arms | Resolves the standing anomaly: relaxed ΔP correlates *better* despite being assumed noisier, which is backwards under attenuation. If relaxed SEM is smaller, the lucky-draw worry shrinks and you can say so quantitatively. |

---

## Part 2 — The k sweep (pre-registration)

**Timestamp this section and do not edit it after the calibration starts.** The entire force of an
intervention comes from having stated the prediction beforehand, and this project has already shown how
easily a test family gets chosen after the winner is known.

### 2.1 — Calibration gate (2 pilots, 1 day, run first)

Generate levels 3 and 4 for **one high-IQR pilot** (tbp, 11.78) and **one low** (runx, 1.29). Measure
cleaned bend IQR at each level.

**Gate:** the levels must separate the axis. k = 1.5 moved mean IQR by 0.33° and the readout by 0.011,
so there is a live possibility k = 5.0 does not move it either.

| result | action |
|---|---|
| IQR separates levels, high-IQR pilot moves more | Proceed to the full sweep. |
| IQR separates levels, both pilots move similarly | Proceed, but note the interaction prediction is already weakened before the endpoint is measured. |
| IQR does not separate levels | **Stop.** The sweep has not manipulated the axis. Go to alternative restraint schemes or to the sequence-design variant. One day spent, two weeks saved. |
| Structures degrade badly at k = ∞ (fnat pass collapses) | Drop level 4, run {0, 1.5, 5.0}, and state why. |

### 2.2 — Design

| level | k | note |
|---|---|---|
| 1 | 0.0 | frozen arm, continuity |
| 2 | 1.5 | existing relaxed arm, continuity |
| 3 | 5.0 | the genuine manipulation |
| 4 | ∞ | ceiling; likely end-fraying |

**Endpoint:** ΔPearson at levels {0, 5.0} — two arms per pilot, not four — **unless §1.2 validates
`r_mean`**, in which case all four levels run inference-only and ΔP is measured only at the two levels
the inference sweep says differ most.

**Unit of analysis:** per-pilot **change** in ΔP between levels, paired. Each pilot is its own control,
which removes every between-pilot confound at once. This is the sweep's real advantage and it should be
stated as such — it does not escape n = 12, but it does escape confounding.

### 2.3 — Stated predictions

1. **Direction.** Softer DNA restraints → more conformational disagreement with a single shared label →
   the augmentation effect **degrades**, most in pilots that already fluctuate.
2. **Interaction (the discriminating prediction).** Δ(ΔP) between k = 0 and k = 5.0 correlates
   **negatively** with intrinsic bend IQR, |ρ| ≳ 0.5 on n ≈ 12.
3. **Main effect alone is not confirmation.** A uniform shift in ΔP with k is consistent with structures
   simply getting worse. Only the interaction is specific to the hypothesis.
4. **Monotonicity**, if four levels are run.

### 2.4 — Mandatory quality controls at every level

fnat pass rate, iRMSD distribution, surviving-state count, and the crystal-referenced range assertion on
`bend_uu`. **If fnat pass rate falls with k, any readout change is confounded with structure quality**
and the dose–response is uninterpretable without it. Report these beside the endpoint, not in an
appendix.

### 2.5 — Reading the outcome

| result | reading |
|---|---|
| Interaction present, QCs stable | The mechanism is supported by an intervention with each pilot as its own control. This is the result that moves the project from correlation to finding. |
| Main effect, no interaction | Softer restraints change the readout uniformly — likely structure quality. Check the QCs before any other interpretation. |
| Interaction present but QCs degrade with k | Uninterpretable as run. Re-run with interface-restricted metrics (§3.2) and a fnat-matched state subsample. |
| Flat | Either the mechanism is wrong or the sweep did not move the axis — **§2.1 already distinguishes these**, which is why it gates. |
| Reversal at k = ∞ only | Almost certainly end-fraying. Report {0, 1.5, 5.0} and say why ∞ was excluded. |
| Interaction in the **wrong** direction | Take it seriously rather than explaining it away. It would mean fluctuating DNA *benefits* from softer restraints, which contradicts the shared-label reading and would need the whole mechanism restated. |

---

## Part 3 — Remaining items

### §3.1 — Size-matched donor control (1 training arm; highest stakes)

**The last live alternative explanation for +0.025.** Donor Y's augmented checkpoint differs from Y's
baseline in two ways: it saw BioEmu frames, *and* it saw ~90 extra training examples.

**Control:** Y's baseline plus ~90 additional **crystal** entries from elsewhere in the benchmark, no
BioEmu. Score targets X ≠ Y as in the cross-pilot design.

| result | reading |
|---|---|
| Size-matched gives ≈ 0 | +0.025 is conformational information. The project's only clean positive generalization effect survives. |
| Size-matched also gives ≈ +0.025 | **Augmentation is just more data.** Conformations contribute nothing beyond volume; Claim 2 becomes purely about invariance, with no accuracy benefit at all. |
| Size-matched gives more | Conformations are *worse* than crystal entries at matched volume. Also a finding, and it sharpens the shared-label argument. |

### §3.2 — Interface restriction (0.5 d)

Per-base-pair bend variance along the duplex; recompute IQR over footprint base pairs only; re-run the
correlation. **Now also a prerequisite for k = ∞**, which will make end-fraying worse.

| result | reading |
|---|---|
| ρ holds or strengthens | Real and mechanically located. |
| ρ collapses | The axis was end-fraying. Claim 4 comes out. |
| Splits by pilot | Check duplex length and footprint fraction first. |

### §3.3 — Open and unexplained

- **HSF's three anomalies** — worst agreement (0.27 → 0.61), sole pilot where across-conformer SD *rises*
  under augmentation (0.386 → 0.519), crystal featurized fresh. Cleared of parse failure. Needs a real
  explanation; the variance exception is the part worth a slide.
- **The +0.59** — per-pilot within-pilot slope vs bend IQR, backwards under every reading. Unaffected by
  cleaning (irf was already excluded). Candidates: range restriction, attenuation. An audience member who
  finds it will ask.
- **Ch. 19½ provenance** — `analysis/analyses/nmr_rmsd/` results still not located. Open since the start
  of the session and the cheapest unresolved item on the list. Until resolved, Ch. 19½ stays out of the
  talk and §7 (apo-fnat × ΔP) cannot be scheduled.
- **Chapter-map reconciliation** — two outlines with incompatible numbering. Fix before reconciling any
  numbers; include the `froz_bend_iqr` → `bend_iqr_frozen` rename, since the arm-ambiguity that caused
  C2 is likely present in other columns.

---

## Part 4 — Order

| # | item | cost |
|---|---|---|
| 1 | §1.1 leakage-free own-family recompute | 1 h |
| 2 | §1.2 `r_mean` endpoint + state-count partial | 2 h |
| 3 | §1.3 per-pilot ΔP SEM | afternoon |
| 4 | **Timestamp Part 2**, then §2.1 calibration | 1 d |
| 5 | §3.1 size-matched donor | 1 arm |
| 6 | §3.2 interface restriction | 0.5 d |
| 7 | Full k sweep, if the gate passes | weeks |
| 8 | Held-out pilots with pre-registered sign prediction | weeks |
| 9 | Sequence-design variant | months |

Items 1–3 total half a day and **all three change how item 7 is run.** Item 4's gate can cancel item 7
outright.

---

## Part 5 — On more pilots

**The argument is selection, not power.** At n = 10–11 with ρ = −0.745 you have adequate power for an
effect that size. The problem is that bend IQR was chosen as the largest of ~16 axis×arm tests, so the
estimate is biased upward by construction, and **no re-analysis of the same pilots can correct a
selection bias in those pilots.**

The target is therefore not a larger regression but a **held-out set with a pre-registered sign
prediction**: fit nothing, predict the sign of ΔP from bend IQR on 10–15 genuinely new pilots, and see.
That is worth more than twenty more points in the existing fit.

For reference (Fisher-z, partial with one covariate): n = 11 gives 0.79 power at true ρ = 0.78 and 0.25
at ρ = 0.45; n = 40 gives 0.83 at ρ = 0.45. Post-selection, ρ ≈ 0.45 is the prudent assumption, which
sets n ≈ 40 **if** the goal is a powered regression — but the held-out sign test is the better goal and
needs far fewer.

Candidate pool: the existing 130-entry benchmark. Per-pilot cost hinges entirely on §1.2.

---

## Part 6 — The talk

**Claims 1–3 unchanged** and none depends on the pilot-level regression.

1. **A crystal is a point estimate of a distribution.** Conformation SD 0.235 vs seed SD 0.16, n = 701.
2. **Shared-label augmentation teaches conformational invariance.** SD 0.235 → 0.180, 11 of 12 pilots,
   p = 0.027, while accuracy on unseen complexes rises only +0.025. The model maps ~90 structurally
   distinct conformations — not near-copies, median feature distance 0.81 — onto one motif. The training
   objective working as specified. *Subject to §3.1.*
3. **The gain decays with distance from the training frames.** +0.246 → +0.03–0.11 → +0.025 → ≈0. Three
   designs, one monotone gradient. Volunteer the in-sample structure rather than letting someone find it.

**Claim 4, final form:**

> Ensemble DNA-bend fluctuation is the largest of roughly sixteen axis×arm correlations tested on ten
> pilots (ρ = −0.745 arm-matched, after dropping irf's broken bend channel; uncorrected p = 0.013). It
> survives every single-covariate control. It **does not survive correction for the family of axes
> tested** — q ≈ 0.07–0.15 depending on how the family is counted, since the tests are not independent
> (the two bend axes correlate at ρ = +0.927). And it **does not operate through individual states**
> (7/11 within-pilot slopes, sign test p = 0.55), which is what a shared-label mechanism predicts and is
> a constraint on the claim rather than support for it. At n = 10 this is underpowered rather than null:
> ρ = −0.745 is a large effect that cannot be distinguished from the best-of-sixteen under a global null.
> It is a hypothesis. The k sweep is its test.

**Three things to volunteer rather than defend.** Step 2's within-pilot null. The failed multiple-testing
correction. The in-sample structure of the augmented arms. Each is what an audience would go looking for,
and finding them first is the single most credibility-buying move available.

**Say out loud:** the binding sample size is the pilot count. And the pilot set is one-per-family across
twelve families, so no pilot-pair design in this project can test family transfer.

---

## References

- Mitra, R. *et al.* **Geometric deep learning of protein–DNA binding specificity.** *Nature Methods*
  **21**, 1674–1683 (2024). doi:10.1038/s41592-024-02372-w
- Jiang, Y., Shewchuk, A. M., Chiu, T.-P., Li, J., Kribelbauer-Swietek, J. F., Gompel, N. & Rohs, R.
  **Readout of intrinsic and induced DNA shape by homeodomain transcription factor complexes.**
  *Biophysical Journal* (2026). doi:10.1016/j.bpj.2026.03.036 — *"minor groove width reflects intrinsic
  DNA structure, whereas minor groove width fluctuations capture protein-induced stabilization and
  reshaping of the DNA."*
- Lewis, S. *et al.* **Scalable emulation of protein equilibrium ensembles with generative deep
  learning.** *Science* (2025). doi:10.1126/science.adv9817 — relevant only if the Ch. 19½ provenance
  question resolves.
