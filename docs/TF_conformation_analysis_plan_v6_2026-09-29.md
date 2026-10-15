# Analysis plan v6 — 2026-09-29

**Supersedes** v1–v5. Written at the end of a day in which roughly a dozen analyses ran: the
all-conditions comparison, the arm-matched bend-IQR addendum, the critique response, the cross-pilot
test, the v4 audit, the v5 gating report, and the Deep DNAshape addendum with its two follow-ups.

**Read this first.** The plan has been growing faster than the corrections have been landing. irf's bend
was diagnosed three separate times before today's run established there was never a base-paired duplex
at all. §4.2 has been "the next thing" for four consecutive versions. **v6 is deliberately shorter than
v5.** Part 3 is six items. Everything that is not one of those six is in Part 4 and is not scheduled.

**The organizing rule for this version:** nothing new gets built on a foundation still under audit.

---

## Part 1 — Where things stand

### What holds

| claim | evidence | depends on |
|---|---|---|
| **1. A crystal is a point estimate of a distribution** | conformation SD 0.235 vs seed SD 0.16, n = 701 states, baseline checkpoints | nothing contested |
| **2. Shared-label augmentation teaches conformational invariance** | across-conformer SD 0.235 → 0.180 / 0.172, 11 of 12 pilots, p = 0.027; accuracy on unseen complexes +0.025; ~90 structurally distinct inputs (median feature distance 0.81), not near-copies | §3.4 for the accuracy half |
| **3. The gain decays with distance from the training frames** | +0.246 own complex → +0.03–0.11 own family → +0.025 unseen → ≈0 cross-benchmark; all three points audited | nothing contested |
| **4. Intrinsic and protein-induced DNA fluctuation are decoupled** | Deep DNAshape predicted `-FL` vs pyCurves measured fluctuation: ρ = +0.030 per duplex, +0.224 (q = 0.295) per position for ProT, null for Roll and Twist | §3.2 as a positive control |

Claim 4 is new today and arrived as a byproduct of a failed gate. It is consistent with the lab's own
published distinction — groove width reports intrinsic structure while groove-width *fluctuation*
reports protein-induced reshaping — so the "failed" positive control was testing an intrinsic predictor
against an induced measurement, which was never going to work. It is a citable negative and it stands
independently of everything below.

### What is in doubt

**The mechanism (formerly Claim 4, now Claim 5).** ρ = −0.745 at n = 10, the largest of ~20 axis×arm
tests, BH q ≈ 0.15–0.19, no within-pilot support (7/11 slopes, sign test p = 0.55), and **69% of the
axis's variance is duplex length** with the length effect running opposite to worm-like-chain
expectation. §3.1 decides whether it is a finding or an artifact.

**Claim 2's accuracy component.** +0.025 is the only positive generalization effect in the project, and
§3.4 could take it to zero.

### Retracted — do not quote

The +0.28 held-out crystal gain (−0.008 cross-pilot). The information-content deficit as a *biological*
claim (reverses under augmentation, 7.50 → 8.67 bits). The SSO7D nonspecific calibration. "The mechanism
axis strengthens in the relaxed arm" (the axis arm is irrelevant; only the effect arm moves). v1's HSF
leakage control (HSF is 5th worst). The frame-count bias argument for reachability (~0.003 Å against a
0.103 Å gap). "The only correlation whose bootstrap CI excludes zero" as a significance claim. The
pooled-arm correlation as an improvement (−0.636, CI [−0.987, +0.059] — it crosses zero because the two
arms' errors are correlated at +0.770 against an independent-error ceiling of 0.701).

### Corrections carried forward

- **The range assertion must be crystal-referenced**, not median-referenced. A median-referenced rule
  cannot catch a pilot whose entire distribution is displaced.
- **irf has no duplex.** 0 of 184 states carry base-pair tables against 2,224 of 2,224 everywhere else,
  while its crystal parses fine. The 170° bend was computed on a structure pyCurves never recognised as
  double-stranded. Dropping it was correct; the reason is now documented.
- **Deep DNAshape N-padding biases `-FL` upward** (ρ = +0.535 vs N-count, p = 3e-16). Layer 4 costs 8 bp
  of every duplex, so usable positions are (L−8)/L — 27% to 62% across the pilots. Any analysis must drop
  padded positions or carry N-count as a covariate.
- **Duplex length must be in the covariate list and the test family** (now ~20 tests).
- **No inference-only endpoint exists.** Dispersion is null and redundant with `r_mean`; `r_mean` itself
  is frozen-arm-only at p = 0.098. The k sweep runs with ΔP.
- **Entry expansion does not fix the noise.** Seed variance dominates in 7 of 10 pilots; infinite eval
  entries buy 0.20 → 0.24 power. hsf only (82% reducible), and for its own sake.
- **pyCurves is a faithful Curves/Curves+ reimplementation**, so "homegrown estimator" was the wrong
  frame. The length dependence is a property of the *measure* — `bend_uu` is the angle between terminal
  axis vectors, and terminal levels are the least-constrained part of a curvilinear axis fit. It would be
  present in Curves+ too.

### The power ceiling (structural, not fixable by scoping)

Seed noise is ⅔ of the between-pilot spread (SEM/between-SD = 0.644 frozen, 0.661 relaxed; no arm
asymmetry, Wilcoxon p = 0.63). At n = 12 pilots the k sweep's discriminating prediction has **~21%
power**, rising to ~31% at four times the seed budget. A held-out sign test on 15 pilots has ~26%. Eighty
percent power needs ~60 pilots.

**No affordable pilot-level design is confirmatory.** This is a property of the experimental system and
belongs in the pre-registration and the talk, not in a limitations paragraph.

---

## Part 2 — The standing anomaly, stated because it cannot be resolved

The two arms are physically near-identical (axis ρ = +0.927, mean IQR 5.46° → 5.79°, readout
Δr = +0.011 p = 0.36), **equally precise** (median SEM 0.0334 vs 0.0346, p = 0.63), and their
measurement errors are **correlated** (+0.770 observed against a 0.701 independent-error ceiling). Yet
relaxed-effect correlates at −0.745 and frozen-effect at −0.564.

Differential noise was the last candidate explanation and it is dead. What remains is a genuine arm
difference too small to detect in the readout test, or chance at n = 10. Both should be stated; neither
can be distinguished with this data. **Report relaxed/relaxed as primary with frozen/frozen beside it**,
and say that the gap is unexplained.

---

## Part 3 — What to do, in order

Six items. Items 1–3 are two days and resolve half of what is contested.

### §3.1 — Interface restriction, fixed-width window ⭐

**This has been deferred for four versions and now fixes four separate problems.** It gates the sweep,
Part 4's sequence experiment, and the mechanism claim. Nothing else should run first.

Per-base-pair bend variance along the duplex; recompute IQR over a **fixed-width central window** — N
base pairs where N is the minimum common footprint — not a footprint-restricted one, since footprints
also vary in width and would reintroduce a length term.

What it fixes:

| problem | how |
|---|---|
| duplex length at ρ = −0.83 against the axis | fixed N removes length by construction |
| Deep DNAshape N-padding bias | interior positions are unpadded |
| terminal restraint / end fraying | terminal levels excluded |
| intrinsic vs induced mixing | interface positions are where protein-induced effects live |

**Cheap inputs.** pyCurves emits per-level tables directly (`--format csv` writes one file per
dataframe), so this is not the 1–1.5 day JSON extraction previously costed. And there is a ten-minute
pre-check while it runs: correlate ΔP against IQR × L and IQR × √L — if a simple scaling kills or
preserves the correlation you will know the regime before the per-level work lands.

**Also run:** both axis conventions (`--axis-convention legacy` vs `curvesplus`) on the same structures.
One flag. If the length dependence differs between conventions it is axis-construction; if both give
ρ ≈ −0.83 it is the measure or the structures.

| result | reading |
|---|---|
| Length-free axis holds | The mechanism claim is real and now measured on a scale-free quantity — substantially stronger than the current version. Re-run the headline, covariates and test family on it. |
| Length-free axis collapses | Claim 5 was duplex length plus terminal-vector noise. It comes out of the talk. Far better found now than after the sweep. |
| Holds but reorders pilots | The sweep's moderator changes; the pre-registration uses the new one. |

### §3.2 — Predicted MGW vs measured MGW (1 h)

The positive control that should have run instead of the `-FL` comparison. Mean against mean, intrinsic
against intrinsic, no protein-induced confound. Both quantities exist per position and the join and
padding filter are already built.

| result | reading |
|---|---|
| Tracks well | The predictor and the sequence-side pipeline are sound, and the `-FL` decoupling is biology. Claim 4 becomes fully supported. |
| Also null | Something is wrong upstream — the join, the position indexing, or the strand convention. **Claim 4 cannot be stated until this is resolved**, because you would not know whether the decoupling is real or a plumbing failure. |

### §3.3 — Pipeline integrity audit (0.5 d)

irf produced 184 states with no base-paired DNA and **nothing downstream noticed** — the analysis simply
returned fewer rows. That went undetected for months and was visible today only because someone counted.

- Assert in the pyCurves batch step that every state yields a `local_base_base` table and the count
  equals the state count.
- Broaden it: row-count assertions at every stage. The Stage-4 helix-guard silent skip is a known
  instance of the same pattern.
- Determine whether irf's duplex genuinely melted under the restraints. If it did, that is a
  restraint-protocol finding directly relevant to the sweep, where restraints get **looser** — and
  duplex integrity goes into §3.6's per-level QC list.

### §3.4 — Size-matched donor control (1 training arm)

The last live alternative explanation for +0.025, the only positive generalization effect in the
project. Donor Y's augmented checkpoint saw BioEmu frames **and** ~90 extra training examples.

Control: Y's baseline plus ~90 additional **crystal** entries from the benchmark, no BioEmu. Score
targets X ≠ Y as in the cross-pilot design. One arm on one or two donors, not twelve.

| result | reading |
|---|---|
| ≈ 0 | +0.025 is conformational information. Claim 2's accuracy half survives. |
| ≈ +0.025 | **Augmentation is just more data.** Claim 2 becomes purely about invariance, with no accuracy benefit. |
| > +0.025 | Conformations are *worse* than crystal entries at matched volume. Sharpens the shared-label argument. |

### §3.5 — Sequence-experiment viability gate (small)

**Replaces the p5/p95 design gate, which tested the wrong thing.** The predictor's role has changed:
predicted intrinsic fluctuation is not the x-variable — **measured induced fluctuation is.** The
predictor is only a design aid for choosing sequences likely to span a range, so the decoupling raises
cost (over-sample and measure) rather than invalidating the design.

But it introduces the risk that actually decides viability: **if induced fluctuation is driven mainly by
the protein rather than the sequence, varying sequence at a fixed protein may produce almost no spread
in the measured axis.** The 12 pilots span 3.88× in measured MGW-FL, but they differ in protein *and*
DNA, so that range may be entirely protein-side. The sequence experiment holds the protein fixed — its
whole advantage and possibly its fatal flaw.

**The gate:** build 4–6 sequence variants on **one pilot (ets1)** using the existing stage-3 relax
machinery, and ask whether measured fluctuation spans a usable range at fixed protein.

| result | reading |
|---|---|
| Usable spread | The experiment is viable and can be costed. It remains the only design that escapes n = pilots. |
| Little or no spread | Induced fluctuation is protein-determined. The sequence experiment is dead and the project's confirmatory route has to come from the protein side or the restraint. |

### §3.6 — The k sweep, descriptive (weeks; only after §3.1)

**Pre-registration, to be timestamped after §3.1 fixes the moderator.**

> The interaction test has approximately **21% power** at 5 seeds and n = 12 pilots for a true ρ of 0.5,
> and ~31% at 20 seeds. A null therefore **does not distinguish "no mechanism" from "underpowered."** The
> pre-registered reading of a null is **"uninformative on the interaction."** The interaction will be
> reported as an effect size with a confidence interval; no hypothesis test will be performed on it.

**Calibration gate first (2 pilots, 1 day).** Levels k ∈ {0, 1.5, 5.0, ∞}; generate 3 and 4 for one
high-IQR and one low-IQR pilot and measure the length-free axis at each. This is a within-pilot geometry
measurement on ~90 states, so it is adequately powered and worth running regardless. **If the levels do
not separate the axis, stop** — one day spent, weeks saved.

**Predictions** (reported as effect sizes, not tested): softer restraints → the augmentation effect
degrades, most in pilots that already fluctuate; the interaction is the only discriminating prediction;
a main effect alone is consistent with structures simply getting worse; **an interaction in the wrong
direction is to be taken seriously, not explained away.**

**Per-level QC, mandatory and reported beside the endpoint:** fnat pass rate, iRMSD, surviving-state
count, crystal-referenced range assertion, and **duplex integrity** (§3.3). If fnat falls with k, any
readout change is confounded with structure quality.

**State-level precursor**, correctly labelled: per-state bend change × per-state readout change across
k, pilot as random effect, n ≈ 700 with x manipulated. The endpoint here cannot be ΔP — ΔP is
irreducibly pilot-level — so it measures *"does restraint strength change what the model reads,"* not
*"does fluctuation drive the augmentation effect."* Same construct substitution that made the within-
pilot analysis not-the-same-question, with a manipulated x. **A precursor, not the mechanism test** —
and it is the only well-powered part of the sweep.

---

## Part 4 — Not scheduled

Listed so they are not forgotten, deliberately not in Part 3.

**Two cheap items deferred every round, both still open:**

- **Ch. 19½ provenance.** `analysis/analyses/nmr_rmsd/` results are not on the cluster. One hour, and a
  chapter's worth of numbers depends on it. Until resolved, Ch. 19½ stays out of the talk and the
  apo-fnat × ΔP join cannot be scheduled.
- **Chapter-map reconciliation.** Two outlines with incompatible numbering (`book/OUTLINE.md`, 21
  chapters; `docs/TF_conformation_book_outline.md`, part-based). **Every number reconciliation will
  generate new conflicts until this is fixed.** Include the `froz_bend_iqr` → `bend_iqr_frozen` rename.

**Open and unexplained:**

- **HSF's three anomalies** — worst agreement (0.27 → 0.61), sole pilot where across-conformer SD *rises*
  under augmentation (0.386 → 0.519), crystal featurized fresh. Cleared of parse failure (11.0°
  deviation). Needs a real explanation; the variance exception is worth a slide.
- **The +0.59** — within-pilot slope vs bend IQR, backwards under every reading.
- **Deep DNAshape layer sweep** — run layers 1–4, not layer 2 alone. Lower layers propagate less context,
  so layer-2 `-FL` approaches a pentamer lookup; a single alternative level is ambiguous between
  "coverage helped" and "long-range context was noise." Check predicted-`-FL` variance at each layer too.
- **hsf entry expansion** — 3 entries, 82% reducible. For its own sake, not as a power fix.

**Dropped:**

- **Held-out sign test on 10–15 pilots.** At true ρ ≈ 0.5, sign accuracy ≈ 0.67, so n = 15 needs 12
  correct for p < 0.05 and P(≥12) ≈ **0.26** — the same power as the sweep. A confirmatory version needs
  30–40 pilots.
- **Pooled-arm correlation.** CI crosses zero; the arms' errors are correlated.
- **Entry expansion as a project-level power fix.** 0.20 → 0.24 at the floor.

**Worth exporting rather than filing:** the N-padding artifact is a real trap for anyone using Deep
DNAshape on short oligos, and the mechanism is clean — a padded position's `-FL` expresses uncertainty
about missing context, not a property of the DNA. Tell Jinsen and Tsu-Pei; it may belong in the tool's
documentation. The first-pass ρ = −0.316 at p = 3e-5 with the wrong sign is also a good teaching example
of two opposite-signed position gradients manufacturing significance.

---

## Part 5 — The talk

**Lead with what is measured cleanly and large-n. Finish with the mechanism as motivation, not result.**

1. **A crystal is a point estimate of a distribution.** Conformation SD 0.235 vs seed SD 0.16, n = 701.
2. **Shared-label augmentation teaches conformational invariance.** SD 0.235 → 0.180, 11 of 12 pilots,
   p = 0.027, while accuracy on unseen complexes rises only +0.025. The model maps ~90 structurally
   distinct conformations — not near-copies, median feature distance 0.81 — onto one motif. The training
   objective working as specified, and why augmentation cannot recover conformational signal.
3. **The gain decays with distance from the training frames.** +0.246 → +0.03–0.11 → +0.025 → ≈0. Three
   designs, one monotone gradient, all three points audited.
4. **Intrinsic and protein-induced DNA fluctuation are decoupled.** Sequence-predicted fluctuation does
   not track ensemble-measured fluctuation, per position or per duplex, for any bend-relevant parameter
   — consistent with the lab's finding that groove-width fluctuations report protein-induced reshaping
   while the mean reports intrinsic structure. **Sequence-based shape prediction cannot substitute for
   ensemble measurement when the question is about a bound complex.**

**Claim 5 — the mechanism, as an observation:**

> Ensemble DNA-bend fluctuation is the largest of ~20 axis×arm correlations on ten pilots (ρ = −0.745;
> uncorrected p = 0.013; BH q ≈ 0.15–0.19). It survives every single-covariate control — and notably,
> **nothing** we measured predicts the augmentation effect, which makes those controls uninformative
> rather than reassuring. It does not operate through individual states. Sixty-nine percent of the
> axis's variance is duplex length, and the length effect runs opposite to worm-like-chain expectation.
> Seed noise is two-thirds of the between-pilot spread. At n = 12 with this noise, ~60 pilots would be
> needed for 80% power on the discriminating prediction: **no affordable pilot-level design can confirm
> this.** It is a hypothesis, and testing it requires moving the unit of analysis off the pilot.

**Volunteer, do not defend:** the within-pilot null; the failed multiple-testing correction; the
in-sample structure of the augmented arms; the power ceiling; the duplex-length confound. Each is what
an audience would go looking for, and finding them first is the most credibility-buying move available.

**Say out loud:** the binding sample size is the pilot count; the pilot set is one-per-family across
twelve families, so no pilot-pair design here can test family transfer; and the next experiment is
designed around escaping that rather than optimizing within it.

**A note on how this week reads.** Roughly a dozen analyses: two clean positives, four solid negatives,
three retractions, and a body of methodology. That looks bad and is not. Every negative was found
internally, and several — the in-sample augmented arms, the axis label mismatch, the padding artifact,
irf having no duplex — would have been fatal if found by a reviewer instead.

---

## References

- Mitra, R. *et al.* **Geometric deep learning of protein–DNA binding specificity.** *Nature Methods*
  **21**, 1674–1683 (2024). doi:10.1038/s41592-024-02372-w
- Jiang, Y., Shewchuk, A. M., Chiu, T.-P., Li, J., Kribelbauer-Swietek, J. F., Gompel, N. & Rohs, R.
  **Readout of intrinsic and induced DNA shape by homeodomain transcription factor complexes.**
  *Biophysical Journal* (2026). doi:10.1016/j.bpj.2026.03.036 — *"minor groove width reflects intrinsic
  DNA structure, whereas minor groove width fluctuations capture protein-induced stabilization and
  reshaping of the DNA"*; *"recognition involves a context-dependent balance between conformational
  selection and induced fit."* The basis for Claim 4.
- Li, J., Chiu, T.-P. & Rohs, R. **Predicting DNA structure using a deep learning method.** *Nature
  Communications* **15**, 1243 (2024). Deep DNAshape; the `-FL` features used in §3.2 and §3.5.
- Li, J. **pyCurves.** https://github.com/JinsenLi/pyCurves — Curves/Curves+ reimplementation; the
  per-level tables that make §3.1 cheap, and the `--axis-convention` flag that tests the length
  dependence.
- Lewis, S. *et al.* **Scalable emulation of protein equilibrium ensembles with generative deep
  learning.** *Science* (2025). doi:10.1126/science.adv9817 — relevant only if Ch. 19½'s provenance
  resolves.
