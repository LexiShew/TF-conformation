# Analysis plan v3 — 2026-09-29

**Supersedes** v1 and v2. Written after four completed runs: the all-conditions comparison, the
arm-matched bend-IQR addendum, the critique response (steps 1–3), and the cross-pilot test.

**What changed since v2.** v2 was a list of analyses to run. Most of them ran, and the picture
inverted: the augmentation-as-invariance story is now the strongest material in the project, measured
with large n and a clean design, while the mechanism story is the most speculative thing in it. The
work left is no longer "check the plan's analyses" — it is one intervention (the k sweep), three
loose ends that could each invalidate something, and a talk reordering.

---

## Part 0 — What is settled

Record these so they are not relitigated. Each has a design that does not depend on n = 11 pilots.

| finding | evidence | status |
|---|---|---|
| A crystal is a point estimate of a distribution | conformation SD 0.235 vs seed SD 0.16, n = 701 states, baseline checkpoints | **holds** |
| Shared-label augmentation teaches conformational invariance | across-conformer SD 0.235 → 0.180 (frozen) / 0.172 (relaxed), 11 of 12 pilots, p = 0.027 | **holds** |
| The augmentation gain decays with structural distance from training frames | +0.246 own complex → +0.03–0.11 own family → +0.025 unseen complex → ≈0 cross-benchmark | **holds**, three points measured independently |
| Memorization is label-side, not input-side | median feature-distance ratio 0.81; ERR 1.07 and NFAT 1.34 are *further* from their own crystal than unrelated training entries | **holds** |
| The bend-IQR axis is fluctuation, not magnitude | mean/median-vs-IQR ρ = +0.35/+0.48; reverse partial flips sign (−0.136 → +0.229); rank regression IQR β = −0.834, p = 0.0065 | **holds** |
| The arms are near-equivalent | axis ρ = +0.927 frozen-vs-relaxed; mean IQR 5.46° → 5.79°; readout Δr = +0.011, p = 0.36 | **holds** |
| Bend fluctuation does not act through individual states | within-pilot slopes 7/11 negative, sign test p = 0.55, ≈ −0.01 r per degree | **holds** — and is what a shared-label mechanism predicts |

**Retracted, do not quote:** the +0.28 held-out crystal gain (−0.008 cross-pilot); the information-
content deficit as a *biological* claim (it reverses under augmentation: 7.50 → 8.67/8.80 bits); the
SSO7D nonspecific-binding calibration, which supported the retracted biological reading; "the
mechanism axis strengthens in the relaxed arm" (the axis arm is irrelevant; only the effect arm moves);
v1's HSF-crystal leakage control (HSF is 5th worst, not worst); the frame-count bias argument for
reachability (~0.003 Å against a 0.103 Å gap).

**The one live alternative explanation for the only positive effect** is in §2.1 below. Read it before
quoting +0.025.

---

## Part 1 — The k sweep (the priority)

Step 2 removed within-pilot support for the mechanism. Everything remaining is observational: eleven
pilots, correlating what they happen to be like against what happened when you trained on them. The
k sweep is the only available design that **manipulates the axis** rather than observing it.

### 1.0 Why it is well-motivated now

Three independent measurements say k = 1.5 barely moves anything:

- bend IQR frozen-vs-relaxed ρ = +0.927, rank shifts ≤ 2 positions for every pilot
- mean IQR 5.46° → 5.79°, 8 of 11 pilots increasing
- readout Δr = +0.011, p = 0.36 over 59 TF × seed pairs

So the current "relaxed" arm is not a manipulation of DNA flexibility — it is the same regime with
slightly different coordinates. **A dose–response requires levels that actually differ.**

### 1.1 Design

| level | k | intent |
|---|---|---|
| 1 | 0.0 (rigid) | the frozen arm, for continuity |
| 2 | 1.5 | the existing relaxed arm, for continuity |
| 3 | 5.0 | genuinely softer — where the regime should start to change |
| 4 | ∞ (unrestrained) | the ceiling; DNA free to relax to force-field preference |

**Run the calibration before the science.** Generate levels 3 and 4 for **two pilots first** — one
high-IQR (tbp or nfat) and one low (runx or ets1) — and measure bend IQR at each level. If IQR does
not separate the levels, the sweep has not manipulated the axis and there is no point running the
other ten pilots. That check costs a day and can save two weeks.

**Pilots:** all 12, subject to the calibration passing. **Seeds:** the sweep is about structure
generation, so the expensive training arms are only needed if ΔP is the endpoint — see §1.3.

### 1.2 Pre-registered predictions

If bend fluctuation is the mechanism:

1. **Monotone dose–response.** The readout should move monotonically with k, not just differ at the
   extremes. Monotonicity across four levels is the claim; a difference between two levels is not.
2. **Interaction with intrinsic IQR.** High-IQR pilots (tbp 11.78, nfat 9.95, irf 8.54 — but see §3.1
   on irf) should move *more* than low-IQR pilots (runx 1.29, csl 2.40, engrailed 2.44). This is the
   discriminating prediction: a main effect of k alone is consistent with any number of boring
   explanations (the structures are simply getting worse), whereas the interaction is specific to the
   hypothesis.
3. **Direction.** Under the shared-label reading, more DNA freedom → more conformational disagreement
   with a single label → the augmentation effect should degrade in the pilots that already fluctuate.

### 1.3 The endpoint problem

ΔPearson costs a full training arm per pilot per level. Four levels × 12 pilots × 5 seeds is the
entire project's training budget again.

**Cheaper design, and probably sufficient:** run the sweep with **inference only** on baseline
checkpoints and measure per-state PWM behaviour across k. Then train arms only at the two levels the
inference sweep says differ most. This costs one training arm per pilot instead of four, and it uses
the inference sweep to choose *where* to spend.

This is contingent on §2.2 finding a usable inference-only endpoint. If §2.2 comes back null, the
sweep runs with ΔP at levels {0, 5.0} only — two arms, not four, losing monotonicity but keeping the
interaction test.

### 1.4 How to read every outcome

| result | reading |
|---|---|
| Monotone in k **and** interaction with intrinsic IQR | The mechanism claim is supported by an intervention. This is the result that turns the project from a correlation into a finding. |
| Monotone in k, **no** interaction | Something about softer restraints changes the readout uniformly — likely structure quality degrading, not fluctuation acting. Check fnat and iRMSD across levels before interpreting. |
| Interaction but not monotone | Suspicious. Probably a level-specific artifact (k = ∞ frays the ends). Check §3.2's interface restriction at each level. |
| Flat in k | Either the mechanism is wrong, or the sweep did not move the axis. **The §1.1 calibration distinguishes these**, which is why it runs first. |
| Non-monotone with a reversal at k = ∞ | Most likely end-fraying dominating at the unrestrained level. Re-run with interface-restricted metrics; if the reversal disappears, report levels {0, 1.5, 5.0} and state why ∞ was excluded. |

### 1.5 Mandatory quality controls at each level

The sweep changes structure generation, so structural quality is no longer held fixed. At each k,
record and report: fnat pass rate, iRMSD distribution, surviving-state count, and the §3.1 range
assertion on `bend_uu`. **If fnat pass rate falls with k, any readout change is confounded with
structure quality** and the dose–response is uninterpretable without it.

---

## Part 2 — Two things that must run before the sweep's endpoint is chosen

### 2.1 Size-matched donor control (the last live alternative explanation)

**The problem.** Donor Y's augmented checkpoint differs from Y's baseline in two ways: it saw BioEmu
frames, **and** it saw ~90 extra training examples. Some or all of the +0.025 cross-pilot gain could
be data volume rather than conformational information.

**The control.** One additional donor arm: Y's baseline plus ~90 additional **crystal** entries drawn
from elsewhere in the benchmark, no BioEmu. Score target pilots X ≠ Y as in the cross-pilot design.

**Cost:** one training arm on one or two donors. Not twelve.

| result | reading |
|---|---|
| Size-matched donor gives ≈ 0 | The +0.025 is conformational information. This is the project's only clean positive effect and it survives. |
| Size-matched donor also gives ≈ +0.025 | **Augmentation is just more data.** The conformations contribute nothing beyond volume, and Claim 2½ becomes purely about invariance with no accuracy benefit at all. |
| Size-matched donor gives more than +0.025 | Conformations are *worse* than crystal entries as training data, volume held constant. Also a finding, and it would sharpen the shared-label argument considerably. |

**This is the highest-stakes single result remaining**, because +0.025 is currently the only positive
generalization effect in the project.

### 2.2 Crystal-free PWM dispersion (the endpoint question, properly tested)

The dispersion test in the critique response used the IQR of per-state **r-to-crystal**, which is
crystal-anchored — so its ρ = −0.776 against `r_mean` was close to inevitable, since it measured the
spread of a quantity whose mean is `r_mean`. That test does not settle the endpoint question.

**The crystal-free version:** mean pairwise distance **among the state PWMs themselves**, no reference.
Mean pairwise (1 − Pearson) over the ensemble, or mean KL of each state's PWM from the ensemble-mean
PWM. `compare_pwms.py` already writes `pwms_<tf>_s<n>.npz` with the aligned stacks, so this is a
distance matrix per pilot per seed. **Half a day.**

Test it against: (a) bend IQR, (b) own-family ΔP, (c) `r_mean` — the third to confirm it is *not*
redundant.

| result | reading |
|---|---|
| Correlates with bend IQR **and** with ΔP, and is not redundant with `r_mean` | Usable inference-only endpoint. Makes the four-level k sweep affordable and makes new pilots cheap. |
| Correlates with ΔP but not bend IQR | Still usable as an endpoint; would mean dispersion and fluctuation are different things, which is itself informative. |
| Null against both, or redundant with `r_mean` | Drop the inference-only endpoint idea properly. The sweep runs at two levels with ΔP, and new pilots stay expensive. |

**The bar, from the critique response and correct:** any substitute endpoint must correlate with ΔP at
n = 11 before it is allowed to replace it.

---

## Part 3 — Loose ends that could invalidate something

### 3.1 irf, and the range assertion

`irf`'s relaxed `bend_uu` is 156–179° across all 93 states against a crystal bend of 16.9°. The same
pilot's MGW-FL channel was already dropped for a groove-parse artifact.

**Do not assume the IQR is clean.** "Quartiles are robust to the offset" assumes the failure *is* an
offset. If `bend_uu` is computed on a mis-assigned helical axis, both location and scale are suspect,
and a 5.5° IQR inside a 156–179° band is unexplained rather than reassuring. The MGW-FL channel in the
same pilot was dropped outright, not partially trusted.

**Actions:**
1. Range assertion in the pyCurves aggregation — flag any state whose `bend_uu` deviates from the pilot
   median by more than ~90°. Also catches the isolated bad frames: err 165°, nfat 170°, runx 146°,
   foxa 144°, dux4 168°.
2. Inspect irf's register directly before using its bend anywhere.
3. **Recompute and report the headline without irf.** It is inside the LOO range so the sign holds, but
   the axis currently contains one pilot whose geometry is known broken in two channels and assumed
   intact in the third. State the without-irf number in the talk; it is not a footnote.
4. Re-run step 2's per-pilot slopes after cleaning. The isolated 144–170° frames are exactly the
   high-deviation points that dominate a per-state slope.

### 3.2 Interface restriction (§4b from v2, still unrun)

Per-base-pair bend variance along the duplex; recompute IQR using only base pairs inside the footprint;
re-run the correlation. If the variance lives in terminal base pairs, bend-IQR is measuring end-fraying.

Still worth doing, and now doubly so because **k = ∞ will make end-fraying worse**, so the sweep needs
an interface-restricted metric available.

| result | reading |
|---|---|
| ρ holds or strengthens | Real and mechanistically located. |
| ρ collapses | The axis was end-fraying. Claim 3 comes out of the talk. |
| Splits by pilot | Check against duplex length and footprint fraction before interpreting. |

### 3.3 HSF1: three anomalies in one pilot

- Worst state-vs-crystal agreement (0.27 baseline, 0.61 augmented) despite the largest ensemble (96)
- The sole pilot where across-conformer SD *rises* under augmentation (0.386 → 0.519), 1 of 12
- Its crystal (5D5U) is absent from the upstream DeepPBS dataset and was featurized fresh

Given the irf precedent — a parse failure undetected in one channel while another channel in the same
pilot had already been dropped — **run the range assertion and a geometry sanity check on HSF before
treating any of this as biology.** If it is clean, HSF's variance exception needs an explanation and is
worth a slide; if it is not, it needs excluding.

### 3.4 The unexplained +0.59

Per-pilot within-pilot slope vs bend IQR gives ρ = +0.59 (relaxed, p = 0.074) — high-IQR pilots showing
*less* negative slopes. That is backwards under every reading, including the shared-label one. Candidate
explanations: range restriction, attenuation, or the §3.1 parse failures contaminating the high-bend
tail. **Re-check after §3.1 cleaning.** It is the one number in the analysis nobody has a story for, and
an audience member who finds it will ask.

### 3.5 Family structure in the cross-pilot deltas

132 ordered donor–target pairs are already computed. Group Δ(X, Y) by whether family(X) = family(Y).

If the +0.025 concentrates in same-family pairs, that is the Ch. 16 own-family story reproduced in a
design with **no leakage at all** — a genuinely strong result, and a groupby on data already on disk.
If it is flat across families, the effect is generic and less interesting. **An hour.**

### 3.6 Still open from v2

- **§0.1 — Ch. 19½ provenance.** Unresolved. The `analysis/analyses/nmr_rmsd/` results are not on the
  cluster. Until located, Ch. 19½ and §7 (apo-fnat × ΔP) cannot be scheduled, and nothing from that
  chapter should go in the talk.
- **§11 — chapter-map reconciliation.** Two outlines with incompatible numbering (`book/OUTLINE.md`,
  21 chapters; `docs/TF_conformation_book_outline.md`, part-based). Fix before reconciling any numbers.
  Include the `froz_bend_iqr` → `bend_iqr_frozen` rename; the arm-ambiguity that caused C2 is likely in
  other columns too.
- **§1 — `baseline × relaxed` condition.** Still absent. Lower priority now: three independent lines say
  k = 1.5 changes little, so this mostly measures how little. Folds naturally into the k sweep as level 2.
- **Multiple testing.** BH over four cells is not the honest family. The real family is every axis tested
  across §17.2–§17.3 — free-state spread, reachability, static crystal bend, interface MGW-FL, bend
  median, bend IQR, in two arms — roughly ten. Bend IQR probably still survives; q = 0.0125 is not the
  number to quote.

---

## Part 4 — Order of work

| priority | item | cost | why here |
|---|---|---|---|
| 1 | §3.1 irf + range assertion; re-run slopes and headline without irf | hours | correctness; contaminates everything downstream |
| 2 | §2.2 crystal-free PWM dispersion | 0.5 d | determines whether the k sweep is affordable at four levels |
| 3 | §3.5 family structure in cross-pilot deltas | 1 h | free; may be a strong result |
| 4 | §1.1 k-sweep calibration on 2 pilots | 1 d | decides whether the sweep is worth running at all |
| 5 | §2.1 size-matched donor control | 1 training arm | the last alternative explanation for +0.025 |
| 6 | §3.2 interface restriction | 0.5 d | needed before k = ∞ is interpretable |
| 7 | §3.3 HSF geometry check | hours | |
| 8 | Full k sweep, if calibration passes | weeks | the intervention |
| 9 | Held-out pilots, pre-registered sign prediction | weeks | see Part 5 |
| 10 | Sequence-design variant | months | the test that would settle it |

Items 1–3 fit in a day and a half. Items 4–7 fit in the rest of the week.

---

## Part 5 — On more pilots

**The argument is selection, not power.** At n = 11 with ρ = −0.78 you have ~79% power; detection was
never the issue. The issue is that bend IQR was chosen as the largest of roughly ten axes tested, at
n = 11. An effect selected that way is biased upward essentially by construction, and **no re-analysis
of the same eleven pilots can correct a selection bias in those eleven.**

So the target is not a larger regression. It is a **held-out set with a pre-registered prediction**:
fit nothing on them, predict the sign of ΔP from bend IQR, and see. Ten to fifteen genuinely new pilots
would do it, and that is worth more than another twenty points added to the existing fit.

Power table, for the record (Fisher-z, partial correlation with one covariate):

| n | power at true ρ = 0.78 | power at true ρ = 0.45 |
|---|---|---|
| 11 | 0.79 | 0.25 |
| 20 | 0.99 | 0.49 |
| 30 | 1.00 | 0.70 |
| 40 | 1.00 | 0.83 |

ρ ≈ 0.45 is the prudent assumption post-selection, which sets the target near n ≈ 40 **if** the goal is
a powered regression. If the goal is a held-out sign test — which is the better goal — 10–15 suffices.

Candidate pool: the existing 130-entry DeepPBS benchmark. Per-pilot cost depends entirely on §2.2: with
an inference-only endpoint, new pilots need BioEmu + Stages 1–4 + a forward pass and no training arms.

---

## Part 6 — The talk

**Reordered.** Lead with what is measured cleanly and large-n; finish with the mechanism as motivation
for the intervention rather than as a finding.

**Claim 1 — a crystal is a point estimate of a distribution.** Conformation SD 0.235 vs seed SD 0.16,
n = 701 states, baseline checkpoints. Independent of everything contested.

**Claim 2 — shared-label augmentation teaches conformational invariance.** Across-conformer SD
0.235 → 0.180, 11 of 12 pilots, p = 0.027, while accuracy on unseen complexes improves by only +0.025.
The model maps ~90 structurally distinct conformations — *not* near-copies, median feature distance 0.81
relative to unrelated entries — onto one motif. That is the training objective working as specified, and
it is why augmentation cannot recover conformational signal. **Subject to §2.1.**

**Claim 3 — the augmentation gain decays with distance from the training frames.** +0.246 own complex →
+0.03–0.11 own family → +0.025 unseen complex → ≈0 cross-benchmark. Three points, three designs, one
monotone gradient. Volunteer the in-sample structure rather than letting someone find it.

**Claim 4 — the mechanism, stated as an observation, not a result.**

> Ensemble DNA-bend fluctuation predicts the augmentation effect across pilots (ρ = −0.78 arm-matched,
> surviving every single-covariate control). It does **not** operate through individual states
> (within-pilot slopes 7/11, n.s.), which is what a shared-label mechanism predicts and is a constraint
> on the claim rather than support for it. n = 11, selected from ~10 axes, so the effect size is an
> upper bound. The k sweep is the test.

**Present step 2's null yourself.** It is the result an audience would try to find, and finding it first
is the single most credibility-buying move available.

**Say out loud:** the binding sample size is the pilot count. Relaxed-arm eval coverage is 12;
DNA-shape-analysis coverage is 7; confirm which applies before quoting either.

---

## References

- Mitra, R. *et al.* **Geometric deep learning of protein–DNA binding specificity.** *Nature Methods*
  **21**, 1674–1683 (2024). doi:10.1038/s41592-024-02372-w — RI definition and within-structure
  normalization; RI vs alanine-scanning ΔΔG PCC 0.605 (Fig. 4h).
- Jiang, Y., Shewchuk, A. M., Chiu, T.-P., Li, J., Kribelbauer-Swietek, J. F., Gompel, N. & Rohs, R.
  **Readout of intrinsic and induced DNA shape by homeodomain transcription factor complexes.**
  *Biophysical Journal* (2026). doi:10.1016/j.bpj.2026.03.036 — *"minor groove width reflects intrinsic
  DNA structure, whereas minor groove width fluctuations capture protein-induced stabilization and
  reshaping of the DNA."* The per-position profile format for §3.2, and the independent convergence with
  the fluctuation axis.
- Lewis, S. *et al.* **Scalable emulation of protein equilibrium ensembles with generative deep
  learning.** *Science* (2025). doi:10.1126/science.adv9817 — *"A hallmark of this may be BioEmu's
  preference to predict holo over apo structures in the cryptic pocket benchmark"*; 85% holo vs 49% apo.
  Relevant only if §3.6's Ch. 19½ provenance question resolves.
