# TF-conformation — state of findings, v7

**2026-09-29.** Supersedes v1–v6. Those were plans; this is mostly results, because the plan resolved.

**What happened today.** Roughly a dozen analyses ran, all from data already on disk apart from one
cross-pilot inference sweep. Four claims survived and were strengthened, one was retracted, three
earlier results were withdrawn, and two pipeline failures were found. The project went in with a
mechanism story and came out with a cleaner set of findings and a sharper negative.

**The one-line version.** Augmentation's effect on DeepPBS is real but small and local; the training
objective teaches conformational invariance rather than conformational sensitivity; and no structural
property we can measure explains how much any given protein benefits.

---

## Part 1 — Findings

### F1. A crystal is a point estimate of a distribution

Across 701 fnat-gated conformations, the standard deviation of predicted-PWM agreement **across
conformations is 0.235** against **0.16 across training seeds**. The variability is a property of the
structures, not of training noise.

Baseline checkpoints only, so unaffected by the in-sample issue in F2. The foundational observation of
the project and the only one that needs no qualification.

### F2. Shared-label augmentation teaches conformational invariance

All 707 states are in-sample for their own pilot's augmented checkpoint (`fracF = 1.000`, verified per
entry) and share **one PWM label per pilot**. Under those checkpoints, state-level agreement with the
crystal rises from 0.504 to 0.746 — +0.242, 11 of 12 pilots.

The inputs are **not** near-copies: median feature distance from a state to its own crystal is **0.81**
relative to unrelated training entries, and ERR (1.07) and NFAT (1.34) are *further* from their own
crystal than unrelated entries are. So this is **label-side memorization with structurally distinct
inputs** — the model learning to map ~90 different conformations onto one motif.

The variance result is the cleaner statement: **across-conformer SD falls from 0.235 to 0.180 (frozen)
and 0.172 (relaxed), in 11 of 12 pilots** (Wilcoxon p = 0.027, 0.021). HSF1 is the sole exception
(0.386 → 0.519).

Augmentation buys consistency across conformers. That is the training objective working as specified,
and it is why augmentation cannot recover conformational signal.

*Contingent on §3.2 (size-matched donor) for the accuracy component only; the invariance result stands
regardless.*

### F3. The augmentation gain decays with distance from the training frames

| test structure's relation to training | Δ*r* | design |
|---|---|---|
| the pilot's own complex (~90 frames in training) | **+0.246** | in-sample, verified per entry |
| genuinely unseen same-family entries | **+0.03 to +0.11** | own-family benchmark; audited clean (0 of 11 pilots' reference PDBs appear in their own-family subsets) |
| unseen complex (cross-pilot donor) | **+0.025** | 12 × 11 ordered pairs × 2 arms × 5 seeds; Wilcoxon vs 0, p = 3.4e-05, positive in 10 of 12 |
| cross-benchmark, other families | **≈ 0** | Ch. 16 |

Three independent designs, one monotone gradient. **+0.025 is the defensible effect size**, and it
agrees with the own-family range measured by a completely different route.

The crystal-readout gain is **indistinguishable from zero on unseen complexes** (−0.008; Wilcoxon
p = 0.24), which retracts the +0.28 reported earlier as clean.

### F4. Intrinsic and protein-induced DNA fluctuation are decoupled

Deep DNAshape predicted fluctuation (`-FL`) versus pyCurves measured ensemble fluctuation:

- per duplex: **ρ = +0.030** (p = 0.93); interface-restricted +0.200 (p = 0.58)
- per position, unpadded windows: ProT **+0.224** (uncorrected p = 0.049, BH q = 0.295, mixed-model
  p = 0.149); **Roll and Twist flat** (|ρ| ≤ 0.06)
- not explained by duplex length (LRT for length terms χ² = 1.92, df = 2, p = 0.383)

**The positive control passes**, which is what makes this a finding rather than a plumbing failure: the
*same* pipeline recovers predicted-vs-measured **mean** MGW at **ρ = +0.386** (p = 3e-4, crystal) and
**+0.384** (p = 5e-4, relaxed) on unpadded positions, with within-pilot correlations of +0.17 to +0.30 —
so the predictor resolves which position has a wide groove, not just pilot-level offsets. The crystal
being the strongest cell is the right ordering, since it is the least protein-perturbed structure.

This is consistent with the lab's published distinction — groove width reports intrinsic DNA structure
while groove-width *fluctuation* reports protein-induced stabilization and reshaping. **Sequence-based
shape prediction cannot substitute for ensemble measurement when the question is about a bound
complex.**

### F5 (negative). No measured structural property predicts the augmentation effect

Tested against own-family ΔP at n = 10–11 pilots: protein reachability, free-state spread, static
crystal bend, interface MGW-FL, whole-duplex bend median, whole-duplex bend IQR (both arms),
length-free windowed bend IQR (4 window widths × 2 arms), `r_mean`, surviving-state count, duplex
length. **None reaches q < 0.05.**

This is the answer to the question the project has been circling since Ch. 16. The augmentation effect
varies across pilots, and the variation is not structural in any way currently detectable.

It also coheres with F2: if augmentation's main action is teaching conformational invariance, its *size*
should be governed by training dynamics — label consistency, ensemble breadth, family coverage — rather
than by the geometry of the bound complex.

**The honest caveat, which travels with it.** ΔP is itself noisy: seed SEM is **two-thirds** of the
between-pilot spread (0.644 frozen, 0.661 relaxed), and entry expansion cannot fix it (seed variance
dominates in 7 of 10 pilots; infinite eval entries move interaction power 0.20 → 0.24). So the claim is
*no measured structural property explains the pilot-to-pilot variation at the precision available* —
weaker than "none exists," and still worth stating. The persuasive force is in the enumeration.

---

## Part 2 — Retracted, and why

### R1. The mechanism claim (bend fluctuation → augmentation effect)

**ρ = −0.745 arm-matched, uncorrected p = 0.013. Retracted.**

Retired by the length-free axis. Whole-duplex `bend_uu` IQR is **69% duplex length** (ρ = −0.827), and
the length effect runs **opposite to worm-like-chain expectation** — WLC says angular variance grows
with contour length. That is the signature of a terminal-axis-vector angle, not a polymer property.

Recomputing on a twist-aware fixed-width central window:

- **length dependence is gone**: |ρ| ≤ 0.28 across all window widths, all p > 0.4
- **the ΔP correlation collapses**: 8 of 8 window × arm cells n.s., minimum p = 0.117
- **the sign reverses**: 7 of 8 cells positive, the eighth at −0.067

Three diagnostics confirm it was not a power failure:

1. **The axes are unrelated.** Window IQR vs `bend_uu` IQR: ρ = −0.055 (p = 0.88). Not two noisy
   measures of one quantity.
2. **R² does not decompose.** Old axis alone 0.556; length alone 0.203; length-free alone 0.050; both
   together 0.228. There is no "real signal diluted by length" story available.
3. **The measurement is sound.** Arm-to-arm reproducibility +0.90 to +0.944, *better* than `bend_uu`
   managed (+0.927). A cleaner instrument failing to see what a dirtier one saw.

**Do not report the +0.527 best cell.** It is the best of eight and would re-enter the same selection
trap that produced R1 in the first place.

**One loose end.** The pre-check found IQR × √L at ρ = −0.867, stronger than the raw axis — unexplained,
and most likely a two-parameter fit finding a combination at n = 10. A permutation test on ΔP labels
across the full set of axis transforms would settle whether it exceeds what the search space produces by
chance. Ten minutes, and it is the one number pointing the other way.

### R2. The +0.28 held-out crystal gain

Reported as "the result worth quoting." Cross-pilot gives **−0.008**, Wilcoxon vs 0 p = 0.24, 4 of 12
targets positive. The state-level number was inflated *and flagged*; the crystal number was inflated and
read as trustworthy, which made it the more dangerous of the two.

### R3. The information-content deficit as a biological claim

The baseline's −1.14-bit deficit reverses under augmentation (7.50 → 8.67 frozen, 8.80 relaxed). The
flattening is a property of the unaugmented model, not of the conformations. The SSO7D nonspecific-
binding calibration, which supported the biological reading, goes with it.

### Also withdrawn

"The mechanism axis strengthens in the relaxed arm" (the axis arm is irrelevant; only the effect arm
moves). The pooled-arm correlation as an improvement (−0.636, CI [−0.987, +0.059]; the arms' errors are
correlated at +0.770 against a 0.701 independent-error ceiling). The HSF crystal-leakage control (HSF is
5th worst, not worst). The frame-count bias argument for reachability (~0.003 Å against a 0.103 Å gap).

---

## Part 3 — Remaining work

Three items. Everything else is in Part 5.

### §3.1 — Pipeline integrity audit (0.5 d) — **first**

Two failures found today were invisible for months and were both caught by **a value looking
physically impossible**, not by any check:

- **irf had no duplex.** 0 of 184 states carry base-pair tables, against 2,224 of 2,224 everywhere else,
  while its crystal parses fine. The 170° bend was computed on a structure pyCurves never recognised as
  double-stranded. Downstream code simply returned fewer rows.
- **The naive window-bend sum** gave TBP 118.8° → 191.8° over 8 steps, because successive base-pair
  steps are rotated by helical twist and wedge angles at opposite phases must partially cancel. Correct
  accumulation carries cumulative twist Ω: `Bx += roll·cosΩ − tilt·sinΩ`, `By += roll·sinΩ + tilt·cosΩ`,
  `bend = √(Bx²+By²)`.

**Build physical assertions, not statistical ones.** Bend in [0, 180]. Every state yields a
`local_base_base` table and the count equals the state count. Row counts match stage to stage. The
Stage-4 helix-guard silent skip is a known instance of the same pattern.

Separately: determine whether irf's duplex genuinely melted under the restraints. If it did, that is a
restraint-protocol finding and duplex integrity belongs in §3.3's QC list.

### §3.2 — Size-matched donor control (1 training arm) — **highest stakes**

**The only remaining live threat to a surviving claim.** Donor Y's augmented checkpoint differs from Y's
baseline in two ways: it saw BioEmu frames **and** ~90 extra training examples. Some of the +0.025 could
be data volume.

Control: Y's baseline plus ~90 additional **crystal** entries from the benchmark, no BioEmu; score
targets X ≠ Y as in the cross-pilot design. One arm on one or two donors.

| result | reading |
|---|---|
| ≈ 0 | +0.025 is conformational information. F3's gradient and F2's accuracy half both survive intact. |
| ≈ +0.025 | **Augmentation is just more data.** F2 becomes purely about invariance, with no accuracy benefit. F3's gradient becomes a statement about memorization alone. |
| > +0.025 | Conformations are *worse* than crystal entries at matched volume. Sharpens F2 considerably. |

### §3.3 — Two small, bounded experiments

**Sequence viability gate.** Build 4–6 DNA sequence variants on **one pilot (ets1)** with the existing
stage-3 relax machinery, and ask whether **measured** fluctuation spans a usable range **at fixed
protein**. The predictor is now a design aid, not the x-variable — F4 established that predicted
intrinsic fluctuation does not track measured induced fluctuation — so the question is whether sequence
causes enough variation in the induced quantity to support an experiment at all. The 12 pilots span
3.88× in measured MGW-FL, but they differ in protein *and* DNA, so that range may be entirely
protein-side.

**Restraint calibration (descoped k sweep).** Two pilots, levels k ∈ {0, 1.5, 5.0, ∞}, measuring whether
restraint strength changes anything. This is a within-pilot geometry measurement on ~90 states, so it is
adequately powered — unlike the 12-point mechanism test it was originally designed to support. **Report
as "does restraint strength change what the model reads," with no mechanism framing**, and do not let
the mechanism claim re-enter through it. Include duplex integrity in the per-level QC, since restraints
are getting looser and irf may show what that does. Three independent measurements already say k = 1.5
changes little (axis ρ = +0.927, mean IQR 5.46° → 5.79°, readout Δr = +0.011, p = 0.36), so a null here
is a real possibility and worth knowing cheaply.

---

## Part 4 — The talk

**Four claims, all large-n or multiply-audited, plus one negative.**

1. **A crystal is a point estimate of a distribution.** Conformation SD 0.235 vs seed SD 0.16, n = 701.
2. **Shared-label augmentation teaches conformational invariance.** SD 0.235 → 0.180, 11 of 12 pilots,
   p = 0.027, while accuracy on unseen complexes rises only +0.025. ~90 structurally distinct
   conformations — not near-copies, median feature distance 0.81 — mapped onto one motif.
3. **The gain decays with distance from the training frames.** +0.246 → +0.03–0.11 → +0.025 → ≈0. Three
   designs, one gradient, all three points audited.
4. **Intrinsic and protein-induced DNA fluctuation are decoupled.** With the mean-MGW positive control
   (ρ = +0.386) showing the pipeline works.
5. **(Negative) No measured structural property predicts the augmentation effect's size.** Enumerate all
   nine axes. The force is in the list.

**Volunteer, do not defend.** The within-pilot null (7/11 slopes, p = 0.55). The failed multiple-testing
correction. The in-sample structure of the augmented arms. The duplex-length confound and the resulting
retraction. The power ceiling. Each is what an audience would go looking for, and finding them first is
the most credibility-buying move available.

**Say out loud.** The binding sample size is the pilot count. The pilot set is one-per-family across
twelve families, so no pilot-pair design here can test family transfer. ΔP's seed noise is two-thirds of
the between-pilot spread, so F5 is bounded by precision as well as by what was measured.

**On how the week reads.** A dozen analyses yielding four surviving claims, one retraction, three
withdrawals and two pipeline bugs looks worse than it is. Every negative was found internally. The
retraction in particular — a mechanism claim that turned out to be duplex length — was caught by a
fixed-width-window recomputation that had been deferred four times and run before any training arms were
spent on it. That is the sequence working.

---

## Part 5 — Not scheduled

**Cheap and deferred every round:**

- **Ch. 19½ provenance.** `analysis/analyses/nmr_rmsd/` results are not on the cluster; the chapter's
  `%across` table, three-way fnat panel, informativeness gate and runx reclassification have no
  traceable source. One hour. Until resolved, Ch. 19½ stays out of the talk.
- **Chapter-map reconciliation.** Two outlines with incompatible numbering (`book/OUTLINE.md`, 21
  chapters; `docs/TF_conformation_book_outline.md`, part-based). Every number reconciliation will
  generate new conflicts until this is fixed. Include the `froz_bend_iqr` → `bend_iqr_frozen` rename.

**Lost their frame with R1** (both were anomalies *relative to* the retracted axis):

- HSF's three anomalies — worst agreement (0.27 → 0.61), sole pilot where across-conformer SD rises
  under augmentation. Cleared of parse failure (11.0° deviation). The variance exception is still worth
  a slide as an F2 exception.
- The +0.59 (within-pilot slope vs bend IQR).

**Methods findings worth exporting rather than filing:**

- **Deep DNAshape N-padding biases `-FL` upward** (ρ = +0.535 vs N-count, p = 3e-16). A padded position's
  `-FL` expresses uncertainty about missing context, not a property of the DNA. Layer 4 costs 8 bp of
  every duplex, so usable positions are (L−8)/L — 27% to 62% across these pilots, 73% even at 30 bp. A
  real trap for anyone using the tool on short oligos; tell Jinsen and Tsu-Pei.
- The first-pass ρ = −0.316 at p = 3e-5 **with the wrong sign** is a clean teaching example: two
  unrelated position gradients with opposite signs manufacturing significance. The sign was the tell.
- A **layer sweep (1–4)**, not layer 2 alone — lower layers propagate less context, so layer-2 `-FL`
  approaches a pentamer lookup and a single alternative level is ambiguous between "coverage helped" and
  "long-range context was noise."

**Dropped with reasons:**

- Held-out sign test on 10–15 pilots (~26% power at true ρ = 0.5; needs 30–40 to be confirmatory).
- Entry expansion as a power fix (0.20 → 0.24 at the floor; hsf only, for its own sake).
- Weighting or filtering training states by predicted-PWM agreement (selects frames the model already
  agrees with; guarantees regression to baseline).
- Crystal-free PWM dispersion as an endpoint (null against the axis and ΔP; −0.982 redundant with
  `r_mean`).
- `r_mean` as an endpoint (frozen-arm only, p = 0.098).

---

## References

- Mitra, R., Li, J., Sagendorf, J. M., Jiang, Y., Cohen, A. S., Chiu, T.-P., Glasscock, C. J. & Rohs, R.
  **Geometric deep learning of protein–DNA binding specificity.** *Nature Methods* **21**, 1674–1683
  (2024). doi:10.1038/s41592-024-02372-w
- Jiang, Y., Shewchuk, A. M., Chiu, T.-P., Li, J., Kribelbauer-Swietek, J. F., Gompel, N. & Rohs, R.
  **Readout of intrinsic and induced DNA shape by homeodomain transcription factor complexes.**
  *Biophysical Journal* (2026). doi:10.1016/j.bpj.2026.03.036 — *"minor groove width reflects intrinsic
  DNA structure, whereas minor groove width fluctuations capture protein-induced stabilization and
  reshaping of the DNA"*; *"recognition involves a context-dependent balance between conformational
  selection and induced fit."* The basis for F4.
- Li, J., Chiu, T.-P. & Rohs, R. **Predicting DNA structure using a deep learning method.** *Nature
  Communications* **15**, 1243 (2024). Deep DNAshape.
- Li, J. **pyCurves.** https://github.com/JinsenLi/pyCurves — Curves/Curves+ reimplementation; the
  per-level tables behind the length-free axis.
- Lewis, S. *et al.* **Scalable emulation of protein equilibrium ensembles with generative deep
  learning.** *Science* (2025). doi:10.1126/science.adv9817 — BioEmu; *"a hallmark of this may be
  BioEmu's preference to predict holo over apo structures in the cryptic pocket benchmark."*
