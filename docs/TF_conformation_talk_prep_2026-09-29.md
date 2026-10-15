# Talk preparation plan — 2026-09-29

**Horizon:** a few days. **Goal:** present what holds, retract what doesn't, and leave with a decision
made about whether to commit 24 training arms.

**The operating rule for this document: stop adding analyses.** Twelve ran today. Four claims hold, one
is retracted, and one conflict is unresolved. That is a talk. Three things still need running — two
because they change what you can say, one because it is a correctness issue — and everything else
freezes.

---

## Part 1 — Run these three

### §1.1 — Terminal-window bend (0.5 d) ⭐ **highest value per hour in the project**

Same twist-aware accumulation as the central-window axis, applied to the duplex **ends** instead. It
decides whether the −0.867 whole-duplex correlation is interpretable, which in turn decides whether 12
new pilots are worth 24 training arms. Half a day gating a multi-week commitment.

**Compute all three combination variants in the same pass** — how the two ends are combined is a choice
that can manufacture or destroy signal:

1. **Each end separately** — two values per state, each correlated against ΔP. Both ends reproducing
   −0.87 is convergent evidence; only one doing so means something asymmetric is happening and needs its
   own explanation.
2. **Mean of the two ends** — the natural summary.
3. **Vector sum across the whole duplex excluding the centre** — *the real comparison.* `bend_uu`
   measures the angle between terminal axis **vectors**, a relative orientation across the whole
   molecule, not a sum of local bending at the ends. It is entirely possible that neither end window
   reproduces −0.87 individually while the end-to-end relative orientation does — which would mean the
   signal lives in the global relative orientation of the two ends, a third physical quantity again.

**Then one cheap follow-up, conditional on the ends carrying signal:** correlate terminal-window IQR
against **interface contact count** and **footprint width**, both already on disk. If terminal mobility
tracks how much of the duplex the protein grips, the axis becomes a protein-side property measured
through the DNA — interpretable, and the branch where new pilots are clearly worth it.

| result | reading | consequence |
|---|---|---|
| Ends reproduce −0.87, and it tracks fraying / restraint edge | The axis is a **protocol property**, not biology. | Mechanism line closes. **No new pilots.** Feeds the restraint calibration directly. |
| Ends reproduce −0.87, and it tracks footprint / interface extent | Interpretable, protein-side. | **Replication at n ≈ 12 is worth 24 arms.** Pre-register and run after the talk. |
| Neither end window nor the end-to-end orientation reproduces it | The signal is in neither the centre nor the ends as constructed. | Present as unresolved; **do not spend.** |
| Still ambiguous | — | Present the conflict; defer the decision. |

**Do not use** the "end contribution" construct (whole-duplex minus rescaled central). It is near-circular
— the two axes are orthogonal, so subtracting one mostly returns the other, and it correlates with
whole-duplex IQR at +0.685. Not independent evidence.

### §1.2 — Pipeline integrity audit (0.5 d)

A correctness issue, and the week produced two reasons to do it now. irf's 184 states with no base-paired
DNA, and the naive window-bend sum giving TBP 191.8° over 8 steps. **Both were caught by a value looking
physically impossible, not by any check.** Neither would have been caught by a statistical assertion.

Build **physical** assertions:

- bend ∈ [0, 180] for every state
- every state yields a `local_base_base` table, and the table count equals the state count
- row counts match stage to stage (the Stage-4 helix-guard silent skip is a known instance)
- duplex integrity per state — relevant to the restraint calibration, where restraints get looser

Separately determine whether irf's duplex genuinely melted under the restraints. If it did, that is a
restraint-protocol finding and belongs in the calibration's QC list.

### §1.3 — Size-matched donor control (1 training arm) — **the only live threat to a surviving claim**

Donor Y's augmented checkpoint saw BioEmu frames **and** ~90 extra training examples. Some of the +0.025
could be volume rather than conformational information.

Control: Y's baseline plus ~90 additional **crystal** entries from the benchmark, no BioEmu; score
targets X ≠ Y as in the cross-pilot design. One arm on one or two donors.

| result | reading |
|---|---|
| ≈ 0 | +0.025 is conformational information. F2's accuracy half and F3's gradient both survive intact. |
| ≈ +0.025 | **Augmentation is just more data.** F2 becomes purely about invariance; F3 becomes a statement about memorization alone. |
| > +0.025 | Conformations are *worse* than crystal entries at matched volume. Sharpens F2. |

**If it does not finish in time**, present F2 as the invariance result and state the +0.025 with the
volume caveat attached. The invariance claim does not depend on it.

---

## Part 2 — Freeze list

Not before the talk. Each has a reason.

| item | why not now |
|---|---|
| Sequence viability gate (ets1 variants) | Scopes a future experiment; cannot change a current claim. §1.1 strictly dominates it for the same half-day. |
| Restraint calibration (descoped k sweep) | Worth a day, but after. Three measurements already say k = 1.5 changes little. |
| Deep DNAshape layer sweep (1–4) | Refines F4, which already holds with a passing positive control. |
| hsf entry expansion | Power fix that isn't one (0.20 → 0.24 at the floor). |
| New pilots | Gated on §1.1. Decide after the talk. |
| Any new axis transform | The permutation family is nine. Adding a tenth re-opens the selection problem the permutation test just closed. |

**Two exceptions, both cheap and both overdue since v2:**

- **Ch. 19½ provenance** (1 h). The `%across` table, three-way fnat panel, informativeness gate and runx
  reclassification have no traceable source on the cluster. **Until resolved, Ch. 19½ stays out of the
  talk** — so this only needs doing if you wanted to use any of it. Decide which, then act.
- **Chapter-map reconciliation** (1 h). Two outlines with incompatible numbering. Not needed for the
  talk; needed before any further writing. Include the `froz_bend_iqr` → `bend_iqr_frozen` rename.

---

## Part 3 — Figures

Most of this exists. Build only what's missing.

### Exists, use as-is

| figure | carries |
|---|---|
| `conf_vs_crystal` panel (a) per-conformation r + ensemble diamond | F1 |
| three-condition comparison, panels (a)–(e) | F2, with the in-sample caveat in the caption |
| leakage-gradient figure, panels (a)–(c) | F3 — this is the single best figure in the project |
| `M2_where_the_signal_lives` | background on where the effect is |
| exemplar motif logos (EGR1 worst, TBP best) | F1, makes the abstraction concrete |

### Build (3)

1. **The variance figure for F2** — across-conformer SD, baseline vs frozen vs relaxed, one point per
   pilot, HSF1 flagged as the sole exception. F2's headline is currently carried by a number in prose.
   **Priority 1.**
2. **The F5 enumeration** — every axis tested against ΔP, ρ with CI, ordered, q-threshold line drawn,
   the whole-duplex axis marked as the named exception. The force of F5 is the list; a figure makes the
   list legible in ten seconds. **Priority 1.**
3. **The two-axes conflict** — whole-duplex IQR × √L vs ΔP beside central-window IQR vs ΔP, with the
   orthogonality (ρ = −0.055) stated. This is the retraction figure and it should be shown, not
   described. **Priority 2**, and it updates with §1.1's result.

### Caption discipline

- Anything from the augmented arms carries the in-sample note.
- F4's figure carries the mean-MGW positive control (ρ = +0.386) in the same panel or caption — without
  it, the null looks like a broken pipeline.
- Do not show the +0.527 best central-window cell anywhere. It is the best of eight and re-enters the
  selection trap.

---

## Part 4 — Slide structure

**Arc:** here is what a single structure can't tell you → here is what training on ensembles actually
does → here is what we can't explain, and what we tried.

| # | slide | content |
|---|---|---|
| 1 | The question | A crystal is one conformation; binding specificity is a distribution. Can ensembles bridge that? |
| 2 | **F1** | Conformation SD 0.235 vs seed SD 0.16, n = 701. Variability is structural, not training noise. |
| 3 | **F2** | SD 0.235 → 0.180, 11/12 pilots, p = 0.027. ~90 structurally distinct conformations (feature distance 0.81 — *not* near-copies) mapped to one label. The objective working as specified. |
| 4 | **F3** | The decay gradient: +0.246 → +0.03–0.11 → +0.025 → ≈0. Three designs. **Volunteer the in-sample structure here**, before anyone asks. |
| 5 | **F4** | Intrinsic vs induced fluctuation decoupled — with the mean-MGW control on the same slide. |
| 6 | **The retraction** | The mechanism claim, what it rested on, and why it's withdrawn. Two orthogonal reproducible axes, opposite answers, unresolved at n = 10. |
| 7 | **F5** | The enumeration. Nine axes, none at q < 0.05 among length-free properties, whole-duplex named as the exception with its interpretability problem. |
| 8 | Limits | Pilot count is the binding n. Seed noise is ⅔ of between-pilot spread. One pilot per family across twelve families — no pilot-pair design here can test family transfer. |
| 9 | Next | Terminal-window result (by then in hand), and the conditional 12-pilot replication. |

**Slides 6 and 7 are the ones to rehearse.** They are the hardest to deliver and the most credible if
delivered well.

---

## Part 5 — Prepared answers

Write these out. Each is a question the talk invites.

**"Isn't the augmentation improvement just memorization?"**
Yes, substantially — and we measured it. All 707 states are in-sample, verified per entry. The clean
number is +0.025 on unseen complexes, from a cross-pilot design. The inputs are not near-copies (feature
distance 0.81), so it's label-side memorization with structurally distinct inputs, which is the more
interesting version.

**"Why did you retract the mechanism claim if p = 0.005?"**
Because we can't say what the quantity is. It's a terminal axis-vector angle, 69% duplex length, and a
length-free measure of central flexing shows nothing and is orthogonal to it. The permutation test rules
out a search artifact, so this is a correlation we can't interpret — not one we think is false. Terminal-
window analysis is the next step.

**"Why not just run more pilots?"**
For the central-window axis, ~85. For the whole-duplex axis, 10–12 would suffice at >90% power — which is
affordable, and is exactly what we'd do if the terminal-window analysis makes the quantity interpretable.
Spending 24 arms to pin down something uninterpretable is the wrong order.

**"Your pilots are all different families — isn't that a problem?"**
It was a deliberate coverage choice and it costs us family-transfer testability: all 132 cross-pilot pairs
are cross-family by construction, so no pilot-pair design here can test it.

**"Could the whole thing be a DeepPBS artifact rather than biology?"**
F1 and F2 are statements about the model, explicitly. F4 is a statement about two measurement methods.
Only the retracted claim was about biology, and it's retracted.

**"What would change your mind about augmentation being useful?"**
A per-frame labelling scheme. The current objective maps every conformation of a protein to one motif,
so it cannot learn conformational sensitivity by construction — that's F2. Shared labels are the
bottleneck, not the frames.

---

## Part 6 — Schedule

| day | work |
|---|---|
| 1 | §1.1 terminal-window, all three variants + the footprint follow-up. Launch §1.3. |
| 2 | §1.2 integrity audit. Build figures 1 and 2. |
| 3 | Figure 3 with §1.1's result. Assemble slides. |
| 4 | Rehearse 6 and 7. Write out Part 5. |

Leave day 4 whole. Slides 6 and 7 carry the talk's credibility and they're the ones that benefit from
being said out loud a few times.

---

## Part 7 — The decision deferred to after

Twelve new pilots, 24 training arms, candidates from the existing 130-entry benchmark (no system-
hunting), pre-registered sign and magnitude. A few weeks of mostly unattended compute, and it would
convert the whole-duplex correlation from a selected observation into a replicated one — which is the
only thing that addresses selection bias, and which no re-analysis of the existing ten can do.

**Gated entirely on §1.1.** Make the call after the talk, with the terminal-window answer in hand and
whatever the audience pushes on.
