# From augmentation to mechanism — a design memo

**Reframe:** stop treating DeepPBS augmentation as an end in itself and treat the *sign* of the
augmentation effect as a **functional probe of recognition mechanism** across TF families.
DeepPBS is a fixed-DNA, protein-structure→PWM model, so it is constitutively blind to DNA
deformation and to induced-fit protein rearrangement. That blindness is what makes it a clean
instrument: the sign of "does a free-state protein ensemble improve prediction" isolates whether
*protein-side* free-state conformers are binding-relevant.

---

## 1 · Apo/holo test — result

The BioEmu ensemble is a **computational apo (free) state**; the crystal is **holo (bound)**.
Two independent axes per pilot (`analysis/data/mechanism_apo_holo.csv`):

| pilot | free-state spread (Å) | reachability d_min (Å) | DNA deformation | aug Δr |
|---|---|---|---|---|
| ETS1 | 1.75 | 0.87 | minimal | **+0.15** |
| TBP | 1.57 | **0.59** (best) | extreme ~80° kink | **−0.26** |
| EGR1 | 2.90 | 0.80 | minimal | +0.10 |
| engrailed | 2.84 | 0.69 | modest | −0.11 |
| FOXA | 5.92 | 1.79 | modest | +0.01 |
| LEF1 | 6.36 | 2.43 | severe ~110° bend | −0.04 |

**No single protein-side axis predicts the sign** (Spearman |ρ|<0.5, all n.s., n=6). The decisive
case is TBP: it has the *best* protein reachability (0.59 Å — its rigid β-saddle reaches the bound
backbone perfectly) yet is hurt most. The protein ensemble is not the problem; the missing
variable is the DNA.

**Reading (M1, panel B):** augmentation helps only in the corner where the protein reaches its
bound pose **and** recognition needs little DNA deformation (ETS1, EGR1). This is the
**conformational-selection** regime — the free protein already samples binding-competent geometry.
Augmentation hurts when recognition is **induced fit on the DNA** (TBP kink, LEF1 bend): free-protein
conformers carry no information about the DNA distortion the model holds fixed, so they add noise.

*Caveats:* n=6, single-seed, correlations not significant; the DNA-deformation axis is qualitative
(from literature, not yet computed). This is a hypothesis with a clean mechanistic story and a
decisive negative control (TBP), not a settled result.

---

## 2 · Incorporating DNA relaxation into the pipeline

The apo/holo result says the pipeline's **rigid-DNA assumption is the ceiling**, not the protein
sampling. Three tiers, increasing cost:

**Tier 1 — relax DNA in the existing minimization (cheapest, closes confound #2-adjacent).**
Stage 3 already runs OpenMM with the DNA restrained to crystal. Release the DNA (or ramp its
restraints down) so protein and DNA co-relax into each docked frame. Gives per-frame DNA geometry
essentially for free. Risk: unrestrained B-DNA drifts/frays under GBSA on this timescale — needs a
**sequence-dependent stiffness prior** to stay physical.

**Tier 2 — cgDNA+/cgNA+ stiffness prior (the hexABC connection).** Don't use hexABC trajectories
directly (they characterize *free* B-DNA; recognition is *protein-induced* deformation). Instead use
the distilled **cgDNA+ ground-state + stiffness matrix** as a harmonic restraint on the DNA
during co-relaxation: the DNA is free to deform toward the protein-preferred geometry but pays a
sequence-dependent elastic cost. This is the principled way to let TBP-like kinks form *if the
protein contacts drive them*, without letting the duplex melt. Deliverable: a per-hexamer stiffness
lookup wired into the Stage-3 energy function.

**Tier 3 — learned DNA deformability (Deep DNAshape) as the prior.** The in-house Deep DNAshape
model predicts sequence→shape (and shape flexibility) directly. Lowest-friction go/no-go: use its
predicted shape as the DNA relaxation target / restraint center, so the augmented frames carry
DNA geometry consistent with the sequence's intrinsic preferences.

**Falsifiable prediction that ties it back to the mechanism result:** a flexible-DNA pipeline should
*recover exactly the families augmentation hurts today* — TBP, LEF1, homeodomain — because those
are the induced-fit-on-DNA cases. If co-relaxing DNA flips TBP from −0.26 toward neutral/positive,
the mechanism hypothesis is confirmed and the method is fixed by the same change.

---

## 3 · Per-frame PWM labels — the decisive control (reiterated)

**The confound (still #1).** Every augmented frame of a complex currently shares **one** crystal-derived
PWM label. So "does the ensemble help" is entangled with "does averaging N frames onto one label just
regularize the model." We cannot tell signal from smoothing.

**The idea.** Give each frame its **own** label instead of the shared crystal PWM:
- *Structure-derived per-frame PWM* — recompute the expected specificity from each frame's own
  protein–DNA contacts (a biophysical readout, e.g. per-frame interface energy → position weights),
  so a frame that has drifted predicts a *different* motif.
- *Even predicted per-frame labels* suffice for the control: the question is whether frame-to-frame
  **label diversity** carries information.

**What it decides.**
- If per-frame labels **help** where shared labels didn't → the ensemble encodes real
  conformation→specificity signal (conformational selection, confirmed) — the strong result.
- If per-frame labels **wash out** the gain → the shared-label benefit was regularization/noise
  averaging, and ETS's edge is mostly its low baseline having headroom. Still publishable, opposite
  conclusion.

**Why it's decisive for the ETS finding specifically.** ETS has the lowest baseline of the helped
families, so "low baseline = most improvable by any regularizer" is the live alternative to
"ETS is genuinely conformationally selected." Per-frame labels separate those two — it is the
experiment that turns the family map from suggestive into mechanistic.

**Sequencing:** per-frame PWM is a bigger build (new labeling stage) but answers the sharper
question; DNA relaxation Tier 1 is cheaper and independently valuable. Run the 5-seed control first
(already in flight) to put error bars on the family effects, then per-frame PWM, then DNA relaxation.

*Figures:* M1 apo/holo mechanism ({artifact:ecf5118a-eb7a-40c2-8a33-6e69268bfcf6}), P3 family effect ({artifact:97e2ae62-5d64-49bb-80d9-e2ca39fcd5ac}).
