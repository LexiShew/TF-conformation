# Presentation structure v4 — ~60 minutes

**Supersedes** v1–v3. Same frame and ordering as v3, expanded from 35 to ~50 minutes of slides plus
10–15 minutes of discussion. The extra time goes to three things: **context from the August deck** (so
the talk stands alone for anyone who wasn't there), **the intrinsic/induced DNA result**, which v3 had
to compress, and **a correctness slide** that earns the room's trust before the negatives land.

---

## Part 0 — The frame

**In August I had a hypothesis and a pipeline. I now have two measurements and a method.** Everything
else — including the retraction — is how I got there.

Three things make that honest rather than spin:

1. **The retraction came from my own stated next steps.** August's slide 105 listed "decompose bend-IQR"
   and "localize fluctuation." I did both, and they're what withdrew the claim.
2. **The nulls are mechanistic.** "Augmentation didn't help" is a null. "Augmentation teaches the model to
   ignore conformation, because every conformation carries one label" is a finding with a fix attached.
3. **The project gained a method it didn't have in August** — the apo-fnat read-out — and a published test
   case to validate it against.

**Title:** _"What conformational ensembles can and can't teach a specificity model."_

**Design rule for the whole deck: never end a slide on a withdrawal.** Every negative is immediately
followed by what it enables.

---

## Part 1 — The arc

| act | question                                                | slides  | min     |
| --- | ------------------------------------------------------- | ------- | ------- |
| 1   | Setup, and where we left off                            | S1–S6   | 8       |
| 2   | **Does the generator reach the bound-competent state?** | S7–S12  | 11      |
| 3   | **What does augmentation actually teach?**              | S13–S19 | 13      |
| 4   | **Can we predict who benefits?**                        | S20–S24 | 10      |
| 5   | What the project can do now                             | S25–S28 | 8       |
| 6   | Open questions                                          | S29     | 10–15   |
|     |                                                         |         | **~60** |

---

## Part 2 — Slide by slide

### Act 1 — Setup and where we left off (8 min)

**S1 — Title** [0.5]
Subtitle: "13 TFs · 701 conformations · 14 apo/holo NMR pairs."

**S2 — The biological question** [1.5]
A co-crystal is one conformation. Binding specificity is a distribution over sequences. DeepPBS reads the
first and predicts the second. Does the answer depend on which conformation you hand it?

Two objects side by side: one structure, one PWM logo.

**S3 — DeepPBS and BioEmu, in one slide each** [2] _(reuse old deck's intro slides)_
What DeepPBS does (structure → PWM, geometric deep learning on a symmetrized DNA helix), and what BioEmu
does (generative equilibrium ensembles from sequence). Thirty seconds each — enough that a newcomer can
follow, no more.

**S4 — The pipeline** [1.5] _(reuse old slide 28 or 40)_
BioEmu → dock → minimize → fnat gate → featurize → DeepPBS. One row of boxes with surviving-state counts.

_For anyone here in August: this is the 55 slides you already sat through, compressed. Ask me and I'll
open any of them._ That line gets a laugh and buys you the compression.

**S5 — August, in my own words** [1.5] _(reuse old slides 8 and 100 verbatim)_

> _"Refined hypothesis: augmentation helps under conformational selection of DNA and hurts under induced
> fit."_
> _"Fluctuating DNA = augmentation hurts. Rigid DNA = augmentation helps or neutral."_

And slide 105's stated next steps: decompose bend-IQR, localize the fluctuation.

**S6 — The verdict table** ⭐ [1]

| August claim                                                | how I tested it                                           | outcome                                             |
| ----------------------------------------------------------- | --------------------------------------------------------- | --------------------------------------------------- |
| Augmentation helps under selection, hurts under induced fit | measured mechanism directly from deposited apo structures | not supported — and the measurement is now a method |
| The effect tracks DNA bend fluctuation                      | slide 105's plan, executed                                | withdrawn — the quantity can't be identified        |
| Augmentation changes what the model reads                   | 701 conformations × three arms                            | **yes — and the direction is the interesting part** |

Put this up by minute eight. It converts the retraction from something that happens _to_ you into the
structure of the talk.

---

### Act 2 — Does the generator reach the bound-competent state? (11 min)

**S7 — The question nobody had asked** [1.5]
The pipeline pairs BioEmu protein frames with crystal DNA. Nobody had checked whether those frames
reconstitute the native interface at all.

**S8 — How you can ask it** [2]
Solution-NMR apo/holo pairs: the same protein deposited with and without DNA. Dock a free-state structure
onto the bound DNA, score fnat. **No DeepPBS, no training, no augmentation.**

Method detail worth two sentences: all RMSDs on one shared atom set — four-way sequence intersection,
expression tags stripped, disordered termini trimmed. Not cosmetic: one pair reads 9.47 Å untrimmed and
0.81 Å on its core.

**S9 — The informativeness gate** [1.5]
A pair earns its place only if apo/holo separation clears the bundles' own heterogeneity — NMR spread is
precision as much as dynamics. 14 of 16 pass.

Include this. It's the slide that tells the room you didn't just take deposited coordinates at face value.

**S10 — Three-way comparison** ⭐ [3]

| source                               | median fnat | pass at 0.5 |
| ------------------------------------ | ----------- | ----------- |
| holo bundle, DNA stripped, re-docked | 0.839       | 100%        |
| **real apo bundle**                  | 0.629       | **92%**     |
| **BioEmu ensemble**                  | 0.396       | **11.5%**   |

Apo beats BioEmu in 13 of 16 pairs. Each row answers a different objection: holo-redocked shows the dock
is lossless and the floor is reachable; real apo mostly clears it, so the gate isn't rejecting free-state
conformations as such.

**S11 — And on our own pilots** ⭐ [2]

> **BioEmu's single best frame is below the deposited apo median in 5 of 5 pilots.**

Across ~90 frames per pilot. Not a sampling-efficiency problem — the conformations aren't in the
distribution.

**S12 — What BioEmu _does_ do, and the open caveat** [1]
Where it lands on the apo→holo axis varies from 2% to 94%, median 55% — **not apo-biased, unsystematic.**
BioEmu's own paper reports a _holo_ preference on cryptic pockets (85% vs 49%); neither bias reproduces
here. What does reproduce: ensembles 1.9–10.2× wider than the deposited bundles, in all 14 pairs.

**Then the caveat, offered before anyone asks:** the free state genuinely contains partially unfolded
members, and NMR determination selects for a determinable fold. Some of the gap may be BioEmu correctly
sampling a real unfolded tail. The folded-subpopulation split is an hour of work and it's next.

---

### Act 3 — What does augmentation actually teach? (13 min)

**S13 — The experiment** [1.5] _(reuse old slide 81 — the paired-training design)_
Baseline vs augmented, paired by seed, differing only in `data_dir`. **One sentence on the shared label:
every conformation of a protein carries that protein's single experimental PWM.** Flag it now.

**S14 — Structural fidelity first** [2] _(reuse old slide 76)_
fnat distributions, gate pass rates, the rigidity ordering. The states that reach training are
geometrically reasonable — iRMSD and fnat correlate at ρ = 0.90.

This matters for the talk's logic: it forecloses "the structures were just bad" as an explanation for what
follows.

**S15 — Conformation matters more than training noise** [2]
Conformation SD **0.235** vs seed SD **0.16**, n = 701. The variation across conformations exceeds the
variation across independently trained models.

**S16 — What that looks like** [1.5]
Exemplar logos: EGR1's worst state losing the poly-G signature; TBP's best indistinguishable from the
crystal.

**S17 — The gain decays with distance from the training frames** ⭐ [3]
+0.246 own complex → +0.03–0.11 own family → +0.025 unseen complex → ≈0 cross-benchmark.

**Volunteer the in-sample structure here.** All 707 states are in-sample for their own pilot's augmented
checkpoint, verified per entry — and the inputs are _not_ near-copies (median feature distance 0.81; two
pilots sit further from their own crystal than unrelated training entries do). **Label-side memorization
with structurally distinct inputs.**

Retroactively explains August's slide 82 ("some evidence of overfitting").

**S18 — Augmentation teaches conformational invariance** ⭐ [2.5]
Across-conformer SD **0.235 → 0.180**, 11 of 12 pilots, p = 0.027. Accuracy on unseen complexes: **+0.025.**

> _The model becomes much more consistent across conformations and barely more accurate. That's exactly
> what we asked for — we told it every conformation of this protein means this one motif, and it learned
> that. Augmentation can't recover conformational signal because the objective removes it._

Pause. Then close forward: _the fix is identifiable — per-frame labels, already on my August slide 106._

HSF1 as the sole exception (0.386 → 0.519).

**S19 — A methods aside worth two minutes** [0.5] _(reuse old slide 86)_
The pseudoreplication point you already made in August — seed is not the experimental unit. One pilot's
pooled entry×seed t-test gives p = 0.001 while its five per-seed deltas are 2/5 negative.

Keep it short, but keep it: it establishes that the statistical care in this talk predates the results.

---

### Act 4 — Can we predict who benefits? (10 min)

**S20 — The question, and why it's the hard one** [1]
Augmentation helps some pilots and hurts others. What explains the difference? That was August's
hypothesis, and it's where most of the last two months went.

**S21 — The retraction, in five beats** ⭐ [4]

1. _The correlation was real and got stronger under scrutiny._ ρ = −0.87, family-corrected permutation
   p = 0.005, leave-one-out stable. Not a search artifact — tested directly.
2. _Then I asked what the quantity was._ `bend_uu` is a terminal axis-vector angle, 69% duplex length,
   with the length effect running **opposite** to worm-like-chain expectation.
3. _So I rebuilt bend from local base-pair-step geometry, twist-aware._ Central window: flat. Each end:
   flat. End-to-end orientation: flat. **Complete accumulation of every step: flat — and correlated with
   `bend_uu` IQR at ρ = +0.05.**
4. _That's the problem._ The complete reconstruction isn't a partial view. `bend_uu` disagrees with the
   geometry it's nominally a summary of.
5. **Withdrawn because the quantity can't be identified** — not because I think it's absent.

Say the last point twice.

**S22 — Can sequence predict it instead?** [2]
A clean positive, and it stands on its own. Deep DNAshape predicts DNA shape _and fluctuation_ from
sequence. Does predicted fluctuation track what we measure in the ensembles?

**No** — ρ = +0.030 per duplex, flat per position for Roll and Twist. **With the positive control on the
same slide:** predicted vs measured _mean_ groove width tracks at ρ = +0.386. The pipeline works; the
quantities are decoupled.

Tie to the lab's own result: groove width reports intrinsic structure, groove-width _fluctuation_ reports
protein-induced reshaping. Two methods, same conclusion. **Sequence-based shape prediction can't
substitute for ensemble measurement when the question is about a bound complex.**

**S23 — The full enumeration** [2]
Every structural axis tested, ordered, significance line drawn: reachability, free-state spread, static
crystal bend, interface MGW-FL, bend median, whole-duplex bend IQR, length-free windowed IQR, complete
reconstruction, predicted fluctuation, `r_mean`, state count, duplex length.

**No reconstructed DNA-geometry axis predicts the augmentation effect.**

**S24 — And that's a finding about the design** [1]
Seed noise is two-thirds of the between-pilot spread. At n = 12 you'd need ~60 pilots for 80% power.
**A pre-registered gate said "don't spend 24 training arms on this," and I didn't.** The question needs a
different unit of analysis, not more of the same.

---

### Act 5 — What the project can do now (8 min)

**S25 — Correctness, briefly** [1.5]
Three bugs found by audit this quarter, each invisible to the code and caught because a value was
physically impossible:

- **irf had no base-paired duplex** — 0 of 184 states against 2,224 of 2,224 elsewhere.
- **TBP's PWM label pointed at NDT80**, a yeast protein in a different Pfam clan. Quarantined.
- **A naive roll+tilt bend sum** gave 192° over 8 steps — successive steps are rotated by helical twist,
  so wedge angles at opposite phases must cancel.

The lesson: **physical assertions catch this class; statistical ones don't.**

Place it here, not earlier — it earns trust for Act 4's negatives _after_ they've landed, and it sets up
Act 5 as a project that knows where its floor is.

**S26 — A model-free mechanism read-out** [2]
Dock a deposited apo structure onto bound DNA, score fnat. No BioEmu, no DeepPBS, no training. It splits
the pair set and agrees with the literature where the literature is clear — Trp repressor and NHP6A fail
(both documented induced-fit), MeCP2, WRKY4, FOXD3 pass.

A structural measurement of **how often the DNA-binding interface is already formed in the free state.**

**S27 — And a way to validate it** ⭐ [3]
_Give this room — it's the most discussable slide in the deck._

Sox2: **a known right answer and a known wrong one.**

- **2LE4** (2011, unpublished): helices significantly reoriented relative to the bound state — this is
  what supported the induced-fit model
- **9QPF** (Orsetti _et al._, _NAR_ 2025): _"a three-helical bundle conformation identical to its
  DNA-bound state"_, 0.8 Å backbone RMSD to the bound structure

Their conclusion, from NOESY, RDCs, ¹⁵N relaxation and MD: **conformational selection, not induced fit.**
2LE4's RDC Q-factor is 0.76 against 0.08 for the refined structure.

Can the read-out tell them apart from interface geometry alone? **Half a day, queued.**

**The pre-registration, which is the discussion hook:** _I expect this may come back "both pass," and I've
written that down in advance. The disorder trim removes half the native contact set, and the helix that
distinguishes the two structures makes only 2 of 72 DNA contacts. If both pass, that's a quantified blind
spot — the read-out scores the ordered half of the interface — and that's worth knowing about a method._

**S28 — Three things running or queued** [1.5]

- **33 more apo/holo pairs**, launching — n from 14 to ~47, turning the mechanism split into a survey.
- **Per-frame labels** — a λ-blend sweep interpolating between the experimental PWM and each
  conformation's own prediction. Tests the invariance result directly.
- **NDT80 as a new pilot** — the silver lining of the TBP label error: label, structure and test anchors
  finally agree.

---

### Act 6 — Open questions (10–15 min)

**S29 — Where I'd value input** ⭐

Five real decision points. Put them up and stop talking.

1. **Is the folded/unfolded split the right control for the BioEmu gap?** The free state really does
   contain partially unfolded members. How would you separate "generator is wrong" from "generator is
   right and the metric measures foldedness"?
2. **The read-out scores about half the native contact set** — blind to the terminal arms, which is often
   where binding-coupled ordering happens. Is "is the DBD fold preformed?" worth answering on its own?
3. **Per-frame labels: whose labels?** Self-labelling with DeepPBS is cheap but circular. Physics-based
   scoring and SELEX-derived labels are both expensive. Better source?
4. **Is the augmentation negative worth writing up on its own?** "Shared-label augmentation teaches
   invariance" generalizes past this model — the field keeps bolting ensembles onto predictors.
5. **Which TF families would make the apo/holo survey most useful?** 33 pairs are launching; the next
   round is choosable.

**This slide stays up through the discussion. Don't advance past it.**

---

## Part 3 — Backup slides

| slide                                             | for                                  |
| ------------------------------------------------- | ------------------------------------ |
| Old deck 25–80 (full methods)                     | any pipeline question                |
| Old deck 9–22 (AF3 vs BioEmu per TF)              | "why BioEmu and not AF3?"            |
| Cross-pilot design schematic                      | "isn't +0.025 leakage?"              |
| Size-matched donor control                        | "isn't it just more training data?"  |
| Within-pilot null (7/11 slopes, p = 0.55)         | "does it work state by state?"       |
| Permutation distribution, best-of-nine            | "didn't you search over transforms?" |
| Per-pilot ΔP with seed error bars                 | "how noisy is the effect?"           |
| TBP label audit detail                            | "what happened to TBP?"              |
| Sox2 interface breakdown (72 contacts by segment) | "why would both pass?"               |
| Twist-aware accumulation formula                  | "how did you reconstruct bend?"      |
| Trim-survival numbers                             | "what does the trimming cost?"       |
| `%across` per-pair table                          | "how do you place the ensembles?"    |

---

## Part 4 — Prepared answers

**"So was August wrong?"**
The pipeline and the augmentation results hold. The mechanism interpretation doesn't — I tested it with
the analyses I proposed on slide 105, and the axis turned out to be a terminal-vector summary that
disagrees with a complete reconstruction of the same geometry.

**"Isn't the augmentation improvement just memorization?"**
Substantially, and I measured it. All 707 states are in-sample, verified per entry. The clean number is
+0.025 on unseen complexes from a cross-pilot design. The inputs aren't near-copies — feature distance
0.81 — so it's label-side memorization with structurally distinct inputs.

**"Why retract a correlation at p = 0.005?"**
Because I can't say what the quantity is. Every reconstruction from local DNA geometry is flat; only the
terminal-vector summary correlates, and it disagrees with a complete reconstruction at ρ = +0.05.

**"Why not run more pilots?"**
For the length-free axis, ~85. For the whole-duplex axis, 10–12 — affordable, and exactly what I'd have
done if the quantity were identifiable. I'd be replicating something I can't name.

**"What would make augmentation work?"**
Per-frame labels. The current objective maps every conformation to one motif, so it can't learn
conformational sensitivity by construction.

**"Is BioEmu just bad at this?"**
The most useful open question I have. Its own paper documents a holo preference on cryptic pockets —
85% vs 49% — and on protein–DNA I see neither bias: 2% to 94%, median 55%. What _is_ consistent is
displacement: 1.9–10.2× wider than the deposited bundles in all 14 pairs.

**"Why not AF3?"**
AF3 gives a single structure, and the published comparison shows it struggles to predict how mutations or
conformational dynamics alter DNA shape. The point of BioEmu is the distribution.

**"What happened to TBP?"**
Its PWM label pointed at an NDT80 matrix, and since that label is a substring filter at the augmentation
stage, the arm trained against the wrong motif and was scored on NDT80 structures. Found by audit,
quarantined. The silver lining: NDT80 is now the strongest new pilot candidate.

---

## Part 5 — Timing and cuts

| act                            | slides  | min     |
| ------------------------------ | ------- | ------- |
| 1 Setup and where we left off  | S1–S6   | 8       |
| 2 Does the generator reach it? | S7–S12  | 11      |
| 3 What augmentation teaches    | S13–S19 | 13      |
| 4 Can we predict who benefits? | S20–S24 | 10      |
| 5 What's possible now          | S25–S28 | 8       |
| 6 Open questions               | S29     | 10–15   |
| **total**                      |         | **~60** |

**If running long, cut in this order:** S19 (pseudoreplication aside), S16 (exemplar logos), S9
(informativeness gate), S25 (correctness) compressed to one line on S24.

**Never cut:** S6 (the verdict table), S11 (best frame below apo median), S18 (the thesis), S27 (Sox2),
S29 (the discussion).

**If the room is interactive and you fall behind:** drop Act 4 to S21 + S24 only. The retraction and the
power argument are the load-bearing parts; S22 and S23 can be told in two sentences.

---

## Part 6 — Delivery

- **Rehearse S6, S18, S21 and S27 aloud.** S6 sets the tone, S18 is the thesis, S21 is where
  reconstructing the reasoning live would cost you, S27 is where the room will engage.
- **Put n on every slide with a statistic.**
- **Effect sizes before p-values in speech.**
- **Never end a slide on a withdrawal.** Every negative is followed by what it enables.
- **Say "withdrawn," not "retracted,"** and give the reason in the same breath.
- **Never say "failed"** about a negative result. Say what it rules out.
- **S29 stays up through the discussion.**
- At an hour you can afford to **take questions during Acts 2 and 3**. Say so at the start — it changes
  the room's posture and it's where the interesting objections live.

---

## Part 7 — Do not include

- Ch. 19½'s runx reclassification until `pilot_mechanism.csv` is regenerated.
- The +0.527 best central-window cell — best of eight.
- The end-contribution construct — near-circular.
- The disattenuated −0.89 — one-sided correction at n = 10.
- The pooled-arm correlation as an improvement — CI crosses zero.
- TBP as a decisive example of anything.
- The pilot-candidate list — invites a scoping discussion that competes with Act 6.
