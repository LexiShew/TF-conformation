# Analysis plan — week of 2026-09-29

**Purpose.** Every analysis on the table between now and the progress presentation, with what it
costs, what question it *decides* (not just informs), and how to read each possible outcome before
you see the result. The last column is the important one: an analysis whose outcomes you haven't
pre-interpreted is an analysis you will rationalize after the fact.

**Scope note.** Eleven of these are planned or proposed. Four are yours (§1–4); seven are additions
(§5–11). Three things are explicitly deferred (§12).

**Standing constraint on everything below.** The binding sample size in this project is the **pilot
count (n = 7–13)**, not the seed count and not the state count. No analysis here changes that. Say it
out loud in the talk before someone says it for you.

---

## Contents

| § | Analysis | Cost | Decides |
|---|---|---|---|
| 1 | Relaxed complexes through DeepPBS | 1–2 d | Does DNA variation add information, or more blur? |
| 2 | RI scores, crystal vs relaxed | 0.5 d (scoped) | Which residues' importance is conformation-dependent |
| 3 | Bend-IQR → per-bp roll/tilt/twist | 0.5 d | What the bend signal is made of, mechanically |
| 4 | Localizing bend fluctuation | 0.5 d | **Whether ρ = −0.76 survives** |
| 5 | De-leak the crystal-vs-experiment headline | 1 h | Whether your most quotable line holds |
| 6 | Per-state fnat × r-to-crystal | 1 h | Whether the fnat gate is a quality filter at all |
| 7 | apo-fnat class × own-family ΔP | 15 min | Whether a measured protein axis predicts the sign |
| 8 | Augmented-checkpoint rerun of 701 states | 0.5 d | Conformer degradation vs distribution shift |
| 9 | Conformer–crystal r × own-family ΔP | 2 h | Denoiser vs genuine-difference reading |
| 10 | Frame-count-matched reachability | 1 h | Whether TBP is still the decisive control |
| 11 | Number reconciliation | 0.5 d | Which per-pilot table is the one true table |
| 12 | *Deferred* | — | NMR pair expansion, at-scale dnarelax eval, PWM-agreement filtering |

---

## 1. Relaxed complexes through DeepPBS

**What.** Repeat the `conf_vs_crystal_pwm` protocol on the `_dnarelax` states: baseline checkpoints,
per-state predicted PWM, scored against (a) the same checkpoint's crystal prediction and (b) the
experimental JASPAR/HOCOMOCO PWM. Pair each relaxed state against its frozen counterpart so the
comparison is within-state, not across arms.

**Why it matters more than the frozen run.** In the frozen arm every state shares crystal DNA
geometry, so the entire r = 0.55 spread you measured is **protein-side blur**. The relaxed arm is
the only place in the project where DNA geometry varies per state. This is the first per-state test
of §17.3, which until now rests on 11 pilot-level points.

**Design decision to make before you run it.** The frozen report treats divergence from the crystal
prediction as degradation. Under the mechanism hypothesis, a relaxed conformer that moves *away*
from the crystal prediction and *toward* the experimental PWM is the **success mode**. Scoring
against the crystal alone will score your own hypothesis as a failure. Pre-register (b) as the
endpoint and (a) as the descriptive comparison.

**Coverage.** 7 pilots have dnarelax data. n = 7.

**How to read the outcomes:**

| relaxed vs frozen | reading |
|---|---|
| spread ↑, agreement with **experiment** ↑ | DNA relaxation adds real information. Strongest possible result; feeds Ch. 18 directly. |
| spread ↑, agreement with experiment flat or ↓ | Relaxation adds noise, not signal. The soft tether is moving DNA somewhere unproductive — check §4 for whether it's the ends. |
| spread ≈ frozen | The k = 1.5 tether is too stiff to matter at the readout level. Argues for a k sweep before any further dnarelax investment. |
| spread ↑, experiment ↑ **only for the high-bend-IQR pilots** | The cleanest possible confirmation of §17.3 — relaxation helps exactly where the fixed-DNA assumption was most wrong. |

---

## 2. Relative-importance scores, crystal vs relaxed

**What.** DeepPBS heavy-atom RI via bipartite edge perturbation, computed for the crystal and for a
handful of relaxed frames. Aggregate to residue level. Compare *which residues rank highest*.

**Scope this down or it eats the week.** Full ensemble × all interface atoms is two forward passes
per heavy atom per state — a different order of magnitude from the inference you've been running.
For the presentation, do **one pilot, crystal + 3–4 relaxed frames**, and show residue rank changes.
The systematic version is a month of work and a worse slide.

**The technical trap.** DeepPBS min-max normalizes RI *within each structure* — the published
definition divides each atom's MAE by the maximum MAE across all atoms in that same structure. Raw
RI is therefore **not comparable across conformations by construction**. Two valid options:

1. Rank correlation (Spearman) of residue orderings between structures — robust, loses magnitude.
2. The un-normalized numerator `MAE(Y_a, Y_~a)` with a denominator shared across the comparison set
   — preserves magnitude, requires you to pick and document the shared denominator.

Decide this before generating anything. It is not fixable after the fact.

**Validation anchor you already have.** The DeepPBS paper reports PCC 0.605 between log-sum
aggregated residue importance and alanine-scanning ΔΔG. That licenses treating RI as approximately
physical. It does *not* license treating small RI differences as meaningful — the paper's own
correlation leaves ~64% of variance unexplained.

**How to read the outcomes:**

| result | reading |
|---|---|
| Residue ranking is stable across frames | Recognition chemistry is conformation-robust; the specificity blur in §1 is diffuse, not residue-localized. |
| A small set of residues changes rank sharply | **The most presentable result.** Those residues are the conformation-dependent contacts — the structural correlate of induced fit. Cross-check against known interface residues from the literature for that pilot. |
| Ranking is unstable everywhere | Likely a normalization artifact. Go back to the denominator choice before believing it. |

---

## 3. Decomposing bend-IQR into per-bp roll, tilt, twist

**What.** Replace the scalar ensemble bend-IQR with its mechanical components and re-run the §17.3
correlation against own-family ΔPearson for each.

**Manage the expectation.** You are splitting one axis into three at n = 11 pilots, in an analysis
where nothing currently survives Benjamini–Hochberg correction. If one component comes out "more
significant" than the composite, that is **selection, not discovery**. Report it as *what the bend
signal is made of*, not as a better predictor. If you want it to be a predictor, it needs pre-
registration and more pilots.

**How to read the outcomes:**

| result | reading |
|---|---|
| Roll dominates | Consistent with groove-width-coupled deformation; connects directly to the MGW-fluctuation framing in the *Biophys J* 2026 paper. |
| Twist dominates | Suggests torsional register variation — a different and arguably more worrying failure mode, because register errors would corrupt the PWM alignment, not just the shape channels. Check `aln_score` in the `.npz` if so. |
| Tilt dominates | Least expected; would want a sanity check that it isn't a 3DNA convention artifact at the duplex ends. |
| No component separates | Bend-IQR is a genuinely composite quantity. Fine — report the composite and stop. |

---

## 4. Localizing the fluctuation along the DNA

**This is the highest-value item on the list, and not for the reason it was proposed.** It is a
validity check on ρ = −0.76, the single correlation in the entire mechanism analysis whose bootstrap
CI excludes zero.

**What.** Per-base-pair variance of bend (and of the §3 components) across the ensemble, as a profile
along the duplex. Then: **interface-restricted bend-IQR** — recompute the IQR using only base pairs
within the protein footprint — and re-run the §17.3 correlation with it.

**Why it gates everything.** If the bend variance lives in the terminal base pairs, bend-IQR is
substantially measuring end-fraying under a soft tether rather than recognition-relevant
deformation. Ch. 18.1 already notes that in the frozen arm 6 of 7 pilots *over*bend relative to the
crystal, which is consistent with ends doing work they shouldn't. This is cheap and it decides
whether §17.3 is a finding or an artifact. **Run it on Day 1.**

**How to read the outcomes:**

| result | reading |
|---|---|
| ρ holds or strengthens on interface-restricted IQR | The result is real *and* now mechanistically located. Best single slide available to you this week. |
| ρ collapses when ends are excluded | You caught an artifact before presenting it. Painful, cheap, and far better than the alternative. Ch. 17 reverts to "no axis identified" and the honest talk is Claims 1–2 only. |
| ρ splits by pilot | Informative either way — likely tracks DNA length or how much duplex sits outside the footprint. Check against per-pilot duplex length before interpreting. |

**Comparable object in the literature.** Jiang, Shewchuk et al. (*Biophys J* 2026) report per-position
MGW-FL profiles partitioned into core-motif / spacer / core-motif regions, with fluctuation reduction
being orientation- and sequence-dependent. Your bend-variance profile is the same kind of object and
can be shown beside theirs.

---

## 5. De-leak the crystal-vs-experiment headline

**What.** Recompute crystal / ensemble-mean / single-conformation agreement with the experimental PWM
on the **five held-out pilots only** (excluding CSL, EGR1, Engrailed, FoxA, IRF1, LEF1, NFAT).

**Why.** "The crystal tracks experiment best, 11 of 12, r = 0.67" is your most quotable line and
currently your most attackable one: 7 of 12 crystals are in the fold-0 training set. The report
already flags this; the talk needs the corrected number as primary and the pooled one as the
optimistic bound.

**Built-in positive control.** HSF1's crystal (5D5U) was featurized fresh and is not in training —
and it is also the pilot where the crystal does worst. Consistent with leakage carrying part of the
crystal's advantage. Worth stating explicitly.

**How to read the outcomes:**

| result | reading |
|---|---|
| Held-out crystal advantage ≈ pooled | Leakage wasn't driving it. Quote freely. |
| Held-out advantage shrinks materially | The crystal's edge is partly memorization. Reframe: "the crystal is the best *available* single representative," not "the crystal is right." |
| Held-out advantage disappears | Conformers and crystal are on equal footing out of sample — which would make the augmentation null much less surprising and is itself worth a slide. |

---

## 6. Per-state fnat × r-to-crystal

**What.** Correlate each state's fnat (already in `<pdb>_fnat.csv`) against its predicted-PWM
agreement with the crystal. Pooled and per pilot.

**Why it's a validity question, not a tuning question.** Ch. 8 sells `FNAT_FLOOR = 0.5` as "the single
structural-quality filter." Ch. 15 hands the resulting pass-rate ranking to Ch. 17 as the rigidity
axis. **HSF passes the gate at ~99% and has the worst crystal agreement in the set (0.27 ± 0.39 on 96
well-sampled states).** If that decoupling is general, both chapters are describing the gate as
something it isn't.

**How to read the outcomes:**

| result | reading |
|---|---|
| Strong negative ρ (fnat ↑ → r ↑) | Gate is doing real work. HSF is an outlier needing its own explanation. |
| ρ ≈ 0 | **The gate filters geometry and nothing else.** Say so plainly in Ch. 8; stop calling it a quality filter; the Ch. 15 rigidity axis becomes "interface-geometry conservation," which is a narrower claim than Ch. 17 currently leans on. |
| Pilot-dependent | Most likely, and most interesting — the gate works for rigid modules and not for floppy ones, which is itself the rigidity ordering showing up as a *property of the metric*. |

---

## 7. apo-fnat mechanism class × own-family ΔPearson

**What.** Cross the Ch. 19½ apo-fnat class against the Ch. 16 own-family effect for the five pilots
that have both. Both tables already exist; this is a fifteen-minute join.

| pilot | apo class (Ch. 19½) | own-family ΔP (Ch. 16) | hypothesis predicts | agrees? |
|---|---|---|---|---|
| ets1 | selection-like (100%) | +0.111 | help | ✓ |
| hsf | selection-like (100%) | +0.037 | help | ✓ |
| nfat | borderline (60%) | +0.015 | help | ✓ (weak) |
| engrailed | selection-like (100%) | −0.015 | help | ✗ |
| runx | **induced-fit** (0%, reclassified) | **+0.084** | hurt | ✗ |

3 of 5 — a coin flip. **runx is the sharpest case:** reclassified induced-fit on the S-switch
evidence, and the third-largest *positive* effect in the set.

**Why this helps rather than hurts.** §17.2 already concluded that protein-side axes don't predict
the sign, using two axes (free-state spread, reachability) that are both derived from the model. This
uses a **measured, model-free** protein axis and reaches the same conclusion. It converts §17.2 from
"our axes failed" to "protein-side axes fail, including one we measured experimentally" — which is
what makes the DNA-side result in §17.3 a discovery rather than a fallback.

**Caveat to state.** n = 5, and Ch. 19½ already warns the read-out measures the DBD fold rather than
the protein (ets1's autoinhibitory module is trimmed away). Don't over-claim; the value is the
direction, not the count.

---

## 8. Augmented-checkpoint rerun of the 701 states

**What.** Run the same 701 frozen states through the `augmented_<tf>_fold0_sN` checkpoints. Inference
only — no training.

**The confound it closes.** The baseline checkpoints have never seen a minimized BioEmu frame. So
"conformers predict worse than the crystal" is currently ambiguous between:

- *conformers carry degraded specificity* (a claim about structures), and
- *the model is out of distribution on them* (a claim about the model).

These are different claims with different fixes, and you do not want to discover which one you have
during Q&A.

**How to read the outcomes:**

| result | reading |
|---|---|
| Per-state r rises materially under augmented checkpoints | A substantial part of the deficit is distribution shift. Headline becomes "conformers are OOD for a crystal-trained model" — and the fix is training-side, not structure-side. |
| r essentially unchanged | The degradation is a property of the structures. Your current headline stands, and stands much more firmly for having been tested. |
| r rises for some pilots only | Cross-tabulate against surviving-state count — pilots with more augmented frames in training should shift more. If they don't, something else is going on. |

---

## 9. Conformer–crystal agreement × own-family ΔPearson

**What.** Per-pilot mean conformer-vs-crystal r against own-family ΔPearson, across all 11 pilots.

**Why.** The eyeballed direction across six pilots is that **low** crystal agreement goes with
**positive** ΔP (HSF 0.27 / +0.037; NFAT 0.33 / +0.015) and **high** agreement with **negative** ΔP
(engrailed 0.84 / −0.015; TBP 0.81 / −0.031). If that holds across all eleven, the gain comes from
conformers **differing** from the crystal, not from repeated observation of it — the opposite of the
denoiser reading in the report's discussion, and the most interesting possible result in the project.

**The mandatory control.** The report shows pilots with more surviving states sit closer to the
crystal (ρ = 0.60, p = 0.039). Surviving-state count **is** the fnat pass rate, which **is** the Ch. 15
rigidity ordering. So conformer–crystal agreement is not independent of the rigidity axis. **Partial
out state count** (or regress on both) before drawing any conclusion, or you will rediscover the
rigidity ordering and label it a mechanism.

**How to read the outcomes:**

| result | reading |
|---|---|
| Negative relationship survives the partial | Augmentation gains come from conformational *difference*. Reframes Ch. 16 and directly motivates per-frame labels (Ch. 21). |
| Relationship disappears after partialling | It was rigidity all along. Report as a negative and move on — cheap to have checked. |
| Positive relationship | The denoiser reading is right: conformers help by being redundant with the crystal. Also publishable, and it makes the regularization horn of confound #1 the answer. |

---

## 10. Frame-count-matched reachability

**What.** Recompute the Ch. 17 reachability `d_min` column with frame counts matched across pilots
(best-k, k = min ensemble size), as Ch. 19½ §19½.3 already does for fnat.

**Why.** TBP's 0.59 Å "best reachability" looks like a min-over-frames statistic, and TBP has 99
frames. Min-over-N grows with N. You already found this exact bias in the NMR analysis and correcting
it moved p from 0.002 to 0.15 — a conclusion-changing amount. The correction has not been applied to
Ch. 17, where the uncorrected number supports the "TBP is the decisive negative control" claim that
Ch. 4, Ch. 17 and Ch. 23 all lean on.

**How to read the outcomes:**

| result | reading |
|---|---|
| TBP still best after matching | The control holds. Report the matched number and note the check was done. |
| TBP's rank drops | The "best reachability yet hurt most" framing dissolves. Not fatal — §17.3's DNA axis doesn't depend on it — but Ch. 4, 17 and 23 all need rewording, so find out now. |

---

## 11. Number reconciliation

**What.** One canonical per-pilot table; delete every competing set.

Conflicts currently on record:

- **Cross-benchmark ΔP signs disagree.** Outline Ch. 16: csl +0.027, nfat +0.015 (net-positive).
  Results Ch. 16: csl −0.050, nfat −0.038. Opposite signs for the same quantity.
- **TBP effect.** Outline Ch. 17: "hurt most (−0.26)". Results table: −0.031. Outline Ch. 16: −0.015.
- **Family-level effects.** Outline: ETS +0.042, IRF +0.039. Results §16.3: ETS +0.008, IRF +0.010.
  4–5×.
- **Interface MGW-FL.** Outline: ρ ≈ −0.55. Results §17.3: ρ = −0.37, n.s.
- **fnat ↔ iRMSD.** ρ ≈ −0.84 (outline) vs −0.90 (results). Already flagged; just pick one.
- **Pilot roster.** Quoted as 6, 11, 12 and 13 in different places. One roster table with an explicit
  inclusion rule (dux4: 0 survivors; lef1: no own-family entries) fixes all of it.
- **NMR pair counts.** Ch. 19½ says 14 informative; §17.4½ says 13. Both are right for different
  things (14 informative, 13 with usable fnat after c03's broken reference). State it once.
- **Un-propagated retraction.** "Relaxed DNA bends toward the crystal" is retracted in Ch. 18.1 and
  the figure plan, but Ch. 20's thesis and Ch. 23's falsifiable through-line still rest on it.

**Unresolved and worth its own paragraph in the book: the fnat scale gap.** Ch. 15 has ets1 at ~100%
pass, engrailed ~91%, hsf ~99%. Ch. 19½'s pilot table has the *same* BioEmu ensembles at 1%, 0%, 18%.
Ch. 19½ attributes it to NMR coordinate precision, but a 100 → 1 gap needs a direct statement. Two
candidate explanations:

1. NMR references are simply much stricter (plausible; testable by re-scoring a pilot against an NMR
   reference where one exists).
2. **The pilot gate is partly circular** — the ensemble is docked onto crystal DNA from the same
   structure it was seeded from, and fnat is scored against that same structure's contact set. If so,
   the pilot pass rates are inflated and the Ch. 15 rigidity ordering inherits the inflation.

Decisive test: score a pilot's ensemble against a **different co-crystal of the same TF**. Cheap, and
it settles which chapter's numbers to trust.

---

## 12. Deferred, with reasons

- **The 21 staged NMR pairs.** Highest data-value-per-GPU-hour in the project (14 → 35 informative
  pairs would turn the mechanism split into a classifier with power, and would make the current
  `%across`-vs-fnat non-correlation, Spearman +0.36 p = 0.23, a real negative rather than an
  underpowered one). But it will not land inside a week. Next month.
- **The at-scale `_dnarelax` training eval.** Nothing survives BH correction across 11 pilots. More
  seeds cannot fix that, because **the n is pilots**. Run it for completeness; report direction only;
  do not expect significance and do not budget time against that hope.
- **Weighting or filtering training states by predicted-PWM agreement.** Do not do this. It selects
  frames the model already agrees with, so the augmented model is guaranteed to regress toward the
  baseline and the result is uninterpretable. Weighting by fnat, iRMSD or bend-IQR is fine — those are
  independent of the model's output. Keep PWM agreement strictly diagnostic.

---

## Week schedule

| day | work |
|---|---|
| 1 | §4 interface-restricted bend-IQR + §17.3 re-run (gating). Launch §1 relaxed inference, both references. |
| 2 | §5 de-leak, §6 fnat × r, §7 apo-fnat cross-table. All three short, all three defensive. |
| 3 | §1 analysis, paired within state. §8 augmented-checkpoint rerun in parallel. |
| 4 | §3 decomposition + §4 per-position profile figure. |
| 5 | §2 single-pilot RI exemplar. |
| 6–7 | §11 reconciliation, figures, talk. |

Leave both of days 6–7. Four converging analyses plus the number reconciliation is more assembly
work than it looks.

---

## What the talk should claim

Three claims, in this order:

1. **A crystal is a point estimate of a distribution.** n = 701 states; conformation SD 0.25 vs seed
   SD 0.16. Seed-independent, large-n, and it frames everything else.
2. **Protein-side conformational variation blurs the specificity readout toward nonspecific; it does
   not redirect it.** The −1.14-bit information deficit, the +0.20 recovery from ensemble averaging,
   58/59 wins. Calibration point: the DeepPBS paper reports that nonspecific complexes (SSO7D) give
   notably lower information content — so the model is reading the average conformer as *more
   nonspecific-like* than the crystal. Since DNA is frozen in this arm, only the protein varies.
3. **The residual signal is a DNA property, and it is localized.** Interface-restricted bend-IQR, the
   per-position fluctuation profile, and — if §1 cooperates — the relaxed arm moving toward experiment.

Claims 1 and 2 hold regardless of how the week goes. **Only claim 3 can fail.** That is the right risk
posture for a progress talk: you have a presentation either way, and the thing that can fail is the
thing worth presenting.

---

## References

- Mitra, R., Li, J., Sagendorf, J. M., Jiang, Y., Cohen, A. S., Chiu, T.-P., Glasscock, C. J. & Rohs,
  R. **Geometric deep learning of protein–DNA binding specificity.** *Nature Methods* **21**,
  1674–1683 (2024). doi:10.1038/s41592-024-02372-w — RI definition and within-structure
  normalization (Methods, "Bipartite edge perturbation…"); RI vs alanine-scanning ΔΔG PCC 0.605
  (Fig. 4h); nonspecific-complex information content (SSO7D, Supplementary Fig. 8); robustness to
  small dynamical fluctuations in MD application.
- Jiang, Y., Shewchuk, A. M., Chiu, T.-P., Li, J., Kribelbauer-Swietek, J. F., Gompel, N. & Rohs, R.
  **Readout of intrinsic and induced DNA shape by homeodomain transcription factor complexes.**
  *Biophysical Journal* (2026). doi:10.1016/j.bpj.2026.03.036 — *"minor groove width reflects
  intrinsic DNA structure, whereas minor groove width fluctuations capture protein-induced
  stabilization and reshaping of the DNA"*; *"recognition involves a context-dependent balance
  between conformational selection and induced fit."* The independent convergence with §17.3, and the
  per-position profile format for §4.
- Lewis, S. *et al.* **Scalable emulation of protein equilibrium ensembles with generative deep
  learning.** *Science* (2025). doi:10.1126/science.adv9817 — *"Presently, oligomer and ligand binding
  state are implicit… A hallmark of this may be BioEmu's preference to predict holo over apo
  structures in the cryptic pocket benchmark"*; 85% holo vs 49% apo success at a 1.5 Å threshold.
  The documented bias that Ch. 19½ measures on protein–DNA systems.
