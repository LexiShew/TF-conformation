# Chapters 15–18 — Part IV Results: Fidelity, Effect, Mechanism, and DNA Relaxation

> **What these four chapters are.** This is the results arc, read as one argument. Chapter 15 establishes
> that the ensembles are structurally faithful in a rigidity-ordered way. Chapter 16 shows the net
> augmentation effect on the cross-benchmark is small and mostly negative — the apparent null that
> motivates the reframe. Chapter 17 resolves that null into a mechanism map — **and this is where the
> data on disk revises the original hypothesis**. Chapter 18 stress-tests the map with the DNA-relaxation
> arm.
>
> **Sourcing note — this arc is data-grounded, not code-grounded.** Unlike the pipeline chapters, these
> report *findings*, so the ground truth is the analysis outputs on endeavour, not engine parameters:
> `analysis/data/{mechanism_apo_holo,reachability,ensemble_diversity,perentry_accuracy,family_annotation}.csv`,
> the `analysis/dna_relax/data/` suite, the 24 `id_benchmark_<tf>.json` eval files, and two authoritative
> memos — `analysis/mechanism/FINDINGS.md` (2026-08-01) and `analysis/mechanism_and_roadmap.md`. Values
> read from these are marked **[data]**; the memos' own conclusions are marked **[FINDINGS]**.
>
> **A correction the successor must know up front.** The book outline's Chapter 17 is written on a
> **static** apo/holo axis (free-state protein spread × qualitative crystal DNA deformation). The
> re-analysis on disk shows that static axis **does not predict the augmentation sign** — it groups ETS1
> with TBP, the exact pair the hypothesis needs to separate. **[FINDINGS]** The mechanism survives, but in
> a *modified, dynamic* form: what matters is how much DNA shape **fluctuates across the ensemble**, not how
> deformed the bound duplex is. Chapter 17 below is written to the revised finding, with the original static
> framing shown as the hypothesis it replaced. Do not present the static apo/holo table as the result.

---

# Chapter 15 — Structural fidelity (the ensembles are trustworthy, unequally)

**The claim of the chapter.** Before any accuracy number can be interpreted, the ensembles have to be shown
structurally faithful — and they are, but *unequally*, in a way that orders the pilots along a single axis.
That axis (recognition-module rigidity, operationally the fnat pass rate) is the one Chapter 17 later splits on.

## 15.1 — Fidelity is rigidity-ordered

**What the data shows.** fnat pass rate, interface RMSD, and Cα-RMSD all co-vary and rank the pilots
monotonically from rigid to floppy: **[data]**

- **Rigid, high fidelity:** ETS1, TBP, HSF (≈99–100% of states clear the 0.5 fnat gate).
- **Intermediate:** EGR1, engrailed (~91%), IRF, FOXA (~75–86%).
- **Floppy, low fidelity:** CSL, ERR (~50%), RUNX, NFAT (~30–37%), LEF1 (19%), DUX4 (0%).

**Fidelity tracks rigidity, not interface size.** The decisive contrast is TBP vs LEF1: near-identical interface
size (≈40 vs 39 interface residues) but opposite fidelity (TBP ~100% pass, LEF1 19%). **[outline+data]** Interface
size does not predict fidelity; recognition-module rigidity does.

**fnat and iRMSD co-vary tightly.** Across all frozen states, fnat is a near-monotone readout of interface
geometry — Spearman ρ ≈ **−0.90** between fnat and global interface RMSD in the current per-state metrics **[data]**
(the outline quotes ρ ≈ −0.84 from the earlier fidelity set; both say the two criteria measure the same thing, so
the 0.5 fnat gate corresponds to a consistent interface-RMSD threshold rather than an arbitrary cut).

## 15.2 — Distortion is localized and minimization is a gentle nudge

- **Localized distortion.** One interface segment always moves more than the rest; the per-residue signature is
  universal — rigid cores, floppy termini. **[outline]** (Figures `I2_iRMSD_seg`, `R2_per_residue_profiles`.)
- **Minimization is a local relaxation.** Global backbone barely moves between the docked (Stage 2) and minimized
  (Stage 3) pose, and per-state changes nudge marginally *toward* the crystal, never away — established
  quantitatively in the Chapter 9 figure suite (backbone RMSD Stage 2 ≈ Stage 3; 70% of states move below the
  no-change diagonal in interface RMSD). **[data, Ch. 9 figs]** This is what licenses treating the minimized
  ensemble as a faithful relaxation of the docked ensemble rather than a new structure.

## 15.3 — Why this matters for the rest of Part IV

The rigidity axis is not decoration — it is the independent variable Chapter 17 correlates the augmentation sign
against. A crucial honesty point carried forward from the "where does the rigidity metric come from" question:
**the rigidity ordering is operationally the fnat pass-rate ranking**, so any figure that plots the augmentation
effect against "rigidity" is really plotting it against a fnat-derived ordering. Chapter 17's genuinely
independent axis is the *DNA* one (measured bend fluctuation), not the protein rigidity ordering — which is
exactly why the mechanism finding is non-circular where a protein-rigidity correlation would not be.

## Figures — Chapter 15
- **Existing:** `I1_iRMSD_distributions`, `I2_iRMSD_seg`, `I4_interface_geometry`, `R1_ca_rmsd_stages`,
  `R2_per_residue_profiles`, `R3_minimization_delta`, `S1_bestworst_{ets1,lef1}`, `S2_stage_progression_ets1`.
- **Reuse from the Ch. 9 ensemble suite:** the fnat violins/pass-rate ladder (fig16/17), backbone-vs-sidechain
  (fig6), interface-tightening (fig7), and fnat-vs-iRMSD (fig19) — all built from the same per-state metrics and
  all rigidity-ordered.

---

# Chapter 16 — The augmentation effect (the productive null)

**The claim of the chapter.** On the 130-entry cross-benchmark, augmentation's net effect is small and mostly
flat-to-negative — and that is *setup, not verdict*. The signal is real but it is not where the cross-benchmark
average looks; it lives in each pilot's **own family**.

## 16.1 — The cross-benchmark average washes the effect away

**What the data shows.** At the seed level (seed as the unit of replication), the 130-entry cross-benchmark
ΔPearson is small for every pilot and has 2–3 of 5 seeds negative for almost all of them — i.e. not
significant. **[data]** The three least-negative-to-positive pilots (ETS, IRF, CSL-adjacent) hint at structure,
but the pooled mean is a null.

**Why the average is the wrong lens.** The cross-benchmark pools ~120 entries from families the pilot's ensemble
says nothing about. A protein's free-state ensemble can only inform predictions for *its own* recognition family;
averaging over unrelated families dilutes any real effect toward zero. **[FINDINGS]**

## 16.2 — The signal is in the own-family subset

**The headline number.** Restricting to each pilot's own motif-level family (seed-paired, seed as unit): **[FINDINGS]**

| pilot | own-family ΔPearson | 95% CI | p | neg. seeds | cross-benchmark ΔP |
|---|---|---|---|---|---|
| **ets1** | **+0.111** | [+0.035, +0.187] | **0.013** | 0/6 | −0.002 |
| csl | +0.103 | [−0.010, +0.215] | 0.065 | 2/6 | −0.050 |
| runx | +0.084 | [−0.048, +0.216] | 0.161 | 2/6 | +0.001 |
| hsf | +0.037 | [−0.055, +0.129] | 0.331 | 2/5 | +0.003 |
| irf | +0.033 | [−0.083, +0.149] | 0.471 | 2/5 | −0.013 |
| nfat | +0.015 | [−0.042, +0.072] | 0.500 | 2/5 | −0.038 |
| egr1 | −0.003 | [−0.058, +0.053] | 0.911 | 4/6 | −0.017 |
| err | −0.003 | [−0.067, +0.062] | 0.924 | 3/6 | −0.023 |
| foxa | −0.006 | [−0.094, +0.082] | 0.877 | 3/6 | −0.004 |
| engrailed | −0.015 | [−0.102, +0.071] | 0.670 | 3/6 | +0.016 |
| tbp | −0.031 | [−0.144, +0.082] | 0.513 | 3/6 | −0.012 |

**ETS1 is the anchor result:** own-family ΔPearson **+0.111 (p = 0.013, 0 of 6 seeds negative)** against a
cross-benchmark value of −0.002 — the effect is entirely invisible in the pooled average and clear in the own
family. **[FINDINGS]** (lef1 is excluded: it is the only HMG-box in the benchmark and has no same-family entries.)

**The within-pilot contrast (the test `fig6` never ran).** Own-minus-other family, at the seed level: 8 of 11
pilots show a positive contrast; csl (+0.145, p = 0.008) and ets1 (+0.119, p = 0.011) reach raw p < 0.05, **but
nothing survives Benjamini-Hochberg correction across the 11 tests.** **[FINDINGS]** The direction is consistent;
the significance is not established.

## 16.3 — Confounds excluded, and the honest ceiling

The own-family effect is **not** explained by the obvious alternatives: **[FINDINGS]**
- **Baseline headroom** (low baseline has more room to improve): own-family effect vs baseline own-family Pearson
  ρ = −0.27, p = 0.42; vs baseline MAE ρ = −0.06, p = 0.85.
- **Subset size** (small n manufactures large deltas): ρ = −0.32, p = 0.34.
- An OLS on DNA-bend IQR + headroom + subset size (n = 11, R² = 0.52) leaves the DNA axis with the largest
  standardized coefficient but no term significant at n = 11.

**Family-level effects.** The three net-positive families are exactly the three pilots with the largest own-family
effects: IRF (+0.010), ETS (+0.008), CSL/RBPJ (+0.006). The only family whose *negative* effect reaches p < 0.05
is **C2H2 zinc finger (−0.016, p = 0.017)**. **[FINDINGS]**

**Framed explicitly as setup.** A flat cross-benchmark mean that resolves into a real own-family effect is the
signature of a hidden variable splitting the families. Chapter 17 identifies that variable.

## Figures — Chapter 16
- **Existing:** `P3_augeffect_by_family` ⭐ (ΔPearson by motif-level family), `P2_baseline_by_family`,
  `P1_family_table`, `figure_scripts/fig1_three_arm_accuracy`, `fig6_within_family_transfer`. **Do not** use the
  earlier per-PDB-mislabeled version (ETS1–RUNX1 co-crystals misfiled as Runt, inflating the ETS median). **[outline]**
- `analysis/mechanism/figures/M2_where_the_signal_lives.png` — own vs other family, contrast, headroom control,
  dynamic-vs-static axis. The single best figure for this chapter's argument. **[FINDINGS]**
- *(create)* **The cross-benchmark-vs-own-family scatter** — each pilot's pooled ΔP (≈0) against its own-family ΔP,
  ETS1 flagged, to show the signal the average hides.
- *(create)* **The pseudoreplication contrast** (shared with Ch. 11/13) — csl pooled p = 0.001 beside its five
  per-seed deltas (2 of 5 negative).

---

# Chapter 17 — The mechanism story (the climax, revised by the data)

**The claim of the chapter.** The hidden variable that splits the families is a property of the **DNA**, not the
protein — and specifically the DNA's *dynamic* behavior across the ensemble, not the static deformation of the
bound crystal. This is the payoff the fidelity and effect chapters build toward, and it is the chapter where the
on-disk re-analysis materially changed the answer.

## 17.1 — The original hypothesis and the apo/holo test

The reframe (Ch. 4): treat the BioEmu ensemble as a **computational apo** state and the crystal as **holo**, and
read the *sign* of the augmentation effect as a probe of recognition mechanism. The prediction: augmentation helps
under **conformational selection** (free protein already samples binding-competent geometry) and hurts under
**induced-fit-on-DNA** (free-protein frames carry no information about the DNA distortion the model holds fixed).

The apo/holo table scores two protein-side axes per pilot — free-state spread and reachability `d_min` — against a
qualitative DNA-deformation label: **[data]**

| pilot | free-state spread (Å) | reachability d_min (Å) | DNA deformation | own-family aug ΔP |
|---|---|---|---|---|
| ETS1 | 1.62 | 0.87 | minimal | **+0.111** |
| TBP | 1.50 | **0.59** (best) | extreme ~80° kink | **−0.031** |
| EGR1 | 1.93 | 0.80 | minimal | −0.003 |
| engrailed | 2.71 | 0.69 | modest | −0.015 |
| FOXA | 5.62 | 1.79 | modest | −0.006 |
| LEF1 | 5.95 | 2.43 | severe ~110° bend | (no own-family entries) |

**The decisive negative control is TBP:** it has the *best* protein reachability (0.59 Å — its rigid β-saddle
reaches the bound backbone perfectly) yet its own family is hurt. **[data]** The protein ensemble is not the
problem; the missing variable is the DNA.

## 17.2 — The static DNA axis fails (the revision)

**What the re-analysis found.** When the qualitative DNA-deformation label is replaced with a *measured* static
index — crystal bend deviation and minor-groove deviation from canonical B-DNA, z-scored — that static axis
**does not predict the augmentation sign**: **[FINDINGS]**
- Spearman ρ = **+0.04** (p = 0.92) between static crystal deformation and the own-family effect — essentially
  zero, and if anything the *wrong* sign.
- A median split on the static index puts **ETS1 and TBP in the same group** — precisely the pair the hypothesis
  uses to separate conformational-selection from induced-fit. **[FINDINGS]**

**Why this matters.** The original Chapter 17 framing — that *how deformed the bound duplex is* predicts the sign —
is not supported. Static crystal deformation is not the discriminating axis. The clean "TBP kink → hurt, ETS
straight → helped" story is true as a pair of anecdotes but does not generalize as a static correlation.

## 17.3 — The dynamic DNA axis is the real discriminator

**The finding that replaces it.** The discriminating variable is **how much DNA shape fluctuates across the
ensemble** — operationalized as the ensemble DNA-bend inter-quartile range (bend IQR). Own-family ΔPearson vs
measured DNA geometry, Spearman with bootstrap 95% CI over pilots (10,000 resamples): **[FINDINGS]**

| axis | arm | n | ρ | boot 95% CI | same sign | p |
|---|---|---|---|---|---|---|
| **Ensemble DNA-bend IQR** | **relaxed** | 11 | **−0.76** | **[−0.96, −0.28]** | **99.6%** | **0.006** |
| Ensemble DNA-bend IQR | frozen | 11 | −0.53 | [−0.91, +0.22] | 92.9% | 0.096 |
| Interface MGW-FL | relaxed | 10 | −0.37 | [−0.92, +0.42] | 83.7% | 0.293 |
| Crystal DNA bend (static) | relaxed | 11 | −0.07 | [−0.80, +0.66] | 57.4% | 0.832 |
| Crystal DNA bend (static) | frozen | 11 | +0.04 | [−0.65, +0.76] | 53.6% | 0.915 |

**Bend IQR in the relaxed arm is the only correlation in the entire mechanism analysis whose bootstrap CI excludes
zero** (ρ = −0.76, CI [−0.96, −0.28], same sign in 99.6% of resamples). **[FINDINGS]** The negative sign is the
mechanism: the more a pilot's DNA shape *wanders across the ensemble*, the more augmentation hurts — because the
model holds DNA geometry fixed, so frames whose DNA shape has drifted carry geometry the shared crystal PWM label
no longer describes.

## 17.4 — The modified hypothesis (state this, not the static one)

> Augmentation helps a pilot's own family when the pilot's DNA is **conformationally stable across the ensemble**,
> and stops helping as ensemble DNA-shape fluctuation grows. The model holds DNA geometry fixed, so frames whose
> DNA shape wanders carry geometry the shared crystal PWM label no longer describes. **[FINDINGS]**

This keeps the mechanistic logic of the original — *the fixed-DNA assumption is the ceiling, not the protein
sampling* — while replacing the axis with a measured, dynamic one. It also re-enters the shared-label confound
(Ch. 12, Ch. 21): fluctuating DNA is exactly the condition under which one crystal PWM label mislabels its frames
most severely. The two problems are the same problem viewed from two sides.

## 17.5 — Honest limits

- **n = 11 pilots.** The frozen-arm bend-IQR CI *crosses zero* ([−0.91, +0.22]); only the relaxed-arm CI excludes
  it. **[FINDINGS]**
- **Nothing survives multiple-comparison correction** across the 11 own-family tests. **[FINDINGS]**
- Own-family subsets are 3–24 entries; the smallest are noisy. Several families are represented by a single pilot,
  so "family effect" and "pilot effect" are partly confounded there. **[FINDINGS]**
- This is a directionally consistent, mechanistically coherent trend with one CI that excludes zero — a strong
  hypothesis with a decisive negative control (TBP) and a measured axis, **not a settled result**.

## 17.6 — A data-integrity note the successor must not skip

The benchmark CSVs label the DNA conditions **`frozen` and `relax`** — not `relaxed`. An earlier
`mechanism_analysis.py` filtered on `"relaxed"`, matched zero rows, and **silently produced an empty relaxed arm**
while 119 relaxed rows sat in the input, with no error raised. **[FINDINGS]** The fix asserts the condition labels
at load time. Two effects were invisible until the relaxed arm was correctly included: **foxa relaxed ΔP −0.036
(p = 0.0004)** — the only pilot×condition effect that survives at the seed level — and **tbp relaxed −0.039
(p = 0.051)**. Anyone reusing these CSVs should check they filtered on the right label.

## Figures — Chapter 17
- **Existing:** `M1_apo_holo_mechanism` ⭐ — but caption it as *the hypothesis*, and pair it with the revision.
- `analysis/mechanism/figures/M3_samefamily_both_arms.png` and `M4_dnaflex_vs_effect.png` — the own-family effect
  in both arms and the DNA-fluctuation-vs-effect correlation with the static control. **These carry the revised
  result.** **[FINDINGS]**
- *(create)* **The static-fails / dynamic-works pair** — own-family effect vs static crystal deformation (ρ = +0.04,
  flat) beside vs ensemble bend IQR (ρ = −0.76, relaxed) — the single figure that makes the revision legible.
- *(create)* **The hypothesis-filled-with-data quadrant** — the Ch. 4 2×2 (reachability × DNA fluctuation) with the
  pilots at measured positions and TBP flagged as the decisive control.

---

# Chapter 18 — DNA relaxation (the stress test and the falsifiable prediction)

**The claim of the chapter.** The `_dnarelax` variant lets protein and DNA co-relax into each frame. Its structural
behavior confirms the DNA moves the right way without fraying; its effect on accuracy is where the mechanism
hypothesis becomes falsifiable — and, at the own-family level, the prediction is met.

## 18.1 — Structural behavior: the DNA moves, without fraying

**TBP is the worked case.** Under relaxation, TBP DNA backbone moves ~1.7× more than frozen — median backbone RMSD
**0.72 Å relaxed vs 0.41 Å frozen** **[data]** — *without* excess fraying: the k = 1.5 soft-tether stiffness floor
holds the duplex batch-wide (Ch. 9.6). Displacement localizes to the central TATA bases (the kink region). **[outline+data]**

**The DNA gains flexibility without collapsing.** The interface minor-groove-width flexibility (MGW-FL) increases
modestly under relaxation for most pilots (mean Δ ≈ +0.2; EGR1 +0.42, HSF +0.38, ETS1 +0.21; two pilots slightly
negative), confirming the DNA gains physical flexibility rather than freezing or fraying. **[data]**

**A correction to the outline: relaxed DNA does *not* measurably converge toward the crystal on the current data.**
The outline's Ch. 18 asserts that relaxed DNA "bends toward the bound conformation." The on-disk convergence data do
not support this. Across the six pyCurves shape features (`bend_uu, shortening, minor_w, major_w, inclin, tip`) for
the 7 covered pilots, relaxation moves *closer* to the crystal than the frozen arm in only **7 of 42** feature×pilot
cases; on the global bend median only **2 of 7** pilots (engrailed, dux4) get closer, and pyCurves flags only **3 of
13** pilots as `closer_to_crystal`. **[data, crystal_convergence + pycurves_bend_summary]** In the frozen arm, 6 of
7 pilots *over*bend relative to the crystal (only TBP underbends), so relaxation is not correcting a systematic
underbend either. The honest statement is: **relaxation adds DNA flexibility and lets the duplex move (TBP DNA
backbone RMSD 0.72 vs 0.41 Å; §18.1), but it does not move the DNA toward the crystal geometry on these metrics.**
Whether any convergence exists at the localized kink region specifically (as opposed to the whole duplex) would need
the per-position `crystal_convergence_bootstrap` CIs checked directly — do not assert whole-duplex convergence.

## 18.2 — The falsifiable prediction, and that it is met at the own-family level

**The prediction** (from `mechanism_and_roadmap.md`): a flexible-DNA pipeline should *recover exactly the families
augmentation hurts today* — TBP, LEF1, homeodomain — because those are the induced-fit-on-DNA cases. **[FINDINGS]**

**On the cross-benchmark average it appears to fail; on the own-family effect it is met.** Comparing frozen vs
relaxed own-family ΔPearson: **[FINDINGS]**

| pilot | frozen ΔP (p, neg/seeds) | relaxed ΔP (p, neg/seeds) |
|---|---|---|
| ets1 | +0.111 (0.013, 0/6) | +0.110 (0.017, 0/6) |
| **csl** | +0.103 (0.065, 2/6) | **+0.158 (0.004, 0/6)** |
| runx | +0.084 (0.161, 2/6) | +0.094 (0.024, 1/6) |
| hsf | +0.037 (0.331, 2/5) | +0.051 (0.228, 1/5) |
| egr1 | −0.003 (0.911, 4/6) | +0.022 (0.291, 1/6) |
| err | −0.003 (0.924, 3/6) | +0.036 (0.629, 2/6) |
| foxa | −0.006 (0.877, 3/6) | +0.030 (0.453, 2/6) |
| engrailed | −0.015 (0.670, 3/6) | +0.014 (0.771, 3/6) |
| **tbp** | **−0.031 (0.513, 3/6)** | **+0.008 (0.888, 2/6)** |

**8 of 11 pilots shift upward under relaxation.** csl and runx reach raw p < 0.05 only in the relaxed arm; **TBP —
the induced-fit exemplar whose ~80° kink motivated the whole hypothesis — crosses from negative to marginally
positive.** **[FINDINGS]** And the dynamic mechanism axis *strengthens* in the relaxed arm (bend IQR ρ −0.53 frozen
→ −0.76 relaxed, the only CI excluding zero; §17.3). The relaxation arm both improves the induced-fit families and
sharpens the mechanism correlation — the two predictions the hypothesis makes.

## 18.3 — Current status, stated honestly

- **Coverage is partial.** The DNA-shape analysis has *not* been run on all 13 pilots: the shape overview spans 7,
  the per-TF panels 11 (no hsf/irf), the bend bootstrap 7, and the AF3/interface comparisons 12. **[outline+data]**
  Every DNA-relaxation figure should be captioned **preliminary / partial-coverage**.
- **Training/eval coverage is thinner still.** The `_dnarelax` arm has data for 7 pilots but was originally trained
  and evaluated for far fewer; the own-family relaxed effects above rest on 4–6 seeds per pilot. **[FINDINGS]**
- **Nothing survives BH correction** across the 11 pilots, and the story rests on small seed counts. **[FINDINGS]**
- The direction is consistent and the mechanism axis sharpens, but this is a **directional confirmation of a
  falsifiable prediction, not a completed at-scale result.** The at-scale `_dnarelax` eval is the Part V to-do.

## 18.4 — Why this ties the arc together

Chapter 15 said the ensembles are faithful; Chapter 16 said the cross-benchmark effect is a null that resolves in
the own family; Chapter 17 said the discriminator is dynamic DNA fluctuation, not static deformation or protein
sampling; Chapter 18 shows that softening the DNA — the one intervention the mechanism predicts should matter —
moves the induced-fit families up and sharpens the very axis Chapter 17 identified. The falsifiable through-line of
the whole book (Ch. 23): *if co-relaxing DNA flips TBP from negative toward neutral/positive, the mechanism
hypothesis is confirmed and the method is fixed by the same change.* The data on disk show TBP's own-family effect
crossing zero — the first, partial-coverage evidence that it does.

## Figures — Chapter 18
- **Existing:** `dna_relax/dna_shape_features` (6-panel violin overview, 7 TFs), `dna_relax/tbp_dna_shape`
  (4-panel TBP: backbone RMSD, per-residue P displacement, P–P gap, Δbend), `dna_relax/crossfamily_bend`
  (Δ axis bend per pilot, bootstrap CIs), `crystal_convergence_bootstrap` (per-position bend vs crystal, with CIs —
  verify its direction before captioning as "toward bound"; the whole-duplex medians do not converge, §18.1), the per-TF
  `<tf>_dna_shape` set (11 pilots, no hsf/irf), and the interactive `pycurves_viz/*.html` viewers.
- `analysis/mechanism/figures/M3_samefamily_both_arms.png` — frozen→relaxed shift, the §18.2 table as a figure.
- *(create)* **The falsifiable-prediction figure** — predicted vs actual aug-sign shift for the induced-fit
  families (TBP, LEF1, homeodomain), frozen → relaxed, with the 8/11-upward result and TBP's zero-crossing flagged.

---

## Reading guide — what changed vs the outline, at a glance

| chapter | outline said | data on disk says | action |
|---|---|---|---|
| 15 | fidelity is rigidity-ordered, ρ≈−0.84 | confirmed; ρ≈−0.90 in current metrics; rigidity = fnat pass-rate ranking (non-independent) | keep; note the ordering's origin |
| 16 | net effect small/mostly negative on general-130 | confirmed as *setup*; the effect is real in the **own family** (ETS1 +0.111, p=0.013) | reframe around own-family |
| 17 | mechanism = static apo/holo (protein spread × crystal deformation) | **static axis fails** (ρ=+0.04, groups ETS1+TBP); real axis is **dynamic ensemble bend IQR** (ρ=−0.76 relaxed) | **rewrite to the dynamic axis**; show static as the hypothesis it replaced |
| 18 | DNA relaxation as consequence + prediction | prediction **met at own-family level** (8/11 up, TBP crosses zero); mechanism axis sharpens under relaxation | keep; state partial coverage plainly |

