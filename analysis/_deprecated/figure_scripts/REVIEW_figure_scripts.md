# Review — analysis/figure_scripts/ (figs 1–9)

Reviewer pass on scientific + statistical soundness. Verified by reloading the three CSVs
and independently re-fitting the key models.

## Overall
Well-engineered suite: clean shared `_common.py`, correct data model (10 pilots × 130-entry
general benchmark × 5 seeds; csl has no relaxed run, handled as n.d.), version-robust plotting.
The **core design decision is right** — because the frozen and relaxed pipelines each retrained
their own baseline, absolute Pearson is not comparable across pipelines, so the cross-treatment
quantity is the within-pipeline seed-paired ΔPearson = aug_sN − base_sN. figs 2 and 5 use exactly this.

**One critical statistical error (fig9) and several over-stated titles.** The figures' own error
bars are honest; the prose on top of them claims more than the data support.

## Critical — fig9 mixed model is pseudoreplication
fig9 fits `metric ~ C(arm)*C(dna)` with crossed random **intercepts** on entry and seed, and
reports CSL p=0.003, RUNX p<0.001 as "power recovered vs the seed-averaged Wilcoxon."
Reproduced and diagnosed:
- The model estimates **seed Var ≈ 0**, so it treats all 50 (csl) / 100 (runx) entry×seed rows as
  independent replicates of the treatment contrast → arm SE shrinks ~2× (csl 0.028 vs seed-level 0.048).
- But augmentation is applied **once per seed** (5 retrainings); the 5–24 same-family entries within a
  seed are correlated. Seed is the experimental unit, n=5.
- Correct analyses agree with each other and contradict fig9:
  - seed-level paired t (n=5): csl p=0.16, runx p=0.20
  - mixed model with random **slope of arm across seeds**: csl p=0.076, runx p=0.128
- The per-seed deltas are a coin-flip: csl [−0.025, −0.043, +0.159, +0.162, +0.164],
  runx [+0.281, −0.020, −0.060, +0.159, +0.108] — 2 of 5 negative. Cannot support p<0.005.

**Fix:** use `re_formula="~arm"` grouped by seed (random slope), or just report the n=5 seed-level
paired test. Drop the "CSL and RUNX now significant" claim.

## Over-stated titles — figs 2 & 5
Titles assert directional family effects ("shifts further negative for the rigid families",
"raises error more for the rigid families"). Verified: at the correct seed level, **only
foxa-relaxed reaches significance** (ΔPearson p=4e-4 raw; ΔMAE similar), and nothing survives
multiple-comparison correction across the 19 pilot×dna tests. Every other CI crosses zero.
The point-estimate *pattern* (relaxation helps lef1/runx, hurts ets1/tbp/egr1/foxa) is a
reasonable narrative but should be stated as a trend, not a result. The error bars themselves are correct.

## Scientific mislabel — figs 2 & 5 "rigid families"
- fig2 names the rigid set "ETS1, TBP, FOXA"; fig5 names it "TBP, EGR1, FOXA" — **two different sets**.
- **TBP is the strongest DNA bender in the panel (~80° kink via Phe wedges)** — the opposite of rigid.
  Calling it a rigid family contradicts the project's own induced-fit framing. Relabel or drop the
  "rigid" descriptor; if a mechanistic split is wanted, define it once in `_common.py` and thread it.

## fig6 visual asymmetry
Jittered per-entry points are overlaid on the **same-family boxes only**, not on the grey
other-family background, making the colored group look denser and more positive than the paired
comparison warrants. The headline "augmentation is family-specific" is **not statistically tested**
in the figure; several same-family boxes (tbp/egr1/foxa frozen; egr1/engrailed relaxed) straddle or
sit below zero. Either add the same jitter to the background, or add a per-pilot same-vs-other test.

## Sound as-is
- **figs 1, 3, 4** — honest. fig1 correctly caveats cross-pipeline non-comparability; figs 3/4 show
  the three arms barely differ against entry-to-entry spread (the correct, deflationary reading).
- lef1 excluded from within-family figs (0 same-family entries) — correct.
- csl n.d. handling throughout — correct.

## Bottom line
The descriptive figures (1–6 as *distributions/CIs*) are scientifically sound and tell an honest,
mostly-null story: augmentation is small, mostly slightly negative, with a non-significant
relaxation-helps-benders / hurts-others trend. **fig9 should not be used** as written — its
significance is a pseudoreplication artifact. Titles on figs 2/5/9 overstate; the underlying error
bars are the trustworthy layer.
