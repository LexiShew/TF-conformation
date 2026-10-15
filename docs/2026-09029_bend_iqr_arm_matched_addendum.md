# Addendum — arm-matched bend-IQR recomputation (C2)

**Run:** 2026-09-29. **Verdict: the §17.3 headline survives. It is not a label artifact.**

## What was wrong, and what it cost

`mechanism_analysis.py:102` computes `froz_bend_iqr` from `pycurves_all_perstructure.csv` filtered on
`cond == "frozen"`, and no relaxed-arm equivalent existed. Both rows of the §17.3 table therefore used
the same frozen x-values. I computed the missing `relax_bend_iqr` with the identical estimator
(`q75 − q25` of `bend_uu`, grouped by pilot) on `cond == "relaxed"` rows — present for all 13 pilots —
and re-ran all four axis×effect combinations.

**Validation:** my frozen recompute reproduces the on-disk `froz_bend_iqr` column to 0.0 absolute
difference on all 11 pilots, so the relaxed column is computed the same way and is directly comparable.

## The 2×2

| axis arm | effect arm | matched | n | Spearman ρ | p | boot 95% CI | CI excl. 0 | same sign | q (BH, 4 tests) | ρ w/o tbp | LOO range |
|---|---|---|---|---|---|---|---|---|---|---|---|
| frozen | frozen | **yes** | 11 | **-0.527** | 0.0956 | [-0.903, +0.204] | no | 0.933 | 0.1173 | -0.370 | [-0.733, -0.370] |
| frozen | relaxed | no | 11 | -0.764 | 0.0062 | [-0.954, -0.271] | yes | 0.995 | 0.0125 | -0.721 | [-0.879, -0.709] |
| relaxed | frozen | no | 11 | -0.500 | 0.1173 | [-0.904, +0.262] | no | 0.916 | 0.1173 | -0.394 | [-0.697, -0.358] |
| relaxed | relaxed | **yes** | 11 | **-0.782** | 0.0045 | [-0.991, -0.280] | yes | 0.997 | 0.0125 | -0.733 | [-0.915, -0.733] |

10,000 bootstrap resamples over pilots, seed 20260929.

**The arm-matched relaxed correlation is ρ = −0.782 (p = 0.0045, CI [−0.991, −0.280]) versus the
published mismatched −0.764.** The result is marginally *stronger* when done correctly. Both relaxed-effect
cells survive BH correction across the four tests (q = 0.0125); both frozen-effect cells do not (q = 0.117).

**So: §17.3, §18.2, and Ch. 18.4 do not need retracting.** The claim "the only correlation whose bootstrap
CI excludes zero" holds, and now holds on an arm-matched axis.

## But the mechanism narrative does need rewording

Three things the recomputation establishes that change how this should be said:

**1. The two axes are nearly the same measurement.** Frozen vs relaxed bend IQR correlate at
Spearman ρ = **+0.927** (p = 4.0e-05; Pearson +0.881). Mean IQR rises only 5.46° → 5.79°, with 8 of 11
pilots increasing. Rank shifts are ≤2 positions for every pilot. **The k = 1.5 soft tether does not
materially change how much the DNA bend fluctuates across the ensemble** — it changes the geometry
(TBP backbone RMSD 0.72 vs 0.41 Å, Ch. 18.1) without changing the fluctuation ordering.

This is independently useful: it is direct evidence for the plan's §1 third outcome row — *"the k = 1.5
tether is too stiff to matter at the readout level, argues for a k sweep."* On the bend-IQR axis
specifically, that is now measured rather than hypothesized.

**2. The "strengthening under relaxation" is entirely an effect-arm phenomenon.** Holding the axis fixed
and swapping only the effect arm:

- frozen axis: vs frozen effect −0.527 → vs relaxed effect −0.764 (Δ −0.236)
- relaxed axis: vs frozen effect −0.500 → vs relaxed effect −0.782 (Δ −0.282)

Holding the effect fixed and swapping the axis moves ρ by ≤0.03. **The axis arm is irrelevant; the effect
arm is everything.** So Ch. 18's "relaxation sharpens the mechanism axis" is the wrong description — the
axis is unchanged. The correct statement is: **DNA-bend fluctuation predicts the relaxed arm's augmentation
effect substantially better than it predicts the frozen arm's.**

That is still a coherent mechanistic claim (the relaxed arm is where DNA geometry can actually respond, so
its effect should be the one DNA fluctuation governs) — but it is a claim about the *effect*, not about the
axis, and Ch. 18.4's tie-together sentence currently says the latter.

**3. The alternative reading is not excluded.** The relaxed-arm ΔP values rest on 4–6 seeds and two pilots
(irf, nfat) flip sign between arms. A correlation that is strong against the noisier of two effect
measurements, on the same x, at n = 11, is consistent with a lucky draw. The bootstrap CI and BH survival
argue against that, but n = 11 pilots remains the binding constraint and the talk should say so.

## Robustness

- **TBP is influential but not load-bearing** for the relaxed arm: dropping it takes ρ = −0.782 → −0.733
  (p 0.0045 → 0.0158). The frozen arm is fragile — −0.527 → −0.370 (p 0.096 → 0.293). Report the relaxed
  arm as primary; the frozen arm does not stand on its own.
- **Leave-one-out, relaxed/relaxed:** ρ ranges [−0.915, −0.733] — sign-stable across all 11 deletions.
- TBP remains the highest-IQR pilot in both arms (11.78 frozen, 8.52 relaxed vs a 2.38–10.23 range).

## Per-pilot values

| pilot | frozen bend IQR (°) | relaxed bend IQR (°) | frozen own-family ΔP | relaxed own-family ΔP | n frozen | n relaxed |
|---|---|---|---|---|---|---|
| csl | 2.40 | 4.27 | +0.1028 | +0.1580 | 84 | 89 |
| egr1 | 5.25 | 6.73 | -0.0025 | +0.0219 | 98 | 96 |
| engrailed | 2.44 | 3.91 | -0.0152 | +0.0135 | 95 | 96 |
| err | 7.93 | 6.06 | -0.0025 | +0.0359 | 90 | 92 |
| ets1 | 2.91 | 3.14 | +0.1109 | +0.1104 | 100 | 98 |
| foxa | 4.07 | 6.66 | -0.0056 | +0.0301 | 91 | 90 |
| hsf | 3.49 | 4.83 | +0.0366 | +0.0508 | 97 | 98 |
| irf | 8.54 | 6.92 | +0.0332 | -0.0525 | 91 | 93 |
| nfat | 9.95 | 10.23 | +0.0151 | -0.0038 | 86 | 88 |
| runx | 1.29 | 2.38 | +0.0841 | +0.0942 | 97 | 96 |
| tbp | 11.78 | 8.52 | -0.0309 | +0.0079 | 94 | 96 |

## Recommended edits

1. **`mechanism_analysis.py`** — add `rel_bend_iqr` to the `rel` aggregation (line ~106) alongside
   `rel_bend_med`; register it in the axis list (line ~231); make `make_dnaflex.py` plot the arm-matched
   pairing. Rename `froz_bend_iqr` → `bend_iqr_frozen` so the arm is in the name and this class of error
   is visible at the call site.
2. **Ch. 17.3 / 18.2 / 18.4** — keep the result, change the wording from "the mechanism axis strengthens
   in the relaxed arm" to "DNA-bend fluctuation predicts the relaxed arm's effect better than the frozen
   arm's." Report the 2×2, not two cells.
3. **Ch. 18.1** — add the ρ = +0.927 axis-invariance finding. It strengthens the "relaxation moves the DNA
   without changing its fluctuation regime" reading and is a measured argument for a k sweep.
4. **§4 proceeds as planned.** The interface-restriction question is unaffected and is still the right
   Day-1 gate — it now tests an axis that is arm-matched and BH-surviving, which makes a collapse under
   end-exclusion more informative, not less.

## Reproduce

```bash
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /project2/rohs_102/shewchuk/TF-conformation
# inputs: analysis/analyses/dna_relax/data/pycurves_all_perstructure.csv  (cond in frozen|relaxed|crystal)
#         analysis/analyses/mechanism/data/selffam_effects.csv            (selffam==True, dna in frozen|relax)
# estimator: groupby(tf)['bend_uu'].quantile(.75) - .quantile(.25), per cond
```
