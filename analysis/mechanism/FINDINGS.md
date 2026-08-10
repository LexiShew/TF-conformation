# Mechanism analysis — induced fit vs conformational selection

**Location:** `analysis/mechanism/` · **Date:** 2026-08-01 · **Status:** provisional

Reanalysis of the project's mechanism hypothesis using the measured DNA-shape data now
on disk, replacing the curated literature-derived deformation labels. Everything here is
computed from pipeline outputs; nothing is carried forward from prose.

---

## Headline

**The augmentation signal is real, but it is not where the benchmark was looking.**

1. The 130-entry cross-benchmark ΔPearson averages the effect away. The signal lives in
   each pilot's **own family**: ETS1 own-family ΔPearson = **+0.111** (p = 0.013,
   0 of 6 seeds negative) against a cross-benchmark value of **−0.002**.
2. It is **not** explained by baseline headroom (ρ = −0.27, p = 0.42) or by own-family
   subset size (ρ = −0.32, p = 0.34).
3. The discriminating DNA variable is **dynamic, not static**. Ensemble DNA-bend IQR
   tracks the own-family effect at ρ = **−0.53**; crystal bend — the static deformation
   the original hypothesis was written on — gives ρ = **+0.04** (p = 0.92).

The mechanism hypothesis survives, but in a modified form: what matters is not how
deformed the bound duplex is, but **how much DNA shape fluctuates across the ensemble**.

---

## What was tested

### The static induced-fit index fails

A measured replacement for the curated `dna_deform` labels was built from crystal DNA
geometry: bend deviation and minor-groove deviation from canonical B-DNA (5.7 Å,
0°), z-scored and averaged.

Against the cross-benchmark ΔPearson it gives Spearman ρ = **+0.25** — the *wrong sign*
relative to the hypothesis. A median split on it groups the pilots as:

| Regime (measured, static) | Pilots | mean ΔPearson |
|---|---|---|
| conformational-selection | csl err foxa hsf irf nfat | −0.021 |
| induced-fit | egr1 engrailed ets1 lef1 runx tbp | −0.002 |

This puts ETS1 and TBP **in the same group**, which is precisely the pair the hypothesis
uses to separate the regimes. Static crystal deformation is not the discriminating axis.
Welch p = 0.085, Mann-Whitney p = 0.13 — and in the direction opposite to the prediction.

### The effect is in the own-family subset

The cross-benchmark average pools ~120 entries from families the pilot's ensemble says
nothing about. Restricting to each pilot's own family (motif-level assignment,
seed-paired, seed as the unit of replication):

| pilot | own-family ΔPearson | 95% CI | p | neg. seeds | cross-benchmark ΔPearson |
|---|---|---|---|---|---|
| ets1 | **+0.111** | [+0.035, +0.187] | **0.013** | 0/6 | −0.002 |
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

lef1 is excluded — it is the only HMG-box in the benchmark and has no same-family entries.

### The within-pilot contrast (the test fig6 never ran)

`REVIEW_figure_scripts.md` flagged that fig6's "augmentation is family-specific" headline
was never statistically tested. Tested here at the seed level as own-minus-other within
each pilot:

| pilot | own | other | contrast | p | q (BH) |
|---|---|---|---|---|---|
| csl | +0.103 | −0.043 | **+0.145** | 0.008 | 0.092 |
| ets1 | +0.111 | −0.008 | **+0.119** | 0.011 | 0.058 |
| runx | +0.084 | −0.003 | +0.087 | 0.148 | 0.406 |
| nfat | +0.015 | −0.033 | +0.048 | 0.257 | 0.565 |
| irf | +0.033 | −0.015 | +0.048 | 0.274 | 0.502 |
| hsf | +0.037 | +0.002 | +0.034 | 0.432 | 0.594 |
| egr1 | −0.003 | −0.033 | +0.030 | 0.053 | 0.194 |
| err | −0.003 | −0.025 | +0.023 | 0.283 | 0.445 |
| foxa | −0.006 | −0.002 | −0.004 | 0.899 | 0.899 |
| engrailed | −0.015 | −0.003 | −0.012 | 0.509 | 0.622 |
| tbp | −0.031 | −0.015 | **−0.017** | 0.681 | 0.749 |

8 of 11 pilots show a positive contrast. csl and ets1 reach raw p < 0.05; **nothing
survives BH correction across 11 tests.** The direction is consistent; the significance
is not established.

### Confounds excluded

| alternative | test | result |
|---|---|---|
| A1 headroom — low baseline has more room | own-family effect vs baseline own-family Pearson | ρ = −0.27, p = 0.42 |
| A1b headroom (MAE) | vs baseline own-family MAE | ρ = −0.06, p = 0.85 |
| A2 subset size — small n manufactures large deltas | vs own-family entry count | ρ = −0.32, p = 0.34 |

An OLS of own-family effect on DNA-bend IQR + baseline headroom + subset size (n = 11,
R² = 0.52) leaves the DNA axis with the largest standardized coefficient
(β = −0.021) but no term significant at n = 11.

### Fragility

Bootstrap over pilots (10,000 resamples) for the headline correlation:

- observed ρ = −0.53, 95% CI **[−0.91, +0.20]**, same sign in **93%** of resamples.

The CI crosses zero. This is a directionally consistent trend at n = 11 pilots, not an
established correlation. Leave-one-out is sign-stable for every axis reported.

---

## Family-level effects (seed-level, frozen pipeline)

| family | ΔPearson | 95% CI | p |
|---|---|---|---|
| IRF | +0.010 | [−0.019, +0.039] | 0.40 |
| ETS | +0.008 | [−0.013, +0.028] | 0.40 |
| CSL/RBPJ | +0.006 | [−0.014, +0.026] | 0.47 |
| Runt | −0.000 | [−0.030, +0.029] | 0.98 |
| Forkhead | −0.003 | [−0.038, +0.032] | 0.83 |
| HSF | −0.006 | [−0.018, +0.005] | 0.21 |
| Rel/NF-κB | −0.011 | [−0.049, +0.028] | 0.51 |
| MADS-box | −0.016 | [−0.033, +0.002] | 0.07 |
| C2H2 zinc finger | −0.016 | [−0.027, −0.004] | **0.017** |
| Nuclear receptor | −0.017 | [−0.033, +0.000] | 0.05 |
| TBP / β-saddle | −0.019 | [−0.048, +0.009] | 0.14 |
| Homeodomain | −0.020 | [−0.050, +0.011] | 0.16 |
| bZIP | −0.028 | [−0.059, +0.003] | 0.07 |
| bHLH | −0.035 | [−0.077, +0.007] | 0.09 |

The three net-positive families (IRF, ETS, CSL/RBPJ) are exactly the three pilots with
the largest own-family effects. C2H2 zinc finger is the only family whose *negative*
effect reaches p < 0.05.

---

## A bug this analysis found

The benchmark CSVs label the DNA conditions **`frozen` and `relax`** — not `relaxed`.
An earlier version of `mechanism_analysis.py` filtered on `"relaxed"`, matched zero rows,
and silently produced an empty relaxed arm while 119 relaxed rows sat in the input. No
error was raised.

Fixed: the condition labels are declared as module constants and **asserted at load
time**, so a mismatch raises rather than dropping an arm. With the relaxed arm present:

- **foxa relaxed: ΔPearson −0.036, p = 0.0004** — the only pilot×condition effect that
  survives at the seed level.
- tbp relaxed: −0.039, p = 0.051.

Both were invisible in the frozen-only analysis.

**Anyone reusing these CSVs should check this.** Worth auditing whether any delivered
figure filtered on the wrong label.

---

## Pilot coverage

`pilot_audit.py` prints the pilot set at every stage. Summary:

- **13** base pilot configs
- **12** benchmarked pilots (dux4 has 0 fnat survivors → no eval JSON)
- **11** pilots in the own-family analysis (lef1 has no same-family benchmark entries)
- MGW-FL tables are 12 (irf drops for a groove-parse artifact)

Any *n* in a mechanism figure traces to this table.

---

## What this means for the hypothesis

The original framing — *augmentation helps under conformational selection, hurts under
induced fit on the DNA* — is **not supported in its static form**. Crystal deformation
does not predict the sign, and it groups ETS1 with TBP.

The modified form that the data support:

> Augmentation helps a pilot's own family when the pilot's DNA is **conformationally
> stable across the ensemble**, and stops helping as ensemble DNA-shape fluctuation
> grows. The model holds DNA geometry fixed, so frames whose DNA shape wanders carry
> geometry the shared crystal PWM label no longer describes.

This keeps the mechanistic logic (the fixed-DNA assumption is the ceiling) while
replacing the axis with a measured, dynamic one. It also re-enters the shared-label
confound: fluctuating DNA is exactly the condition under which one crystal PWM label
mislabels its frames most severely.

---

## Caveats

- n = 11 pilots. The headline bootstrap CI crosses zero.
- Nothing survives multiple-comparison correction.
- Own-family subsets are 3–24 entries; the smallest are noisy.
- Family effects use 6 seeds; several families are represented by a single pilot, so
  "family effect" and "pilot effect" are partly confounded for those rows.
- The dynamic axis (bend IQR) is measured on the frozen pipeline. Whether it predicts the
  *relaxed* arm's effect is untested here.

---

## Reproduce

```bash
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /project2/rohs_102/shewchuk/TF-conformation
python analysis/mechanism/mechanism_analysis.py   # table + axis correlations + family effects
python analysis/mechanism/mechanism_probe.py      # static vs dynamic, own-family, family shape
python analysis/mechanism/mechanism_confound.py   # headroom, subset size, bootstrap, contrast
python analysis/mechanism/pilot_audit.py          # pilot coverage at every stage
python analysis/mechanism/make_mechanism_fig.py   # figures/M2_where_the_signal_lives.png
```

Scripts run in order; each writes to `analysis/mechanism/data/` and later scripts read
earlier outputs.


---

# Update — relaxed DNA as a co-equal condition (2026-08-01, later)

The analysis above treated the frozen-DNA pipeline as primary. Re-running with the
relaxed-DNA arm as a co-equal condition changes the conclusion.

## The mechanism prediction is met at the own-family level

`mechanism_and_roadmap.md` §2 stated a falsifiable prediction: *a flexible-DNA pipeline
should recover exactly the families augmentation hurts today.* On the 130-entry
cross-benchmark average it appeared to fail. On the own-family effect it is met.

| pilot | frozen ΔP | p | neg/seeds | relaxed ΔP | p | neg/seeds |
|---|---|---|---|---|---|---|
| ets1 | +0.111 | **0.013** | 0/6 | +0.110 | **0.017** | 0/6 |
| csl | +0.103 | 0.065 | 2/6 | **+0.158** | **0.004** | **0/6** |
| runx | +0.084 | 0.161 | 2/6 | +0.094 | **0.024** | 1/6 |
| hsf | +0.037 | 0.331 | 2/5 | +0.051 | 0.228 | 1/5 |
| irf | +0.033 | 0.471 | 2/5 | −0.053 | 0.376 | 3/5 |
| nfat | +0.015 | 0.500 | 2/5 | −0.004 | 0.915 | 3/5 |
| egr1 | −0.003 | 0.911 | 4/6 | +0.022 | 0.291 | 1/6 |
| err | −0.003 | 0.924 | 3/6 | +0.036 | 0.629 | 2/6 |
| foxa | −0.006 | 0.877 | 3/6 | +0.030 | 0.453 | 2/6 |
| engrailed | −0.015 | 0.670 | 3/6 | +0.014 | 0.771 | 3/6 |
| **tbp** | **−0.031** | 0.513 | 3/6 | **+0.008** | 0.888 | 2/6 |

**8 of 11 pilots shift upward under relaxation.** csl and runx reach raw p < 0.05 only in
the relaxed arm. TBP — the induced-fit exemplar whose ~80° kink motivated the whole
hypothesis — crosses from negative to marginally positive.

Caveat: nothing survives BH correction across 11 pilots, and this rests on 4–6 seeds.

## The dynamic DNA axis strengthens in the relaxed arm

Own-family ΔPearson vs measured DNA geometry, Spearman with bootstrap 95% CI over pilots
(10,000 resamples):

| axis | arm | n | ρ | boot 95% CI | same sign | p |
|---|---|---|---|---|---|---|
| **Ensemble DNA-bend IQR** | **relaxed** | 11 | **−0.76** | **[−0.96, −0.28]** | **99.6%** | **0.006** |
| Ensemble DNA-bend IQR | frozen | 11 | −0.53 | [−0.91, +0.22] | 92.9% | 0.096 |
| Interface MGW-FL | relaxed | 10 | −0.37 | [−0.92, +0.42] | 83.7% | 0.293 |
| Interface MGW-FL | frozen | 10 | −0.16 | [−0.71, +0.62] | 68.9% | 0.651 |
| Crystal DNA bend (static) | relaxed | 11 | −0.07 | [−0.80, +0.66] | 57.4% | 0.832 |
| Crystal DNA bend (static) | frozen | 11 | +0.04 | [−0.65, +0.76] | 53.6% | 0.915 |

**Bend IQR × relaxed is the only correlation in this analysis whose bootstrap CI excludes
zero.** The static crystal-deformation axis is flat in both arms — confirming that the
discriminating variable is DNA-shape *fluctuation across the ensemble*, not how deformed
the bound duplex is.

## Corrections to the section above

- The earlier claim that relaxation "does not rescue" the effect was an artifact of
  reading the cross-benchmark average and of the `relaxed`/`relax` label bug.
- An earlier draft caption asserted every bootstrap CI crosses zero. That is false once
  the relaxed arm is included (bend IQR × relaxed does not). Corrected in `make_dnaflex.py`.

## Provenance gap this work closed

`analysis/dna_relax/figures/samefamily_augmentation_effect.png` and
`crossfamily_augmentation_effect.png` had **no generating script anywhere in the repo**
(verified by grep across `analysis/` and `scripts/`). They covered 5 pilots and could not
be regenerated or extended. `analysis/mechanism/make_samefamily.py` supersedes both:
11 pilots, both arms co-equal, palette-threaded, seed-level statistics.

`analysis/dna_relax/iface_mgwfl.py` (which does have a script) is **extended, not
duplicated**, by `make_dnaflex.py`: it added the relaxed arm and the own-family effect to
the cross-benchmark frozen-only correlation it already computed.

## Colour

All three figures now import `palette.py` at the repo root — GREY = baseline/reference,
TEAL = augmented·frozen DNA, GREEN = augmented·relaxed DNA, ALARM = annotation only.
An earlier version of `make_mechanism_fig.py` invented its own hues and used a pink as a
data-series colour for negative bars, which both broke per-entity threading and misused
the reserved alarm hue. Fixed; no literal hex remains in the mechanism scripts.

## Figures

| file | what |
|---|---|
| `figures/M2_where_the_signal_lives.png` | own vs other family, contrast, headroom control, dynamic-vs-static axis |
| `figures/M3_samefamily_both_arms.png` | own-family effect both arms, family specificity, frozen→relaxed shift |
| `figures/M4_dnaflex_vs_effect.png` | DNA fluctuation vs effect, both arms, with static control |

## Reproduce (order matters)

```bash
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /project2/rohs_102/shewchuk/TF-conformation
python analysis/mechanism/mechanism_analysis.py
python analysis/mechanism/mechanism_probe.py
python analysis/mechanism/mechanism_confound.py
python analysis/mechanism/pilot_audit.py
python analysis/mechanism/make_samefamily.py    # writes data/samefamily_both_arms.csv
python analysis/mechanism/make_dnaflex.py       # reads it
python analysis/mechanism/make_mechanism_fig.py
```
