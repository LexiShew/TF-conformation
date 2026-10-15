# v3 execution report — priority items 1–3 plus §2.2 and the multiple-testing family

**Run:** 2026-09-29, all from data already on disk. No new training, no new structure generation.

**Bottom line.** §3.1's challenge to my addendum was **correct and I was wrong**; the headline nevertheless
survives cleaning. §3.5 is **not computable** as written. §2.2 comes back in its **third outcome row**
(drop the idea). And the honest multiple-testing family **removes q < 0.05 from the mechanism claim
entirely** — which changes what Claim 4 can say.

---

## 1. §3.1 — irf. The challenge was right; my justification was not.

v3 wrote: *"'Quartiles are robust to the offset' assumes the failure IS an offset. If `bend_uu` is computed
on a mis-assigned helical axis, both location and scale are suspect."* That is the correct objection, and
my addendum asserted robustness without testing it. Here is the test that settles it.

**The frozen arm is the control.** In the frozen condition the DNA backbone is restrained to the reference
crystal (`stage3_minimize.py` restrains P and C1'), so median frozen-arm bend **must** reproduce the
crystal's bend. Deviations measure parse error, not biology:

| pilot | median frozen bend | crystal bend | deviation |
|---|---|---|---|
| ets1 | 19.3° | 19.2° | **0.1°** |
| runx | 7.3° | 8.0° | 0.7° |
| nfat | 12.1° | 10.9° | 1.2° |
| … | | | … |
| dux4 | 28.6° | 9.5° | 19.1° |
| **irf** | **170.7°** | **16.9°** | **153.9°** |

Median deviation excluding irf: **3.8°**. irf is **41× that**. So irf's bend is broken, confirmed against a
physical constraint rather than assumed.

**Is it a recoverable offset?** No. `bend_pp` — the independent point-to-point estimator in the same file —
also reads **171.4°** for irf against a crystal `bend_pp` of 19.1°. Both estimators fail together, in both
arms, across **all 184 states**. A sign flip or a 180° wrap would show up in one channel or as a
bimodal distribution; this is a uniformly wrong helical-axis assignment. **irf's bend must be dropped, not
corrected** — matching the treatment its MGW-FL channel already received.

**Your range assertion needs rewriting.** A ">90° from the pilot median" rule **does not catch irf** — its
entire distribution is shifted, so the median moves with it. It flags 34 isolated states in
dux4/err/foxa/nfat/runx and misses the one pilot that matters. **Assert against the crystal bend instead**
(|state − crystal| > 90°): that flags all 184 irf states plus the same 34 isolated frames.

### The headline, recomputed on the cleaned axis (irf dropped, 34 outlier states removed)

| axis arm | effect arm | matched | ρ | p | boot 95% CI | excl. 0 | LOO range |
|---|---|---|---|---|---|---|---|
| frozen | frozen | **yes** | -0.564 | 0.0897 | [-0.925, +0.205] | no | [-0.783, -0.400] |
| frozen | relaxed | no | -0.721 | 0.0186 | [-0.985, -0.132] | yes | [-0.883, -0.633] |
| relaxed | frozen | no | -0.527 | 0.1173 | [-0.925, +0.231] | no | [-0.733, -0.383] |
| relaxed | relaxed | **yes** | -0.745 | 0.0133 | [-1.000, -0.170] | yes | [-0.933, -0.650] |

**Arm-matched relaxed: ρ = −0.745 (p = 0.0133), CI [−1.000, −0.170], still excludes zero**, versus −0.782
raw. Cleaning also changes nfat's frozen IQR by −1.59° (its 3 outlier frames were inflating it) and err's
relaxed by −0.30°; every other pilot moves < 0.07°.

**So the claim survives, but v3 was right that it had to be checked and right that the without-irf number
belongs in the talk rather than a footnote.**

---

## 2. §3.6 multiple testing — this is the finding that changes Claim 4

v3: *"BH over four cells is not the honest family… roughly ten [tests]. Bend IQR probably still survives;
q = 0.0125 is not the number to quote."* Correct on the first two points, **wrong on the third.**

Assembling the real family — every axis in §17.2–§17.3 (free-state spread, protein reachability, static
crystal bend, interface MGW-FL in both arms, bend median, bend IQR in both arms) × both effect arms, on the
cleaned n = 10 panel — gives **16 tests**:

| axis | effect | ρ | p | BH q |
|---|---|---|---|---|
| bend IQR (relax) | relaxed | -0.745 | 0.0133 | 0.1486 |
| bend IQR (froz) | relaxed | -0.721 | 0.0186 | 0.1486 |
| bend IQR (froz) | frozen | -0.564 | 0.0897 | 0.4692 |
| bend IQR (relax) | frozen | -0.527 | 0.1173 | 0.4692 |
| protein reachability | frozen | +0.285 | 0.4250 | 0.9602 |
| free-state spread | frozen | -0.176 | 0.6272 | 0.9602 |
| free-state spread | relaxed | -0.176 | 0.6272 | 0.9602 |
| interface MGW-FL (froz) | relaxed | -0.164 | 0.6515 | 0.9602 |

**Nothing survives q < 0.05. Bend IQR (relaxed/relaxed) reaches q = 0.149.**

This is not a small revision. The addendum's "q = 0.0125, both relaxed cells survive BH" was computed over
a 4-test family that was chosen *after* seeing which axis won — which is the selection problem v3's Part 5
identifies, applied to the correction itself. On the honest family, bend IQR is the **best of sixteen
tests** at an uncorrected p = 0.013, and that is exactly what you would expect from sixteen tests on eleven
points if nothing were there.

**Claim 4 should therefore drop "the only correlation whose bootstrap CI excludes zero" as a
significance claim.** What remains defensible: the axis is the largest effect among the axes tested, it
survives every single-covariate control, and it does not survive correction for having tested sixteen
things. That is a hypothesis-generating observation, and it makes the k sweep the load-bearing experiment
rather than a confirmation of something already established.

---

## 3. §3.5 — not computable. All 12 pilots are in distinct families.

v3 costs this at "an hour" and calls it possibly "a genuinely strong result." The design cannot produce it:
per `analysis/data/mechanism_apo_holo.csv`, the twelve pilots occupy **twelve different families** (CSL/RBPJ,
C2H2, Homeodomain, Nuclear receptor, ETS, Forkhead, HSF, IRF, HMG-box, Rel/NFAT, Runt, TBP/β-saddle). Of
the **132 ordered donor–target pairs, exactly 0 are same-family.**

(Note `analysis/data/family_annotation.csv` looks like it answers this but does not — it maps the 130
benchmark *entries* per pilot, not pilot identity, and only 6 pilots appear in it at all. Using it yields a
3-pilot family map and a spurious answer; I hit this first.)

The nearest computable proxy, grouping by broad structural class (HTH-like / Zn / β-scaffold / HMG):

| metric | same class | other class | MWU p |
|---|---|---|---|
| Δ vs experiment | +0.0017 (n = 34) | +0.0183 (n = 98) | 0.75 |
| Δ vs crystal | +0.0231 (n = 34) | +0.0250 (n = 98) | 0.87 |

Flat. So the +0.025 cross-pilot gain is **generic, not family-mediated** — which is the "less interesting"
branch of v3's own outcome table, but it is also a real answer: it means the own-family effect in Ch. 16
is *not* reproduced in a leakage-free design. Worth one line in the talk.

---

## 4. §2.2 — crystal-free dispersion. Third outcome row: drop it.

Implemented as specified: mean pairwise (1 − Pearson) among the state PWMs themselves, from the
`pwms_<tf>_s<n>.npz` stacks, no crystal reference.

Two implementation notes. The stacks are **NaN-padded to the experimental motif width** — only positions
where a state's prediction overlaps the motif are populated (tbp_s1: 54% NaN, 0 of 94 rows fully finite).
Requiring positions covered in *every* state drops 5 pilots; computing each pairwise distance on that
pair's own overlap keeps all 12 (median overlap 5–10 positions, 93–100% of pairs usable). **The strict
version gives ρ = +0.94 against ΔP on n = 6 — and that is a small-sample artifact.** On the full panel:

| test | ρ | p | n |
|---|---|---|---|
| dispersion vs bend IQR (relaxed) | −0.079 | 0.83 | 10 |
| dispersion vs own-family ΔP (relaxed) | +0.309 | 0.38 | 10 |
| dispersion vs own-family ΔP (frozen) | +0.552 | 0.098 | 10 |
| **dispersion vs mean r-to-crystal** | **−0.982** | **8.4e-08** | **11** |

| pilot | dispersion | mean r to crystal | states |
|---|---|---|---|
| tbp | 0.214 | 0.812 | 94 |
| engrailed | 0.248 | 0.835 | 85 |
| egr1 | 0.436 | 0.600 | 90 |
| foxa | 0.448 | 0.574 | 70 |
| ets1 | 0.502 | 0.543 | 99 |
| csl | 0.643 | 0.399 | 57 |
| runx | 0.659 | 0.433 | 35 |
| err | 0.692 | 0.375 | 41 |
| lef1 | 0.730 | 0.349 | 12 |
| nfat | 0.747 | 0.325 | 17 |
| hsf | 0.771 | 0.272 | 96 |

**It is null against both targets and almost perfectly redundant with `r_mean`** — a quantity already on
disk. Partialling `r_mean` out leaves dispersion vs ΔP at ρ = −0.021.

v3 correctly noted my earlier test was crystal-anchored and therefore circular. The crystal-free version is
not circular and **still lands in the same place**: a pilot whose ensemble disagrees internally is a pilot
whose states agree poorly with the crystal, and the crystal reference was never what made these two
redundant — ensemble spread and mean accuracy are coupled through the same underlying quantity.

**Consequence for Part 1: the k sweep cannot use an inference-only endpoint.** Per v3 §1.3's own
contingency, it runs **with ΔP at levels {0, 5.0} — two arms, not four** — losing monotonicity and keeping
the interaction test. That is the discriminating prediction anyway (§1.2 point 2), so the design survives;
it just costs two arms per pilot and cannot deliver a dose–response curve.

---

## 5. Revised priorities

| was | item | now |
|---|---|---|
| 1 | §3.1 irf | **done** — drop irf's bend; crystal-referenced assertion, not median-referenced |
| 2 | §2.2 dispersion | **done — negative.** k sweep goes to 2 levels with ΔP |
| 3 | §3.5 family structure | **done — not computable** (0 same-family pairs); proxy is flat |
| — | §3.6 multiple testing | **done — and it demotes Claim 4.** Do this before writing the talk |
| 4 | §1.1 k-sweep calibration | **now priority 1** — unchanged in design, and now the only route to a causal claim |
| 5 | §2.1 size-matched donor | **priority 2** — unchanged; still the highest-stakes single result |
| 6 | §3.2 interface restriction | **priority 3** — and now it is also the last chance to find the axis is end-fraying before spending training arms |
| 7 | §3.3 HSF geometry check | **cheap and warranted** — see below |

**On §3.3 (HSF):** the frozen-vs-crystal control clears HSF of irf's failure mode — deviation 11.0°, sixth
largest and an order of magnitude below irf. Its bend geometry is not broken. So HSF's three anomalies are
not a parse artifact and need a real explanation; the variance exception (0.386 → 0.519) is the one worth a
slide. Still worth running the range assertion, but do not expect it to fire.

**Re-running step 2's slopes after cleaning** (v3 §3.1 action 4) is worth doing but should not be expected
to rescue anything: irf contributed 5 of 701 states, and the 34 outlier frames are spread across five
pilots whose slopes were already mixed in sign.

**§3.4's unexplained +0.59** is unaffected by cleaning — irf was already excluded from that analysis for
having too few states. It stays open.

---

## 6. What the talk can now claim

Claims 1–3 are untouched; all three rest on n = 701 states or on three independently-measured points, not
on the pilot-level regression.

**Claim 4 needs rewriting one step further than v3 proposed.** v3's draft already says "n = 11, selected
from ~10 axes, so the effect size is an upper bound." Add the correction result:

> Ensemble DNA-bend fluctuation is the largest of sixteen axis×arm correlations tested on eleven pilots
> (ρ = −0.745 arm-matched, cleaned; uncorrected p = 0.013). It survives every single-covariate control but
> **does not survive correction for the family of axes tested (BH q = 0.149)**, and it does not operate
> through individual states (7/11 slopes, n.s.). It is a hypothesis, and the k sweep is its test.

That is a weaker claim than v3's draft, and it is the one the data supports. Volunteering the failed
correction is the same move as volunteering step 2's null — and for the same reason.

---

## Reproduce

All inputs on disk; no new compute.
- irf diagnosis: `pycurves_all_perstructure.csv`, compare `cond=="frozen"` median `bend_uu`/`bend_pp` per
  pilot against the `cond=="crystal"` row.
- cleaned axis: drop `tf=="irf"`; drop states with |`bend_uu` − crystal `bend_uu`| > 90°; then
  `groupby(tf,cond)['bend_uu'].quantile(.75)-quantile(.25)`.
- family structure: `crosspilot_deltas.csv` × `mechanism_apo_holo.csv` `family` column (**not**
  `family_annotation.csv`).
- dispersion: `pwms_<tf>_s<n>.npz` key `states`, shape (n_states, motif_len, 4), NaN-padded; compute each
  pairwise correlation on that pair's mutually-finite positions.
