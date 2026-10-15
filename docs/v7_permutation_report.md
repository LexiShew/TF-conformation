# v7 loose end — the IQR × √L anomaly, resolved to a sharper problem

**Run:** 2026-09-29, from artifacts. **The permutation test does not dismiss the anomaly — it confirms
it (p = 0.0054 after correcting for the transform search). R1's retraction should stand, but its stated
reason is wrong, and the correct reason is more interesting.**

v7 costed this at ten minutes as a cleanup. It is not cleanup.

---

## 1. The permutation test

v7's concern was that a two-parameter search at n = 10 would find a strong combination by chance. Tested
directly: nine axis transforms, 20,000 shuffles of the ΔP labels, recording the **best-of-family** ρ
under each shuffle.

| transform | ρ vs relaxed ΔP | p | ρ vs duplex length |
|---|---|---|---|
| IQR x sqrt(L) | -0.867 | 0.0012 | -0.654 |
| IQR x L | -0.855 | 0.0016 | -0.377 |
| IQR x L^1.5 | -0.794 | 0.0061 | -0.154 |
| IQR | -0.745 | 0.0133 | -0.827 |
| IQR x L^2 | -0.745 | 0.0133 | -0.012 |
| IQR / sqrt(L) | -0.745 | 0.0133 | -0.827 |
| IQR + z(L) | -0.697 | 0.0251 | +0.105 |
| IQR / L | -0.673 | 0.0330 | -0.858 |
| L | +0.451 | 0.1911 | +1.000 |

- naive p for IQR × √L, ignoring the search: **0.0008**
- **family-corrected p (best of 9): 0.0054**
- by chance, searching nine transforms at n = 10 typically yields ρ ≈ **−0.31** (5th percentile −0.709)

**So the search explanation fails.** −0.867 is well outside what the search space produces from
permuted labels. Leave-one-out is also tight and sign-stable: **[−0.917, −0.817]** across all ten
pilots — no single pilot carries it.

And the anomaly is not length in disguise. **IQR × √L is *less* length-correlated than raw IQR**
(ρ = −0.654 vs −0.827), yet correlates *more* strongly with ΔP. Multiplying by √L moves the axis away
from length and toward the effect. A pure length artifact cannot do that.

## 2. Which creates a real conflict

Two axes, both measured on the same ten pilots, both reproducible across arms, giving opposite answers:

| axis | arm-arm reliability | frozen ρ vs ΔP | relaxed ρ vs ΔP | rank-R² |
|---|---|---|---|---|
| whole-duplex `bend_uu` IQR × √L | +0.903 | −0.564 (p = 0.090) | **−0.867 (p = 0.001)** | **0.751** |
| twist-aware length-free window IQR | +0.903 | +0.418 (p = 0.229) | +0.224 (p = 0.533) | 0.050 |

Both are **sign-consistent across arms** — the whole-duplex axis negative in both, the window axis
positive in both. Neither is behaving like noise. And they are **orthogonal to each other**
(ρ = −0.055, p = 0.83), so they are not two estimates of one quantity with one of them noisier.

Partialling one out of the other does not help: IQR × √L vs ΔP controlling the window axis is
**−0.874 (p = 0.002)**. The window axis contributes nothing (R² 0.751 → 0.776 when added).

## 3. What R1 actually established, and what it did not

The v6 report concluded the mechanism correlation "does not survive" length removal. Re-reading the
evidence against this test, that conclusion overreached in one specific way.

**What is solid:** the central 4–10 base-pair steps of these duplexes carry **no** signal about ΔP.
That is measured cleanly, with a reliable instrument, and it does retire the claim *as stated* — the
mechanism story was about local DNA flexibility, and local flexibility is where the signal is not.

**What does not follow:** that the whole-duplex signal was therefore length plus noise. The two axes
measure different physical quantities. `bend_uu` is a **terminal axis-vector** angle — it sees
end-to-end orientation, dominated by motion at the duplex termini. The central window is **blind to the
ends by construction**. So the correct reading is not "a cleaner instrument saw nothing" but **"two
instruments looking at different parts of the molecule disagree,"** and the one looking at the ends is
the one that correlates.

I tested whether the signal localises to the ends and **could not settle it**:

| W steps | arm | mean fraction of duplex spanned | ρ vs ΔP | p |
|---|---|---|---|---|
| 4 | frozen | 0.32 | +0.309 | 0.385 |
| 4 | relaxed | 0.32 | -0.067 | 0.855 |
| 6 | frozen | 0.48 | +0.527 | 0.117 |
| 6 | relaxed | 0.48 | +0.273 | 0.446 |
| 8 | frozen | 0.64 | +0.418 | 0.229 |
| 8 | relaxed | 0.64 | +0.224 | 0.533 |
| 10 | frozen | 0.78 | +0.233 | 0.546 |
| 10 | relaxed | 0.78 | +0.317 | 0.406 |

If the signal lived at the termini, widening the window toward them should drive ρ negative. The trend
is **−0.400 in the frozen arm and +0.800 in the relaxed arm** — contradictory, and neither significant.
An explicit "end contribution" (whole-duplex minus rescaled central) gives ρ = −0.830 (p = 0.003), but
**that test is near-circular**: because the two axes are orthogonal, subtracting one from the other
mostly returns the other, and the construct correlates with whole-duplex IQR at +0.685. It is not
independent evidence and should not be quoted.

## 4. Recommendation

**Keep R1 retracted. Change the reason, and add the conflict to the "volunteer" list.**

- **Do not say** "the mechanism correlation was duplex length." The permutation test contradicts it, and
  an audience member who recomputes will find −0.867 at a corrected p = 0.005.
- **Do say**: the claim rested on a whole-duplex terminal-angle measure that is 69% duplex length; a
  length-free measure of *local central* flexing shows nothing; the two are orthogonal, and **which of
  them is the biologically meaningful quantity is unresolved at n = 10**. Retract on *interpretability*,
  not on absence of signal.
- **Still do not report the +0.527 best cell** — that recommendation was right and is unaffected.

This is a better retraction than the one in v6. It is also more defensible: "we found a correlation we
cannot interpret and will not build on" survives questioning in a way that "it was an artifact" does
not, when the artifact explanation is testable and fails.

**One consequence for F5.** The negative claim enumerates axes that fail to predict ΔP. The
whole-duplex axis **does** predict ΔP at a corrected p = 0.005, so F5's enumeration must say
*length-free* structural properties, and must list the whole-duplex axis as the explicit exception with
its interpretability problem attached. As written, F5 is contradicted by a number in the same document.

### What would settle it (not scheduled, but cheap and specific)

**Terminal-window bend.** The same twist-aware accumulation applied to the *first and last* 4 steps
rather than the central 8. If end-window IQR reproduces the −0.87 while central-window IQR stays flat,
the signal is terminal motion and becomes interpretable — possibly as fraying or restraint edge effects
rather than biology, which would be a restraint-protocol finding, not a mechanism. One extraction pass,
reusing `extract_v6_twist.py` with a different window selector. **This is a better use of half a day
than §3.3's sequence-viability gate**, because it can retire or rehabilitate a claim rather than scope
a future experiment.

## Reproduce

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

Inputs are artifacts, not cluster files: `v6_prescreen_scaling.csv` (per-pilot IQR, ΔP, duplex length)
and `v6_window_bend_twistaware.csv` (twist-aware window IQR). Permutation: shuffle `dP_relax` labels,
recompute Spearman ρ for all nine transforms, keep the minimum; the family-corrected p is the fraction
of shuffles whose best-of-family ρ is ≤ the observed −0.867. Seed 20260929, 20,000 shuffles.
