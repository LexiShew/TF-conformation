# v6 report — §3.1 and §3.2

**Run:** 2026-09-29, from data on disk. **§3.1 lands in v6's second outcome row: the length-free axis
collapses, and Claim 5 should come out of the talk. §3.2 passes, which makes Claim 4 fully supported.**

Both items also turned up method problems that had to be fixed before the results meant anything — one
in my own construction of the window axis.

---

## §3.1 — Interface restriction / fixed-width window

### The ten-minute pre-check, run first as v6 suggested

| axis form | ρ vs relaxed ΔP | p | controlling L |
|---|---|---|---|
| IQR (raw) | −0.745 | 0.013 | −0.743 |
| IQR × L | −0.855 | 0.0016 | −0.828 |
| **IQR × √L** | **−0.867** | **0.0012** | **−0.847** |
| IQR / L | −0.673 | 0.033 | −0.624 |
| L alone | +0.451 | 0.191 | — |

Undoing the length dilution **strengthens** the correlation. Read alone, that argues the raw axis was
length-diluted and the true quantity is extensive — encouraging for the mechanism. It is also exactly
what you would see if the axis were mostly a length artifact and multiplying by L sharpened the length
signal, so the pre-check is suggestive and not decisive. It made the per-level work worth doing rather
than answering it.

Also confirmed: the length effect runs **opposite to worm-like-chain expectation**. WLC says bend
variance grows with contour length, so IQR should rise with L; observed ρ(L, IQR) = **−0.827**. That is
the signature of a terminal-axis-vector angle, not a polymer property — consistent with v6's own
diagnosis of `bend_uu`.

### A correction to my construction before the result

My first fixed-window implementation summed per-step roll and tilt over the central window. **That is
wrong**, and the data said so: it gave TBP a median bend of 192° over 8 steps. Successive base-pair steps
are rotated by the helical twist, so wedge angles at opposite helical phases must partially **cancel**,
not add. The correct accumulation carries cumulative twist Ω:

    Bx += roll·cos Ω − tilt·sin Ω ;  By += roll·sin Ω + tilt·cos Ω ;  bend = √(Bx² + By²)

| pilot | naive median (°) | twist-aware (°) | naive IQR | twist-aware IQR |
|---|---|---|---|---|
| tbp | 191.8 | 118.8 | 6.03 | 10.87 |
| lef1 | 101.3 | 84.5 | 26.29 | 31.45 |
| ets1 | 58.1 | 32.4 | 8.85 | 13.02 |
| dux4 | 48.5 | 37.5 | 30.59 | 26.09 |
| runx | 45.8 | 23.9 | 11.72 | 23.49 |
| egr1 | 37.5 | 18.6 | 11.51 | 8.91 |
| csl | 31.1 | 30.0 | 12.66 | 17.01 |
| foxa | 25.7 | 20.0 | 18.57 | 21.03 |
| nfat | 24.2 | 30.0 | 14.80 | 22.35 |
| engrailed | 23.3 | 17.8 | 8.19 | 7.90 |
| err | 19.9 | 25.7 | 16.56 | 17.25 |
| hsf | 19.0 | 13.6 | 8.90 | 11.74 |

Twist-aware medians span 13.6–118.8°, which is physical. Both values are in the saved table so the
correction is auditable.

**The axis is sound as a measurement.** Arm-to-arm reproducibility is ρ = +0.90 to **+0.944** across
window widths — better than whole-duplex `bend_uu` managed (+0.927) — so what follows is not a noise
result.

### The result

| W (steps) | arm | n | ρ vs **length** | p | ρ vs **ΔP** | p | ctrl L |
|---|---|---|---|---|---|---|---|
| 4 | frozen | 10 | -0.086 | 0.812 | +0.309 | 0.385 | +0.337 |
| 4 | relaxed | 10 | +0.043 | 0.906 | -0.067 | 0.855 | -0.097 |
| 6 | frozen | 10 | +0.191 | 0.596 | +0.527 | 0.117 | +0.507 |
| 6 | relaxed | 10 | +0.204 | 0.572 | +0.273 | 0.446 | +0.207 |
| 8 | frozen | 10 | +0.198 | 0.584 | +0.418 | 0.229 | +0.392 |
| 8 | relaxed | 10 | +0.154 | 0.670 | +0.224 | 0.533 | +0.175 |
| 10 | frozen | 9 | -0.077 | 0.844 | +0.233 | 0.546 | +0.263 |
| 10 | relaxed | 9 | +0.282 | 0.462 | +0.317 | 0.406 | +0.219 |

**Length dependence is gone** (|ρ| ≤ 0.28, all p > 0.4, against −0.827 for the raw axis) — the
construction worked. **And the correlation with ΔP collapses.** None of the eight window × arm cells
reaches p < 0.05 (minimum p = 0.117), and **seven of eight are positive** — the opposite sign to the
−0.745 headline. The sole negative cell is W = 4 relaxed at ρ = −0.067, i.e. indistinguishable from
zero; the sign reversal is therefore near-unanimous across window widths and arms, not marginal.

Three diagnostics confirm what happened:

1. **The two axes are unrelated.** Twist-aware window IQR vs whole-duplex `bend_uu` IQR:
   **ρ = −0.055 (p = 0.88)**; controlling length, +0.132. They are not two measures of one quantity.
   `bend_uu` IQR was almost entirely duplex length.
2. **Rank-regression R² against ΔP:** old axis alone **0.556**; length alone 0.203; length-free axis
   alone **0.050**; both together 0.228. The old axis's explanatory power does not decompose into
   length-plus-flexibility — it was a composite that the length-free version does not contain.
3. **The sign reverses.** Best length-free cell is W = 6 frozen at **+0.527 (p = 0.117)**.

### Verdict: v6's second outcome row

> **Claim 5 was duplex length plus terminal-vector noise. It comes out of the talk.**

This is v6's own pre-registered reading, and the plan's judgment that finding it now rather than after
the sweep was worth four versions of deferral is vindicated. One nuance worth stating precisely: the
length-free axis is not a *null* — it is a weak **positive** trend, n.s. at n = 10, pointing the
opposite way from the retracted claim. Do not replace one under-powered directional claim with another.

**Downstream consequences:**

- **§3.6's calibration gate loses its moderator.** The gate asked whether k separates the axis; the axis
  it was to separate no longer predicts anything. The 2-pilot calibration is still worth running as a
  within-pilot geometry measurement — it is adequately powered and tells you whether the restraint knob
  does anything at all — but it is no longer gating a mechanism test.
- **The ~20-test family shrinks in relevance, not size.** The axis that headed it is retracted; the
  remaining members were never significant.
- **Part 4's HSF and +0.59 items lose their frame.** Both were anomalies *relative to the axis*.
- **§3.5's premise is unaffected** and arguably strengthened: whether measured induced fluctuation varies
  at fixed protein is now the open question rather than a supporting one.

---

## §3.2 — Predicted MGW vs measured MGW. **Passes.**

| condition | positions | n | ρ | p | Pearson | within-pilot ρ | p |
|---|---|---|---|---|---|---|---|
| crystal | all positions | 124 | +0.409 | 2.4e-06 | +0.444 | +0.300 | 7.1e-04 |
| crystal | unpadded only | 84 | +0.386 | 2.9e-04 | +0.394 | +0.171 | 1.2e-01 |
| frozen | all positions | 136 | +0.276 | 1.2e-03 | +0.314 | +0.165 | 5.5e-02 |
| frozen | unpadded only | 77 | +0.332 | 3.2e-03 | +0.340 | +0.091 | 4.3e-01 |
| relaxed | all positions | 141 | +0.314 | 1.5e-04 | +0.360 | +0.248 | 3.1e-03 |
| relaxed | unpadded only | 78 | +0.384 | 5.1e-04 | +0.378 | +0.214 | 6.0e-02 |

**Predicted intrinsic MGW tracks measured mean MGW at ρ = +0.386 (p = 3e-4) in the crystal and +0.384
(p = 5e-4) in the relaxed ensemble**, on unpadded positions. The crystal being the strongest cell is the
right ordering — it is the least protein-perturbed structure — and within-pilot correlations stay
positive (+0.17 to +0.30), so the predictor resolves *which position* has a wide groove, not just
pilot-level offsets.

**So the plumbing is sound**: the join, position indexing, strand convention and padding filter all
work, and the predictor performs as published on this data.

> **Claim 4 can be stated.** Intrinsic and protein-induced DNA fluctuation are decoupled, and the
> decoupling is biology rather than a pipeline failure — the same pipeline recovers the intrinsic
> mean-vs-mean relationship cleanly while showing nothing for fluctuation.

This is the cleanest positive of the day and it cost an hour, exactly as v6 costed it.

---

## Where this leaves the talk

Claims 1–3 untouched. **Claim 4 is now fully supported** and is the fourth clean result. **Claim 5 is
retracted.** That is a better position than it sounds: the talk's four surviving claims are all
large-n or multiply-audited, and the one that rested on ten pilots and a confounded axis is gone before
an audience found it.

Part 5's framing needs one change. The line *"nothing we measured predicts the augmentation effect,
which makes those controls uninformative rather than reassuring"* is now the **main** mechanism
statement rather than a caveat inside Claim 5. With the length-free axis included, nothing measured —
protein reachability, free-state spread, static crystal bend, interface MGW-FL, whole-duplex or
length-free bend fluctuation — predicts ΔP at n = 10 with q < 0.05. That is a defensible and interesting
negative: **the augmentation effect's size is not explained by any structural property we can measure**,
which is itself an argument for moving the unit of analysis off the pilot.

### Suggested Part 3 reorder

| was | item | now |
|---|---|---|
| §3.1 | interface restriction | **done — Claim 5 retracted** |
| §3.2 | predicted vs measured MGW | **done — passes; Claim 4 supported** |
| §3.3 | pipeline integrity audit | **now item 1.** irf's silent failure is the same class of bug as my own twist error: both were caught by a value looking physically wrong, not by an assertion. |
| §3.4 | size-matched donor control | **item 2.** Now the highest-stakes live question — it is the only remaining threat to the surviving claims. |
| §3.5 | sequence viability gate | item 3, premise unchanged |
| §3.6 | k sweep | **descope to the 2-pilot calibration.** Report as "does restraint strength change what the model reads," with no mechanism framing. |

---

## Reproduce

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
cd /scratch1/shewchuk/deepdnashape_seqshape
python extract_v6.py          # groove MGW per position + naive window bend
python extract_v6_twist.py     # twist-aware window bend (the correct one)
```

§3.1 axis: per state, take the W central steps of `local_inter_base_pair`, accumulate roll/tilt with
cumulative twist as above, then IQR across states. Use `iqr_ta`, never `iqr_naive`.
§3.2: predicted `MGW` (layer 4) joined to `groove.data[level].sub_levels["0"].minor_width` on
`(tf, position)`, filtered to `nN == 0`.
