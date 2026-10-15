# Addendum — Deep DNAshape run, and what it settles about Part 3

**Run:** 2026-09-29. Answers v5 Part 3 question 4 with the real tool rather than the composition
proxies. **The positive control fails again — but for a different and more actionable reason than the
proxy version suggested.**

## Availability: it is already on this host

Deep DNAshape (Li, Sagendorf, Chiu, Rohs — the webserver at deepdnashape.usc.edu) is installed in the
lab space at `/project2/rohs_102/jinsenli/deepDNAshape-webserver`, the author's own copy. Critically it
ships **precomputed k-mer query tables** (`querytables/`, ~8 GB across 26 features plus `odd.txt` /
`even.txt` k-mer indices), so prediction is a table lookup — **no TensorFlow, no GPU, no model load.**
That closes the tooling gap identified in the main report: DNAshapeR does not need installing at all,
and Deep DNAshape is its successor from the same lab.

**Two operational notes.**
1. The shipped `fDeepDNAshape.py` hardcodes `/srv/www/deepdnashape/querytables/` (the webserver's
   deployment path), so it cannot be imported as-is here. The lookup contract is simple and stable
   though: line *i* of `<feature>.txt` corresponds to line *i* of `odd.txt` (base-pair features) or
   `even.txt` (base-pair-step features); layer *L* means k-mer width 2L+1 (bp) or 2L+2 (bp-step), with
   `N`-padding at the ends. My script reimplements exactly that by streaming, and touches nothing in
   the lab directory.
2. **`querytables/` is group-writable.** Any script pointed at it must open read-only. Everything from
   this run was written under `/scratch1/shewchuk/deepdnashape_seqshape/`.

## The features that matter

Deep DNAshape predicts both mean shape and **shape fluctuation** — the `-FL` features. That is the
sequence-side analogue of the pyCurves-measured `mgwfl` / `iface_mgwfl` columns already in
`selffam_vs_shape.csv`, which makes a **like-for-like positive control** possible: predicted fluctuation
against measured fluctuation, same quantity, two independent routes.

Predicted at layer 4 for all 13 pilot reference duplexes (longest DNA chain of each docked reference):

| pilot | PDB | bp | MGW-FL | Roll-FL | HelT-FL |
|---|---|---|---|---|---|
| egr1 | 1aay | 11 | 0.939 | 4.982 | 3.743 |
| lef1 | 2lef | 15 | 0.986 | 5.041 | 3.713 |
| err | 1lo1 | 13 | 0.992 | 5.105 | 3.758 |
| hsf | 5d5u | 12 | 0.998 | 4.918 | 3.777 |
| runx | 1hjc | 16 | 1.012 | 5.108 | 3.782 |
| dux4 | 5z6z | 19 | 1.019 | 5.088 | 3.734 |
| tbp | 1tgh | 12 | 1.045 | 5.408 | 3.633 |
| foxa | 1vtn | 13 | 1.047 | 5.117 | 3.732 |
| ets1 | 1k79 | 15 | 1.048 | 5.264 | 3.772 |
| csl | 3brg | 15 | 1.083 | 5.422 | 3.817 |
| irf | 1if1 | 26 | 1.087 | 5.232 | 3.801 |
| engrailed | 3hdd | 21 | 1.091 | 5.298 | 3.755 |
| nfat | 1a66 | 12 | 1.096 | 5.498 | 3.815 |

## Result 1 — the like-for-like control fails

**Deep DNAshape predicted MGW-FL vs pyCurves measured MGW-FL: ρ = +0.030 (p = 0.93)** whole-duplex,
**+0.200 (p = 0.58)** interface-restricted, Pearson +0.16 / +0.24. No relationship.

Against measured ensemble bend IQR, every predicted fluctuation feature is null:

| predicted feature | vs relaxed bend IQR | p | vs frozen bend IQR | p |
|---|---|---|---|---|
| MGW-FL | -0.042 | 0.907 | -0.164 | 0.651 |
| ProT-FL | +0.212 | 0.556 | +0.103 | 0.777 |
| Opening-FL | -0.055 | 0.881 | -0.139 | 0.701 |
| Roll-FL | +0.200 | 0.580 | +0.091 | 0.803 |
| HelT-FL | -0.345 | 0.328 | -0.479 | 0.162 |
| Tilt-FL | +0.212 | 0.556 | +0.224 | 0.533 |
| Slide-FL | -0.273 | 0.446 | -0.321 | 0.365 |

**0 of 7 reach p < 0.05.** Best is HelT-FL at −0.345 (p = 0.33). And none predicts the endpoint either
(own-family relaxed ΔP: HelT-FL +0.479 p = 0.16, everything else |ρ| < 0.27).

So the negative result from the composition proxies was **not** an artifact of the proxies being crude.
The real predictor gives the same answer.

## Result 2 — the reason, and why it changes the recommendation

The proxy run left three readings open. This run discriminates between them:

| quantity | range across the 10 pilots | fold |
|---|---|---|
| Deep DNAshape predicted MGW-FL | 0.939 – 1.096 Å | **1.17×** |
| pyCurves measured MGW-FL (relaxed) | 0.393 – 1.525 Å | 3.88× |
| pyCurves measured MGW-FL (frozen) | 0.197 – 1.487 Å | 7.56× |

**Predicted fluctuation is almost constant across these duplexes** (CV 0.049) while measured
fluctuation spans nearly fourfold (CV 0.52). A predictor with essentially no variance cannot correlate
with anything, however good it is. This is **range restriction, not predictor failure** — which is
reading 2 from the main report, not reading 3.

Scoping it against achievable space: streaming all **262,144** real-base 9-mers gives a predicted MGW-FL
distribution with mean 1.012 Å, SD 0.104 Å, and full range 0.766 – 1.580 Å. The 13 pilots span
percentile **25 to 79** — about **19% of the achievable range**. The pilot duplexes are compositionally
unremarkable, because they were selected for protein-side reasons (family coverage, monomeric DNA-reading
unit, clean crystal) with no attention to DNA sequence at all.

This is a real finding about the pilot set, and it is reassuring rather than alarming: **the sequence
axis was never varied, so its null result carries little information.**

## Consequence for Part 3

Part 3 was designed to *manipulate* predicted flexibility. The achievable contrast is large relative to
what the pilots happen to sample:

- p5 → p95: **0.340 Å**, about **2.2×** the pilot span
- p1 → p99: **0.506 Å**, about **3.2×** the pilot span

Applying a range-restriction correction to the best observed pilot-range coefficient (ρ = 0.345, itself
not significant, so treat this as an **optimistic ceiling rather than a forecast**):

| design | effective ρ | power n = 10 | n = 16 | n = 24 |
|---|---|---|---|---|
| pilot-like range (current) | 0.345 | 0.16 | 0.25 | 0.38 |
| p5 – p95 designed contrast | 0.63 | 0.50 | 0.76 | 0.92 |
| p1 – p99 designed contrast | 0.76 | 0.75 | 0.95 | 1.00 |

**Part 3 is the only experiment in the project that reaches adequate power at an affordable n**, and the
reason is that it manipulates its x instead of accepting the range nature supplied. That is the same
argument that promoted the state-level design in the v4 report, applied on the DNA side.

**But the prerequisite has not been met, only re-scoped.** Question 4 asked whether the predictor tracks
the measurement; the answer is that at the pilot duplexes it does not, and we cannot yet distinguish
"range restriction" from "predictor and measurement are different quantities." **Design power is
conditional on the first.**

## Recommended gate before costing Part 3

A cheap decisive test, one pilot, no training:

1. Predict MGW-FL for a library of 12-mers; select a **p5 pair and a p95 pair** (≈0.34 Å predicted
   contrast, ~5× the pilot-to-pilot spread).
2. Build those duplexes onto **one** pilot's existing relaxed-DNA states (the stage-3 relax machinery
   already runs per-pilot).
3. Run pyCurves and ask a single question: **does measured MGW-FL separate the p5 from the p95 group?**

If yes, the predictor tracks the measurement once range is adequate, and Part 3's power table above is
usable. If no, predicted and measured fluctuation are genuinely different quantities — reading 3 — and
Part 3's manipulation and readout are disconnected, which must be known **before** designing sequences
and buying training arms. Either way this is hours of compute on one pilot and it converts Part 3 from
a speculative experiment into a costed one.

**One caveat that survives either outcome.** Measured ensemble bend IQR is the fluctuation of a duplex
that is backbone-restrained toward a crystal while the protein varies. Deep DNAshape predicts the
fluctuation of free DNA in solution. Those need not be the same thing even with adequate range, and the
interface-restricted metric (§4.2) is the version most likely to bridge them — another reason §4.2
should move up.

## Reproduce

Read-only inputs: the lab `querytables/`, and `output/stage2_docked/<tf>/<pdb>_state_001.pdb` for the
DNA sequences. All outputs under `/scratch1/shewchuk/deepdnashape_seqshape/`.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
cd /scratch1/shewchuk/deepdnashape_seqshape && python run_deepdnashape.py
# -> deepdnashape_pilots.json : 13 pilots x 12 features, layer 4
```

Lookup contract: line *i* of `<feature>.txt` ↔ line *i* of `odd.txt` (bp) / `even.txt` (bp-step);
layer *L* → k = 2L+1 (bp) or 2L+2 (bp-step); pad both ends with `N` × L. Never open `querytables/`
for writing — it is group-writable.


---

# Gate result — per-position, not per-pilot

**Run:** 2026-09-29. **The gate does not pass.** It also came out better designed than the protocol I
proposed, because the data needed was already on disk, and the redesign is itself the main finding.

## The gate got redesigned before it ran

My proposed protocol was to design p5/p95 duplexes, rebuild them on one pilot, and re-run pyCurves —
hours of new compute. Checking the inputs first showed that unnecessary and, more importantly, that my
**pilot-level test had been the wrong test**.

The pilot-level comparison averaged predicted MGW-FL over positions, and averaging is what destroyed the
dynamic range (1.17× across pilots). But Deep DNAshape predicts **per position**, and pyCurves reports
**per position** too — `dataframes.local_base_base.propel` (ProT per bp) and
`dataframes.local_inter_base_pair.roll` / `.twist` (per step), each carrying a `sequence_index`. Taking
the IQR across states at each position gives measured fluctuation at that position, which is the exact
like-for-like partner of the `-FL` prediction.

That converts the gate from 10 pilot averages into **960 paired positional observations across all 13
pilots** (1,474 extracted rows), needs no sequence design, and asks the sharper question: *does the
predictor know which position fluctuates?*

## First pass — a significant result that is entirely an artifact

| feature | cond | n | pooled ρ | p | within-pilot ρ | p | pilots positive |
|---|---|---|---|---|---|---|---|
| ProT | frozen | 168 | -0.223 | 3.7e-03 | -0.291 | 1.3e-04 | 2/12 |
| ProT | relaxed | 168 | -0.213 | 5.6e-03 | -0.316 | 2.9e-05 | 2/12 |
| Roll | frozen | 156 | -0.082 | 3.1e-01 | -0.124 | 1.2e-01 | 4/12 |
| Roll | relaxed | 156 | -0.095 | 2.4e-01 | -0.151 | 5.9e-02 | 5/12 |
| Twist | frozen | 156 | -0.024 | 7.7e-01 | -0.014 | 8.7e-01 | 6/12 |
| Twist | relaxed | 156 | +0.010 | 9.0e-01 | +0.017 | 8.4e-01 | 7/12 |

ProT looked significant — within-pilot ρ = −0.316, p = 3e-5, 10 of 12 pilots negative (sign p = 0.039).
But the **sign is wrong**: predicted high fluctuation where measured fluctuation is low. A correctly
working predictor cannot do that, so the sign was the tell.

**Cause: position padding, with opposite-signed end trends in the two channels.**

- predicted ProT-FL vs distance-from-end: **ρ = −0.465** (rises toward the ends)
- measured ProT IQR vs distance-from-end: **ρ = +0.436** (rises toward the interior)
- predicted vs measured, controlling distance-from-end: **ρ = −0.013 (p = 0.87)** — gone
- interior positions only (≥2 bp from an end): ρ = −0.037 (p = 0.69)

And the mechanism is explicit. A position whose k-mer window runs off the duplex end gets `N`-padded,
and Deep DNAshape reports **higher** fluctuation for `N`-containing context:

| N's in window | positions | mean predicted ProT-FL |
|---|---|---|
| 0 | 96 | 5.721 |
| 1 | 26 | 5.792 |
| 2 | 26 | 5.777 |
| 3 | 26 | 6.225 |
| 4 | 26 | 8.479 |

Predicted ProT-FL vs N-count: **ρ = +0.535, p = 3e-16**. The `-FL` value at a padded position is the
predictor **expressing uncertainty about missing context**, not a claim about the DNA. Meanwhile
measured fluctuation genuinely *falls* at duplex termini here, because the frozen/relaxed restraints
anchor the ends. Two unrelated position effects with opposite signs manufactured a p = 3e-5 correlation
pointing the wrong way.

**This is a reusable gotcha, and it applies to Part 3 directly.** Only **96 of 200 pilot positions
(48%)** have a fully real 9-mer window at layer 4. Any analysis using Deep DNAshape on short duplexes
must either drop padded positions or carry N-count as a covariate — otherwise the padding gradient will
correlate with anything else that varies along the duplex.

## Second pass — unpadded windows only

Restricting to positions with a fully real 9-mer window (n = 78 per cell):

| feature | cond | ρ | p | BH q |
|---|---|---|---|---|
| ProT | frozen | +0.186 | 0.102 | 0.306 |
| ProT | relaxed | +0.224 | 0.049 | 0.295 |
| Roll | frozen | +0.035 | 0.759 | 0.915 |
| Roll | relaxed | +0.013 | 0.907 | 0.915 |
| Twist | frozen | -0.058 | 0.611 | 0.915 |
| Twist | relaxed | -0.012 | 0.915 | 0.915 |

**The sign flips positive, which is the direction a working predictor should give** — the first
encouraging signal in this whole line of work. But it does not clear the bar:

- ProT relaxed: pooled ρ = +0.224, uncorrected **p = 0.049**
- BH over the six tests: **q = 0.295** — nothing survives
- within-pilot (each pilot's level removed): ρ = +0.198, **p = 0.082**
- mixed model with pilot random intercept: β = +4.51, se = 3.12, **p = 0.149**
- per-pilot: **5 of 8 positive** (sign p = 0.73), median ρ = +0.229, but range −0.39 (ets1) to +0.69
  (engrailed)
- Roll and Twist: flat (|ρ| ≤ 0.06)

So: right sign, right order of magnitude, wrong significance — at n = 78 positions, which is not a
small-sample excuse.

## Verdict, and what it means for Part 3

**The gate fails, but the failure is now specific rather than global.** Three things changed:

1. **The pilot-level null was partly an artifact of averaging.** Per position, ProT shows a weak
   correctly-signed relationship that pilot averages destroyed. Deep DNAshape carries *some* positional
   information about ProT fluctuation in this system.
2. **Roll and Twist carry none.** That matters because the project's axis is a *bend* measure, and bend
   is built from roll and twist, not from propeller twist. The one feature showing signal is not the one
   the mechanism claim runs on.
3. **The padding artifact must be handled explicitly** in any Part 3 analysis, and it halves the usable
   positions on 11–21 bp duplexes.

**Recommendation: do not cost Part 3 as a mechanism test yet.** The design-power table earlier in this
addendum assumed the predictor tracks the measurement once range is adequate. Per position — where range
is *not* restricted — it tracks weakly for ProT and not at all for Roll/Twist. That is reading 3 for the
bend-relevant features: predicted free-DNA fluctuation and measured restrained-ensemble fluctuation are
substantially different quantities.

Two consequences worth acting on:

- **The cheap follow-up is now about longer duplexes, not designed sequences.** The padding constraint
  and the interior-vs-terminus effect both shrink with duplex length. A pilot with a 25–30 bp reference
  would have ~80% unpadded positions and a genuine interior. If the ProT signal strengthens there, Part
  3 becomes fundable; if it does not, the disconnect is real. `irf` (26 bp) and `engrailed` (21 bp) are
  the existing pilots closest to this, and `engrailed` already has the strongest per-pilot ρ (+0.69) —
  which is consistent with the length explanation and worth checking before anything is designed.
- **§4.2 (interface restriction) rises again.** Interface positions are interior by construction, so a
  footprint-restricted metric sidesteps both the padding artifact and the terminus effect. It is now the
  fix for three separate problems: duplex-length confounding, padding, and end restraint.

## Reproduce

All inputs on disk; outputs under `/scratch1/shewchuk/deepdnashape_seqshape/`.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
cd /scratch1/shewchuk/deepdnashape_seqshape
python extract_perposition.py     # -> measured_perposition_fluct.csv (1,474 rows, 13 pilots)
```

Measured fluctuation at a position = `q75 − q25` of that feature across the arm's states, from
`pycurves/<tf>/<cond>_state_*_curvesplus.json`. Join to the predictions on `(tf, feature, position)`
with ProT ← `local_base_base.propel`, Roll ← `local_inter_base_pair.roll`,
Twist ← `local_inter_base_pair.twist`, all keyed by `sequence_index`. **Filter to positions whose
k-mer window contains no `N`** (`nN = max(0, 4−(i−1)) + max(0, 4−(L−i)) == 0` at layer 4) before
interpreting any correlation.


---

# Longer-duplex check — the length explanation fails

**Run:** 2026-09-29. Tests the one hypothesis left standing after the gate: that the weak ProT signal
and the null Roll/Twist result are consequences of short duplexes, and would strengthen on longer ones.
**It does not hold.** The check also produced a hard structural finding about irf and a constraint on
any future Deep DNAshape work here.

## irf is dropped, and now for a documented reason

irf was the longest reference (26 bp) and therefore the natural test case. It cannot be used, and the
reason is more serious than the earlier register diagnosis implied.

Counting states whose pyCurves output contains base-pair-level tables at all:

| | states | with bp tables |
|---|---|---|
| irf (frozen + relaxed) | 184 | **0 (0.0%)** |
| every other pilot, both arms | 2,224 | **2,224 (100.0%)** |

irf's frozen/relaxed JSONs are missing `curvesplus_base_pair_axis`, `local_base_base`,
`local_inter_base_pair`, and `groove` entirely — they retain only `local_inter_base` (single-strand),
`backbone`, and annotations. **pyCurves could not identify base pairs in any irf state**, so no duplex
was assigned and no base-pair-level geometry exists. Its crystal file is fine (13 bp rows, sane ProT of
+5.6°), so this is a failure on the generated states, not on the reference.

This independently confirms and explains the v3 finding. The 170° bend was not a mis-assigned helical
axis in an otherwise-valid duplex — **there was no duplex.** The bend numbers were computed on a
structure pyCurves did not recognise as double-stranded. That closes the question of whether irf's
values were recoverable: they were not, and dropping the pilot's DNA geometry was correct.

**Worth investigating separately** (it is a pipeline finding, not an analysis one): 184 states of
stage-3 output for one pilot produced no base-paired DNA. Whether that is a strand-pairing problem in
the relax setup, a chain-ID issue, or a genuine melting of the duplex under the restraints is unknown,
and the same failure mode could affect future pilots silently — the extract simply returns fewer rows.
**Recommend adding an assertion to the pyCurves batch step: every state must yield a
`local_base_base` table, and the count must equal the state count.** This run is the only reason the
failure was visible at all.

## The length test, on the 11 remaining pilots with bp geometry

Per-pilot ProT signal (relaxed arm, unpadded windows only):

| pilot | bp | unpadded positions | ρ |
|---|---|---|---|
| hsf | 12 | 4 | +0.00 |
| nfat | 12 | 4 | -1.00 |
| tbp | 12 | 4 | +0.00 |
| err | 13 | 5 | +0.00 |
| foxa | 13 | 5 | +0.30 |
| lef1 | 15 | 7 | +0.29 |
| ets1 | 15 | 7 | -0.39 |
| csl | 15 | 7 | -0.25 |
| runx | 16 | 8 | +0.57 |
| dux4 | 19 | 11 | +0.17 |
| engrailed | 21 | 13 | +0.69 |

**Per-pilot ρ vs duplex length: +0.571 (p = 0.066)**, and for ProT frozen +0.682 (p = 0.021). Read
alone that supports the hypothesis. It does not survive two checks:

1. **It does not firm up when noise is removed.** Restricting to the 6 pilots with ≥6 unpadded
   positions gives ρ = +0.698 but **p = 0.123** — and several per-pilot values in the full set are
   exactly ±1.00 on 4 positions, which is not an estimate. The trend is carried by the pilots least
   able to support one.
2. **Pooling positions reverses it.** The direct test — pool all unpadded positions from short vs long
   duplexes — gives ProT **+0.293 (n = 46, ≤15 bp, p = 0.048)** versus **+0.229 (n = 32, >15 bp,
   p = 0.206)**. The signal is if anything *weaker* in long duplexes. Roll: +0.122 short vs −0.059
   long. Twist: −0.016 vs +0.027.

The mixed model settles it. On all unpadded positions with pilot random intercepts:

| term | β | se | p |
|---|---|---|---|
| predicted ProT-FL (z) | +1.165 | 0.719 | 0.106 |
| duplex length (z) | −0.336 | 1.339 | 0.802 |
| **interaction** | **+0.905** | **0.654** | **0.167** |

LRT for adding both length terms: χ² = 1.92, df = 2, **p = 0.383**. **No moderation by length.**

So engrailed's ρ = +0.69 — the observation that motivated this check — is the best of eleven per-pilot
coefficients on 13 positions, and the length story it suggested does not generalise.

## A constraint that matters more than the result

| pilot | bp | unpadded positions | fraction |
|---|---|---|---|
| egr1 | 11 | 3 | 27% |
| hsf | 12 | 4 | 33% |
| tbp | 12 | 4 | 33% |
| nfat | 12 | 4 | 33% |
| err | 13 | 5 | 38% |
| foxa | 13 | 5 | 38% |
| csl | 15 | 7 | 47% |
| ets1 | 15 | 7 | 47% |
| lef1 | 15 | 7 | 47% |
| runx | 16 | 8 | 50% |
| dux4 | 19 | 11 | 58% |
| engrailed | 21 | 13 | 62% |

The fraction of usable positions follows **(L − 8)/L** exactly: layer 4 costs 8 bp of every duplex, 4 at
each end. Even a 30-mer gives **73%**, not 100%, and the pilots here range 27–62%. Two implications:

- **Any Part 3 design must budget for this.** A 20 bp designed duplex yields 12 usable positions; the
  padded 8 are not merely noisy, they carry a systematic upward bias in predicted `-FL` (ρ = +0.535
  against N-count) that will correlate with anything varying along the duplex.
- **Lowering the layer trades context for coverage.** Layer 2 (k = 5) costs only 4 bp and would give
  ~80% coverage on a 20-mer, at the price of shorter sequence context. Worth testing whether the ProT
  signal survives at layer 2 — that is a free re-run of the same tables and a better next step than
  either designing sequences or finding longer pilots.

## Where this leaves Part 3

The gate failed; the length explanation for that failure now also fails. Honest summary of the
sequence→shape line:

- Deep DNAshape's predicted fluctuation shows a **weak, correctly-signed, non-significant** relation to
  measured per-position ProT fluctuation (pooled ρ = +0.224, q = 0.295, mixed-model p = 0.149).
- It shows **nothing** for Roll or Twist, which are the parameters the project's bend axis is built
  from.
- Neither result is explained by duplex length, so the pilot set's short references are not the
  obstacle.

**Recommendation unchanged and now firmer: do not cost Part 3 as a mechanism test.** The remaining cheap
probes, in order, are (i) re-run the per-position gate at layer 2 for coverage, (ii) restrict to
interface positions per §4.2, which are interior by construction. If neither moves the Roll/Twist
result, the conclusion is that predicted free-DNA fluctuation and measured restrained-ensemble
fluctuation are different quantities in this system — and that belongs in the talk as a finding about
the method, not as a gap in the plan.

## Reproduce

```bash
cd /scratch1/shewchuk/deepdnashape_seqshape
python bpchk.py                  # base-pair table coverage per pilot/cond -> isolates irf at 0/184
```

Length test: join `measured_perposition_fluct.csv` to the Deep DNAshape predictions on
`(tf, feature, position)`, drop irf, keep positions with `nN == 0`, then regress `iqr ~ z_pred * z_bp`
with `groups = tf` (`statsmodels.formula.api.mixedlm`, ML) and compare to the main-effect model by LRT.
