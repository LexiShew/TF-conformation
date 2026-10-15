# v4 audit report — Part 1 items §1.1, §1.2, §1.3

**Run:** 2026-09-29, from data on disk. All three audits complete. **Each changes how the sweep is run,
and §1.3 changes whether it should be run in its current form at all.**

**Headline.** §1.1 comes back clean (first outcome row). §1.2 comes back in a state v4 did not
anticipate — it survives the mandatory control but was never significant to begin with, and the number
v4 inverted belongs to the *frozen* arm. §1.3 comes back in v4's **second** outcome row, and the
consequence is more severe than "spend on more seeds": **the sweep's discriminating prediction has ~21%
power as designed, and cannot be brought to 80% by adding seeds at any affordable count.**

---

## §1.1 — Own-family leakage audit. **Clean. First outcome row.**

Checked each pilot's own-family benchmark subset for its own reference PDB:

| pilot | own-family entries | distinct PDBs in subset |
|---|---|---|
| egr1 | 24 | 1a1g,1a1i,1a1j,2lt7,2wbu,4f6n,4r2d,4r2p,4x9j,5ke6,5ke7,5ke8,5kl3,5kl4,5und,5vmu,5vmv,5yef,6b0q,6dfc |
| nfat | 15 | 1pzu,1ram,1vkx,2i9t,3gut,7cli,7vuq |
| engrailed | 11 | 1le8,4xrm,5bng,5lty,5z2t,6a8r,6dfy,6u81 |
| ets1 | 10 | 3wts,3wtu,3wtw,4l0z,6vg2,6vg8,6vge,6vgg |
| err | 10 | 1r0n,4cn3,5emp,5emq,6fbq,6fbr,6x6e |
| tbp | 6 | 2etw,2euz,2evf,2evi,2evj,7oh9 |
| csl | 5 | 3iag,3v79,6dks,6wqu |
| runx | 5 | 1hjb,3wts,3wtu,3wtw,4l0z |
| foxa | 5 | 3coa,6lbi |
| irf | 4 | 7jm4 |
| hsf | 3 | 5hdn,7dcj |

**No pilot's reference PDB appears in its own-family subset** (csl/3brg, egr1/1aay, engrailed/3hdd,
err/1lo1, ets1/1k79, foxa/1vtn, hsf/5d5u, irf/2irf, nfat/1a66, runx/1hjc, tbp/1tgh — 0 of 11).
The middle point of Claim 3's gradient is clean as measured, and the gradient stands fully audited.
No recompute needed.

**One thing the audit surfaced that is worth recording.** ets1 and runx share four PDBs in their
own-family subsets — 3wts, 3wtu, 3wtw, 4l0z, the RUNX–ETS composite crystals. This is **correct
handling, not leakage**: each pilot scores a different chain with a different motif label
(`3wts_C_ETS1_HUMAN…` for ets1 vs `3wts_F_RUNX1_MOUSE…` for runx), so the two pilots read opposite
halves of the same co-crystal. Worth a methods sentence, because it looks like double-counting until
you check the chain IDs. It also means 4 of runx's 5 own-family entries and 4 of ets1's 8 sit on shared
structures, so those two pilots' own-family effects are not fully independent of each other.

**On v4's secondary point:** confirmed — the pilot set is one-per-family across twelve families, so no
pilot-pair design in this project can test family transfer. This is the same fact that made v3's §3.5
uncomputable, and it is a design property rather than a defect: the coverage breadth was bought at the
cost of family-transfer testability.

---

## §1.2 — `r_mean` as an endpoint. **Survives the control; was never significant.**

v4 read the +0.552 out of §2.2's dispersion table and inverted it. Two corrections:

**1. The number is frozen-arm, not relaxed.** The dispersion table's +0.552 was against
`dP_frozen`. Computing `r_mean` against both arms directly:

| relation | unadjusted | controlling state count | controlling bend IQR |
|---|---|---|---|
| `r_mean` vs ΔP **frozen** | −0.552 (p = 0.098) | **−0.584 (p = 0.098)** | −0.668 (p = 0.049) |
| `r_mean` vs ΔP **relaxed** | −0.321 (p = 0.366) | −0.395 (p = 0.293) | −0.516 (p = 0.155) |

**2. v4's stated confound is weaker than assumed.** v4: *"Surviving-state count correlates with crystal
agreement at ρ = 0.60, and state count **is** the fnat pass rate **is** the Ch. 15 rigidity ordering."*
On the cleaned panel, state count vs `r_mean` is **ρ = +0.321**, and state count vs ΔP is
**+0.006 (frozen, p = 0.99)** and +0.152 (relaxed, p = 0.68). State count does not predict the endpoint
at all, so it cannot be mediating this. The partial is therefore slightly *larger* than the raw
correlation, not smaller.

**So the result is v4's first outcome row on the control question and its third on the redundancy
question — but neither matters, because the effect never reaches p < 0.05 in the arm that carries
the project's claims.** `r_mean` vs bend IQR is ρ = −0.030, so the two axes are **not** redundant; they
are independent, which means adding `r_mean` genuinely expands the test family rather than collapsing it.

**Consequences:**
- **Do not promote `r_mean` to the sweep's endpoint.** A frozen-arm-only, p = 0.098 relation is not a
  validated substitute for ΔP. v4's §1.3-contingency stands: the sweep runs with ΔP.
- **The test family grows, it does not shrink.** `r_mean` × 2 arms are 2 more independent tests, taking
  the family from 16 to 18. Bend IQR's q moves from 0.149 to ≈0.17. Claim 4's "q ≈ 0.07–0.15 depending
  on how the family is counted" should become **q ≈ 0.15–0.17**.
- The frozen-arm relation is still interesting as a *descriptive* observation — pilots whose states
  agree poorly with their crystal gain more from augmentation, which is what a shared-label mechanism
  would predict — but it is one more axis in a family of eighteen, not an endpoint.

---

## §1.3 — Per-pilot ΔP standard error. **Second outcome row, with a harder consequence.**

Computed from `perseed_perentry.csv`, own-family entries, per pilot per arm, seed-level ΔP:

| pilot | relaxed ΔP | per-seed SD | seeds | SEM |
|---|---|---|---|---|
| egr1 | +0.0219 | 0.0403 | 5 | 0.0180 |
| csl | +0.1580 | 0.0581 | 5 | 0.0260 |
| runx | +0.0942 | 0.0597 | 5 | 0.0267 |
| ets1 | +0.1104 | 0.0631 | 5 | 0.0282 |
| nfat | -0.0038 | 0.0747 | 5 | 0.0334 |
| hsf | +0.0508 | 0.0799 | 5 | 0.0357 |
| foxa | +0.0301 | 0.0811 | 5 | 0.0362 |
| engrailed | +0.0135 | 0.0969 | 5 | 0.0433 |
| tbp | +0.0079 | 0.1178 | 5 | 0.0527 |
| err | +0.0359 | 0.1341 | 4 | 0.0670 |

| arm | between-pilot SD of ΔP | median per-pilot SEM | SEM / between-SD | attenuation factor |
|---|---|---|---|---|
| frozen | 0.0518 | 0.0334 | **0.644** | 0.841 |
| relaxed | 0.0523 | 0.0346 | **0.661** | 0.834 |

**Seed noise is about two-thirds of the between-pilot spread in both arms.** This is v4's second outcome
row: the observed ρ is attenuated. Disattenuating, the relaxed-arm −0.745 corresponds to a true
correlation near **−0.89** — which is also a warning, since a disattenuated estimate above 0.85 at
n = 10 is more likely to indicate the attenuation model is too simple than that the true effect is
that large.

**On v4's third outcome row (arm asymmetry): no.** Median SEM is 0.0334 frozen vs 0.0346 relaxed,
Wilcoxon p = 0.63, with only 4 of 10 pilots showing smaller relaxed SEM. **The arms are equally noisy.**
So the standing anomaly — relaxed ΔP correlating better despite being assumed noisier — is not resolved
by differential precision. The relaxed arm is not more precise; it simply correlates better. That keeps
the lucky-draw worry alive rather than shrinking it, and it cannot be quantified away.

### The part that changes the sweep

v4's reading of this row was *"spend on more seeds at two levels rather than more levels."* Directionally
right, but the arithmetic does not deliver a usable design.

The sweep's endpoint is a **paired difference** Δ(ΔP) between two k levels. Seed noise is independent per
pilot (seed main-effect LRT χ² = 1.84, df = 4, **p = 0.77** — there is no consistent seed offset to
cancel), so the paired difference has SD ≈ √2 × SEM:

| seeds per level | SEM | SD of paired Δ(ΔP) | ratio to between-pilot SD | attenuation |
|---|---|---|---|---|
| 5 (current) | 0.0346 | 0.0489 | 0.935 | **0.731** |
| 10 | 0.0244 | 0.0346 | 0.661 | 0.834 |
| 15 | 0.0200 | 0.0282 | 0.540 | 0.880 |
| 20 | 0.0173 | 0.0244 | 0.467 | 0.906 |

Power for v4's **prediction 2** — the interaction, |ρ| ≥ 0.5 at n ≈ 12, which §2.3 point 3 correctly
identifies as the only discriminating prediction:

| seeds per level | true ρ = 0.5 | true ρ = 0.7 |
|---|---|---|
| 5 (current) | **0.21** | 0.40 |
| 10 | 0.27 | 0.52 |
| 15 | 0.29 | 0.58 |
| 20 | 0.31 | 0.61 |
| 40 | 0.34 | 0.67 |

**The sweep as pre-registered has roughly a one-in-five chance of detecting its own discriminating
prediction, and quadrupling the seed budget takes that to one in three.** The ceiling is n = 12 pilots,
not seeds: at 5 seeds, reaching 80% power for a true ρ = 0.5 would need ≈ 60 pilots.

This does not say do not run the sweep. It says **the sweep cannot be the confirmatory test v4 casts it
as**, and the pre-registration should say so before it is timestamped. Three defensible responses:

1. **Run it as a calibration-and-direction experiment, not a hypothesis test.** §2.1's gate (does k move
   the axis?) is well-powered because it is a within-pilot geometry measurement on ~90 states, not a
   12-point correlation. Report the interaction as an effect size with a CI and explicitly decline to
   test it. Honest, cheap, and still informative.
2. **Change the unit of analysis to states, not pilots.** Δ(ΔP) per pilot is 12 points, but per-state
   bend change × per-state PWM change across k is ~700 points with pilot as a random effect. This is
   the step-2 design applied to an *intervention* rather than to observational variation — and unlike
   step 2, here the x is manipulated. It is the only version of the sweep that escapes n = 12.
   **Recommend this.**
3. **Spend the training budget on the held-out pilots instead** (v4 Part 5). A pre-registered sign test
   on 10–15 new pilots addresses selection bias, which is the actual problem; the sweep addresses
   confounding, which §1.2 just showed is not where the weakness is.

Option 2 and option 3 are complementary and together cost less than the four-level sweep.

---

## Consolidated effect on v4

| v4 item | audit result | change to the plan |
|---|---|---|
| §1.1 | clean, 0/11 self-leak | close; add the ets1/runx shared-PDB methods note |
| §1.2 | survives control, but frozen-arm only and p = 0.098 | **do not** promote to endpoint; family grows 16 → 18; Claim 4 q ≈ 0.15–0.17 |
| §1.3 | SEM/SD = 0.66 both arms; no arm asymmetry | **sweep is underpowered for its discriminating prediction**; re-scope per options above |
| §2.2 endpoint | unchanged — ΔP at {0, 5.0} | but unit of analysis should move to states |
| §2.1 gate | unaffected and still correct | it is the well-powered part of the sweep; run it |
| Part 5 | strengthened | held-out sign test now competes with the sweep on cost-effectiveness, not just on rigour |

**One addition to Part 2 before it is timestamped.** Add a stated power line: *"the interaction test has
~21% power at 5 seeds and n = 12 for a true ρ of 0.5; a null result therefore does not distinguish
'no mechanism' from 'underpowered', and the pre-registered reading of a null is 'uninformative on the
interaction' rather than 'mechanism wrong'."* Pre-registering a test while knowing it is underpowered is
defensible; pre-registering it without recording that is not.

---

## Reproduce

- §1.1: `analysis/data/perseed_perentry.csv`, filter `selffam == True`, group by `tf`, take
  `entry.str.split("_").str[0]` for PDB; compare against `PDB_ID` in `config/pilots/<tf>.sh`.
- §1.2: per-pilot `r_mean` from `per_state_all_conditions.csv` (`condition == "baseline"`,
  `kind == "state"`); state count = distinct states per pilot; partials by rank-residualizing both
  variables on the covariate.
- §1.3: seed-level ΔP = mean own-family `m_pearsonr` (augmented) − (baseline) per (tf, dna, seed);
  SEM = SD/√n_seeds; attenuation = 1/√(1 + (SEM/between-pilot SD)²); paired-difference SD = √2 × SEM,
  justified by the seed main-effect LRT (p = 0.77).
