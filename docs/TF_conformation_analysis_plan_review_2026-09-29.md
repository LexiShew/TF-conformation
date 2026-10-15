# Review of `TF_conformation_analysis_plan_2026-09-29.md`

**Reviewed against:** `/project2/rohs_102/shewchuk/TF-conformation` on endeavour, read-only, 2026-09-29.
**Method:** every numeric claim in the plan was traced to the file that produces it. Where the file
existed, the number was recomputed. Where it did not, that is stated as a blocker rather than assumed.

The plan is unusually good. Pre-registering outcome readings, naming the RI normalization trap before
generating anything, and demanding the partial in §9 are all correct and rarely done. Most of what
follows is not disagreement — it is that **four of the eleven analyses cannot run as written**, and two
of them are already done and came out differently than the plan assumes.

---

## Part 1 — Blockers (these will not run, or will produce a number that means something else)

### B1. §8 cannot be run at all: the augmented checkpoints were trained on all 701 states

This is the most important finding in the review, and it kills §8 in its stated form.

`check_membership.py` already computed this and its output is on disk. Reading
`train_membership.json`:

| pilot | frozen states scored | of those, in that pilot's augmented training set |
|---|---|---|
| ets1 | 99 | 99 |
| tbp | 94 | 94 |
| hsf | 96 | 96 |
| egr1 | 90 | 90 |
| engrailed | 85 | 85 |
| foxa | 73 | 73 |
| csl | 57 | 57 |
| err | 41 | 41 |
| runx | 35 | 35 |
| nfat | 19 | 19 |
| lef1 | 13 | 13 |
| irf | 5 | 5 |
| **total** | **707** | **707** |

**Every single state is in-sample for the `augmented_<tf>_fold0_sN` checkpoint.** `fracF = 1.000` for
all twelve pilots; `frozen_states_heldout_from_frozen` is empty for every pilot. The same holds for the
relaxed arm (671/671).

§8's premise is that comparing baseline-checkpoint and augmented-checkpoint scores on the same 701
states separates "conformers carry degraded specificity" from "the model is out of distribution on
them." It cannot, because the augmented checkpoint **memorized these exact structures**. The comparison
measures training-set recall, not distribution shift.

And it is not a subtle effect. That analysis has effectively already been run — `compare_all_conditions.py`
scored all three conditions on 2026-09-29 and `per_tf_all_conditions.csv` is on disk:

- mean state-vs-own-crystal agreement: **baseline r = 0.504 → augmented/frozen r = 0.746 (Δ +0.242)**
- 11 of 12 pilots improve; csl +0.459, ets1 +0.329, hsf +0.343

If §8 is run as written, it will return exactly this, and its own outcome table says to read a material
rise as *"a substantial part of the deficit is distribution shift → the fix is training-side."* That
reading would be wrong. The number is the in-sample fit of a model that saw every one of those frames.

**What to do instead — one of:**

1. **Cross-pilot checkpoints (recommended, and cheap).** Score pilot X's states under
   `augmented_Y_fold0_sN` for Y ≠ X. Those checkpoints have seen *minimized BioEmu frames* (so the OOD
   question is answered) but never *pilot X's* frames (so it is out of sample). This is the clean version
   of the experiment §8 wants, and it needs only inference on existing checkpoints — 11 × 12 pilot pairs
   against existing `.npz`, no training.
2. **Leave-states-out retraining.** Correct but expensive; a fresh augmented arm per pilot holding out a
   state subset. Not a week-scale task.
3. **Reframe §8 as what it actually measures** and keep it: "the augmented model reproduces its training
   frames' PWMs at r = 0.75 vs 0.50 for a model that never saw them." That is a real result about
   label-sharing and memorization capacity — it belongs in the shared-label / regularization discussion
   (Ch. 12/17/21 in `book/OUTLINE.md`), not in the conformer-degradation argument.

Note the plan calls §8 a half-day. The inference already exists; what is left is *not running it* and
rewriting the claim.

### B2. §5's positive control is factually wrong — HSF's crystal is not the worst

§5 says: *"HSF1's crystal (5D5U) was featurized fresh and is not in training — and it is also the pilot
where the crystal does worst. Consistent with leakage carrying part of the crystal's advantage."*

The first half is right (`crystal_leakage.json`: `exact_in_baseline = false`, no same-structure entries).
The second half is not. From `per_tf_all_conditions.csv`, baseline condition, `crystal_r_exp` worst-first:

| rank | pilot | crystal r vs experiment | leakage class |
|---|---|---|---|
| 1 | tbp | 0.128 | **held out** |
| 2 | lef1 | 0.324 | exact |
| 3 | err | 0.405 | **held out** |
| 4 | nfat | 0.624 | exact |
| 5 | **hsf** | **0.708** | **held out** |
| … | | | |
| 12 | csl | 0.924 | exact |

HSF is 5th worst, not worst. The argument also doesn't survive being fixed: of the four pilots whose
crystal does worse than HSF, two (lef1, nfat) have **exact** leakage. A leaked crystal sitting at the
bottom of the ranking is the opposite of the claimed pattern.

**Worse, the headline §5 is trying to protect does not need protecting.** I ran the de-leaked comparison.
Crystal advantage over the ensemble mean, in agreement with the experimental PWM:

| subset | n | crystal r | ensemble-mean r | crystal − ensemble | crystal wins |
|---|---|---|---|---|---|
| pooled | 12 | 0.660 | 0.559 | **+0.102** | 11/12 |
| exact-leakage only | 6 | — | — | **+0.106** | 6/6 |
| held out (exact) | 5 | 0.548 | 0.453 | **+0.095** | 4/5 |
| held out + structure-only | 6 | 0.598 | 0.501 | **+0.097** | 5/6 |

**The crystal's edge is +0.106 on leaked pilots and +0.095 on held-out pilots.** Leakage is not carrying
it. This is §5's best outcome ("leakage wasn't driving it — quote freely") and it is already in hand;
the analysis is ~0 hours, not 1, because the numbers are computed. What needs an hour is deciding how
to say it.

Two caveats to state when you do:
- The absolute level drops a lot on held-out pilots (0.660 → 0.548). The *advantage* is stable; the
  *ceiling* is leakage-inflated. Those are different claims and the talk should make only the first.
- n = 5. The advantage's stability across subsets is reassuring but it is five points.

Also: §5 says "eleven of twelve" and excludes seven pilots to leave five. The held-out set is 5 (exact)
or 6 (including foxa, which is structure-only: same PDB+chain under a different PWM label). Pick one and
define it — and note the plan's exclusion list (CSL, EGR1, Engrailed, FoxA, IRF1, LEF1, NFAT) matches the
leakage classes on disk, so the roster is right even though the HSF claim isn't.

### B3. §3 has no input data: `roll`, `tilt`, and `twist` are not in any distilled table

§3 proposes decomposing bend-IQR into per-bp roll/tilt/twist. The tables that look like they'd hold this
do not.

- `analysis/analyses/dna_relax/data/<tf>_dna_perbp.csv` → columns are `bp_Bres, cond, mean_c1c1, sd_c1c1, n`.
  That is C1′–C1′ distance only. No roll, no tilt, no twist.
- `<tf>_dna_perstate.csv` → `cond, state, dna_bb_rmsd, bend_docked, bend, dbend, maxppgap`. Scalar bend only.
- `pycurves_all_perstructure.csv` → `bend_uu, bend_pp, shortening, minor_w, major_w, inclin, tip, xdisp, ydisp, minor_w_cp, major_w_cp`. Per-structure, not per-bp, and still no roll/tilt/twist.

The data **does** exist, but only in the raw per-frame JSON: `pycurves/<tf>/{frozen,relaxed}_state_NNN_curvesplus.json`,
under `dataframes.local_inter_base_pair` (per-step: roll, tilt, twist) and `dataframes.local_inter_base`
(shift, slide, rise). That is ~326 MB across 13 pilots and needs a new aggregation pass — the existing
`pycurves_batch/aggregate_perposition.py` is the natural place to extend.

So §3 is **not 0.5 d**. It is a new per-step extraction over 2,400+ JSON files plus the correlation
re-run. Call it 1–1.5 d. And note it shares its entire extraction step with §4 — see S1 below, because
that changes the schedule.

### B4. §7's input table does not exist in the repo

§7 is described as *"Both tables already exist; this is a fifteen-minute join."* The Ch. 16 own-family
column does exist (`analysis/analyses/mechanism/data/selffam_effects.csv` — I verified all five values:
ets1 +0.111, hsf +0.037, nfat +0.015, engrailed −0.015, runx +0.084 all match the plan exactly).

The **Ch. 19½ apo-fnat class column does not.** Searching the repo for `apo_fnat`, `apo-fnat`,
`selection-like`, `S-switch`, and the signature values (`100%`, `60%`, `0%, reclassified`) returns nothing
in `docs/`, `book/`, or `analysis/`. There is no `mechanism_apo_holo`-style table carrying an apo-fnat
classification, no Ch. 19½ document, and no NMR-pair results table. What exists is
`structures/nmr/` (71 dirs of BioEmu ensembles for NMR PDB IDs) and `tools/run_new_nmr_examples.sh`
(the generation script) — i.e. **inputs, not results**.

Either that analysis lives outside this repo (a laptop, a notebook, another host) or it hasn't been
written down. §7 cannot be a fifteen-minute join against a table that isn't there. Before scheduling it,
locate the source. If it only exists as prose, §7's first cost is transcribing it into a CSV with
provenance — and then it is a 15-minute join.

The plan's own §7 table is also internally at odds with itself: it states the hypothesis predicts "help"
for engrailed (selection-like) and scores it ✗, and predicts "hurt" for runx (induced-fit) and scores it
✗. That's right, but then §7 is 3/5 *by construction of a hypothesis that §17.2 already rejected*. Its
value is exactly what the plan says — a model-free axis reaching the same negative conclusion — so it
should be framed as a **confirmation of a null**, and at n = 5 it cannot be more.

---

## Part 2 — Corrections to specific claims

### C1. §10's premise is quantitatively wrong: the frame-count bias is ~5× too small to matter

§10 argues TBP's 0.59 Å `d_min` is a min-over-N artifact inflated by TBP's frame count, citing the NMR
precedent where the correction moved p from 0.002 to 0.15.

Two problems.

**First, TBP does not have an unusually large ensemble.** From `analysis/data/reachability.csv`, TBP has
**n = 91**, which is 6th of 13. Frame counts span 77–100 — a 1.3× range, not the kind of spread that
drives min-over-N. (The plan says "TBP has 99 frames"; 99 is TBP's *pre-gate* frame count from
`RESULTS_INVENTORY.md`. The reachability calculation runs on n = 91.)

**Second, I simulated the bias.** Moment-matching a gamma to each pilot's (`d_med`, `spread`) and taking
min-over-n versus min-over-77:

- mean expected upward shift when matching to k = 77: **+0.018 Å**
- maximum over all pilots: +0.139 Å (dux4, n = 90 — and dux4 is excluded from the benchmark anyway)
- TBP's own expected shift: **+0.003 Å**
- the gap §10 needs to explain (TBP 0.587 → engrailed 0.690): **0.103 Å**

TBP's ranking is robust to the correction by a factor of ~30. There is a real `n`-vs-`d_min` correlation
(Spearman ρ = −0.66, p = 0.015), but it is confounded: floppy pilots lose more states at the fnat gate
*and* reach less well, so n and reachability share a common cause (rigidity). Matching frame counts
does not break that confound — it just removes the ~0.02 Å mechanical part.

**Recommendation:** demote §10 from "1 h, decides whether TBP is still the decisive control" to a
20-minute robustness footnote. Run best-k anyway (it is cheap and the check should be on record), but
pre-register that it will not move the ranking, and do not budget a Ch. 4/17/23 rewrite against it.
The `d_min` column also isn't recomputable from the distilled tables — `compute_reachability.py` writes
only the min/median/spread summary, so best-k means re-running it over stage-3 frames. Still cheap, but
it is a pipeline job, not a CSV join.

### C2. §4's headline ρ = −0.76 is an artifact of a label mismatch — and this is the real Day-1 finding

§4 is correctly identified as the highest-value item. The reason is different and more serious than the
plan thinks.

The axis is `froz_bend_iqr`. From `mechanism_analysis.py:102`, it is computed as
`groupby(tf).agg(froz_bend_iqr=("bend_uu", lambda s: s.quantile(.75)-s.quantile(.25)))` filtered on
**`cond == "frozen"`**. There is no relaxed-arm equivalent anywhere in the repo — I grepped every
`.py` for `bend_iqr`, and `selffam_vs_shape.csv`'s only IQR column is `froz_bend_iqr`.

So in the headline row of §17.3:

| axis | arm | n | ρ |
|---|---|---|---|
| Ensemble DNA-bend IQR | **relaxed** | 11 | **−0.76** |
| Ensemble DNA-bend IQR | frozen | 11 | −0.53 |

**both rows use the same frozen-arm x-values.** Only the y changes: `dna == "frozen"` vs `dna == "relax"`
own-family ΔPearson. I reproduced both from `selffam_vs_shape.csv` + `selffam_effects.csv`:

- froz_bend_iqr vs frozen ΔP: ρ = **−0.527**, p = 0.096 (matches −0.53)
- froz_bend_iqr vs relaxed ΔP: ρ = **−0.764**, p = 0.006 (matches −0.76)

This is confirmed by `dnaflex_correlations.csv`, whose own `axis` column reads `froz_bend_iqr` for
**both** the `frozen` and `relax` rows — the file is self-documenting about it.

**Why this matters.** "The mechanism axis *strengthens* in the relaxed arm" (§17.3, §18.2, and Ch. 18.4's
tie-together) is the project's single strongest claim — the only bootstrap CI excluding zero. But the
axis did not change between arms. What the two rows show is that **frozen-arm DNA-bend fluctuation
predicts the relaxed arm's augmentation effect better than it predicts the frozen arm's.** That is a
much stranger statement than "relaxation sharpens the axis," and it has a plain alternative reading: the
relaxed-arm ΔP values are the noisier of the two (4–6 seeds, and `irf` and `nfat` flip sign between
arms), so this may be a lucky draw at n = 11, not a sharpening.

Note the naming doesn't help: `FINDINGS.md`'s §17.6 flags a *previous* bug of exactly this shape (filtering
`"relaxed"` where the data says `"relax"`, silently producing an empty arm). This is the same class of
error surviving in a different column.

**Do this before anything else (~2 h):**
1. Compute `relax_bend_iqr` from `pycurves_all_perstructure.csv` — `cond == "relaxed"` rows are present
   for **all 13 pilots** — and re-run the correlation with the arm-matched axis.
2. Report the 2×2 honestly: {frozen, relaxed} axis × {frozen, relaxed} effect. Four numbers, not two.
3. Only then ask §4's interface-restriction question.

Also worth knowing before you interpret §4: the correlation is fragile in the way that matters for the
plan's own outcome table. Dropping TBP alone takes the frozen-arm ρ from −0.527 to −0.370 (and p from
0.096 to 0.293). The relaxed-arm ρ is sturdier (−0.764 → −0.721 without TBP; leave-one-out range
−0.709 to −0.879). TBP is both the highest-IQR pilot (11.78 vs a 1.29–9.95 range for everyone else) and
the induced-fit exemplar, so it is simultaneously the most influential point and the one the narrative
leans on hardest. State the LOO range in the talk.

### C3. §6 is a good question but the gate has already truncated its x-axis

§6 wants per-state fnat × r-to-crystal to decide whether `FNAT_FLOOR = 0.5` is a quality filter. The
gate makes this partly unanswerable: PWM predictions only exist for states that **passed** the gate,
so every point available has fnat ∈ [0.5, 1.0]. Counting from the `*_fnat.csv` files against the states
that have predictions:

| pilot | stage-3 states | fnat < 0.5 (gated out) | states with predictions |
|---|---|---|---|
| lef1 | 83 | 67 | 13 |
| runx | 97 | 61 | 35 |
| nfat | 86 | 60 | 19 |
| err | 90 | 44 | 41 |
| csl | 84 | 27 | 57 |
| irf | 91 | 13 | 5 |

For the floppy pilots, 70–80% of the fnat range is gone. A ρ ≈ 0 on the surviving half **cannot**
distinguish "the gate filters geometry and nothing else" (§6's stated reading) from "the gate worked, and
above the floor fnat stops mattering" — which is what a well-chosen floor should look like. §6's middle
outcome row, which tells Ch. 8 to stop calling it a quality filter and narrows the Ch. 15 rigidity axis,
is therefore not licensed by the truncated data.

To actually decide it, featurize and score a sample of **failing** states (say 20 per pilot at
fnat < 0.5) and extend the correlation across the floor. Stage 4 + inference on ~240 extra states; that's
a real job, not an hour, but it's the only version that answers the question. Alternatively, keep §6
as a within-passing-range descriptive and delete the Ch. 8 / Ch. 15 consequences from its outcome table.

The plan's HSF observation is right, though, and worth keeping: HSF passes at ~99% and has the worst
state-vs-crystal agreement in the set (r = 0.272 — I confirmed it is rank 1 of 12 worst on
`r_mean`, baseline). That decoupling is real and is a single-pilot anecdote strong enough for a slide.

### C4. §1 needs a fourth condition that doesn't exist yet

§1 asks to pair each relaxed state against "its frozen counterpart" and score both against the
experimental PWM. The scoring harness (`compare_all_conditions.py`) has three conditions:
`baseline`, `aug_frozen`, `aug_relax`. Note what `aug_relax` is: **augmented-arm checkpoints** scoring
relaxed states (`KIND["aug_relax"] = "augmented"`, features from `stage5_aug_dnarelax`).

There is **no `baseline × relaxed` condition** — no baseline checkpoint has scored the relaxed states.
So the available relaxed-vs-frozen comparison is contaminated by B1: both arms are in-sample. Adding the
fourth condition is a small change to `COND` (point a `baseline` kind at the relax id-file and relaxed
feature dir) plus an inference run, and it is what makes §1 interpretable. Do that first.

Two smaller things on §1:
- **Coverage is better than the plan says.** §1 claims "7 pilots have dnarelax data. n = 7." On disk,
  `output/stage3_min_dnarelax/`, `stage5_aug_dnarelax/`, and the relaxed rows of
  `pycurves_all_perstructure.csv` all cover **12–13 pilots** (12 in the aug assembly, 13 in pyCurves).
  `per_tf_all_conditions.csv` already has `aug_relax` for all 12. The "7 pilots" figure comes from
  `RESULTS_INVENTORY.md`, which is stale (it also says "only tbp_dnarelax has a completed eval"). **n = 12,
  not 7** — this is the single biggest upward revision available to the plan, and it affects §1, §11,
  and the Ch. 18 "partial coverage" caption.
- The plan's design decision — pre-register agreement with the **experimental** PWM as the endpoint and
  crystal agreement as descriptive — is correct and important. Keep it exactly as written.

### C5. §2's normalization analysis is right; one addition

The min-max-within-structure trap is correctly identified and the two options are the right two. Add a
third that is cheaper than option 2 and more interpretable than option 1: report the **rank correlation
plus the un-normalized MAE distribution per structure** side by side. If the raw MAE scale shifts between
conformations, the rank correlation alone will hide it, and that scale shift is itself a result (it says
the model's sensitivity to the interface changed).

Also: on-disk RI exists for **ets1 only** (`output/interpret_results/` holds
`{baseline,augmented}_ets1_fold0_importance.npz` plus `augmented_ets1_fold0_dnarelax_s1_importance.npz`).
So the "one pilot, crystal + 3–4 relaxed frames" scoping is not just advisable, it matches what is
already computed — and the dnarelax importance file means part of §2 may already be done. Check that file
before generating anything.

### C6. §11's conflict list is right, and there is a bigger one it misses

I verified the specific conflicts. They are real, and the plan's diagnosis (one canonical table) is
correct. Two notes:

- **The "Outline Ch. 16 / Ch. 17" references point at a numbering scheme that no longer exists.** There
  are two outlines on disk with **incompatible chapter numbers**: `book/OUTLINE.md` (21 chapters —
  family structure is Ch. 13, mechanism is Ch. 14, DNA relaxation is Ch. 15) and
  `docs/TF_conformation_book_outline.md` (Part-based, v2 — where the mechanism material is the Ch. 16/17
  the plan cites). The `docs/TF_conformation_Ch15-18_results.md` results doc uses the *second* scheme.
  So "Outline Ch. 16 says csl +0.027" is `docs/TF_conformation_book_outline.md:226`, while `book/OUTLINE.md`
  has that content at Ch. 13. **Reconciling numbers without first reconciling the chapter map will
  generate new conflicts.** Fix the numbering first; it is the cheapest item on the whole list.
- The specific values check out: `csl +0.027 / nfat +0.015` at `docs/TF_conformation_book_outline.md:226`
  vs `selffam_effects.csv`'s csl −0.043 / nfat −0.033 for the *other-family* rows — these are
  **different quantities** (cross-benchmark pooled vs own-family vs other-family), not a contradiction.
  Same for "TBP −0.26": that value appears exactly once, at
  `docs/TF_conformation_book_outline.md:238`, and matches nothing on disk — it looks like a typo for
  −0.026 or a stale figure. The ETS +0.042 / IRF +0.039 pair is a **median over entries**
  (`book/OUTLINE.md:384`) vs the +0.008 / +0.010 **seed-level family means** in
  `family_effects_seedlevel.csv` — again different estimators, so the "4–5×" discrepancy is a units
  problem, not an error. Your canonical table needs an **estimator column**, not just one number per
  pilot.

On the **fnat scale gap** (100% vs 1%), the plan's explanation 2 is correct and I can strengthen it:
the pilot gate is structurally circular. `stage2_redock.py` docks each frame onto the reference
crystal's DNA by the protein's Kabsch transform, and fnat is then scored against that same reference's
contact set. The states are placed to match the reference and then graded on matching it. The plan's
proposed test (score against a *different* co-crystal of the same TF) is exactly right and is the single
highest-value item in §11.

---

## Part 3 — Schedule and scope issues

### S1. §3 and §4 share their expensive step; sequence them together

Both need per-bp/per-step geometry extracted from the raw pyCurves JSON (see B3). Doing §4 on Day 1 and
§3 on Day 4 means writing and running that extraction twice. Merge them: one extraction pass that emits
per-step roll/tilt/twist *and* per-bp bend variance, then two cheap correlation analyses on top. Saves
roughly half a day and guarantees the two analyses use the same register convention — which matters,
because §3's own outcome table flags a 3DNA end-convention artifact as a risk.

### S2. The week is over-committed by roughly 2.5 days

Re-costing against what's on disk:

| § | plan | actual | why |
|---|---|---|---|
| 1 | 1–2 d | 1 d + inference | needs the 4th condition (C4); scoring harness exists |
| 2 | 0.5 d | 0.5 d | ok — and may be partly done (C5) |
| 3 | 0.5 d | **1–1.5 d** | no input data; new JSON extraction (B3) |
| 4 | 0.5 d | **1 d** | +2 h for the arm-matched axis (C2); shares extraction with §3 |
| 5 | 1 h | **~0 h** | already computed (B2) |
| 6 | 1 h | 1 h, or **+1 d** | 1 h for the truncated version; a day if you score failing states (C3) |
| 7 | 15 min | **blocked** | input table absent (B4) |
| 8 | 0.5 d | **blocked / reframe** | in-sample (B1) |
| 9 | 2 h | 2 h | ok |
| 10 | 1 h | **20 min** | bias is negligible (C1) |
| 11 | 0.5 d | **1 d** | chapter-map reconciliation first (C6) |

**Suggested reshuffle:**

| day | work |
|---|---|
| 1 | **C2 arm-matched bend-IQR** (2 h — this gates §17.3 and may retract the headline). Then start the shared §3+§4 per-step extraction. Launch the §1 fourth-condition inference. |
| 2 | §3 + §4 on the extracted per-step data, including interface restriction. |
| 3 | §1 paired analysis. §9 partial (2 h). §10 footnote (20 min). |
| 4 | §5 write-up (numbers in hand), §6 truncated version + HSF anecdote, §2 single-pilot RI. |
| 5 | §11 — chapter map first, then the canonical table with an estimator column. Locate §7's source table; join if found. |
| 6–7 | figures, talk, slack. |

Two analyses move off the week: §8 (reframed, no new compute) and §7 (until its input is located).

### S3. §12's three deferrals are correctly reasoned — keep all three

Particularly the third. Weighting training states by predicted-PWM agreement is exactly the circularity
the plan describes, and the reasoning (fnat/iRMSD/bend-IQR are model-independent, PWM agreement is not)
is the right test to apply. No notes.

---

## Part 4 — On the three claims the talk should make

Claim 1 (crystal as point estimate; conformation SD 0.25 vs seed SD 0.16, n = 701) and Claim 2
(blur toward nonspecific, not redirection) are as robust as the plan says, and both are independent of
everything above. Agreed: these hold regardless of the week.

**Claim 3 is more exposed than the plan accounts for.** The plan says only Claim 3 can fail, and treats
§4's interface-restriction as the risk. The bigger risk is C2: if the "relaxed arm sharpens the axis"
result is a frozen-axis-vs-relaxed-effect mismatch, then the "residual signal is a DNA property" claim
loses the one bootstrap CI that excludes zero — *before* §4's end-fraying question is even asked. The
plan's own §4 outcome table ("ρ collapses when ends are excluded → Ch. 17 reverts to 'no axis identified'")
is the right posture; it just needs to be triggered by C2 as well.

One addition to the talk that the plan doesn't include, and should: **state the n = 707-states /
in-sample structure of the augmented arms explicitly** (B1). It is the kind of thing an audience member
finds in the methods and asks about, and it is much better volunteered. It also makes Claim 1 stronger,
not weaker — the 701-state conformation-vs-seed SD comparison is a baseline-checkpoint result and is
unaffected.

The standing note about pilot count (n = 7–13) being the binding sample size is correct and should be
said out loud, as the plan says. Add the correction from C4: **for the relaxed arm it is n = 12, not 7.**

---

## Summary table

| item | status | action |
|---|---|---|
| §8 augmented rerun | **blocked** — all 707 states in-sample | cross-pilot checkpoints, or reframe as memorization result |
| §7 apo-fnat join | **blocked** — input table not in repo | locate source before scheduling |
| §3 roll/tilt/twist | **blocked** — no distilled per-step data | extract from raw pyCurves JSON; merge with §4 |
| §4 bend-IQR | **runs, but axis is mislabeled** | compute relaxed-arm IQR first; report the 2×2 |
| §5 de-leak | **already computed** | crystal edge +0.106 leaked vs +0.095 held out; fix the HSF claim |
| §6 fnat × r | **runs, cannot decide the question** | keep descriptive, or score failing states |
| §10 frame-matching | **runs, premise wrong** | 20-min footnote; bias is ~0.02 Å vs a 0.103 Å gap |
| §11 reconciliation | **runs, underscoped** | reconcile the two chapter-numbering schemes first |
| §1, §2, §9, §12 | **sound** | §1 needs a baseline×relaxed condition; §1 coverage is n = 12 |

---

## Reproducing the checks in this review

```bash
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /project2/rohs_102/shewchuk/TF-conformation

# B1 leakage / in-sample membership (already on disk, re-run to confirm)
python analysis/analyses/conf_vs_crystal_pwm/check_membership.py
python analysis/analyses/conf_vs_crystal_pwm/check_leakage.py

# B2 / C4 / §5 numbers  -> per_tf_all_conditions.csv
python analysis/analyses/conf_vs_crystal_pwm/compare_all_conditions.py

# C2 axis provenance
grep -n "bend_iqr" analysis/analyses/mechanism/*.py
head -1 analysis/analyses/mechanism/data/selffam_vs_shape.csv   # only froz_bend_iqr
cat analysis/analyses/mechanism/data/dnaflex_correlations.csv   # 'axis' is froz_bend_iqr in BOTH arms

# C1 reachability frame counts
column -s, -t analysis/data/reachability.csv

# B3 confirm roll/tilt/twist live only in raw JSON
head -1 analysis/analyses/dna_relax/data/tbp_dna_perbp.csv
python -c "import json;d=json.load(open('analysis/analyses/dna_relax/pycurves/tbp/crystal_curvesplus.json'));print(list(d['dataframes']));print(list(d['dataframes']['local_inter_base_pair'][0]))"
```
