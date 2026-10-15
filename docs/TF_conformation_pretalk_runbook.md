# Pre-talk runbook — three jobs

**2026-09-29.** Everything to run before the presentation, in order, with how to read each outcome and
what each one changes on a slide. Nothing here is exploratory — each job either lands in the deck or
protects something already in it.

**Total hands-on time: ~2.5 hours.** Job 1 runs unattended for days.

| # | job | hands-on | changes |
|---|---|---|---|
| 1 | Launch 33 `needs_ensembles` NMR pairs | 30 min | next-steps slide; bonus if partial results land |
| 2 | Regenerate `pilot_mechanism.csv` | 1 h | makes the Act-2 pilot table fully defensible |
| 3 | General-130 relaxed deltas, 11 pilots | 1 h | possible new slide, or one sentence |

**Not before the talk:** single-seed figure audit (correctness check on written material — do it after),
synthetic-duplex curvature check, size-matched donor control, terminal-window follow-ups, 2euv pilot,
λ sweep. None changes what you say on Tuesday.

---

## Job 1 — Launch the 33 staged NMR pairs

**Why first:** highest-value compute in the project independent of the talk, and launching it converts a
next-steps bullet from "we could expand this" into "running now, n goes from 14 to ~47." BioEmu only —
no docking, no DeepPBS, no training.

### Steps

1. **Inspect the queue.** `pairs.csv` has 5 ready and 33 at `status: needs_ensembles`. Confirm the
   status column and that each row names both an apo and a holo PDB.
2. **Screen before launching.** Three checks, every one learned from a failure this week:
   - **Both DNA chains present in the holo member.** The irf failure was Type B — deposition had both
     strands, local prep dropped one. RCSB metadata will not catch this. Check
     `structures/source_chains/<pdb>_chains/` after extraction, not before.
   - **Protein:duplex stoichiometry is 1:1** in the holo member. Two protein copies on one duplex is
     what broke irf.
   - **Construct sequences differ between apo and holo members** — if byte-identical, you have a
     noise-floor control pair rather than an informative pair (nhp6a and vnd are the existing examples,
     and they're how the 2–3% floor was measured). Keep them, label them.
3. **Launch** via `tools/run_new_nmr_examples.sh`, batched so a single failure doesn't take the array
   down. Both members of each pair, separately — the apo and holo constructs are run independently, and
   that separation is itself the control.
4. **Log the launch time and the pair list** so partial results can be scoped by what finished.

### Interpretation — only if a subset lands before the talk

Do not block the talk on this. If 8+ new pairs complete, recompute the three headline numbers:

| result | reading | slide |
|---|---|---|
| `%across` spread stays wide (roughly 0–100%, median near 55%) | Finding 1 holds at larger n. The unsystematic-landing result is robust. | update the number, same slide |
| `%across` narrows toward one end | The n = 14 spread was partly noise, and BioEmu *does* have a directional bias on protein–DNA. **A more interesting result than the current one** — report it as the headline. | rewrite Act 2 |
| Residual ratio stays above 1 in all pairs | The displacement signature is a measured constant, not a small-n artifact. Strengthens the correction to old slide 7. | strengthen the claim |
| Apo-beats-BioEmu count holds proportionally (~13/16) | Finding 2 holds. Quote the larger n. | update the number |
| Apo-beats-BioEmu weakens | Report both n's honestly. The pilot result (Finding 3) is independent and unaffected. | caveat on the slide |

**If nothing lands:** say it's running and give the expected n. That's a stronger next-steps slide than
any speculation.

---

## Job 2 — Regenerate `pilot_mechanism.csv`

**Why:** it is the one table in the NMR tree with no identifiable writer. `fnat_on_pilots.py` writes only
`fnat_pilots_summary.csv`, and the string `nmr_pass` appears in no script in the tree. The file was
written 5 s after the summary, so it was almost certainly an unsaved post-processing step.

**What depends on it:** the runx reclassification, and part of the five-pilot table that is going into
Act 2 of the talk. You want that table defensible end to end, not four-fifths of it.

### Steps

1. **Read `fnat_pilots_summary.csv`** and establish what columns it has versus what
   `pilot_mechanism.csv` claims.
2. **Reconstruct the `apo_pass` / `nmr_pass` distinction.** This is the substantive part. The
   reclassification rests on runx's eight apo entries disagreeing — six crystals passing at 100%, both
   NMR bundles (`1CMO` 43 models, `1CO1` 10) at 0%, fnat 0.25–0.38. So the two columns are "pass rate
   over crystal apo entries" and "pass rate over NMR apo entries," and the distinction only matters for
   pilots that have both.
3. **Write the regeneration as a script**, committed, so this gap doesn't recur.
4. **Verify** the regenerated table reproduces the five-pilot numbers in the findings doc: nfat
   0.509/60%, runx 0.729/100%, ets1 0.739/100%, hsf 0.793/100%, engrailed 0.807/100%.
5. **Fix the README's stale "Known limits" paragraph** while you're in the directory — it still says
   "Five pairs, three informative" against a data state of 16 pairs / 14 informative.

### Interpretation

| result | reading | slide |
|---|---|---|
| Regenerated table matches the findings doc | Fully sourced. **Use the whole five-pilot table including runx**, and the reclassification is citable. | Act 2 slide B as written |
| Numbers match but the crystal/NMR split can't be reconstructed | Use the table; **leave runx's reclassification out.** The pilot fnat numbers don't depend on it. | drop one sentence |
| Numbers don't match | Something upstream differs from what was reported. **Use only `fnat_pilots_summary.csv` values** and say the table is being re-derived. | conservative version |

**Note for the talk either way:** Finding 3's headline — BioEmu's best frame below the apo median in 5 of
5 — comes from columns that are in `fnat_pilots_summary.csv`, not from the disputed ones. **That line is
safe regardless of how Job 2 goes.**

---

## Job 3 — General-130 relaxed deltas across the other 11 pilots

**Why:** the TBP audit found the `_dnarelax` arm robustly harmful on the general-130 benchmark — 0/5 or
1/5 seeds improving, four metrics significant, harm 2–4× the frozen arm's. **That result is
label-independent**, since general-130 scores each entry against its own PWM, so the NDT80 mislabel never
enters.

Nobody has asked whether it generalizes. The evaluations exist — this is a read, not a run.

### Steps

1. **Pull per-seed general-130 metrics** for `augmented_<tf>_fold0_dnarelax_s1..s5` versus the paired
   `baseline_<tf>_fold0_s1..s5`, for all 11 non-TBP pilots.
2. **Compute paired deltas per seed**, sign-corrected so **+ = augmentation better**, on the same six
   metrics the TBP audit used: pearson r, spearman ρ, AUROC, IC-weighted PCC, Brier, MAE.
3. **Report seeds-improving counts**, not just means. The TBP result is persuasive because it is 0/5 on
   four metrics, not because the means are large.
4. **Compare against the frozen arm's** general-130 deltas for the same pilots — the TBP finding is that
   relaxation is 2–4× worse than frozen, which is a within-pilot comparison and more robust than either
   absolute.
5. **Exclude TBP** from the pooled statistic; report it separately as the pilot that prompted the check.

### Interpretation

| result | reading | slide |
|---|---|---|
| Degradation holds across most pilots, relaxed consistently worse than frozen | **A clean negative about DNA relaxation**, out-of-sample and paired. Its own slide. Note it does *not* contradict the all-conditions coin flip (Δr = +0.011, p = 0.36) — that was measured within **in-sample augmented** arms; this is out-of-sample against baseline. Different comparison, different answer. | new slide in Act 5 |
| TBP-only | One sentence, and a candidate explanation: the mislabel may interact with relaxation. | sentence in the TBP slide |
| Mixed by pilot | Report the split. Correlate the degradation against fnat pass rate and ensemble size — **not** bend IQR, which is retired. | sentence + backup slide |
| Relaxed *better* on some pilots | Report honestly. Would mean TBP is the outlier and relaxation is pilot-dependent, which is its own finding. | sentence |

**Why this matters beyond the slide.** If relaxation is generically harmful, the restraint calibration
(the descoped k sweep) changes from "does restraint strength do anything" to "characterizing a known
harm," and the Chapter 18 DNA-relaxation material needs reframing. Worth knowing before someone asks why
you relaxed the DNA at all.

---

## Order and timing

| when | job |
|---|---|
| **Tonight** | Job 1 launch — screening checks, then let it run |
| **Tomorrow AM** | Job 3 (1 h) — it's a read, fastest to resolve |
| **Tomorrow AM** | Job 2 (1 h) — then the README fix |
| **Tomorrow PM onward** | Slides only |

Jobs 2 and 3 are independent; do them in whichever order suits. Job 1 must go first because it runs for
days.

---

## What each job changes in the deck

**Job 1 → next-steps slide.** "Running now; n goes from 14 to ~47." If partial results land, update
Act 2's numbers.

**Job 2 → Act 2, slide B.** Determines whether the five-pilot table is fully sourced and whether runx's
reclassification is citable. The headline line is safe either way.

**Job 3 → Act 5, or one sentence.** Possible new slide on DNA relaxation being generically harmful.

**None of them changes the talk's structure**, which is the point of stopping here:

| act | content | status |
|---|---|---|
| 1 | Where we left off — Aug 5 claims, verdict table | fixed |
| 2 | **The generator doesn't produce binding-competent conformations** | Jobs 1–2 refine |
| 3 | What augmentation actually does — F1, F2, F3 | fixed |
| 4 | The retraction | fixed |
| 5 | What else — F4, audit findings, limits | Job 3 may add |
| 6 | Where this goes — three papers, 33 pairs running | Job 1 feeds |

---

## The line to protect

Whatever happens with these three jobs, this survives and should be the largest text on any slide:

> **BioEmu's single best frame is below the deposited apo median in 5 of 5 pilots.**

Across ~90 frames per pilot. It comes from `fnat_pilots_summary.csv`, it needs no DeepPBS, and it is
unaffected by the shared-label and in-sample issues that qualify everything downstream.

It is also the sentence that makes the whole project's arc coherent: the generator doesn't produce the
right conformations, and the training objective couldn't use them if it did. Both measured, both this
week.
