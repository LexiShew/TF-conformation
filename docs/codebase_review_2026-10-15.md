# TF-conformation — codebase review

**Date:** 2026-10-15 · **Scope:** full repo (1,920 tracked files; ~25k lines of project Python/shell, excluding vendored `lib/`)
**Reviewer note:** every claim below is labelled **CONFIRMED** (verified by running a check against files in this repo) or **INFERRED** (reasoned from code, not executed). Commands used are given so each check is re-runnable.

> ## Remediation status (updated 2026-10-15, after the review)
>
> This document records the repo **as reviewed**. Work done since is tracked here so it does
> not read as an open to-do list. See `docs/pilot_operations.md`, `docs/retired_pilots.md`,
> and `docs/tbp_label_error_downstream.md`.
>
> | finding | status |
> |---|---|
> | **F1** TBP label live in config | **FIXED** — `PWM_LABEL=TBP_HUMAN.H11MO.0.A`, filter `MA0386\|TBP_`, `ON_TARGET_USABLE=no` |
> | **F2** mislabel in analysis layer | **FIXED** — 30 rows in `family_annotation.csv` + 1,285 in `perseed_perentry.csv` relabelled to `NDT80`; downstream artifacts inventoried as needing regeneration |
> | **F3** HAP1 misfiled | **FIXED** — relabelled to `Zn2Cys6 cluster` (6 + 257 rows) |
> | **F4** no invariant enforcing labels | **FIXED** — `scripts/pipeline/validate_pilot.py` (9 checks, C1 is the TBP check); `config/pilots.tsv` + `gen_pilot_configs.py` make the error class structurally impossible; both wired as `run_pilot.sh` pre-flight gates |
> | **F5** inconsistent on-target subsets | **PARTLY FIXED** — foxa's 0-match filter corrected to n=5; engrailed/lef1/tbp now declare `ON_TARGET_USABLE=no`; err/nfat renamed "(family transfer)". The validator warns on pseudoreplication. `irf` retired. |
> | **F6** runs never seeded | **FIXED in code** — `driver.py` now seeds python/numpy/torch/cuda from `C["random_seed"]`; `--random_seed` CLI added; unseeded runs warn loudly; check **C10** guards both halves against regression (negative-tested). **The 264 existing checkpoints are still unseeded and unpaired — a genuinely paired result requires retraining.** |
> | **F7** silent stage-4 attrition | **PARTLY FIXED** — `pilot_queue.py --status` now surfaces it; `irf` (94% loss) retired. Stage 4 itself still does not fail on attrition. |
> | **F8** evaluator prints 1 of 6 metrics | **FIXED** — print moved inside the loop; all six metrics now print, plus q-values, Cohen's dz, and a pseudoreplication banner. Verified by executing the real print block. |
> | **F9** normal-approx p-values | **FIXED** — `compute_paired_stats` now delegates to `analysis/common/stats.py` (t-distribution, df=n−1). Also fixed in the same pass: the bootstrap RNG was re-seeded to 42 *inside* the per-metric loop, so all six CIs shared one resample draw; each metric now gets a distinct seed (verified 6/6 distinct). Output schema is a superset — nothing downstream breaks. |
> | **F10** silent frame fallback | **OPEN** |
> | **F12** BH-FDR missing step-up | **FIXED in 3 places** — `mechanism_analysis.py`, `make_samefamily.py`, and `mechanism_confound.py` (the last of these wrote the published `q_bh` column and I had missed it on the first pass). AST-verified that no hand-rolled BH remains in executable code. |
> | **F13** headline fails FDR | **FIXED (documented)** — and the picture is better than this review first stated. `FINDINGS.md` *does* report FDR for the own-vs-other contrast test; only the Headline's own-family p was raw. With BH corrected and the two disqualified pilots dropped (`irf` retired, `tbp` contaminated), csl and ets1 reach **q = 0.0495** — borderline but below 0.05. Correction boxes added to both sections with the caveats. |
> | **F14** uncorrected axis correlations | **FIXED** — `correlate_axes()` now emits `spearman_q`/`pearson_q`. Confirmed still latent: 0 of 34 tests reach p<0.05 uncorrected (min p=0.118). |
> | **F15** doc drift | **PARTLY FIXED** — `docs/README.md` index added naming the canonical doc per topic; `analysis/README.md` still stale |
> | **stale derived data** | **ACTION NEEDED** — the statistics code is fixed but the CSVs on disk were written by the old code. `analysis/analyses/mechanism/data/same_vs_other_contrast.csv` still holds the non-monotonic `q_bh`; `axis_correlations.csv` has no `q` columns at all; every `output/stage7_eval/id_benchmark_*.json` still holds normal-approximation p-values. Re-run the analysis scripts to regenerate. Numbers quoted from those files before regeneration are the old ones. |
> | Part 4 cleanup | **DONE** — 11 duplicate files, 30 Zone.Identifier, 29 `__pycache__`, 3 scratch `view_test*`, 4 `_deprecated/` trees removed; 11 `_dnarelax` configs collapsed into `ARM=dnarelax` |

---

## Summary

The pipeline is in better shape than the analysis layer, and the analysis layer is in better shape than its own documentation. Three findings are serious enough to change published numbers:

| # | Finding | Severity |
|---|---|---|
| **F1** | The TBP→NDT80 label error is **still live** in `config/pilots/tbp.sh` | Critical |
| **F2** | That same mislabel **propagated into the analysis layer** and contaminates the mechanism analysis's central ETS1/TBP contrast | Critical |
| **F6** | **None of the 264 training runs was ever seeded.** The "paired seed" design is nominal — `random_seed` is a dead config key | Critical |
| **F7** | Stage 4 silently drops states. **IRF trains on 5 conformers out of 78**; frozen/relaxed arms get unequal training sets | Critical |
| **F13** | The mechanism headline (ETS1 own-family, p=0.013) **does not survive** multiplicity correction across the 11 pilots in its own table (q=0.143) | High |

The good news on your specific worry: the data/label mismatch class is **mechanically checkable**, and when I ran that check across all 14 pilots, **TBP is the only one wrong**. Details and the check itself are in F4.

---

# Part 1 — Initial data and labels

## F1 · The TBP mislabel is documented but never fixed — CRITICAL

`docs/tbp_label_audit.md` (2026-10-05) correctly diagnoses that `config/pilots/tbp.sh` labels a TBP structure (1tgh) with NDT80's motif. **The config was never corrected.** It still reads:

```bash
export PWM_LABEL="MA0343.1.jaspar"   # NDT80's motif, on a TBP structure
export TEST_PWM_FILTER="MA0343"      # scores TBP against NDT80 test entries
```

`config/pilots/tbp_dnarelax.sh` sources `tbp.sh`, so it inherits both unchanged.

A new `config/pilots/ndt80.sh` (2euv, `MA0343.1`) was added — that is recommendation #2 from the audit, and it is correct. But `tbp.sh` and `ndt80.sh` now resolve to **identical** `TEST_PWM_FILTER` values, so they evaluate on the same five NDT80 entries.

**CONFIRMED.** `tbp` and `ndt80` both match exactly `2evi, 2evf, 2evj, 2etw, 2euz`.

**Risk:** a known-poisoned config sitting in `config/pilots/` is one `run_pilot.sh tbp` away from regenerating contaminated results. Either fix it to `PWM_LABEL="TBP_HUMAN.H11MO.0.A"` / `TEST_PWM_FILTER="MA0386|TBP_"` (accepting that this leaves n=1, as the audit explains), or delete `tbp.sh`/`tbp_dnarelax.sh` and let `ndt80.sh` carry that structure family. Leaving it as-is is the one option with no upside.

## F2 · The mislabel propagated into the analysis layer — CRITICAL

This is the direct answer to "a mismatch that persisted and messed up the results." It is **not** confined to the config.

`analysis/data/family_annotation.csv` assigns the five NDT80 structures to family **`TBP / β-saddle`** — while, in the same rows, correctly recording `tf_name="NDT80 protein"` and `pfam=PF05224` (the NDT80/PhoG DNA-binding domain):

```
2euz_A_MA0343.1   tf_name=NDT80 protein   pfam=PF05224   family=TBP / β-saddle   family_src=curated
2evj_A_MA0343.1   tf_name=NDT80 protein   pfam=PF05224   family=TBP / β-saddle   family_src=curated
2etw_A_MA0343.1   tf_name=NDT80 protein   pfam=PF05224   family=TBP / β-saddle   family_src=curated
2evf_A_MA0343.1   tf_name=NDT80 protein   pfam=PF05224   family=TBP / β-saddle   family_src=curated
2evi_A_MA0343.1   tf_name=NDT80 protein   pfam=PF05224   family=TBP / β-saddle   family_src=curated
```

The annotation pipeline got it right; a hand-curation (`family_src=curated`) overrode it.

The sixth and only other member of that bucket is `7oh9_K_MA0386.1` — `tf_name="Histone H3.2"`, `pfam=PF00125` (core histone). 7oh9 is a nucleosome/TFIIA/TBP assembly and chain K is a histone, not TBP.

**So the `TBP / β-saddle` family bucket contains 5 NDT80 chains and 1 histone chain, and zero TBP protein chains.**

**CONFIRMED,** and it propagates all the way into the mechanism analysis: `analysis/data/perseed_perentry.csv` — the direct input to `analysis/analyses/mechanism/mechanism_analysis.py` — carries `family="TBP / β-saddle"` on all five NDT80 entries across **1,285 rows**.

Consequences:

- The `tbp` row of the own-family table in `analysis/analyses/mechanism/FINDINGS.md` (ΔPearson −0.031, p=0.513) is computed over 5 NDT80 entries + 1 histone chain.
- `FINDINGS.md` argues the static induced-fit index fails because *"This puts ETS1 and TBP in the same group, which is precisely the pair the hypothesis uses to separate the regimes."* **That contrast pair is not usable**: the TBP arm was trained against NDT80's motif (F1) *and* its own-family readout is NDT80 (F2).
- Anything grouping by `family` and reporting a "TBP" row inherits this — including the P2/P3 family figures in `analysis/README.md`, whose *"TBP and bHLH worst"* claim is about a bucket with no TBP in it.

**Recommended:** relabel those five rows to a distinct `NDT80` family, drop or re-derive the `tbp` row everywhere it appears, and regenerate `perseed_perentry.csv` and everything downstream of it.

## F3 · A second, independent curation error: HAP1 — Low

`1le8_B_MA0328.2` is filed as family `Homeodomain`. It is HAP1 — a yeast Zn(2)-Cys(6) binuclear-cluster TF, `pfam=PF00172`, and `MA0328.2` is HAP1's own motif. Both the structure and the motif say Zn2Cys6; only the hand-set `family` column says Homeodomain (`family_src=curated`).

**CONFIRMED.** Effect is modest — it sits among DUX4/MEIS2/CDX2 in an 11-member Homeodomain bucket and should move to `Zn2Cys6 cluster` (currently 2 members).

## F4 · The invariant that catches this class of error — and the all-pilots result

Two mechanical checks catch every instance of this bug class. I ran both.

**Check A — config label vs. the master PWM table.** Each pilot's `PWM_LABEL` must appear in `all_pwms` for its `PDB_ID` in `config/deeppbs_pdb_pwms.csv`:

| pilot | PDB | PWM_LABEL | in CSV? |
|---|---|---|---|
| csl, dux4, egr1, engrailed, err, ets1, foxa, hsf, irf, lef1, ndt80, nfat, runx | — | — | **YES** (13/13) |
| **tbp** | 1tgh | `MA0343.1.jaspar` | **NO** — CSV says `TBP_HUMAN.H11MO.0.A` |

**CONFIRMED: TBP is the only config-level mismatch among all 14 pilots.** Your worry that others might be lurking does not bear out at the config layer. Note the implication: `deeppbs_pdb_pwms.csv` was right all along, so the root cause is that pilot configs are **hand-written rather than generated from the CSV**, with no check tying them together. That check is ~10 lines and worth adding to `run_pilot.sh` as a pre-flight assertion.

**Check B — annotation consistency.** Where an entry's own motif and its Pfam agree on a family, the `family` column must agree too. Of 130 entries, 13 show a `pfam`/`family` mismatch — and they split cleanly:

- **7 are the convention working correctly** (motif-level assignment, as `analysis/README.md` documents): the five ETS1–RUNX1 co-crystal entries (`4l0z_B`, `3wtu_C`, `3wtw_H`, `3wtw_C`, `3wts_C`) carry a Runt `pfam` from the crystal's first chain but an ETS1 motif, so `family=ETS` is right; likewise `1hjb_F` (bZIP chain, RUNX1 motif → `family=Runt`).
- **6 are genuine errors**: the 5 NDT80 rows (F2) and HAP1 (F3).

**CONFIRMED.** This is a good invariant to encode as a test, because it has no false positives on the current data once the motif-level convention is accounted for.

**Latent risk in the same file:** the motif-level fix was applied to the `family` column **only**. `tf_name` and `pfam` are still crystal-level, so `pfam=PF00853` splits 30 ETS / 24 Runt. Live consumers (`make_P.py`, `extract_metrics.py`) correctly read `family` — and `make_P.py:32` even asserts one family per entry, which is exactly the right guard — so this is **not** currently producing wrong numbers. But any future script that groups by `pfam` or `tf_name` silently reproduces the bug `analysis/README.md` says was fixed. Either drop those two columns or rename them `pfam_crystal` / `tf_name_crystal`.

## F5 · `TEST_PWM_FILTER` means four different things across pilots — High

The "on-target subset" is not a consistent construct. **CONFIRMED** against `data/folds/id.txt` (130 entries):

| pilot | n matched | what the subset actually is |
|---|---|---|
| egr1, ets1, dux4, csl, runx | 5 each | same TF / close paralogs — reasonable |
| hsf | 3 | 2 unique PDBs |
| **irf** | 4 | **all 4 chains of a single PDB (7jm4)** |
| **err** | 10 | 7 PDBs spanning much of the nuclear-receptor superfamily (RXRA, GCR, NR5A, ESRR) |
| **nfat** | 15 | 7 PDBs, **13 of 15 are NF-κB** (TF65/RelA, NFKB1/p50, NFKB2/p52), not NFAT |
| **engrailed, foxa, lef1** | **0** | no on-target readout at all |
| **tbp** | 5 | the wrong family entirely (F1) |

Three distinct problems:

1. **Pseudoreplication.** IRF's "n=4" is four chains of one crystal; HSF's "n=3" is three chains of two. `compute_paired_stats` treats entries as independent, so these n's are inflated and the p-values correspondingly deflated. The effective structural n is 1–2.
2. **Over-broad subsets.** Calling NF-κB entries NFAT's "on-target" set tests family transfer at best; `nfat` is a 1a66 NFAT augmentation scored predominantly on a different protein family. Same pattern for `err` across the nuclear-receptor superfamily. If these come out positive it is not evidence of specificity transfer.
3. **Recoverable n=0 subsets.** `foxa`'s filter `MA1683|FOXA` matches nothing — but the benchmark holds **5 Forkhead entries**, all named `FOXO1_HUMAN`. The filter simply misses them; `FOXO|FOXA|MA1683` recovers n=5. `engrailed`'s `MA0220` matches nothing while the benchmark holds 11 Homeodomain entries. `lef1` is genuinely unrecoverable — it is the only HMG-box in the benchmark.

Also worth stating plainly: **most pilots run at n=5**, and `docs/tbp_label_audit.md` already concluded that *"a 5-entry subset at these effect sizes is not a reliable readout."* That lesson generalizes to every pilot, not just TBP.

---

# Part 2 — Pipeline

The pipeline architecture is genuinely solid: self-locating paths, `require_var` guards, one wrapper per stage, a sequence-match guard in Stage 2, filename-derived state labels in the fnat gate (an explicit anti-drift measure — good), and a metal-coordination cage in Stage 3. The issues below are specific, not structural.

## F6 · No training run was ever seeded; the "paired seed" design is nominal — CRITICAL

`stage5_build_aug/build_training_configs.py:58-62` writes:

```python
# CRITICAL: paired seed so baseline and augmented share parameter init
# and shuffle order. The augmentation effect is then isolated from
# init/shuffle noise.
"random_seed": random_seed,
"no_random": False,
```

Two independent failures:

1. **`random_seed` is a dead key.** `grep -n random_seed stage6_train/driver.py` returns **nothing**. The driver never reads it. Its only seeding is hardcoded and ignores the config value entirely (`driver.py:144-147`):
   ```python
   if C["no_random"] or C["debug"]:
       np.random.seed(8)
       torch.manual_seed(0)
   ```
2. **That block never executes.** It is gated on `no_random`, which the generator writes as `False`, and `debug` defaults to `False`.

**CONFIRMED across all 264 real training configs** in `output/stage6_train/`: `random_seed` ∈ {1,2,3,4,5,42}, `no_random=False` in **264/264**, `debug` absent in **264/264**. `no_random` is never set `True` anywhere in the repo.

**Therefore `torch.manual_seed` was never called for any training run in this project.**

What this does and does not invalidate — this distinction matters:

- **Unaffected:** the data split. Folds come from fixed files (`train0.txt`/`valid0.txt`), so splits are identical across runs.
- **Unaffected:** per-entry pairing inside a run. `compute_paired_stats` pairs by *benchmark entry* — both models score the same entries — and that pairing is real.
- **Invalid:** the baseline↔augmented pairing at a given seed. `baseline_ets1_s1` and `augmented_ets1_s1` share **nothing** — not parameter init, not shuffle order. The suffix `s1` is a filename, not a controlled variable. So Δ at each "seed" carries the full init+shuffle noise the comment claims to have removed.
- **Invalid:** reproducibility. Re-running any config produces different numbers.
- **Salvageable:** five unseeded replicates *are* a legitimate run-to-run variance estimate. The cross-seed SD you report is real. It is simply not a paired design, and per-seed Δ is noisier than documented.

This is item #5 of `docs/concerns.txt`, which reads as fixed (a `random_seed` key was added) but is not. Fix: set `"no_random": True` **and** change `driver.py:145-147` to seed from `C["random_seed"]` rather than the constants 8/0.

## F7 · Stage 4 silently drops states; IRF loses 94% — CRITICAL

`stage4_preprocess/stage4_preprocess.sh` compares `CONTACT COUNT` log lines against npz count, but **never compares gated inputs against npz outputs**. `process_co_crystal.py` does `continue` on failure, so attrition is silent. **CONFIRMED** from paired `Linked …` / `DONE …` lines in `slurm_output/`:

| pilot | arm | gated in | npz out | % lost |
|---|---|---|---|---|
| **irf** | frozen | 78 | **5** | **94%** |
| **irf** | dnarelax | 74 | **7** | **91%** |
| nfat | dnarelax | 13 | 8 | 38% |
| nfat | frozen | 26 | 19 | 27% |
| err | dnarelax | 37 | 29 | 22% |
| lef1 | frozen | 16 | 13 | 19% |
| err | frozen | 46 | 41 | 11% |
| foxa | both | 75 / 66 | 73 / 64 | 3% |
| csl, egr1, engrailed, ets1, hsf, runx, tbp | both | — | — | 0–3% |

**The IRF pilot's augmentation is built from 5 conformers.** `analysis/README.md` reports IRF as one of only two net-positive families (*"IRF (n=4, +0.039)"*) — that is 5 training conformers evaluated on 4 chains of one PDB (F5). That number should not be load-bearing.

**Root cause CONFIRMED:** `grep "ERROR: helix count problem" slurm_output/preprocess_*.out` shows 139 occurrences for `1if1` (IRF), all with count `2`. 3DNA resolves IRF's duplex into two helical segments, so `process_co_crystal.py` fails its `assert len(dna_helices) == 1` and skips the state. Other pilots show the same failure at far lower rates (1lo1: 11, 1a66: 11, 1vtn: 6, 2lef: 4, 1k79: 2, 3hdd: 1).

**This is a selection bias, not random loss** — INFERRED but strongly supported. The failure triggers on states whose DNA deviates enough from a single ideal helix to be split in two. It therefore preferentially removes the most-deformed conformations, which is exactly the conformational diversity the project is about. The surviving 5 IRF states are a biased sample of the most crystal-like frames.

**Secondary: the two arms get unequal training sets.** Since frozen-vs-relaxed is a headline comparison (the audit's §D "the DNA-relaxed arm is robustly harmful"), differential set size is a confound — the relaxed arm is not only relaxed but smaller and differently filtered:

| pilot | frozen npz | relax npz | ratio |
|---|---|---|---|
| **nfat** | 19 | 8 | **0.42** |
| **err** | 41 | 29 | **0.71** |
| **foxa** | 73 | 64 | **0.88** |
| **irf** | 5 | 7 | **1.40** |
| others | — | — | 0.92–1.02 |

Fixes: (a) make Stage 4 fail or loudly warn when `n_outputs < n_inputs`; (b) record per-state drop reasons to a CSV so attrition is auditable; (c) for the frozen/relaxed comparison, restrict both arms to the intersection of states that survived in both.

*(One clarification so the numbers aren't misread: `grep` also shows 95 `helix count problem 0` failures for `5d5u`/HSF, but those come from an aborted 2026-07-21 run that never reached `DONE`. HSF's current runs lose 0 states.)*

## F8 · The evaluator prints only 1 of 6 metrics — Medium

In `stage7_eval/evaluate_id_benchmark.py:362-369`, the `print` is dedented out of its loop:

```python
for k in KEYS:
    if k not in stats: continue
    s = stats[k]
    sig = "*" if s['p_value'] < 0.05 else " "
print(f"    Δ{k:18s}  mean={s['mean_delta']:+.4f}  ...")   # <-- outside the loop
```

Only the final key in `KEYS` is ever printed — `brier_multi`. **CONFIRMED** in real logs: `slurm_output/eval_5163028.out`, `eval_4310682.out`, `eval_4606462.out` each show exactly one `Δbrier_multi` line per pair and nothing else.

The JSON payload contains all six metrics, so analyses reading the JSON are unaffected. But every human who read a Stage-7 console log saw one metric and may have read it as the summary. Fix is one indent level.

## F9 · p-values use a normal approximation at n=5, understating them 2.6–15× — High

`evaluate_id_benchmark.py:196-203` computes the t-statistic correctly, then converts it with the **normal** CDF:

```python
p_two_sided = 2 * (1 - 0.5 * (1 + erf(abs(t_stat) / sqrt(2))))
```

The comment concedes this (*"approximate via normal for simplicity (small n is unfortunate…)"*), but the consequence is large at the n the pilots actually run. **CONFIRMED** numerically:

| n | t | p as coded | correct p (t, df=n−1) | understated by |
|---|---|---|---|---|
| 5 | 2.0 | 0.0455 | 0.1161 | 2.6× |
| 5 | 2.5 | 0.0124 | 0.0668 | 5.4× |
| 5 | 3.0 | 0.0027 | **0.0399** | 14.8× |
| 10 | 2.5 | 0.0124 | 0.0339 | 2.7× |
| 130 | 2.5 | 0.0124 | 0.0137 | 1.1× |

Most on-target subsets are n=5, so **every "significant" subset p-value in the Stage-7 logs and in reports derived from them is overstated**, some by an order of magnitude. The general-130 numbers are essentially unaffected. Fix: `scipy.stats.t.sf(abs(t_stat), n-1) * 2`.

Two related notes on the same function:
- `docs/tbp_label_audit.md` states that its table-A p-values are *"the evaluator's own (10k paired bootstrap)."* They are not — they are this normal-approximation t-test. The bootstrap produces only the CIs. The audit's own caveat mis-describes its provenance.
- `rng = np.random.default_rng(seed=42)` is re-created **inside** the per-metric loop, so all six metrics bootstrap on identical resample indices. Their CIs are perfectly correlated rather than independent. Minor, but it makes "4 of 6 metrics agree" weaker evidence than it looks — and those six metrics (pearson, spearman, IC-weighted PCC, AUROC, MAE, Brier) are themselves highly correlated measures of the same thing.

## F10 · `interface_rmsd.score` silently falls back to the wrong coordinate frame — Medium

`fnat_gate/interface_rmsd.py:174`:

```python
dna = md if (use_model_dna and len(md) == len(R['rd'])) else R['rd']
```

If `use_model_dna=True` but the nucleotide counts differ, it silently uses the **reference** DNA — computing contacts against crystal DNA coordinates for a model whose DNA has moved. That yields a spuriously low fnat and no warning.

**INFERRED** (I did not catch it firing). It matters most for the `dnarelax` arm, where model DNA genuinely differs from reference, and it is a candidate explanation for part of the arm asymmetry already visible at the gate (nfat 26 frozen vs 13 relaxed states passing; err 46 vs 37; foxa 75 vs 66). Fix: raise on the mismatch instead of falling back, or at minimum log it per state.

## F11 · Status of the 11 items in `docs/concerns.txt`

That file is a prior review and is now partly stale. Verified against current code:

| # | Concern | Status |
|---|---|---|
| 1 | DNA carried over rigidly from crystal | **Stands by design** — the `_dnarelax` arm was added to test exactly this. Fair framing caveat, not a bug. |
| 2 | `removeHeterogens` silently strips structural metals | **FIXED** — `stage3_minimize.py` records metal positions before removal and applies a coordination cage (`--metal-cage-k`, `--metal-coord-cutoff`), with `--ignore-metals` to reproduce old behaviour. |
| 3 | Implicit solvent (gbn2) over-compacts a charged interface | **Stands** — `ForceField("amber14-all.xml", "implicit/gbn2.xml")`. Method choice; worth stating in the methods text. |
| 4 | All augmented frames share one PWM target | **Stands by design.** Combined with F7's bias toward crystal-like survivors, the up-weighting concern is sharper than when written. |
| 5 | No fixed random seed | **NOT FIXED** — see **F6**. Appears fixed; isn't. |
| 6 | Best-state checkpoint selected on general validation, not the pilot's TF | **Stands** — all 264 configs use `best_state_metric=mae` on `validation`. Fair across arms, as the note says; absolute numbers are not TF-optimal. |
| 7 | Per-condition denominators can differ silently | **FIXED** — `check_denominators()` warns and restricts to the intersection. |
| 8 | Effect-size pairing by regex drops odd names silently | **FIXED** — unpaired and incomplete pairs are now reported. |
| 10 | Cα count check too loose | **FIXED** — Stage 2 now uses `find_common_subsequence` + sequence identity rather than a bare count. |
| 11 | Zero-mass freeze brittle if dynamics added | **Stands** — still `setParticleMass(i, 0.0)`. Fine for L-BFGS; a trap if anyone adds Langevin. |

(Item 9 is missing from the file — the text jumps from 8 to a stray "be." then 10.)

---

# Part 3 — Analysis

First, credit where due: `analysis/analyses/mechanism/mechanism_analysis.py` is the most statistically careful code in the repo. It uses the proper t-distribution (`stats.ttest_1samp`, `stats.t.ppf`), it applies BH-FDR, and it explicitly averages across pilots within a seed *"so pilots do not inflate n"* — a deliberate anti-pseudoreplication step that F5 shows the Stage-7 evaluator lacks. `analysis/README.md` also documents a previously-caught family-assignment bug and the figure it corrected, which is exactly the right practice.

## F12 · The BH-FDR implementation is missing its monotonicity step — Medium

`mechanism_analysis.py:299-304`:

```python
ranked = p[ok].rank(method="first")
out.loc[p[ok].index, "dP_q"] = (p[ok] * ok.sum() / ranked).clip(upper=1)
```

Benjamini–Hochberg requires a step-up cumulative minimum over p-ordered values; this omits it, so q-values can be non-monotonic in p. **CONFIRMED** with a worked example:

| p | q as coded | correct BH q |
|---|---|---|
| 0.001 | 0.0050 | 0.0050 |
| 0.030 | **0.0750** | 0.0667 |
| 0.040 | **0.0667** | 0.0667 |

p=0.030 receives a *worse* q than p=0.040. At a 0.05 threshold a more-significant test can be rejected while a weaker one passes. Fix:

```python
o = np.argsort(p.values)
q = np.minimum.accumulate((p.values[o] * m / np.arange(1, m+1))[::-1])[::-1]
```

## F13 · The mechanism headline does not survive correction across its own table — High

`FINDINGS.md` opens with *"The augmentation signal is real, but it is not where the benchmark was looking,"* resting on ETS1 own-family ΔPearson **+0.111, p=0.013**. That p is one of **11 pilot tests presented in the same table**, reported raw. Applying BH-FDR across them (**CONFIRMED**):

| pilot | ΔPearson | raw p | BH q |
|---|---|---|---|
| ets1 | +0.111 | 0.013 | **0.143** |
| csl | +0.103 | 0.065 | 0.358 |
| runx | +0.084 | 0.161 | 0.590 |
| hsf | +0.037 | 0.331 | 0.806 |
| irf | +0.033 | 0.471 | 0.806 |
| …remaining 6 | ≤ +0.015 | ≥0.500 | ≥0.806 |

**Nothing survives q<0.05.** Note the irony: the script applies FDR to its *family*-level table but not to this *pilot*-level one.

Two further caveats compound it:

- **The unit of replication is a "seed" that isn't seeded.** The test is over 5–6 unseeded replicates (F6) of a single structure per pilot. So p=0.013 asks "is this model's improvement consistently positive across reruns," not "does augmentation help ETS-family TFs." With n=1 structure, the inference to a general claim is a bigger leap than the p-value suggests.
- **The `tbp` row is contaminated** (F2), and the `irf` row rests on 5 training conformers (F7).

**Recommended wording:** "ETS1 shows the largest positive own-family effect (+0.111, raw p=0.013), which does not survive multiplicity correction across the 11 pilots (q=0.14); it is the strongest lead in the set, not an established result." That is still a genuinely interesting finding — the own-family-vs-cross-benchmark contrast (+0.111 vs −0.002) is a real and useful observation about where to look. It just isn't significance.

## F14 · `axis_correlations.csv` gets no multiplicity correction — Low (latent)

`correlate_axes()` tests **17 axes × 2 DNA conditions × 2 correlation types = 68 tests**, and `main()` writes them straight to CSV with no FDR — the correction is applied only inside `family_effects_seedlevel()`.

**Important and reassuring: this is currently harmless.** I checked the actual output — **0 of 34 Spearman tests reach p<0.05 even uncorrected**, at n=11–12 pilots. Nothing is being over-claimed, and `FINDINGS.md` honestly reports the static index as failing. Flagging it because the moment one of those correlations does come in at p=0.04, it will look like a result and ~1.7 false positives are expected by chance at that test count.

## F15 · Documentation drift in the analysis layer — Medium

Not bugs, but they make the analysis hard to trust at a glance:

- **`analysis/README.md` is stale.** It describes a **six-pilot** study (*"Pilots span six DNA-binding-domain families"*) and says *"Single-seed results (5-seed re-run in flight)."* There are now 14 pilot configs and the 5-seed runs are long done. Its headline P3 claim is single-seed and includes the contaminated TBP bucket.
- **Inconsistent pilot coverage across `analysis/data/`.** `family_annotation.csv` covers **6** pilots; `perentry_accuracy.csv` and `perseed_perentry.csv` cover **12**; `config/pilots/` defines **14**. Nothing records which is current. **CONFIRMED.**
- **`docs/` has 46 markdown files and no index.** Six versions of the analysis plan (`_2026-09-29`, `_v3_`, `_v4_`, `_v5_`, `_v6_`, `findings_v7_`), all carrying the same `2026-09-29` date, all introduced in a single commit (`e54c5fa`) — so neither the filename nor git history identifies the canonical one. This is the "so many versions" problem directly: the ambiguity is the risk, because a reader cannot tell which numbers are current.
- **`INTERPRET_SUITE_README.md`** cites `interpret_tfconf_all.py` as "(210 lines)" and `compare_importance_all.py` as "(100 lines)"; actual counts are 149 and ~100.

---

# Part 4 — Cleanup, refactoring, deletions

## Safe deletions

**1. Ten sets of byte-identical duplicate files — an incomplete migration.** Root-level `analysis/*.py` and `analysis/*.md` were copied into `analysis/analyses/importance/` and `analysis/docs/` and the originals were never removed. **CONFIRMED by md5.** All executable runners (`submit_interpret.sh`, `run_all_interpretations.sh`, `submit_figs.sh`) already reference the canonical `analysis/analyses/importance/` paths — **only stale docs point at the root copies**, so the root copies are safe to delete:

```
analysis/interpret_tfconf.py              analysis/interpret_tfconf_all.py
analysis/compare_importance.py            analysis/compare_importance_all.py
analysis/make_importance_figures.py
analysis/RESULTS_INVENTORY.md             analysis/IMPORTANCE_ANALYSIS_RESULTS.md
analysis/ATOM_IMPORTANCE_DELIVERABLES.md  analysis/mechanism_and_roadmap.md
analysis/conditions_inventory.md          (= analysis/docs/CONDITIONS_INVENTORY.md)
```
`ATOM_IMPORTANCE_DELIVERABLES.md` exists in **three** places. `analysis/INTERPRET_SUITE_README.md` vs `analysis/docs/` differ by one line (a path) — keep the `docs/` one. After deleting, update the path references inside the surviving docs.

**2. Untracked junk (none of it in git, so no history to lose):**
- **30 `*:Zone.Identifier` files** — WSL/NTFS download artifacts, pure noise. (Correction: these *are* already covered by `*:Zone.Identifier` at `.gitignore:231`, verified with `git check-ignore`; an earlier draft of this review said otherwise. They just needed deleting, which has been done.)
- **29 `__pycache__/` directories, 248 `.pyc` files**, including a stray `__pycache__/` at repo root. Untracked; `.gitignore` should cover them.

**3. `slurm_output/` — 34 MB, 9,370 files, untracked.** It is the evidence base for F7/F8, so do not simply delete it. Worth archiving per-pilot tarballs and pruning the aborted runs (e.g. the 2026-07-21 HSF attempt), so the directory listing stops costing 389 KB of inode metadata.

**4. `_deprecated/` trees** (34 tracked files across `analysis/_deprecated/figure_scripts`, `analysis/figures/_deprecated`, `data/_deprecated`, `scripts/_deprecated`). Everything is in git history. Deleting them is safe and removes the main source of "which version is real" confusion. If you want a safety margin, tag the current commit first.

## Refactoring worth doing

1. **Generate pilot configs from `deeppbs_pdb_pwms.csv` instead of hand-writing them** — this is the root cause of F1. At minimum add Check A from F4 as a pre-flight assertion in `run_pilot.sh`.
2. **Extract the paired-stats code into one shared module.** `evaluate_id_benchmark.py` (normal-approx, no FDR, no monotonicity) and `mechanism_analysis.py` (correct t, BH but non-monotonic) implement overlapping statistics with different bugs. One tested `analysis/common/stats.py` with `paired_delta()` and `bh_fdr()` fixes F9, F12, and F14 at once, and prevents the next divergence.
3. **Add the two invariants from F4 as a test file.** Both are cheap, both have zero false positives on current data, and both catch precisely the bug class that cost you the TBP pilot.
4. **Make Stage 4 account for every input state.** Emit `state,status,reason` per state; fail the stage if the loss exceeds a threshold. This would have surfaced F7 the day IRF ran.
5. **Write `docs/README.md`** naming the canonical plan/findings doc and marking the superseded versions. Cheapest high-value fix in the repo.
6. **`palette.py` at repo root duplicates the role of `analysis/common/fig_common.py`.** Consolidate into `analysis/common/`.

---

# Prioritized actions

**Before any further results are quoted:**

1. **F6** — fix the seeding (`no_random: True` + have `driver.py` read `C["random_seed"]`). Until then, no result is reproducible and no "paired seed" claim is supportable. *This is the most consequential item in this review.*
2. **F1 + F2** — fix or delete `tbp.sh`/`tbp_dnarelax.sh`; relabel the 5 NDT80 rows to their own family; regenerate `perseed_perentry.csv` and everything downstream; retire the ETS1/TBP contrast pair from `FINDINGS.md` and the book outline until it can be rebuilt on `ndt80`.
3. **F7** — add the input/output accounting; treat IRF (5 conformers) and NFAT (19 vs 8 across arms) as unusable until rerun; restrict frozen-vs-relaxed comparisons to the intersection of surviving states.
4. **F9 + F13** — switch to the t-distribution; restate the mechanism headline with FDR applied. Re-check any claim currently resting on an n=5 subset p-value.

**Next:**

5. **F5** — settle on one definition of "on-target," fix `foxa`'s filter, and either drop or explicitly relabel the superfamily subsets (`nfat`, `err`) and the single-PDB ones (`irf`, `hsf`).
6. **F8, F10, F12** — three small, contained code fixes.
7. **F3, F4-latent, F15** — HAP1 relabel; rename or drop the stale `pfam`/`tf_name` columns; refresh `analysis/README.md`; add `docs/README.md`.
8. Cleanup from Part 4.

---

## What I did not review

- **Vendored code:** `lib/deeppbs/` and the 3DNA toolchain were treated as external dependencies and not audited, except where pipeline code calls into them (the `assert len(dna_helices)==1` in F7 is in `process_co_crystal.py`, which *is* project code).
- **Nothing was executed on the cluster.** No training, evaluation, or structure processing was run. All quantitative claims come from reading committed files, stored CSV/JSON outputs, and archived `slurm_output/` logs.
- **Not independently verified:** the structural biology of individual PDB chain assignments (beyond the Pfam/motif cross-checks in F4), and the correctness of BioEmu/HPacker sampling.
- **`analysis/analyses/nmr_rmsd/`** (44 scripts, 230 plots — the largest single analysis) got only a survey pass. It is the biggest unreviewed surface in the repo and the natural next target.
- **One flagged-but-unconfirmed item:** F10's silent frame fallback is a code-path risk I did not catch firing.
