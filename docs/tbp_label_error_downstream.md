# TBP label error — downstream inventory and remediation status

**Error:** `config/pilots/tbp.sh` labelled a TBP structure (1TGH) with NDT80's motif
(`MA0343.1.jaspar`) and scored it against NDT80's test entries (`TEST_PWM_FILTER="MA0343"`).
**Diagnosed:** 2026-10-05 (`docs/tbp_label_audit.md`). **Config corrected:** 2026-10-15.

This document exists because the audit diagnosed the error but nothing downstream was
corrected at the time. It is the inventory of every place the error reached, what has been
fixed, and what still has to be regenerated.

> **The error's signature is greppable.** `PWM_LABEL` is embedded in every npz filename, so
> contaminated artifacts are self-identifying:
> ```bash
> grep -rl '1tgh_state_[0-9]*_MA0343' .          # any file listing contaminated npz
> ls output/stage4_npz*/tbp/output/ | grep MA0343 # the npz themselves
> ```
> A correctly labelled TBP npz is `1tgh_state_NNN_TBP_HUMAN.H11MO.0.A.npz`.

---

## Status summary

| Layer | Artifact | Count | Status |
|---|---|---|---|
| Config | `config/pilots/tbp.sh` | 1 | **FIXED** 2026-10-15 |
| Config | `config/pilots/tbp_dnarelax.sh` | 1 | **FIXED** (inherits, no change needed) |
| Analysis | `analysis/data/family_annotation.csv` | 30 rows | **FIXED** 2026-10-15 |
| Analysis | `analysis/data/perseed_perentry.csv` | 1,285 rows | **FIXED** 2026-10-15 |
| Stage 4 | `output/stage4_npz/tbp/output/*.npz` | 94 | **INVALID — regenerate** |
| Stage 4 | `output/stage4_npz_dnarelax/tbp/output/*.npz` | 96 | **INVALID — regenerate** |
| Stage 5 | `stage5_build_aug/folds_aug/train0_aug_tbp.txt` | 97 lines | **INVALID — regenerate** |
| Stage 5 | `stage5_build_aug/folds_aug/train0_aug_legacy_tbp.txt` | — | **INVALID — regenerate** |
| Stage 5 | `output/stage5_build_aug/combined_assembly_tbp` | 1 tree | **INVALID — regenerate** |
| Stage 6 | `output/stage6_train/*tbp*` run dirs | 23 | **INVALID — retrain** |
| Stage 6 | `output/stage6_train/config_*tbp*.json` | 22 | **INVALID — regenerate** |
| Stage 7 | `output/stage7_eval/id_benchmark_tbp.json` | 1 | **INVALID — re-evaluate** |
| Stage 7 | `output/stage7_eval/id_benchmark_tbp_dnarelax.json` | 1 | **INVALID — re-evaluate** |
| Stage 7 | `output/stage7_eval/tbp_label_audit/` | 2 json | **KEEP** — deliberately measures the error |
| Docs | 6 files with contaminated TBP claims | 6 | **ANNOTATE — see below** |

Nothing from Stage 4 onward was corrected in place. The 190 npz files carry NDT80's `Y_pwm`
inside them (written by `computeYAndMask` in `process_co_crystal.py`), so no amount of
relabelling fixes them — they have to be rebuilt.

---

## What was fixed, and how

### 1 · `config/pilots/tbp.sh` — corrected

```diff
-export PWM_LABEL="MA0343.1.jaspar"      # NDT80's motif
-export TEST_PWM_FILTER="MA0343"         # NDT80's test entries
+export PWM_LABEL="TBP_HUMAN.H11MO.0.A"  # correct for 1tgh
+export TEST_PWM_FILTER="MA0386|TBP_"    # genuine TBP-family entries
+export ON_TARGET_USABLE="no"            # n=1: no paired statistic possible
+export ON_TARGET_N_EXPECTED=1
```

Cross-checked against `config/deeppbs_pdb_pwms.csv`, which lists exactly
`TBP_HUMAN.H11MO.0.A` for `pdb_id=1tgh` — i.e. the master table was right all along and the
hand-written config diverged from it. That is the root cause, and it is what
`scripts/pipeline/validate_pilot.py` now checks.

`tbp_dnarelax.sh` sources `tbp.sh`, so it picks up the correction with no edit.

### 2 · Analysis annotation tables — corrected

The mislabel had propagated past the config into the analysis layer, where it was *more*
damaging because it silently redefined a family bucket. `analysis/data/family_annotation.csv`
assigned the five NDT80 structures to family `TBP / β-saddle` while, in the same rows,
correctly recording `tf_name="NDT80 protein"` and `pfam=PF05224`. The annotation pipeline got
it right; a hand-curation (`family_src=curated`) overrode it.

Relabelled `TBP / β-saddle` → `NDT80` for `2etw, 2euz, 2evf, 2evi, 2evj`:

| file | rows relabelled | `TBP / β-saddle` bucket |
|---|---|---|
| `analysis/data/family_annotation.csv` | 30 | 36 → 6 rows |
| `analysis/data/perseed_perentry.csv` | 1,285 | 1,542 → 257 rows |

A second, unrelated curation error was fixed in the same pass: `1le8_B_MA0328.2` (HAP1, a
fungal Zn(2)-Cys(6) TF, `PF00172`) was filed as `Homeodomain`. Relabelled to
`Zn2Cys6 cluster` — 6 rows in `family_annotation.csv`, 257 in `perseed_perentry.csv`.

Resulting family bucket sizes (per-pilot view, n=130 entries):

```
C2H2 zinc finger  24     bZIP               8     bHLH              4
Rel/NF-κB (RHD)   15     CSL/RBPJ           5     IRF               4
MADS-box          12     Runt               5     Zn2Cys6 cluster   3  (was 2)
ETS               10     NDT80              5  (new)
Homeodomain       10     p53-like           5     HSF               3
  (was 11)              Forkhead           5     Homeodomain (POU) 1
Nuclear receptor  10                              TBP / β-saddle    1  (was 6)
```

Pre-fix copies are kept as `*.pre_tbp_fix.bak` next to each file.

**`TBP / β-saddle` is now n=1**, and that one member is `7oh9_K_MA0386.1` — `tf_name="Histone
H3.2"`, `pfam=PF00125`, a histone chain in a nucleosome/TFIIA/TBP assembly. It carries TBP's
motif, so under the documented motif-level convention `family=TBP` is correct, but **no TBP
protein chain exists anywhere in the benchmark**. Any future "TBP family" result is a
statement about one histone chain.

---

## What still has to be regenerated

### Stage 4 — 190 npz files

```
output/stage4_npz/tbp/output/*.npz            94 files   (frozen arm)
output/stage4_npz_dnarelax/tbp/output/*.npz   96 files   (relaxed arm)
```

All named `1tgh_state_NNN_MA0343.1.jaspar.npz` and all carrying NDT80's `Y_pwm`. Rerun
stage 4 after the config fix; the new files will be named `..._TBP_HUMAN.H11MO.0.A.npz`, so
the old and new sets cannot be confused.

### Stage 5 — augmented folds and combined assembly

`stage5_build_aug/folds_aug/train0_aug_tbp.txt` (tracked, 516 lines, 97 of them contaminated
`1tgh_*_MA0343.1` entries) and `train0_aug_legacy_tbp.txt`. Plus the
`output/stage5_build_aug/combined_assembly_tbp` tree, whose augmented entries carry the
NDT80 `Y_pwm`.

### Stage 6 — 23 run directories, 22 configs

Every `output/stage6_train/{baseline,augmented}_tbp_fold0*` directory. Both arms must be
retrained: the baselines are unaffected by the label in principle, but they were paired
against contaminated augmented runs and — per finding F6 of
`docs/codebase_review_2026-10-15.md` — were never seeded, so they cannot be reused as
matched controls anyway. **Fix the seeding before retraining**, or the rerun inherits the
same non-reproducibility.

One directory is a dead partial run and can simply be deleted:
`baseline_tbp_fold0_dnarelax_s3.incomplete_1784044503`.

### Stage 7 — 2 evaluation JSONs

`id_benchmark_tbp.json` and `id_benchmark_tbp_dnarelax.json`.

**Keep `output/stage7_eval/tbp_label_audit/`.** Those two JSONs exist specifically to measure
what the mislabelled pilot did, so they are evidence, not contamination.

---

## Documentation carrying contaminated TBP claims

These need an explicit retraction note rather than a silent edit, since the numbers were
published:

| file | claim to retract |
|---|---|
| `docs/pilot_assessment_memo_07-08-2026.md` | TBP "uniformly worse, all significant"; on-target Δ of −0.263 (r), −0.180 (ρ), −0.299 (MAE). Single-seed, and the arm was NDT80-labelled. |
| `docs/TF_conformation_book_outline.md` | TBP as the rigid β-saddle protagonist that "reaches its bound pose yet is hurt worst" — used with ETS1 to falsify protein-side explanations. **The control pair does not hold.** |
| `docs/TF_conformation_Ch15-18_results.md` | TBP on-target results. |
| `docs/TF_conformation_analysis_plan_2026-09-29.md` · `_review_` | TBP pilot framing. |
| `docs/pilot_candidate_recommendations.md` | TBP/MA0343 references. |
| `analysis/analyses/mechanism/FINDINGS.md` | the `tbp` own-family row (ΔPearson −0.031, p=0.513) was computed over 5 NDT80 entries + 1 histone chain; and the argument that the static index fails *because* "this puts ETS1 and TBP in the same group" rests on that contaminated row. |
| `analysis/README.md` | P3 family figure: "TBP and bHLH worst" — that bucket had no TBP in it. |

`docs/tbp_label_audit.md` itself stays as-is: it is the diagnosis. One correction worth
noting inside it, though — its caveat says the table-A p-values are "the evaluator's own (10k
paired bootstrap)". They are not; they are the evaluator's normal-approximation t-test, which
at n=5 understates p by 2.6–15× (finding F9).

---

## Why it happened, and what now prevents it

**Root cause:** pilot configs were hand-written, and `PWM_LABEL` was never checked against
`config/deeppbs_pdb_pwms.csv`. The master table had the right answer for 1tgh the whole time.

**Contributing cause:** the same wrong label was then hand-curated into
`analysis/data/family_annotation.csv`, overriding correct automated annotation — so the
analysis layer agreed with the config and the error became self-consistent and invisible.

**Now enforced** by `scripts/pipeline/validate_pilot.py`, run automatically as a pre-flight
step in `scripts/pipeline/run_pilot.sh` (and runnable standalone over all pilots):

1. `PWM_LABEL` must appear in `all_pwms` for the pilot's `PDB_ID` in `deeppbs_pdb_pwms.csv`.
2. `TEST_PWM_FILTER` must match ≥1 `id.txt` entry, and the matched entries' families must be
   consistent with the pilot's own family — a filter matching a *different* family fails.
3. If `ON_TARGET_USABLE=no`, no on-target subset statistic may be reported for that pilot.
4. The pilot's source PDB must not appear in `id.txt` (test-set leakage).
5. Annotation invariant: where an entry's own motif and its Pfam agree on a family, the
   `family` column must agree too.

Check 1 alone would have caught this error on the day `tbp.sh` was written.
