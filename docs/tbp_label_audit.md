# TBP pilot — PWM label audit
**Job** 6172207 (4m33s, 1 GPU) · 10 existing checkpoints (baseline/augmented × seeds 1–5) · **no retraining, nothing overwritten**
**Outputs** `output/stage7_eval/tbp_label_audit/{audit_ndt80_subset.json, audit_tbp_true.json}`

## The error, restated precisely
`config/pilots/tbp.sh` sets both of these to NDT80's motif:

```
export PWM_LABEL="MA0343.1.jaspar"   # -> stage4 Y_pwm target + stage5 augmentation filter
export TEST_PWM_FILTER="MA0343"      # -> stage7 on-target subset
```

MA0343.1 is **NDT80** (*S. cerevisiae*, JASPAR, UniProt P38830), not TBP. The comment in the config
("All 5 TBP test entries in id.txt use this PWM") is wrong: those five entries are
`2etw/2euz/2evf/2evi/2evj`, all **NDT80** structures (Pfam PF05224).

Because `PWM_LABEL` sets the **`Y_pwm` target written inside every npz** (`process_co_crystal.py`
→ `computeYAndMask`), the mislabel is not cosmetic — the 94 TBP conformers were trained against
NDT80's motif. And because `TEST_PWM_FILTER` selects the on-target subset, they were also *scored*
on NDT80 entries. **Both ends of the experiment were NDT80; only the structures were TBP.**

Direct evidence the target was wrong for this structure:

| | consensus | total IC |
|---|---|---|
| MA0343.1 (NDT80, as used) | `GGCCACAAAAA` | 12.98 |
| TBP_HUMAN.H11MO.0.A (correct) | `CTATAAAAGG` | 8.51 |
| **1TGH duplex (actual DNA)** | **`CGTATATATACG`** | — |

The crystal DNA is a canonical TATA box. It matches TBP's consensus and not NDT80's.

## What the audit shows
Re-scoring the **existing** checkpoints, sign-corrected so **+ = augmentation better**:

### A · NDT80 subset — the pairing the pilot actually ran (n=5 entries × 5 seeds)
| metric | mean Δ | seeds improving | median p |
|---|---|---|---|
| pearson r | +0.0062 | 2/5 | 0.15 |
| spearman ρ | −0.0177 | 1/5 | 0.056 |
| AUROC | +0.0221 | 4/5 | 0.024 |
| IC-weighted PCC | −0.0101 | 3/5 | 0.15 |
| Brier (multi) | +0.0141 | 2/5 | 0.072 |
| MAE | +0.0050 | 2/5 | 0.022 |

### B · General-130 (all id.txt entries, n=130 × 5 seeds)
| metric | mean Δ | seeds improving | median p |
|---|---|---|---|
| pearson r | −0.0116 | 2/5 | 0.084 |
| spearman ρ | −0.0171 | 2/5 | 0.33 |
| AUROC | −0.0065 | 1/5 | 0.55 |
| IC-weighted PCC | −0.0193 | 1/5 | 0.30 |
| Brier (multi) | −0.0146 | 2/5 | 0.30 |
| MAE | −0.0204 | 2/5 | 0.11 |

### C · True TBP-family subset
`TEST_PWM_FILTER` corrected to `MA0386|TBP_` matches exactly **one** id.txt entry
(`7oh9_K_MA0386.1`, yeast SPT15 in a nucleosome+TFIIA complex). The evaluator reports
*"too few paired entries; skipped"* for all five seeds. **n=1 — no paired statistic is possible.**

## D · The `_dnarelax` arm (job 6174398, 10 checkpoints, 4 min)
Same mislabel applies: `tbp_dnarelax.sh` sources `tbp.sh`, so it inherits
`PWM_LABEL=MA0343.1.jaspar` and `TEST_PWM_FILTER=MA0343` unchanged.

**NDT80 subset (n=5 × 5 seeds):**

| metric | mean Δ | seeds improving | median p |
|---|---|---|---|
| pearson r | +0.0218 | 2/5 | 0.027 |
| spearman ρ | +0.0532 | 4/5 | 0.017 |
| AUROC | −0.0082 | 2/5 | 0.014 |
| IC-weighted PCC | −0.0247 | 1/5 | 0.14 |
| Brier (multi) | +0.0264 | 3/5 | 0.090 |
| MAE | +0.0247 | 2/5 | 0.47 |

**General-130 (n=130 × 5 seeds):**

| metric | mean Δ | seeds improving | median p |
|---|---|---|---|
| pearson r | −0.0392 | **0/5** | 0.082 |
| spearman ρ | −0.0201 | 1/5 | 0.30 |
| AUROC | −0.0203 | 2/5 | 5.6e-04 |
| IC-weighted PCC | −0.0308 | **0/5** | 4.4e-03 |
| Brier (multi) | −0.0304 | **0/5** | 0.023 |
| MAE | −0.0472 | **0/5** | 2.7e-04 |

This is the one genuinely robust effect in the whole TBP pilot, and it is **negative and
general**: on the full 130-entry benchmark the DNA-relaxed augmented model is worse than its
paired baseline on every metric, in 0/5 or 1/5 seeds improving, with four metrics significant.
The harm is roughly 2–4× the frozen arm's (compare §B: −0.0116 r, −0.0204 MAE).

Note the sign split: on the 5-entry NDT80 subset the relaxed arm looks mildly *positive*
(ρ +0.053, 4/5 seeds), while on the general-130 it is uniformly negative. A 5-entry subset at
these effect sizes is not a reliable readout — which is the same lesson as §1 below.

**Consequence for the attribution analysis.** The `interpret_results_all` occlusion study scored a
`augmented_<tf>_fold0_dnarelax_s1` arm per pilot and concluded relaxed-DNA training moves the model
*away* from real contact atoms. For TBP that conclusion rests on a model trained against NDT80's
motif, so the TBP row of `all_pilots_importance_summary.csv` should be dropped or recomputed. The
general finding may well survive on the other 11 pilots — only TBP's contribution is contaminated.

## Reading
0. **The DNA-relaxed arm is robustly harmful on the general benchmark** (§D) — 0/5 seeds improving
   on four metrics, several significant. This is the most solid TBP-pilot result, but it is a
   statement about DNA relaxation in general, not about TBP specificity (the label was NDT80's).
1. **The headline TBP result does not survive at multiseed.** `pilot_assessment_memo` reports TBP as
   "uniformly worse, all significant" with on-target deltas of −0.263 (r), −0.180 (ρ), −0.299 (MAE).
   Across 5 paired seeds on the same subset those deltas are an order of magnitude smaller, mostly
   not significant, and **not uniform in sign** (AUROC and MAE favour augmentation, 4/5 and 2/5 seeds).
   The memo's figure is a single-seed result; the effect it describes is not reproducible across seeds.
2. **Whatever that effect was, it was not about TBP.** Both arms were NDT80-labelled, so the
   comparison never tested TBP specificity. The general-130 column is the only part of the TBP
   pilot that remains interpretable, and it is mildly negative and non-significant.
3. **The corrected experiment cannot be run as designed.** Properly labelled, TBP has one on-target
   test entry. There is no paired on-target test to recover.

## Consequences for the write-up
TBP is one of the two protagonists in the book outline — the rigid β-saddle that reaches its bound
pose yet is "hurt worst", used (with ETS1) to falsify protein-side explanations and force the
DNA-deformation mechanism. **That control pair does not hold up:** the harm is single-seed, and the
arm was mislabelled besides. The ETS1 side is unaffected by this audit and still stands on its own.

Recommended, in order:
1. Retire the TBP on-target claim from Ch. 15–18 and the memo. The general-130 result can stay if
   it is reported as mildly negative and non-significant.
2. **Run 2euv as the NDT80 pilot** (rank 1 in the candidate memo). MA0343.1 is genuinely NDT80's
   motif and the five test anchors are genuinely NDT80 — label, structure, and anchors finally agree.
   This converts the error into the strongest new pilot.
3. If TBP is still wanted as a pilot, it needs a different evaluation design: n=1 on-target means
   either general-130 only, or a within-TBP-family test set built from the train/valid pool
   (1cdw, 1c9b, 1jfi, 1ngm, 1ytf, … all carry TBP_HUMAN / MA0386.1).

## Caveats
- Both the frozen-DNA and `_dnarelax` arms were re-scored (seeds 1–5 each); neither was retrained.
- The general-130 p-values here are normal-approximation one-sample t-tests computed from the stored
  per-entry metrics, not the evaluator's own bootstrap; treat them as indicative. Subset p-values in
  table A are the evaluator's own (10k paired bootstrap).
- Both audits read the existing `combined_assembly_tbp`, whose augmented entries carry the NDT80
  `Y_pwm`. That is deliberate — the point was to measure what the pilot as run actually did.
