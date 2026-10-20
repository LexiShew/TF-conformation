# docs/ index

46 markdown files accumulated here with no index, including six same-dated versions of the
analysis plan. This names which document is current for each topic so a reader does not have
to guess. **Superseded files are kept, not deleted** — they record what was actually run.

Last updated 2026-10-15.

---

## Start here

| topic | current document |
|---|---|
| **Code review / known defects** | [`codebase_review_2026-10-15.md`](codebase_review_2026-10-15.md) — 15 findings across data labels, pipeline, analysis |
| **Pipeline overview & how to run** | [`../README.md`](../README.md) |
| **Adding / running pilots** | [`pilot_operations.md`](pilot_operations.md) — the registry, generator, and continuous queue |
| **Retired pilots** | [`retired_pilots.md`](retired_pilots.md) — dux4 and irf, why, and what was kept |
| **The TBP label error** | [`tbp_label_audit.md`](tbp_label_audit.md) (diagnosis) → [`tbp_label_error_downstream.md`](tbp_label_error_downstream.md) (full footprint + remediation status) |
| **Pipeline fix spec / invariants** | [`PIPELINE_FIXES.md`](PIPELINE_FIXES.md) |

## Analysis plans — the version chain

Six files share the date `2026-09-29` and were all committed together, so neither filename
nor git history identifies the latest. The chain, oldest to newest:

```
TF_conformation_analysis_plan_2026-09-29.md          v1
  └─ _review_2026-09-29.md                           critique of v1
TF_conformation_analysis_plan_v3_2026-09-29.md       v3
TF_conformation_analysis_plan_v4_2026-09-29.md       v4   → v4_audit_report.md
TF_conformation_analysis_plan_v5_2026-09-29.md       v5   → v5_gating_report.md, v5_filled_in.md
TF_conformation_analysis_plan_v6_2026-09-29.md       v6   → v6_report.md
TF_conformation_findings_v7_2026-09-29.md            v7   ← CURRENT  → v7_permutation_report.md
```

**Current: `TF_conformation_findings_v7_2026-09-29.md`.** Treat v1–v6 as historical.

## Results & findings

| document | status |
|---|---|
[`TF_conformation_findings_v7_2026-09-29.md`](TF_conformation_findings_v7_2026-09-29.md) | **current** findings
[`../analysis/analyses/mechanism/FINDINGS.md`](../analysis/analyses/mechanism/FINDINGS.md) | mechanism analysis — **headline needs restating**, see review F13 (ETS1 p=0.013 → BH q=0.143); `tbp` row contaminated, see F2
[`../analysis/README.md`](../analysis/README.md) | figure catalogue — **stale**: describes 6 pilots and "single-seed", and its P3 "TBP worst" claim is about a bucket with no TBP in it
[`TF_conformation_Ch15-18_results.md`](TF_conformation_Ch15-18_results.md) | book chapter draft — TBP on-target claims withdrawn
[`pilot_assessment_memo_07-08-2026.md`](pilot_assessment_memo_07-08-2026.md) | **withdrawn** TBP numbers (single-seed, NDT80-labelled)
[`2026-09-29_all_conditions_report.md`](2026-09-29_all_conditions_report.md), [`2026-09-29-crosspilot_report.md`](2026-09-29-crosspilot_report.md), [`conf_vs_crystal_pwm_report.md`](conf_vs_crystal_pwm_report.md) | condition/cross-pilot reports
[`nmr_rmsd_findings.md`](nmr_rmsd_findings.md) | NMR apo/holo validation; see also [`../analysis/analyses/nmr_rmsd/README.md`](../analysis/analyses/nmr_rmsd/README.md)

## Book / talk drafts

`TF_conformation_book_outline.md` · `TF_conformation_Ch9_minimization.md` ·
`TF_conformation_Ch10-11_featurization_training.md` · `TF_conformation_Ch10-18_figure_plan.md` ·
`TF_conformation_Ch15-18_results.md` · `TF_conformation_Ch19h_apo_validation.md`

Talk: `TF_conformation_presentation_structure_v5.md` (**current**, v4 superseded) ·
`TF_conformation_talk_prep_2026-09-29.md` · `TF_conformation_pretalk_runbook.md`

> The book outline's ETS1/TBP control pair **does not hold** — the TBP arm was trained against
> NDT80's motif and its own-family readout was NDT80. See `tbp_label_error_downstream.md`.

## SOX2 validation guides

`sox2_validation_guide.md` → `_v2` → `_v3` → **`_v4`** (current), plus addenda
`sox2_guide_addendum.md`, `sox2_guide_v2_addendum.md`, `sox2_guide_v3_addendum.md`.

## Background / method notes

`dna_shape_recognition_story.md` · `deepdnashape_addendum.md` ·
`idr_binding_conformation_review.md` · `2026-09029_bend_iqr_arm_matched_addendum.md` ·
`pilot_candidate_recommendations.md` · `status_august_2026.md` ·
`ALL_PILOTS_SUBMISSION_STATUS.md` · `terminal_window_report.md` ·
`2026-09-29-critique_response.md` · `2026-09-29-v3_execution_report.md` · `concerns.txt`

`concerns.txt` is a prior code review; its item-by-item status is tracked as finding **F11**
of `codebase_review_2026-10-15.md` (items 2, 7, 8, 10 fixed; 5 **not** fixed; rest stand).
