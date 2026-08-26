# Figure-generation scripts — DeepPBS conformational-augmentation analysis

Regenerates figures 1–9 from the three-arm comparison (baseline DeepPBS vs
augmented-with-frozen-DNA vs augmented-with-relaxed-DNA), broken out by Pfam
family. Each `make_figN_*.py` writes one PNG into the current directory.

## Layout
    _common.py                     shared data load, colours, labels, helpers  (imported by every script)
    make_fig1_three_arm_accuracy.py    per-family mean Pearson, 3 arms  (mean ± 95% CI over seeds)
    make_fig2_augmentation_delta.py    seed-paired ΔPearson = aug − base, frozen vs relaxed   [PRIMARY]
    make_fig3_box_pearson.py           per-entry Pearson box-and-whisker (n=130), 3 arms
    make_fig4_box_mae.py               per-entry MAE box-and-whisker (n=130), 3 arms
    make_fig5_mae_delta.py             seed-paired ΔMAE, frozen vs relaxed
    make_fig6_within_family_transfer.py  same-family vs other-family ΔPearson, per pilot
    make_fig7_abs_pearson_groups.py    six-group per-entry Pearson (base/aug × all/other/same)
    make_fig8_abs_mae_groups.py        six-group per-entry MAE
    make_fig9_mixedmodel_effects.py    mixed-effects augmentation effect on own family (forest plot)

## Data inputs (place in this directory, or edit the paths at the top of _common.py)
    perentry_condition.csv    seed-mean metrics, one row per (tf, arm, dna, entry).   == artifact table4
    perseed_summary.csv       per-seed benchmark-mean, one row per (tf, arm, dna, seed).  == artifact table2
    perseed_perentry.csv      per-seed × per-entry metrics (family columns populated).  == the perseed checkpoint

All three are extracted from the Stage-7 evaluation JSONs in
`output/stage7_eval/` (the 10 standard `id_benchmark_<tf>.json` files plus the 9
`id_benchmark_<tf>_dnarelax.json` files; csl has no relaxed run). Each JSON's
`results` block maps checkpoint name → `<entry>.npz` → 8 metrics; checkpoints are
named `{baseline,augmented}_<tf>_fold0[_dnarelax][_sN]`, where `_sN` is training
seed N (bare `fold0` = seed 0, excluded from these figures to match the delivered
set). Family per benchmark entry comes from
`analysis/data/family_annotation.csv` (assigned at motif level).

## Running
    cd <this directory>
    PYTHONPATH=. python make_fig2_augmentation_delta.py      # or any other script

(`PYTHONPATH=.` lets the script find `_common.py`; alternatively run from a shell
where this directory is on sys.path.) Requires: numpy, pandas, scipy, matplotlib,
and — for fig9 only — statsmodels.

On endeavour, the `bioemu` and `deeppbs` conda envs have numpy/pandas/scipy/
matplotlib (figs 1–8 run as-is; the box plots are version-robust across
matplotlib 3.5–3.10 via the `hbox` helper). Neither env ships statsmodels, so
fig9 exits with an install hint until you add it, e.g.
`pip install statsmodels` (or `conda install -c conda-forge statsmodels`).

## Key conventions (defined once in _common.py)
- Two independent seed-paired experiments: the frozen and relaxed pipelines each
  retrained their OWN baseline, so absolute Pearson is not comparable across
  pipelines. The cross-treatment quantity is the within-pipeline, seed-matched
  augmentation effect ΔPearson = aug_sN − base_sN.
- Colours: grey = baseline, teal = augmented·frozen DNA, green = augmented·relaxed DNA.
- lef1 (HMG-box) has no other same-family benchmark entries, so it is excluded
  from the within-family figures (fig6–fig9); those iterate over `order9`.
- csl has no relaxed-DNA run — its relaxed cells are drawn as "n.d.".

## Notes on fig9
fig9 fits, per pilot, a linear mixed model on the pilot's same-family entries:
`metric ~ C(arm)*C(dna)` with crossed random intercepts on entry and seed
(baseline and frozen-DNA set as the reference levels, so the arm coefficient is
the augmentation effect). The relaxed effect is the arm main effect + the
arm:dna interaction (SE from the coefficient covariance). p-values are BH-FDR
corrected within each metric. This is the analysis that recovers power at n=5
(CSL, RUNX significant) where the seed-averaged Wilcoxon test is floor-locked.
