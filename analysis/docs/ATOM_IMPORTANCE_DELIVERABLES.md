# Atom-Level Importance Analysis: Final Deliverables

## Overview
Analysis of DeepPBS atom-importance scores across baseline, augmented-frozen, and augmented-relaxed DNA augmentation arms. Method: occlusion-based attribution measuring MAE shifts in predicted PWM when protein atoms are masked.

## Key Scripts
- **interpret_tfconf.py** — Computes per-atom occlusion MAE for all three arms on shared 99-complex ETS1 benchmark set
- **compare_importance.py** — Statistical analysis of importance shifts (paired t-tests, effect sizes)
- **make_importance_figures.py** — Publication-quality figures

## Output Files

### Analysis Results
- **importance_stats.json** — Comprehensive statistics (n_complexes, n_contact_atoms, mean shifts, p-values, effect sizes)
- **importance_comparison_table.csv** — Per-complex shifts (99 rows × 5 columns: baseline/frozen/relaxed importance, frozen/relaxed shifts)

### Figures
- **importance_shift_distributions.png** — Histograms of frozen and relaxed shift distributions
- **importance_boxplot_comparison.png** — Boxplot comparing baseline, frozen, relaxed importance distributions
- **importance_scatter_comparison.png** — Scatterplots showing relationship between baseline and augmented importance

### Documentation
- **IMPORTANCE_ANALYSIS_RESULTS.md** — Detailed interpretation, key findings, biochemical significance, recommendations

## Key Findings

| Metric | Frozen DNA | Relaxed DNA |
|--------|-----------|------------|
| Mean Importance Shift | +0.000489 | -0.000198 |
| p-value | 0.266 (NS) | 0.026 * |
| Effect Size (Cohen's d) | +0.113 | -0.228 |
| Biochemical Alignment | Marginal ✓ | Counterproductive ✗ |

**Conclusion**: Frozen-DNA augmentation shows weak directional benefit toward real biochemistry. Relaxed-DNA augmentation is harmful, training the model to ignore actual protein-DNA contact atoms.

## Methods Summary
- **Model**: DeepPBS with occlusion-based attribution (model_interpret variant)
- **Test Set**: 99 ETS1 benchmark complexes (1k79 X-ray, 12 DNA sequences)
- **Contact Atoms**: Identified via protein→DNA edges (2 atoms with strong contacts to DNA)
- **Metric**: Occlusion MAE — change in predicted DNA specificity (PWM) when atom is masked

## Location
All files in: `/project2/rohs_102/shewchuk/TF-conformation/analysis/`

## Generated
2026-07-30
