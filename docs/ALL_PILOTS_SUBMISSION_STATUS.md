# All-Pilots Attribution Analysis: Submission Summary

## Status: SUBMITTED ✓

**Date**: 2026-07-30
**All 12 pilots submitted** with batch_interpret_all.sh

## Jobs Submitted

| Pilot | Job ID | Status |
|-------|--------|--------|
| csl | 5303457 | QUEUED/RUNNING |
| egr1 | 5303458 | QUEUED/RUNNING |
| engrailed | 5303459 | QUEUED/RUNNING |
| err | 5303460 | QUEUED/RUNNING |
| ets1 | 5303461 | QUEUED/RUNNING |
| foxa | 5303462 | QUEUED/RUNNING |
| lef1 | 5303463 | QUEUED/RUNNING |
| nfat | 5303464 | QUEUED/RUNNING |
| runx | 5303465 | QUEUED/RUNNING |
| tbp | 5303466 | RUNNING ✓ |
| hsf | 5303467 | QUEUED |
| irf | 5303468 | QUEUED |

## What's Running

Each job computes **occlusion-based attribution** (per-atom importance via masking) for:
- Baseline model
- Augmented-frozen DNA model
- Augmented-relaxed DNA model (seed 1)

Per pilot:
- ~50-200 benchmark complexes (varies by TF)
- ~2-20 protein-DNA contact atoms (varies by structure)
- Estimated time: 30-60 minutes per pilot on GPU

## Monitoring

```bash
# Check queue
squeue -u shewchuk

# Track specific job
squeue -j 5303457

# View log (once running)
tail -f /project2/rohs_102/shewchuk/TF-conformation/slurm_output/interpret_csl_*.log

# Expected output
ls /project2/rohs_102/shewchuk/TF-conformation/output/interpret_results_all/*/
```

## Next Steps (Automatic with Master Workflow)

Once all jobs complete:

1. **Comparison**: `python analyze/compare_importance_all.py`
   - Computes shifts, p-values, effect sizes per pilot
   - Outputs: `interpret_results_all/all_pilots_importance_summary.csv`

2. **Master Workflow**: `bash run_all_interpretations.sh`
   - Waits for all jobs to finish
   - Runs comparison automatically
   - One-command solution for full analysis

## Scripts in TF-conformation/analysis

### Interpretation (GPU)
- `interpret_tfconf_all.py` — Generalized occlusion attribution
  - Auto-detects available arms per pilot
  - Usage: `python interpret_tfconf_all.py --pilot ets1 csl`

### Comparison
- `compare_importance_all.py` — Cross-pilot statistics
  - Loads all .npz importance files
  - Computes paired t-tests, effect sizes
  - Outputs summary CSV

### Submission
- `batch_interpret_all.sh` — Submit 12 GPU jobs
- `run_all_interpretations.sh` — Full workflow (submit + wait + compare)

## Output Structure

```
output/interpret_results_all/
├── csl/
│   ├── baseline_csl_fold0_importance.npz
│   ├── augmented_csl_fold0_importance.npz
│   └── augmented_csl_fold0_dnarelax_s1_importance.npz
├── egr1/
│   └── ... (same structure)
├── ... (10 more pilots)
└── all_pilots_importance_summary.csv  <-- FINAL OUTPUT
```

Each `.npz` contains per-atom importance scores (occlusion MAE) for all complexes.

## Key Metrics

**all_pilots_importance_summary.csv** will show:
- `mean_shift`: How much augmentation changes atom importance (vs baseline)
- `p_value`: Statistical significance of shift
- `effect_size`: Magnitude of shift (Cohen's d)

**Interpretation**:
- shift > 0 AND p < 0.05: Augmentation moves model toward real biochemistry ✓
- shift < 0 AND p < 0.05: Augmentation moves model away from biochemistry ✗
- p > 0.05: No significant effect (noise)

## Quick Stats (From ETS1 Single-Pilot Run)

For reference, the ETS1-only pilot showed:
- **Frozen DNA**: shift = +0.000489, p = 0.266 (not significant, marginal benefit)
- **Relaxed DNA**: shift = -0.000198, p = 0.026 (significant, counterproductive)

All-pilot analysis will reveal whether this pattern holds across all TFs.

---

**Contact**: shewchuk@rohs
**Project**: TF-conformation augmentation analysis
