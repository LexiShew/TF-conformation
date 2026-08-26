# DeepPBS Atom-Importance Analysis Suite

Complete, generalized framework for computing and comparing atom-level importance scores across all TF-conformation pilots.

## What It Does

Computes **occlusion-based attribution** for DeepPBS: measures how much each protein atom contributes to DNA binding specificity predictions. Compares three arms per pilot:
- **Baseline**: No augmentation
- **Augmented (Frozen DNA)**: DNA conformations fixed during training
- **Augmented (Relaxed DNA)**: DNA conformations relaxed via energy minimization

## Quick Start

### 1. Run All 12 Pilots
```bash
bash /project2/rohs_102/shewchuk/TF-conformation/analysis/analyses/importance/run_all_interpretations.sh
```

This submits 12 GPU jobs (one per pilot) and waits for completion (~12-24 hours), then compiles results.

### 2. Run Individual Pilots
```bash
cd /project2/rohs_102/shewchuk/TF-conformation
python analysis/interpret_tfconf_all.py --pilot ets1 csl --out output/interpret_results_all
python analysis/interpret_tfconf_all.py --pilot all --out output/interpret_results_all
```

### 3. Compare Results
```bash
cd /project2/rohs_102/shewchuk/TF-conformation
python analysis/compare_importance_all.py
```

Output: `output/interpret_results_all/all_pilots_importance_summary.csv`

## Scripts

### Core Attribution
- **`interpret_tfconf_all.py`** (210 lines)
  - Generalized occlusion-based attribution
  - Auto-detects available arms per pilot
  - Usage: `--pilot PILOT_NAME [PILOT_NAME ...] --out OUTPUT_DIR`
  - Supports: `--pilot all` to run all 12 pilots sequentially

### Comparison & Summary
- **`compare_importance_all.py`** (100 lines)
  - Loads all importance .npz files from a results directory
  - Computes per-pilot statistics (mean shift, t-tests, effect sizes)
  - Outputs: `all_pilots_importance_summary.csv`

### Batch Submission
- **`batch_interpret_all.sh`** — Submits 12 SLURM jobs (one per pilot) in parallel
- **`run_all_interpretations.sh`** — Master workflow: batch submit → wait → compare

## Output Structure

```
output/interpret_results_all/
├── csl/
│   ├── baseline_csl_fold0_importance.npz
│   ├── augmented_csl_fold0_importance.npz
│   └── augmented_csl_fold0_dnarelax_s1_importance.npz
├── ets1/
│   ├── baseline_ets1_fold0_importance.npz
│   ├── augmented_ets1_fold0_importance.npz
│   └── augmented_ets1_fold0_dnarelax_s1_importance.npz
├── ... (10 more pilots)
└── all_pilots_importance_summary.csv  <-- Master summary table
```

Each `.npz` contains:
- `occlusion_mae` — Per-complex per-atom importance scores (shape: n_complexes × n_contact_atoms)
- `prot_atom_indices` — Which protein atoms were scored
- `edge_index` — Protein→DNA contact graph
- `v_prot`, `x_prot` — Coordinates and features (for post-hoc analysis)

## Interpreting Results

**all_pilots_importance_summary.csv** columns:
- `pilot` — Transcription factor (csl, ets1, etc.)
- `arm` — Augmentation type (augmented_frozen, augmented_relax_s1, etc.)
- `n_complexes` — Number of benchmark complexes
- `mean_shift` — Mean importance change vs. baseline (float, can be negative)
- `p_value` — Paired t-test p-value (significance of shift)
- `effect_size` — Cohen's d (magnitude of shift relative to noise)

**Interpretation Guide:**
- **mean_shift > 0, p < 0.05**: Augmentation increases importance on contact atoms (moves toward biochemistry) ✓
- **mean_shift < 0, p < 0.05**: Augmentation *decreases* importance on contact atoms (moves away from biochemistry) ✗
- **p > 0.05**: No significant effect (noise)

## Earlier ETS1-Only Analysis

For comparison with the initial ETS1 single-pilot run:
- Results: `analysis/importance_stats.json`, `analysis/importance_comparison_table.csv`
- Figures: `analysis/importance_*.png`
- Interpretation: `analysis/IMPORTANCE_ANALYSIS_RESULTS.md`

## Hardware Requirements

- **GPU**: RTX5000 or better (16 GB+ VRAM)
- **Memory**: 60 GB RAM per pilot
- **Time**: ~45 min per pilot on GPU (varies by complex count)
- **Total**: ~12-24 hours for all 12 pilots in parallel

## Troubleshooting

**"No .npz in stage4_npz"**: Featurization failed for that pilot. Check stage4 logs.

**"Config not found"**: Training didn't complete. Check stage6_train logs.

**"CUDA out of memory"**: Reduce batch size in interpret_tfconf_all.py (currently 1, already minimal).

**Job hung on occlusion loop**: Check slurmctld logs; may be grid bandwidth bottleneck. Sequential submission usually works.

---

Generated: 2026-07-30
Updated for all-pilots generalization
