#!/usr/bin/env python3
"""
compare_importance_all.py — Compare importance shifts across all pilots/arms.
Generates per-pilot statistics and combined multi-pilot summary.
"""
import numpy as np, pandas as pd, json
from scipy.stats import ttest_rel
from pathlib import Path

R = Path('/project2/rohs_102/shewchuk/TF-conformation')
results_root = R / 'output' / 'interpret_results_all'

pilots = sorted([d.name for d in results_root.iterdir() if d.is_dir()])
print(f"Analyzing {len(pilots)} pilots: {', '.join(pilots)}")

summary_rows = []

for pilot in pilots:
    pilot_dir = results_root / pilot
    npz_files = sorted(list(pilot_dir.glob('*_importance.npz')))
    
    if len(npz_files) < 2:
        print(f"  {pilot}: Skipping (< 2 arms)")
        continue
    
    print(f"\n{pilot}:")
    
    # Load all arms
    data = {}
    for npz in npz_files:
        arm_name = npz.stem.replace('_importance', '')
        data[arm_name] = np.load(npz, allow_pickle=True)
    
    # Identify arm types
    baseline_key = [k for k in data.keys() if k.startswith('baseline')][0]
    augmented_keys = [k for k in data.keys() if 'augmented' in k and 'dnarelax' not in k]
    relax_keys = [k for k in data.keys() if 'dnarelax' in k]
    
    baseline_mae = data[baseline_key]['occlusion_mae']
    baseline_imp = baseline_mae.mean(axis=1)
    
    # Compare each augmented arm
    for aug_key in augmented_keys:
        aug_mae = data[aug_key]['occlusion_mae']
        aug_imp = aug_mae.mean(axis=1)
        shift = aug_imp - baseline_imp
        
        t, p = ttest_rel(aug_imp, baseline_imp)
        
        row = {
            'pilot': pilot,
            'arm': 'augmented_frozen',
            'n_complexes': len(shift),
            'mean_shift': float(shift.mean()),
            'std_shift': float(shift.std()),
            't_stat': float(t),
            'p_value': float(p),
            'effect_size': float(shift.mean() / (shift.std() + 1e-8)) if shift.std() > 0 else 0,
        }
        summary_rows.append(row)
        print(f"  frozen: shift={shift.mean():.6f}, p={p:.4f}")
    
    # Compare each relaxed arm
    for relax_key in relax_keys:
        relax_mae = data[relax_key]['occlusion_mae']
        relax_imp = relax_mae.mean(axis=1)
        shift = relax_imp - baseline_imp
        
        t, p = ttest_rel(relax_imp, baseline_imp)
        seed = relax_key.split('_s')[-1] if '_s' in relax_key else 'unknown'
        
        row = {
            'pilot': pilot,
            'arm': f'augmented_relax_s{seed}',
            'n_complexes': len(shift),
            'mean_shift': float(shift.mean()),
            'std_shift': float(shift.std()),
            't_stat': float(t),
            'p_value': float(p),
            'effect_size': float(shift.mean() / (shift.std() + 1e-8)) if shift.std() > 0 else 0,
        }
        summary_rows.append(row)
        print(f"  relax_s{seed}: shift={shift.mean():.6f}, p={p:.4f}")

# Save summary table
if summary_rows:
    df = pd.DataFrame(summary_rows)
    df.to_csv(results_root / 'all_pilots_importance_summary.csv', index=False)
    print(f"\nSaved: all_pilots_importance_summary.csv")
    print(df[['pilot','arm','mean_shift','p_value','effect_size']].to_string())
