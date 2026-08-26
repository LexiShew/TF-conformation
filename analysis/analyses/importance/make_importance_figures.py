#!/usr/bin/env python3
import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from pathlib import Path

R = Path("/project2/rohs_102/shewchuk/TF-conformation")
results_dir = R / 'output' / 'interpret_results'

# Load data
baseline = np.load(results_dir / 'baseline_ets1_fold0_importance.npz', allow_pickle=True)
augmented = np.load(results_dir / 'augmented_ets1_fold0_importance.npz', allow_pickle=True)
augmented_relax = np.load(results_dir / 'augmented_ets1_fold0_dnarelax_s1_importance.npz', allow_pickle=True)

baseline_mae = baseline['occlusion_mae']
augmented_mae = augmented['occlusion_mae']
augmented_relax_mae = augmented_relax['occlusion_mae']

# Mean importance per complex
baseline_imp = baseline_mae.mean(axis=1)
augmented_imp = augmented_mae.mean(axis=1)
augmented_relax_imp = augmented_relax_mae.mean(axis=1)

# Shifts
frozen_shift = augmented_imp - baseline_imp
relax_shift = augmented_relax_imp - baseline_imp

# Figure 1: Distributions
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
ax.hist(frozen_shift, bins=20, alpha=0.7, color='steelblue', edgecolor='black')
ax.axvline(frozen_shift.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {frozen_shift.mean():.4f}')
ax.set_xlabel('Importance Shift (Augmented - Baseline)')
ax.set_ylabel('Frequency')
ax.set_title('Augmented-Frozen vs Baseline (99 complexes)')
ax.legend()
ax.grid(alpha=0.3)

ax = axes[1]
ax.hist(relax_shift, bins=20, alpha=0.7, color='coral', edgecolor='black')
ax.axvline(relax_shift.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {relax_shift.mean():.4f}')
ax.set_xlabel('Importance Shift (Augmented - Baseline)')
ax.set_ylabel('Frequency')
ax.set_title('Augmented-Relaxed vs Baseline (99 complexes)')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(results_dir / 'importance_shift_distributions.png', dpi=300, bbox_inches='tight')
print("✓ importance_shift_distributions.png")

# Figure 2: Boxplot
fig, ax = plt.subplots(figsize=(10, 6))
data = [baseline_imp, augmented_imp, augmented_relax_imp]
labels = ['Baseline', 'Augmented\n(Frozen DNA)', 'Augmented\n(Relaxed DNA)']
colors = ['lightgray', 'steelblue', 'coral']
bp = ax.boxplot(data, labels=labels, patch_artist=True, widths=0.6)
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
ax.set_ylabel('Mean Occlusion MAE per Complex')
ax.set_title('DeepPBS Atom Importance: ETS1 Benchmark Set')
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(results_dir / 'importance_boxplot_comparison.png', dpi=300, bbox_inches='tight')
print("✓ importance_boxplot_comparison.png")

# Figure 3: Scatter
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
lim_min = min(baseline_imp.min(), augmented_imp.min(), augmented_relax_imp.min())
lim_max = max(baseline_imp.max(), augmented_imp.max(), augmented_relax_imp.max())

ax = axes[0]
ax.scatter(baseline_imp, augmented_imp, alpha=0.6, s=50, color='steelblue')
ax.plot([lim_min, lim_max], [lim_min, lim_max], 'k--', alpha=0.3, label='No change')
ax.set_xlabel('Baseline Importance')
ax.set_ylabel('Augmented (Frozen) Importance')
ax.set_title('Frozen DNA Augmentation Effect')
ax.legend()
ax.grid(alpha=0.3)

ax = axes[1]
ax.scatter(baseline_imp, augmented_relax_imp, alpha=0.6, s=50, color='coral')
ax.plot([lim_min, lim_max], [lim_min, lim_max], 'k--', alpha=0.3, label='No change')
ax.set_xlabel('Baseline Importance')
ax.set_ylabel('Augmented (Relaxed) Importance')
ax.set_title('Relaxed DNA Augmentation Effect')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(results_dir / 'importance_scatter_comparison.png', dpi=300, bbox_inches='tight')
print("✓ importance_scatter_comparison.png")

print("\nAll figures saved")
