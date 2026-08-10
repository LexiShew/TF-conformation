#!/usr/bin/env python3
import numpy as np, pandas as pd, json
from scipy.stats import ttest_rel
from pathlib import Path

R = Path('/project2/rohs_102/shewchuk/TF-conformation')
results_dir = R / 'output' / 'interpret_results'

baseline = np.load(results_dir / 'baseline_ets1_fold0_importance.npz', allow_pickle=True)
augmented = np.load(results_dir / 'augmented_ets1_fold0_importance.npz', allow_pickle=True)
augmented_relax = np.load(results_dir / 'augmented_ets1_fold0_dnarelax_s1_importance.npz', allow_pickle=True)

baseline_mae = baseline['occlusion_mae']
augmented_mae = augmented['occlusion_mae']
augmented_relax_mae = augmented_relax['occlusion_mae']

prot_atom_indices = baseline['prot_atom_indices']
ei = baseline['edge_index']  # protein->DNA contact edges

# All atoms in prot_atom_indices are by definition contact atoms (they have edges to DNA)
is_contact = np.ones(len(prot_atom_indices), dtype=bool)

# Compute shifts
frozen_shift = augmented_mae - baseline_mae
relax_shift = augmented_relax_mae - baseline_mae

# Mean importance per complex
baseline_imp = baseline_mae.mean(axis=1)  # mean over atoms
augmented_imp = augmented_mae.mean(axis=1)
augmented_relax_imp = augmented_relax_mae.mean(axis=1)

# Shifts in mean importance
frozen_effect = augmented_imp - baseline_imp
relax_effect = augmented_relax_imp - baseline_imp

t_frozen, p_frozen = ttest_rel(augmented_imp, baseline_imp)
t_relax, p_relax = ttest_rel(augmented_relax_imp, baseline_imp)

# Top contact atoms (highest baseline importance)
top_atoms = np.argsort(baseline_mae.mean(axis=0))[-5:]
print("Top 5 important contact atoms (by baseline):", prot_atom_indices[top_atoms])

results = {
    'n_complexes': int(len(baseline_mae)),
    'n_contact_atoms': int(len(prot_atom_indices)),
    'baseline_mean_imp': float(baseline_imp.mean()),
    'augmented_mean_imp': float(augmented_imp.mean()),
    'augmented_relax_mean_imp': float(augmented_relax_imp.mean()),
    'frozen_mean_shift': float(frozen_effect.mean()),
    'frozen_std_shift': float(frozen_effect.std()),
    'frozen_t': float(t_frozen),
    'frozen_p': float(p_frozen),
    'frozen_effect_size': float(frozen_effect.mean() / frozen_effect.std() + 1e-8) if frozen_effect.std() > 0 else 0,
    'relax_mean_shift': float(relax_effect.mean()),
    'relax_std_shift': float(relax_effect.std()),
    'relax_t': float(t_relax),
    'relax_p': float(p_relax),
    'relax_effect_size': float(relax_effect.mean() / relax_effect.std() + 1e-8) if relax_effect.std() > 0 else 0,
}

print(json.dumps(results, indent=2))

# Save comparison table
comparison_table = pd.DataFrame({
    'complex_idx': range(len(baseline_mae)),
    'baseline_imp': baseline_imp,
    'augmented_frozen_imp': augmented_imp,
    'augmented_relax_imp': augmented_relax_imp,
    'frozen_shift': frozen_effect,
    'relax_shift': relax_effect,
})

comparison_table.to_csv(results_dir / 'importance_comparison_table.csv', index=False)
print(f"\nSaved: {results_dir / 'importance_comparison_table.csv'}")

# Save results
with open(results_dir / 'importance_stats.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"Saved: {results_dir / 'importance_stats.json'}")
