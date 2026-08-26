# scripts/structure_viz — PyMOL structure rendering toolkit

The single home for PyMOL rendering: a shared primitive library plus the general,
pipeline-level render drivers. (Absorbs the former `scripts/pymol/` and `scripts/viz/`.)

## Shared library
| file | what |
|---|---|
| `pymol_lib.py` | Reusable colouring primitives — `hex_to_rgb`, per-state gradient colouring (`gradient_states`), and named palettes (`RAINBOW9`, `BLUE_PURPLE_PINK`). Back-compat command names `gradient_protein_states` / `gradient_protein_split` / `spectrum_states` are preserved. `run` it inside a PyMOL session. |

## Render drivers (general, any pilot)
| file | renders |
|---|---|
| `render_ensemble.sh` → `_render_ensemble.py` | Crystal bound pose vs a docked + minimized BioEmu ensemble for one pilot, on a common DNA frame. |
| `render_state_stages.sh` → `_render_stages.py` | One TF–DNA state across pipeline stages (apo → docked → minimized), 4-panel. |
| `render_all.sbatch` | Batch `render_ensemble.sh` across pilots. |

The `.sh` wrappers resolve their own directory (`SCRIPT_DIR`) and `cd` to the repo root,
so they work regardless of where the repo lives. Renders are written under
`analysis/figures/pymol/`.

## Related (not moved — they live with the analysis that owns them)
Theme-specific renderers stay in their theme dir and may reuse `pymol_lib.py`:
- `analysis/analyses/diversity/{render_ensembles,montage_ensembles}.py` — ensemble stacks + montage for the D1 diversity figure. *(Note: `render_ensemble` here vs `render_ensembles` in diversity are different outputs — a single-pilot overlay vs the diversity stack montage.)*
- `analysis/analyses/conformation/{render_S.sbatch,render_prepost_min.py,view_test*.py}` — S-figure and pre/post-minimization renders.

These still carry their own inline colouring; migrating them to import `pymol_lib.py`
is a follow-up (needs PyMOL on the cluster to verify).
