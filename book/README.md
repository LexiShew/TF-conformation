# book/ — long-form write-up of the TF-conformation project

## Contents
| File | What it is |
|---|---|
| `OUTLINE.md` | Full chapter-by-chapter outline (6 parts, 21 chapters, 7 appendices). Each chapter lists its argument, its subpoints with the concrete numbers, the figures that ALREADY EXIST (with repo paths), and the figures that still need to be created (IDs `N<chapter>.<n>`). |

## How the outline was built
A read-only survey of the whole repo on endeavour:
- all 20 project markdown docs (READMEs, `RESULTS_INVENTORY.md`, `CONDITIONS_INVENTORY.md`,
  `mechanism_and_roadmap.md`, `INTERPRET_SUITE_README.md`, `REVIEW_figure_scripts.md`,
  `docs/concerns.txt`, `analysis/dna_relax/README.md`, `stiffness_prior/README.md`)
- the figure inventories: `analysis/figures/` (F/I/R/S/P/M/D, 20 PNGs),
  `analysis/figure_scripts/` (fig1-9), `analysis/dna_relax/figures/` (24 PNGs),
  `rmsd_analysis/plots/`
- the result tables in `analysis/data/` and `analysis/dna_relax/data/`
- `deck/deck_spec.json` (29 slide titles - the existing narrative skeleton)
- `git log` 2026-04-24 -> 2026-07-30 (project chronology)

Numbers quoted in the outline were re-derived from the CSVs on disk where possible
(`perentry_accuracy.csv`, `mechanism_apo_holo.csv`, `perstate_metrics.csv`,
`ensemble_diversity.csv`) rather than copied from prose, so a few disagree with older
markdown summaries. Those disagreements are called out in Ch. 12.2 rather than silently
reconciled.

## Reproducing the numbers quoted in the outline
```bash
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /project2/rohs_102/shewchuk/TF-conformation
python - <<'EOF'
import pandas as pd
d = pd.read_csv('analysis/data/perentry_accuracy.csv')
d['dP'] = d.aug_pearsonr - d.base_pearsonr
print(d.groupby('pilot').dP.mean().round(4).to_string())
EOF
```

## Status
v0.1 draft outline. Open questions for the author are tracked at the top of `OUTLINE.md`
section 0 and in the "Highest-value new work" ranking at the end.
