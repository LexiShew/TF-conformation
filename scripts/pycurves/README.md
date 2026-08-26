# scripts/pycurves — pyCurves DNA-shape toolkit

Shared, parameterized tools for running [pyCurves](https://github.com/…) (Curves+/legacy
DNA conformational analysis) on any protein–DNA structure and reducing its output to tables.
Reusable across analyses — currently `analysis/analyses/dna_relax/` and `af3/`, and available
to future DNA-shape work (e.g. `dna_mismatch/`).

pyCurves' JAX/XLA thread pool exceeds the login-node `RLIMIT_NPROC`, so the runners **must**
execute on a compute node. CPU-only JAX (the `pycurves` conda env) is sufficient.

| file | what | interface |
|---|---|---|
| `run_pycurves.sh` | Run pyCurves on one PDB (both `curves` + `json` formats). | `bash run_pycurves.sh <input.pdb> <out_prefix>` |
| `pycurves_array.sh` | SLURM array batch: runs `legacy` + `curvesplus` conventions over a worklist. | `sbatch pycurves_array.sh <worklist.tsv>` (tab cols: `TF COND STATE INPUT_PDB OUT_PREFIX`; `CHUNK` env = rows/task, default 20) |
| `aggregate_pycurves.py` | Reduce a tree of `*_legacy.json`/`*_curvesplus.json` to per-structure + per-(tf,cond) summary tables. | `python aggregate_pycurves.py <pycurves_root> <out_prefix>` |
| `parse_pycurves_grooves.py` | Parse minor/major groove widths + overall bend from pyCurves `.txt`. | `python parse_pycurves_grooves.py <txt_dir>` |

## Callers provide their own worklist + outputs
The tools are worklist-driven and write wherever the worklist / args point — they hold **no
analysis-specific paths**. Each consumer builds its own worklist and keeps its own outputs:

- **dna_relax** — `analysis/analyses/dna_relax/pycurves_batch/build_worklist.sh` writes
  `pycurves_worklist.tsv` (crystal / frozen / relaxed ensembles); then
  `sbatch scripts/pycurves/pycurves_array.sh analysis/analyses/dna_relax/pycurves_worklist.tsv`.
  Its dna_relax-specific reduction (`aggregate_perposition.py`) and representative-viz decks
  (`pycurves_viz_reps.sh`) stay under `pycurves_batch/`.
- **af3** — `af3/build_af3_worklist.sh` writes `af3/af3_pycurves_worklist.tsv` from the AF3
  prediction dirs, then `CONVENTIONS=legacy sbatch --array=… scripts/pycurves/pycurves_array.sh
  af3/af3_pycurves_worklist.tsv` (see the builder's header for the exact submit line).
  `CONVENTIONS=legacy` preserves af3's legacy-only output; the outputs land next to each PDB,
  unchanged. (This replaced the old standalone `af3_pycurves.sh`.)
