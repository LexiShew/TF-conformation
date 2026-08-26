# analysis/ layout & how to add a new analysis

The `analysis/` tree separates **shared foundation**, **canonical data**, **cross-cutting
figures**, and **self-contained analysis themes**. New analyses (e.g. DNA mismatches,
intrinsically disordered regions) are added as drop-in theme dirs — no restructuring.

```
analysis/
├── README.md          # the F/I/R/S/P/M conformation-suite science report
├── common/            # SHARED foundation, imported by every theme
│   └── fig_common.py  #   repo-root-anchored paths, pilot discovery, savefig, family maps
│                      #   (colour palette is palette.py at the REPO ROOT, re-exported here)
├── data/              # canonical cross-cutting tables (per-state/-residue metrics, family map,
│                      #   perseed_*, augmentation/conditions summaries, inventories)
├── figures/           # cross-cutting paper figures (F/I/R/S/P/D) + pymol/ renders
├── docs/              # write-ups, inventories, roadmap, this file
├── analyses/          # ONE self-contained dir per theme
│   ├── conformation/  #   the F/I/R/S/P/M figure suite (make_*.py)
│   ├── importance/    #   atom-importance attribution (scripts + data/ + figures/)
│   ├── mechanism/     #   apo/holo mechanism (code + its own data/ + figures/)
│   ├── dna_relax/     #   DNA relaxation / pyCurves (bulk per-state JSON is gitignored)
│   ├── align_compare/ #   interface-vs-global alignment (data + sbatch)
│   ├── stage3/        #   stage-3 minimization figures
│   └── diversity/     #   ensemble diversity + PyMOL render utilities
└── _deprecated/       # superseded three-arm figure_scripts suite (kept for provenance)
```

## Path conventions (important — scripts are run from the repo root)

- **Shared code:** import the foundation with a one-line bootstrap, then use its paths:
  ```python
  import os, sys
  sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "common"))
  from fig_common import *          # DATA_DIR, FIG_DIR, discover_pilots, savefig, palette colours, ...
  ```
  `fig_common` anchors `TFCONF` to the repo root and exposes `DATA_DIR` (`analysis/data`) and
  `FIG_DIR` (`analysis/figures`) — never hardcode those.
- **A theme's own data/figures:** anchor file-relative so they move with the code:
  `HERE = Path(__file__).resolve().parent; DATA = HERE/"data"; FIGS = HERE/"figures"`.
- **Reading another theme's data:** go through the tree explicitly, e.g.
  `os.path.join(os.path.dirname(DATA_DIR), "analyses", "dna_relax", "data", ...)`.
- **Repo-root anchor by depth:** a script at `analysis/analyses/<theme>/x.py` reaches the repo
  root with `Path(__file__).resolve().parents[3]`; one more level deep (`.../<theme>/scripts/x.py`)
  uses `parents[4]`.

## Adding a new theme (e.g. `dna_mismatch/`, `idr/`)

1. `mkdir -p analysis/analyses/<theme>/{data,figures}` and add a short `README.md`.
2. Put scripts in `analysis/analyses/<theme>/`; import shared code via the bootstrap above.
3. Write theme figures to `HERE/"figures"` (or `savefig(..., subdir="<theme>")` into the
   top-level `figures/`), and derived tables to `HERE/"data"`.
4. Large/regenerable outputs (raw pyCurves JSON, bulk renders) — add a `.gitignore` rule; keep
   only the distilled summary tables. See how `analyses/dna_relax/pycurves/` is ignored.
5. If the theme reuses DNA-shape/pyCurves machinery, factor the shared helpers into
   `common/dna_shape.py` at that point (deferred until a second DNA theme needs it).
