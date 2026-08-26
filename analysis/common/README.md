# analysis/common — shared foundation for every analysis theme

- `fig_common.py` — repo-root-anchored paths, pilot auto-discovery, `savefig`, family/label maps.
  Import from a theme with:
  ```python
  import os, sys
  sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "common"))
  from fig_common import *
  ```
- The colour palette is `palette.py` at the **repo root** (re-exported through `fig_common`).

New themes (dna_mismatch, idr, ...) live in `analysis/analyses/<theme>/` and import from here.
