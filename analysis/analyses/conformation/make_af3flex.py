#!/usr/bin/env python3
"""
make_af3flex.py — AF3-vs-ensemble DNA-flexibility figure (pilots on x).

Rebuilds analysis/dna_relax/figures/af3_vs_ensemble_mgwfl.png from
analysis/dna_relax/data/af3_vs_ensemble_mgwfl.csv (cols: pilot, af3, frozen,
relaxed, ens_mean, ...). Pilots on the x-axis (ordered rigid→flexible by
ens_mean), mean minor-groove-width fluctuation on y; AF3 (diamond) sits near the
floor while the physical frozen/relaxed ensembles spread wide — the DNA-side
version of the "why BioEmu, not AF3" motivation (Ch. 3), reproducing the paper's
Fig 3 flexibility claim on our pilots.

Was an uncommitted ad-hoc figure; this is the committed, palette-threaded,
pilots-on-x version. Env: deeppbs.  Usage: python make_af3flex.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "common"))
from fig_common import DATA_DIR, TEAL, GREEN, AF3, short_of, savefig
try:
    from fig_common import apply_style; apply_style()
except Exception:
    pass
import pandas as pd

# dna_relax data dir is a sibling of the standard analysis/data
DR_DATA = os.path.join(os.path.dirname(DATA_DIR), "analyses", "dna_relax", "data")
a = pd.read_csv(os.path.join(DR_DATA, "af3_vs_ensemble_mgwfl.csv")).sort_values("ens_mean")
x = np.arange(len(a))

fig, ax = plt.subplots(figsize=(8.6, 4.3))
for i, (_, row) in enumerate(a.iterrows()):
    lo, hi = sorted([row.frozen, row.relaxed])
    ax.plot([i, i], [lo, hi], color="0.75", lw=3, alpha=0.6, zorder=1,
            solid_capstyle="round")
ax.scatter(x, a.af3,     marker="D", s=46, color=AF3,   zorder=3, edgecolor="white",
           linewidth=0.5, label="AF3 (10 samples)")
ax.scatter(x, a.frozen,  marker="o", s=40, color=TEAL,  zorder=3, edgecolor="white",
           linewidth=0.5, label="frozen ensemble")
ax.scatter(x, a.relaxed, marker="o", s=40, color=GREEN, zorder=3, edgecolor="white",
           linewidth=0.5, label="relaxed ensemble")
ax.set_xticks(x)
ax.set_xticklabels([short_of(t) for t in a.pilot], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("mean MGW fluctuation (Å)")
ax.set_xlabel("pilot (ordered by ensemble flexibility ➜)")
ax.set_title("AF3 underestimates DNA conformational flexibility vs physical ensembles")
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.margins(x=0.02)
fig.text(0.01, -0.04,
         "12 pilots. AF3 (10 samples) collapses minor-groove-width fluctuation toward zero; the "
         "physical frozen/relaxed ensembles spread far wider. Reproduces the paper's Fig 3 flexibility claim.",
         fontsize=6, ha="left", va="top", wrap=True)
savefig(fig, "af3_vs_ensemble_mgwfl.png", subdir=os.path.join("..", "dna_relax", "figures"))
plt.close(fig)
print("wrote af3_vs_ensemble_mgwfl.png")
