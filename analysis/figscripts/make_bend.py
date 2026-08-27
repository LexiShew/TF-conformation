#!/usr/bin/env python3
"""
make_bend.py — DNA-bend change under relaxation, per pilot (pilots on x).

Rebuilds analysis/dna_relax/figures/crossfamily_bend.png from
analysis/dna_relax/data/bend_delta_bootstrap.csv (cols: tf, family, delta, lo,
hi, abs_delta). Pilots on the x-axis (ordered by Δbend), Δ overall axis bend
(relaxed − frozen, °) on y, bootstrap 95% CI error bars, bars coloured by
induced-fit vs rigid family. Most CIs cross zero — the honest reading is that
bend change under relaxation is small and mostly not significant.

Partial coverage: only the 7 pilots with a bootstrap bend estimate (the subset
with DNA-shape analysis run). Was an uncommitted ad-hoc figure; this is the
committed, palette-threaded, pilots-on-x version. Env: deeppbs.
Usage: python make_bend.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from fig_common import DATA_DIR, TEAL, GREY, ALARM, short_of, savefig
try:
    from fig_common import apply_style; apply_style()
except Exception:
    pass
import pandas as pd

DR_DATA = os.path.join(os.path.dirname(DATA_DIR), "dna_relax", "data")
b = pd.read_csv(os.path.join(DR_DATA, "bend_delta_bootstrap.csv")).sort_values("delta")
b = b.reset_index(drop=True)
x = np.arange(len(b))
colmap = {"induced-fit": TEAL, "rigid": GREY}
cols = [colmap.get(f, GREY) for f in b.family]

fig, ax = plt.subplots(figsize=(7.2, 4.2))
yerr = np.array([b.delta - b.lo, b.hi - b.delta])
ax.bar(x, b.delta, color=cols, width=0.66, edgecolor="0.3", linewidth=0.5, zorder=2)
ax.errorbar(x, b.delta, yerr=yerr, fmt="none", ecolor="0.35", lw=1.0, capsize=3, zorder=3)
ax.axhline(0, color=ALARM, lw=1.0, ls="--", zorder=1)
ax.set_xticks(x)
ax.set_xticklabels([short_of(t) for t in b.tf], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Δ overall axis bend, relaxed − frozen (°)")
ax.set_xlabel("pilot")
ax.set_title("DNA bend change under relaxation, per pilot")
ax.legend(handles=[Patch(facecolor=TEAL, label="induced-fit family"),
                   Patch(facecolor=GREY, label="rigid family")],
          frameon=False, fontsize=8, loc="upper left")
fig.text(0.01, -0.05,
         "7 pilots with a bootstrap bend estimate (subset with DNA-shape analysis run). Error bars = "
         "bootstrap 95% CI; most CIs cross zero — bend change under relaxation is small and mostly not significant.",
         fontsize=6, ha="left", va="top", wrap=True)
savefig(fig, "crossfamily_bend.png", subdir=os.path.join("..", "dna_relax", "figures"))
plt.close(fig)
print("wrote crossfamily_bend.png")
