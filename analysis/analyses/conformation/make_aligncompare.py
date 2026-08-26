#!/usr/bin/env python3
"""
make_aligncompare.py — interface vs global (all-Cα) alignment: DNA misplacement,
per pilot. Justifies the pipeline's default --align-mode interface.

WHY THIS METRIC (not fnat): the pipeline scores fnat with --use_model_dna
(protein vs the DNA carried in the same docked file). That DNA moves rigidly with
the Kabsch transform, so post-hoc fnat is INVARIANT to alignment mode — it cannot
distinguish interface from global. The alignment-sensitive quantity is where the
carried DNA LANDS. Because both modes keep the BioEmu protein in the same frame
(only the carried reference DNA is transformed), the per-state DNA-P RMSD between
the two docked outputs IS the alignment-induced misplacement — no crystal
superposition needed, so it is robust to chain-mapping ambiguity.

Empirically the modes differ only for flexible folds: interface-Cα and all-Cα
Kabsch coincide for a compact domain (ETS1, 0.00 Å) but diverge sharply for a
mobile one (LEF1 HMG box, ~3 Å median, up to 25 Å) — global alignment averages
placement error across the whole chain, swinging the DNA off the binding site.
That is the justification for interface alignment as the safe default.

Data: analysis/analyses/align_compare/dna_displacement_interface_vs_global.csv
      (cols: pilot, state, dna_disp) — built by compute_align_displacement.py,
      which docks each Stage-1 frame both ways (stage2_redock.py --align-mode
      interface|all) and measures the per-state DNA-P displacement between modes.
Pilots on x (rigid→flexible), displacement on y, per-state strip + median tick.
Env: deeppbs.  Usage: python make_aligncompare.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "common"))
from fig_common import DATA_DIR, TEAL, GREY, ALARM, savefig
import pandas as pd

CSV = os.path.join(os.path.dirname(DATA_DIR), "analyses", "align_compare",
                   "dna_displacement_interface_vs_global.csv")
df = pd.read_csv(CSV)

# rigid → flexible ordering, with a fold descriptor for the tick label
FOLD = {"ets1": "ETS1\n(compact)", "tbp": "TBP\n(β-saddle)", "lef1": "LEF1\n(HMG box)"}
order = [p for p in ["ets1", "tbp", "lef1"] if p in df.pilot.unique()]
order += [p for p in sorted(df.pilot.unique()) if p not in order]
x = np.arange(len(order))

fig, ax = plt.subplots(figsize=(2.1 * len(order) + 0.5, 4.4))
rng = np.random.default_rng(0)
for xi, p in zip(x, order):
    d = df[df.pilot == p].dna_disp.values
    med = float(np.median(d))
    col = GREY if med < 1.0 else ALARM          # near-zero = benign, large = alarm
    ax.scatter(xi + rng.uniform(-0.14, 0.14, size=len(d)), d, s=14, color=col,
               alpha=0.45, edgecolor="none", zorder=3)
    ax.plot([xi - 0.22, xi + 0.22], [med, med], color=TEAL, lw=2.6, zorder=5,
            solid_capstyle="round")
    ax.annotate(f"med {med:.2f} Å", (xi + 0.26, med), fontsize=8, va="center",
                color=TEAL, fontweight="bold")

ax.axhline(0, color="0.6", lw=0.8, zorder=1)
ax.set_xticks(x)
ax.set_xticklabels([FOLD.get(p, p) for p in order], fontsize=10)
ax.set_ylabel("DNA misplacement under global alignment\n"
              "(per-state RMSD vs interface alignment, Å)")
ax.set_xlabel("pilot  ·  rigid → flexible fold")
ax.set_title("Global alignment misplaces the DNA where the fold is flexible")
ax.margins(x=0.10)
ax.set_ylim(-1.0, df.dna_disp.max() * 1.08)
fig.text(0.01, -0.04,
         "Each Stage-1 conformer docked both ways (interface-Cα vs all-Cα Kabsch). Both share the protein "
         "frame, so the per-state DNA displacement between modes is the alignment-induced error. For a compact "
         "fold (ETS1) the two are identical; for a flexible one (LEF1, HMG box) global alignment swings the "
         "carried DNA up to 25 Å off — which is why interface alignment is the pipeline default.",
         fontsize=6.2, ha="left", va="top", wrap=True)
savefig(fig, "F5_align_dna_displacement.png")
plt.close(fig)
print("wrote F5_align_dna_displacement.png")
