#!/usr/bin/env python3
"""D1b - AF3 vs BioEmu vs BioEmu2 conformational diversity.

The three-arm version of plot_diversity.py (D1). Same quantity, same axis, same
ordering, one extra arm: per-pilot within-ensemble pairwise Calpha-RMSD, with
BioEmu2's complex-sampled ensembles alongside the original BioEmu-1 docked frames
and AF3's diffusion samples.

Reads
    analysis/data/ensemble_diversity_pairwise.csv            af3 + bioemu (canonical)
    analyses/diversity/data/bioemu2_diversity_pairwise.csv   bioemu2 (this theme)
    plus the two matching summary tables, for medians and n.

Writes
    analyses/diversity/figures/D1b_diversity_bioemu2.png

Two things the figure states explicitly rather than hiding
---------------------------------------------------------
1. n differs by an order of magnitude between arms (AF3 10 models, BioEmu ~90-100
   frames, BioEmu2 100 frames). A violin's width is a density, so unequal n does
   not bias the median, but the tails of a 10-sample ensemble are not comparable
   to those of a 100-sample one. n is printed under every arm.
2. Pilots whose BioEmu2 ensemble has not finished are drawn with the arm absent
   and marked, instead of being dropped from the axis - so a missing arm reads as
   missing rather than as zero diversity.

Run in any env with pandas + matplotlib (e.g. deeppbs); cap BLAS threads on a
login node.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
DATA_CANON = os.path.join(BASE, "analysis", "data")
DATA_THEME = os.path.join(HERE, "data")
FIGS = os.path.join(HERE, "figures")

sys.path.insert(0, BASE)
from palette import TEAL as C_BIO, AF3 as C_AF3                   # noqa: E402

# BioEmu2 needs a third cool hue distinct from teal (BioEmu) and lavender (AF3).
# Kept local rather than edited into palette.py unilaterally; promote it there if
# the three-arm comparison becomes a standing figure.
C_BE2 = "#C98BB8"   # soft rose - BioEmu2 (complex sampling)

# fnat pass-rate order, same as plot_diversity.py
ORDER = ["ets1", "tbp", "egr1", "engrailed", "foxa", "lef1",
         "csl", "err", "nfat", "runx", "hsf"]
LABELS = {"ets1": "ETS1", "tbp": "TBP", "egr1": "EGR1", "engrailed": "engrailed",
          "foxa": "FOXA", "lef1": "LEF1", "csl": "CSL", "err": "ERR",
          "nfat": "NFAT", "runx": "RUNX", "hsf": "HSF"}

ARMS = [("af3", C_AF3, "AlphaFold3"),
        ("bioemu", C_BIO, "BioEmu 1 (docked chain)"),
        ("bioemu2", C_BE2, "BioEmu2 (complex)")]


def load():
    pw = pd.read_csv(os.path.join(DATA_CANON, "ensemble_diversity_pairwise.csv"))
    summ = pd.read_csv(os.path.join(DATA_CANON, "ensemble_diversity.csv"))
    p2 = os.path.join(DATA_THEME, "bioemu2_diversity_pairwise.csv")
    s2 = os.path.join(DATA_THEME, "bioemu2_diversity.csv")
    if os.path.exists(p2) and os.path.getsize(p2) > 0:
        pw = pd.concat([pw, pd.read_csv(p2)], ignore_index=True)
    if os.path.exists(s2) and os.path.getsize(s2) > 0:
        summ = pd.concat([summ, pd.read_csv(s2)], ignore_index=True)
    return pw, summ


def main():
    pw, summ = load()
    med = {(r.pilot, r.source): r.median_pairwise for r in summ.itertuples()}
    nmap = {(r.pilot, r.source): int(r.n) for r in summ.itertuples()}

    fig, ax = plt.subplots(figsize=(11.5, 6.8))
    xpos = np.arange(len(ORDER))
    offs = {"af3": -0.26, "bioemu": 0.0, "bioemu2": +0.26}
    half = 0.085
    YMAX = 13.0

    missing = []
    for i, pilot in enumerate(ORDER):
        for src, col, _lab in ARMS:
            vals = pw[(pw.pilot == pilot) & (pw.source == src)]["pairwise_ca_rmsd"].values
            x = xpos[i] + offs[src]
            if len(vals) == 0:
                if src == "bioemu2":
                    missing.append(pilot)
                    ax.text(x, 0.22, "no\nensemble", ha="center", va="bottom",
                            fontsize=5.6, color="#8A8A8A", style="italic")
                continue
            parts = ax.violinplot([np.clip(vals, None, YMAX)], positions=[x],
                                  widths=0.24, showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(col)
                b.set_edgecolor("black")
                b.set_linewidth(0.5)
                b.set_alpha(0.5)
            m = med.get((pilot, src), float(np.median(vals)))
            ax.plot([x - half, x + half], [m, m], color="black", lw=1.3, zorder=5)
            ax.text(x, m + 0.12, f"{m:.2f}", ha="center", va="bottom",
                    fontsize=6.2, color=col, zorder=6, fontweight="bold")
            n = nmap.get((pilot, src))
            if n:
                ax.text(x, -0.62, f"{n}", ha="center", va="top",
                        fontsize=5.4, color=col)

    ax.set_xticks(xpos)
    ax.set_xticklabels([LABELS[p] for p in ORDER], fontsize=9)
    ax.set_ylabel(r"within-ensemble pairwise C$\alpha$-RMSD ($\AA$)", fontsize=10)
    ax.set_ylim(-1.5, YMAX)
    ax.set_xlim(-0.6, len(ORDER) - 0.4)
    ax.axhline(0, color="#CCCCCC", lw=0.6, zorder=0)
    ax.text(-0.55, -0.62, "n =", ha="right", va="top", fontsize=5.6, color="#555555")
    ax.set_title("Conformational diversity of protein-DNA ensembles: "
                 "AlphaFold3 vs BioEmu 1 vs BioEmu2", fontsize=11.5, pad=12)
    ax.grid(axis="y", color="#E8E8E8", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    handles = [mpatches.Patch(facecolor=c, edgecolor="black", lw=0.5,
                              alpha=0.5, label=lab) for _s, c, lab in ARMS]
    ax.legend(handles=handles, loc="upper right", fontsize=8.5, frameon=False)

    note = ("Pilots ordered by fnat pass-rate (repo convention). Median marked; "
            r"violins clipped at %.0f $\AA$. BioEmu2 frames restricted to one "
            "protein chain (the deposited assembly can hold two copies)." % YMAX)
    if missing:
        note += "  BioEmu2 pending: " + ", ".join(LABELS[p] for p in missing) + "."
    fig.text(0.012, 0.012, note, fontsize=7.0, color="#444444")

    fig.tight_layout(rect=[0, 0.035, 1, 1])
    os.makedirs(FIGS, exist_ok=True)
    out = os.path.join(FIGS, "D1b_diversity_bioemu2.png")
    fig.savefig(out, dpi=200)
    print("wrote", out)

    cols = ["pilot"] + [f"{s}_median" for s, _c, _l in ARMS] + \
           [f"{s}_n" for s, _c, _l in ARMS]
    rows = []
    for p in ORDER:
        row = {"pilot": p}
        for s, _c, _l in ARMS:
            row[f"{s}_median"] = med.get((p, s), "")
            row[f"{s}_n"] = nmap.get((p, s), "")
        rows.append(row)
    tbl = pd.DataFrame(rows)[cols]
    tpath = os.path.join(DATA_THEME, "diversity_three_arm_summary.csv")
    tbl.to_csv(tpath, index=False)
    print("wrote", tpath)
    print(tbl.to_string(index=False))


if __name__ == "__main__":
    main()
