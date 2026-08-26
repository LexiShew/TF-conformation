#!/usr/bin/env python3
"""
make_M.py — mechanism (apo/holo reachability) climax figure.

M1_apo_holo_mechanism.png : two panels.
  A (main)  : protein-reachability (x = d_min, the minimum RMSD the apo ensemble
              achieves toward the bound pose — LOWER = better reach) vs
              conformational spread (y = spread), for the pilots whose
              coordinate pass has been run (reach_status == 'computed').
              Marker colour = sign of aug_dP (TEAL = augmentation helped DeepPBS,
              ALARM = hurt, GREY = no benchmark / aug_dP NaN); marker area
              scaled by |aug_dP|; each point annotated with its short label and
              literature DNA-deformation call. TBP is flagged explicitly — it
              has the BEST reachability (lowest d_min ≈ 0.59 Å) yet augmentation
              hurt it most (most negative ΔPCC), the central paradox of the
              mechanism story.
  B (companion): conformational spread for ALL pilots, so the reach-pending
              pilots are still shown; bars coloured by reach_status
              (computed vs 'reachability pending').

Honest about n: panel A carries only the pilots with a computed reachability
coordinate (6 today); the rest live in panel B labelled pending.

Why d_min (not reach_ratio) on the x-axis: d_min is the direct reachability
observable — the closest the free (apo) ensemble gets to the bound pose — and is
the axis the original mechanism memo used ("reachability d_min (Å)"). reach_ratio
is a derived normalisation of it and is retained in the CSV for reference. The
choice is the established one, independent of any single pilot's outcome.

DATA  analysis/data/mechanism_apo_holo.csv
      cols: pilot,family,spread,d_min,reach_ratio,rmsf_mean,aug_dP,dna_deform,
            reach_status

Usage:  cd analysis/figscripts && python make_M.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "common"))
from fig_common import (
    DATA_DIR, short_of, dna_deform_of, savefig, apply_style,
    TEAL, ALARM, GREY,
)

MECH_CSV = os.path.join(DATA_DIR, "mechanism_apo_holo.csv")


def _sizes(vals, lo=60, hi=340):
    """Scale |aug_dP| to marker area. NaN (no benchmark) -> smallest size so
    the point is still drawn (just visibly the least-weighted)."""
    a = np.abs(np.asarray(vals, float))
    finite = a[np.isfinite(a)]
    if finite.size == 0 or finite.max() == finite.min():
        out = np.full_like(a, (lo + hi) / 2)
    else:
        out = lo + (hi - lo) * (a - finite.min()) / (finite.max() - finite.min())
    out[~np.isfinite(out)] = lo   # NaN aug_dP -> smallest marker
    return out


def fig_M1(mech):
    computed = mech[mech["reach_status"] == "computed"].copy()
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(12.4, 5.2), gridspec_kw=dict(width_ratios=[1.35, 1.0]))

    # ---- Panel A: reachability (d_min) vs spread (computed pilots only) ----
    def col_of(d):
        if not np.isfinite(d):
            return GREY          # no benchmark (aug_dP NaN)
        return TEAL if d >= 0 else ALARM
    # Off-scale handling: a far-unreachable pilot (e.g. the dux4 dimer at
    # d_min≈10 Å) would compress the informative 0.5–4.5 Å range. Cap the axis
    # and draw such points at the cap with an arrow + value, so the resolved
    # cluster stays legible without hiding the outlier.
    XCAP = 5.0
    computed = computed.copy()
    computed["x_plot"] = computed["d_min"].clip(upper=XCAP)
    computed["offscale"] = computed["d_min"] > XCAP
    xr = computed["x_plot"].to_numpy()
    yr = computed["spread"].to_numpy()
    dP = computed["aug_dP"].to_numpy()
    cols = [col_of(d) for d in dP]
    sizes = _sizes(dP)
    axA.scatter(xr, yr, s=sizes, c=cols, alpha=0.85, edgecolor="0.25",
                linewidth=0.8, zorder=3)

    # Deterministic label de-collision: sort by (x,y), alternate the vertical
    # offset direction, and drop a thin leader line so a nudged label still
    # points at its marker (§6.9). Bold TBP; annotate off-scale pilots specially.
    rows_sorted = computed.sort_values(["x_plot", "spread"]).reset_index(drop=True)
    for i, row in rows_sorted.iterrows():
        tf = row["pilot"]
        x, y = row["x_plot"], row["spread"]
        if row["offscale"]:
            axA.annotate(f"{short_of(tf)}  (d$_{{min}}$≈{row['d_min']:.0f} Å, off-scale →)",
                         (x, y), textcoords="offset points", xytext=(-6, 8),
                         fontsize=7.5, ha="right", va="bottom", color="0.35")
            continue
        dy = 12 if (i % 2 == 0) else -14
        va = "bottom" if dy > 0 else "top"
        txt = f"{short_of(tf)}"
        fw = "bold" if tf == "tbp" else "normal"
        axA.annotate(txt, (x, y), textcoords="offset points", xytext=(7, dy),
                     fontsize=7.5, fontweight=fw, ha="left", va=va,
                     arrowprops=dict(arrowstyle="-", lw=0.4, color="0.6",
                                     shrinkA=0, shrinkB=2))
        if tf == "tbp":
            axA.annotate("best reach, hurt most", (x, y),
                         textcoords="offset points", xytext=(7, -26),
                         fontsize=7, color=ALARM, style="italic")
    axA.set_xlim(0.2, XCAP + 0.3)
    axA.set_xlabel("d$_{min}$ to bound pose (Å)  ·  ← better reachability")
    axA.set_ylabel("conformational spread (Å)")
    axA.set_title("A · does the apo ensemble reach the bound pose?")

    from matplotlib.lines import Line2D
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=TEAL,
                  ms=9, label="augmentation helped (Δ ≥ 0)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=ALARM,
                  ms=9, label="augmentation hurt (Δ < 0)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=GREY,
                  ms=9, label="no benchmark (ΔPCC n/a)")]
    axA.legend(handles=leg, frameon=False, loc="best", fontsize=8,
               title="marker area ∝ |ΔPCC|")
    axA.get_legend().get_title().set_fontsize(8)

    # ---- Panel B: spread for ALL pilots, coloured by reach_status ---------
    # Pilots on x, metric (spread) on y — consistent with the F/I/R panels.
    allm = mech.sort_values("spread", ascending=False)   # tallest first, L→R
    xpos = np.arange(len(allm))
    barcol = [GREY if s == "computed" else "0.80" for s in allm["reach_status"]]
    axB.bar(xpos, allm["spread"], color=barcol, edgecolor="0.3", linewidth=0.6,
            width=0.72)
    # mark computed ones with a teal/alarm dot at the bar top by aug_dP sign
    for x, (_, row) in zip(xpos, allm.iterrows()):
        if row["reach_status"] == "computed":
            d = row["aug_dP"]
            dot = GREY if not np.isfinite(d) else (TEAL if d >= 0 else ALARM)
            axB.plot(x, row["spread"], marker="o", ms=6, color=dot, zorder=4)
    axB.set_xticks(xpos)
    axB.set_xticklabels([short_of(t) for t in allm["pilot"]], rotation=45,
                        ha="right", fontsize=8)
    axB.set_ylabel("conformational spread (Å)")
    axB.set_title("B · spread across all pilots")
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    n_pending = int((mech["reach_status"] != "computed").sum())
    legB = [Line2D([0], [0], marker="o", color="w", markerfacecolor=TEAL, ms=8,
                   label="aug helped (Δ ≥ 0)"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=ALARM, ms=8,
                   label="aug hurt (Δ < 0)"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=GREY, ms=8,
                   label="no benchmark")]
    if n_pending:   # only advertise the pending state when some pilot is pending
        legB.append(Patch(facecolor="0.80", label="reachability pending"))
    axB.legend(handles=legB, frameon=False, loc="lower right", fontsize=7.5)

    fig.suptitle("M1 · apo-ensemble reachability vs augmentation outcome",
                 fontsize=13, y=1.00)
    fig.tight_layout()
    return savefig(fig, "M1_apo_holo_mechanism.png")


def main():
    apply_style()
    mech = pd.read_csv(MECH_CSV)
    print("mechanism rows:", len(mech), "| reach_status counts:",
          mech["reach_status"].value_counts().to_dict())
    print(mech[["pilot", "spread", "reach_ratio", "d_min", "aug_dP",
                "reach_status"]].to_string(index=False))
    p = fig_M1(mech)
    print("wrote", p)


if __name__ == "__main__":
    main()
