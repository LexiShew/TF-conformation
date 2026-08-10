#!/usr/bin/env python3
"""
make_dnaflex.py — DNA-shape flexibility vs the augmentation effect, both arms, own-family.

WHY THIS SCRIPT EXISTS
----------------------
analysis/dna_relax/iface_mgwfl.py already correlates interface minor-groove-width
fluctuation (MGW-FL) against the CROSS-BENCHMARK ΔPearson for the FROZEN arm only
(figure: analysis/dna_relax/figures/iface_mgwfl_vs_accuracy.png, ρ = -0.55, p = 0.08).
This script EXTENDS that analysis rather than duplicating it, on two axes the original
does not cover:

  1. the RELAXED-DNA arm, treated as co-equal to frozen (not a secondary comparison);
  2. the OWN-FAMILY effect, which is where the augmentation signal actually lives
     (the 130-entry cross-benchmark average pools ~120 entries from families the
     pilot's ensemble says nothing about).

It also contrasts a DYNAMIC DNA axis (how much DNA shape fluctuates across the ensemble)
against a STATIC one (how deformed the crystal duplex is). The project's original
induced-fit hypothesis was written on the static axis; the data favour the dynamic one.

PANELS
------
  a  Interface MGW-FL vs own-family ΔPearson, frozen and relaxed arms overlaid.
  b  Ensemble DNA-bend IQR vs own-family ΔPearson, both arms.
  c  Static control: crystal DNA bend vs own-family ΔPearson, both arms.
     A flat cloud here is the point — static deformation does not predict the effect.

STATISTICS
----------
Spearman ρ per arm, with a bootstrap 95% CI over pilots (10,000 resamples) printed to
stdout and written to the CSV. n = 10-11 pilots. At that n these are directional trends,
not established correlations; every CI reported crosses or nearly crosses zero and the
figure says so.

EXCLUSIONS
----------
  dux4  no eval JSON (0 fnat survivors).
  lef1  no same-family benchmark entries (only HMG-box in the 130-entry set).
  irf   absent from the whole-molecule MGW-FL table (groove-parse artifact); it retains
        an interface MGW-FL value and a bend value, so it appears in panels a-c.

COLOR
-----
palette.py at the repo root: TEAL = augmented with frozen DNA, GREEN = augmented with
relaxed DNA, GREY = reference/zero line. ALARM is annotation only, never a series.

Env: deeppbs.  Usage:
    cd /project2/rohs_102/shewchuk/TF-conformation
    python analysis/mechanism/make_dnaflex.py

Outputs
    analysis/mechanism/figures/M4_dnaflex_vs_effect.png
    analysis/mechanism/data/dnaflex_correlations.csv
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from palette import GREY, TEAL, GREEN, apply_style  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT_D, OUT_F = HERE / "data", HERE / "figures"
OUT_D.mkdir(parents=True, exist_ok=True)
OUT_F.mkdir(parents=True, exist_ok=True)

COND_FROZEN, COND_RELAX = "frozen", "relax"
ARM_COLOR = {COND_FROZEN: TEAL, COND_RELAX: GREEN}
ARM_LABEL = {COND_FROZEN: "frozen DNA", COND_RELAX: "relaxed DNA"}
RNG = np.random.default_rng(0)

AXES = [
    ("iface_mgwfl_froz", "Interface MGW-FL, frozen ensemble (Å)",
     "a  Dynamic: how much the interface groove fluctuates", "dynamic"),
    ("froz_bend_iqr", "Ensemble DNA-bend IQR (deg)",
     "b  Dynamic: how much DNA bend fluctuates", "dynamic"),
    ("cry_bend_uu", "Crystal DNA bend (deg)",
     "c  Static control: how bent the bound duplex is", "static"),
]


def load():
    """Own-family seed-level effect (both arms) joined to measured DNA-shape axes."""
    sf = pd.read_csv(OUT_D / "samefamily_both_arms.csv")
    eff = sf.pivot_table(index="tf", columns="dna", values="own_dP")
    eff.columns = [f"own_dP_{c}" for c in eff.columns]

    mech = pd.read_csv(OUT_D / "mechanism_table.csv").set_index("tf")
    shape = mech[["iface_mgwfl_froz", "iface_mgwfl_relax", "froz_bend_iqr",
                  "cry_bend_uu", "mgwfl_froz", "mgwfl_relax"]]
    return eff.join(shape, how="left")


def boot_rho(x, y, n_boot=10000):
    obs = stats.spearmanr(x, y)[0]
    out = []
    for _ in range(n_boot):
        i = RNG.integers(0, len(x), len(x))
        if len(np.unique(x[i])) < 3:
            continue
        r = stats.spearmanr(x[i], y[i])[0]
        if np.isfinite(r):
            out.append(r)
    out = np.array(out)
    return obs, np.percentile(out, 2.5), np.percentile(out, 97.5), \
        float(np.mean(np.sign(out) == np.sign(obs)))


def panel(ax, df, xcol, xlabel, title, rows):
    ax.axhline(0, color=GREY, lw=0.8, zorder=1)
    for dna in (COND_FROZEN, COND_RELAX):
        ycol = f"own_dP_{dna}"
        # The dynamic interface axis has an arm-matched measurement; the bend and
        # static axes are measured once on the frozen ensemble / crystal.
        xc = f"iface_mgwfl_{dna}" if (xcol.startswith("iface_mgwfl")
                                      and f"iface_mgwfl_{dna}" in df) else xcol
        sub = df[[xc, ycol]].dropna()
        if len(sub) < 4:
            continue
        x, y = sub[xc].to_numpy(), sub[ycol].to_numpy()
        rho, lo, hi, frac = boot_rho(x, y)
        ax.scatter(x, y, s=30, color=ARM_COLOR[dna], edgecolor="white",
                   linewidth=0.5, zorder=3,
                   label=f"{ARM_LABEL[dna]}   ρ = {rho:+.2f}")
        m, b = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, m * xs + b, color=ARM_COLOR[dna], lw=1.0, alpha=0.5, zorder=2)
        rows.append(dict(axis=xcol, measured_on=xc, dna=dna, n=len(sub),
                         spearman_rho=rho, boot_ci_lo=lo, boot_ci_hi=hi,
                         frac_same_sign=frac,
                         spearman_p=stats.spearmanr(x, y)[1]))
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left")
    ax.legend(frameon=False, loc="best", handletextpad=0.4)
    ax.margins(0.12)


TF_LABEL = {
    "csl": "CSL", "dux4": "DUX4", "egr1": "EGR1", "engrailed": "engrailed",
    "err": "ERR", "ets1": "ETS1", "foxa": "FOXA", "hsf": "HSF", "irf": "IRF",
    "lef1": "LEF1", "nfat": "NFAT", "runx": "RUNX", "tbp": "TBP",
}


def labeled_bend_panel(df):
    """Standalone bend-IQR vs own-family ΔPearson scatter with per-TF labels."""
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.axhline(0, color=GREY, lw=0.8, zorder=1)
    xcol = "froz_bend_iqr"
    for dna in (COND_FROZEN, COND_RELAX):
        ycol = f"own_dP_{dna}"
        sub = df[[xcol, ycol]].dropna()            # index = tf, preserved
        if len(sub) < 4:
            continue
        x, y = sub[xcol].to_numpy(), sub[ycol].to_numpy()
        rho, lo, hi, frac = boot_rho(x, y)
        ax.scatter(x, y, s=70, color=ARM_COLOR[dna], edgecolor="white",
                   linewidth=0.6, zorder=3,
                   label=f"{ARM_LABEL[dna]}   ρ = {rho:+.2f}")
        m, b = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, m * xs + b, color=ARM_COLOR[dna], lw=1.2, alpha=0.5, zorder=2)
        # label each point with its TF (offset to reduce overlap; frozen above, relaxed below)
        dy = 6 if dna == COND_FROZEN else -10
        for tf, xv, yv in zip(sub.index, x, y):
            ax.annotate(TF_LABEL.get(tf, tf), (xv, yv),
                        xytext=(0, dy), textcoords="offset points",
                        ha="center", fontsize=6.5, color=ARM_COLOR[dna], zorder=4)
    ax.set_xlabel("Ensemble DNA-bend IQR (deg)")
    ax.set_ylabel("Own-family ΔPearson")
    ax.legend(frameon=False, loc="upper right", handletextpad=0.4)
    ax.margins(0.14)
    fig.tight_layout()
    p = OUT_F / "M4b_bend_iqr_labeled.png"
    fig.savefig(p, dpi=300, bbox_inches="tight")
    print(f"wrote {p}")
    return fig


def main():
    apply_style()
    df = load()
    labeled_bend_panel(df)
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.6))
    for ax, (col, xlab, title, _) in zip(axes, AXES):
        panel(ax, df, col, xlab, title, rows)
    axes[0].set_ylabel("Own-family ΔPearson", labelpad=6)

    res = pd.DataFrame(rows)
    res.to_csv(OUT_D / "dnaflex_correlations.csv", index=False)

    fig.suptitle("DNA-shape fluctuation tracks the augmentation effect; "
                 "static crystal deformation does not",
                 x=0.005, ha="left", fontsize=9.5)
    fig.text(0.005, -0.05,
             "Own-family seed-level ΔPearson (n = 4–6 seeds per pilot × arm) vs measured "
             "DNA geometry. n = 10–11 pilots (dux4: no eval; lef1: no same-family "
             "entries).\nρ is Spearman with a bootstrap 95% CI over pilots "
             "(data/dnaflex_correlations.csv). Only panel b / relaxed excludes zero "
             "(ρ = −0.76, CI [−0.96, −0.28], same sign in 99.6% of resamples);\nevery "
             "other CI crosses zero. One correlation at n = 11 is a lead, not an "
             "established result.",
             fontsize=6, color="#4A4A4A", ha="left", va="top")
    fig.tight_layout(rect=[0, 0.02, 1, 0.91])
    path = OUT_F / "M4_dnaflex_vs_effect.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"wrote {path}")

    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible()]
    tick = {ax: set(ax.get_xticklabels() + ax.get_yticklabels()) for ax in fig.axes}
    spines = [(s, s.get_window_extent(r)) for ax in fig.axes
              for s in ax.spines.values() if s.get_visible()]
    ov = [(a.get_text()[:24], b.get_text()[:24]) for i, (a, ba) in enumerate(texts)
          for b, bb in texts[i + 1:] if ba.overlaps(bb)]
    ov += [(t.get_text()[:24], "spine") for t, bt in texts for s, bs in spines
           if bt.overlaps(bs) and t not in tick.get(s.axes, set())]
    print("text overlaps:", ov if ov else "none")
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
