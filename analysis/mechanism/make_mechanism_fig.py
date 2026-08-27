#!/usr/bin/env python
"""
make_mechanism_fig.py -- the mechanism figure: where the augmentation signal actually lives.

One figure, four panels, one sentence:
  The augmentation effect is invisible in the 130-entry benchmark average, lives in each
  pilot's OWN family, is not explained by baseline headroom, and tracks how much DNA
  shape FLUCTUATES across the ensemble rather than how deformed the crystal duplex is.

Panels
  a  Own-family vs other-family augmentation effect, per pilot (seed-paired).
     Shows the averaging-away: the own-family effect is large where the cross-benchmark
     effect is ~0.
  b  Within-pilot same-vs-other contrast, seed level, with 95% CI.
     The test REVIEW_figure_scripts.md asked for and fig6 never performed.
  c  Headroom control: own-family effect vs baseline own-family accuracy.
     A flat/negative-free relationship is what excludes "low baseline has more room".
  d  Dynamic vs static DNA axis: own-family effect vs ensemble DNA-bend IQR (dynamic)
     and vs crystal bend (static), same y-axis.

Unit of replication is the SEED throughout; CIs are seed-level t intervals.
n = 11 pilots (dux4 has no eval; lef1 has no same-family benchmark entries).

Usage
-----
    source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    cd /project2/rohs_102/shewchuk/TF-conformation
    python analysis/mechanism/make_mechanism_fig.py

Output
------
    analysis/mechanism/figures/M2_where_the_signal_lives.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DATA = Path(__file__).resolve().parent / "data"
FIGS = Path(__file__).resolve().parent / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

# Cool-pastel palette, consistent with palette.py conventions in this repo.
C_SELF = "#4C9BC0"      # own-family  (blue)
C_OTHER = "#B0B7BE"     # other-family (neutral grey)
C_POS = "#5FB89A"       # positive contrast (green)
C_NEG = "#C77FA6"       # negative contrast (pink)
C_DYN = "#8E7CC3"       # dynamic DNA axis (purple)
C_STAT = "#C9C2DE"      # static DNA axis (pale purple)
GREY = "#4A4A4A"

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "figure.dpi": 300,
})


def load():
    sf = pd.read_csv(DATA / "selffam_effects.csv")
    sf = sf[(sf.dna == "frozen")]
    self_e = (sf[sf.selffam].set_index("tf")
              [["dP_mean", "dP_lo", "dP_hi", "dP_p", "n_seeds", "n_entries"]]
              .add_prefix("self_"))
    other_e = (sf[~sf.selffam].set_index("tf")[["dP_mean"]].add_prefix("other_"))
    ct = pd.read_csv(DATA / "same_vs_other_contrast.csv").set_index("tf")
    mech = pd.read_csv(DATA / "mechanism_table.csv").set_index("tf")
    conf = pd.read_csv(DATA / "confound_tests.csv")

    pe = pd.read_csv(ROOT / "analysis/figure_scripts/perseed_perentry.csv")
    base = (pe[(pe.arm == "baseline") & (pe.dna == "frozen") & pe.selffam]
            .groupby("tf").m_pearsonr.mean().rename("base_self_P"))

    df = (self_e.join(other_e).join(ct[["contrast", "p", "q_bh", "n_neg_seeds"]])
          .join(mech[["froz_bend_iqr", "cry_bend_uu", "froz_dP_mean"]])
          .join(base).dropna(subset=["self_dP_mean"]))
    return df.sort_values("self_dP_mean", ascending=False), conf


def panel_a(ax, df):
    y = np.arange(len(df))[::-1]
    ax.axvline(0, color=GREY, lw=0.8, zorder=1)
    ax.hlines(y, df.other_dP_mean, df.self_dP_mean,
              color="#D8DDE1", lw=1.6, zorder=2)
    ax.scatter(df.other_dP_mean, y, s=22, color=C_OTHER, zorder=3,
               label="other families (~120 entries)")
    ax.scatter(df.self_dP_mean, y, s=34, color=C_SELF, zorder=4,
               edgecolor="white", linewidth=0.5, label="own family")
    ax.set_yticks(y)
    ax.set_yticklabels(df.index)
    ax.set_xlabel("Augmentation effect  (ΔPearson, augmented − baseline)")
    ax.set_title("a  The signal is in the pilot's own family", loc="left",
                 fontweight="bold")
    ax.legend(frameon=False, loc="lower right", handletextpad=0.4)
    ax.margins(x=0.10, y=0.04)


def panel_b(ax, df):
    d = df.dropna(subset=["contrast"]).sort_values("contrast")
    y = np.arange(len(d))
    colors = [C_POS if v > 0 else C_NEG for v in d.contrast]
    ax.axvline(0, color=GREY, lw=0.8, zorder=1)
    ax.barh(y, d.contrast, color=colors, height=0.62, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(d.index)
    for i, (tf, row) in enumerate(d.iterrows()):
        if row.p < 0.05:
            ax.text(row.contrast + 0.006, i, f"p={row.p:.3f}",
                    va="center", ha="left", fontsize=6, color=GREY)
    ax.set_xlabel("Own-family minus other-family ΔPearson")
    ax.set_title("b  Within-pilot contrast, tested at the seed level", loc="left",
                 fontweight="bold")
    ax.margins(x=0.22, y=0.04)


def panel_c(ax, df):
    ax.axhline(0, color=GREY, lw=0.8, zorder=1)
    ax.scatter(df.base_self_P, df.self_dP_mean, s=34, color=C_SELF,
               edgecolor="white", linewidth=0.5, zorder=3)
    rho, p = stats.spearmanr(df.base_self_P, df.self_dP_mean)
    for tf in ("ets1", "csl", "tbp", "egr1"):
        if tf in df.index:
            r = df.loc[tf]
            ax.annotate(tf, (r.base_self_P, r.self_dP_mean),
                        textcoords="offset points", xytext=(5, 4),
                        fontsize=6, color=GREY)
    ax.set_xlabel("Baseline accuracy on own family (Pearson)")
    ax.set_ylabel("Own-family ΔPearson")
    ax.set_title(f"c  Not headroom  (ρ = {rho:+.2f}, p = {p:.2f})", loc="left",
                 fontweight="bold")
    ax.margins(0.10)


def panel_d(ax, df):
    ax.axhline(0, color=GREY, lw=0.8, zorder=1)
    d1 = df.dropna(subset=["froz_bend_iqr"])
    rho_d, p_d = stats.spearmanr(d1.froz_bend_iqr, d1.self_dP_mean)
    ax.scatter(d1.froz_bend_iqr, d1.self_dP_mean, s=34, color=C_DYN,
               edgecolor="white", linewidth=0.5, zorder=3,
               label=f"ensemble DNA-bend IQR   ρ = {rho_d:+.2f}")
    m, b = np.polyfit(d1.froz_bend_iqr, d1.self_dP_mean, 1)
    xs = np.linspace(d1.froz_bend_iqr.min(), d1.froz_bend_iqr.max(), 50)
    ax.plot(xs, m * xs + b, color=C_DYN, lw=1.0, alpha=0.55, zorder=2)
    for tf in ("ets1", "tbp", "lef1", "csl"):
        if tf in d1.index:
            r = d1.loc[tf]
            ax.annotate(tf, (r.froz_bend_iqr, r.self_dP_mean),
                        textcoords="offset points", xytext=(5, 4),
                        fontsize=6, color=GREY)
    ax.set_xlabel("Ensemble DNA-bend IQR, frozen pipeline (deg)")
    ax.set_ylabel("Own-family ΔPearson")

    d2 = df.dropna(subset=["cry_bend_uu"])
    rho_s, p_s = stats.spearmanr(d2.cry_bend_uu, d2.self_dP_mean)
    ax.set_title(f"d  Dynamic DNA shape, not static deformation\n"
                 f"    (static crystal bend: ρ = {rho_s:+.2f}, p = {p_s:.2f})",
                 loc="left", fontweight="bold")
    # Legend goes lower-left: the upper-right corner holds the ets1/csl point labels.
    ax.legend(frameon=False, loc="lower left", handletextpad=0.4)
    ax.margins(0.12)


def main():
    df, conf = load()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4))
    panel_a(axes[0, 0], df)
    panel_b(axes[0, 1], df)
    panel_c(axes[1, 0], df)
    panel_d(axes[1, 1], df)
    fig.suptitle("Where the conformational-augmentation signal lives",
                 fontsize=9.5, x=0.01, ha="left", fontweight="bold")
    fig.text(0.01, 0.005,
             "n = 11 pilots (dux4: no eval; lef1: no same-family benchmark entries). "
             "Seed-paired effects, frozen-DNA pipeline; seed is the unit of replication.",
             fontsize=6, color=GREY, ha="left")
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])
    out = FIGS / "M2_where_the_signal_lives.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")

    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible()]
    ticklabels = {ax: set(ax.get_xticklabels() + ax.get_yticklabels()) for ax in fig.axes}
    spines = [(s, s.get_window_extent(r)) for ax in fig.axes
              for s in ax.spines.values() if s.get_visible()]
    ov = [(a.get_text(), b.get_text()) for i, (a, ba) in enumerate(texts)
          for b, bb in texts[i + 1:] if ba.overlaps(bb)]
    ov += [(t.get_text(), "spine") for t, bt in texts for s, bs in spines
           if bt.overlaps(bs) and t not in ticklabels.get(s.axes, set())]
    print("text overlaps:", ov if ov else "none")


if __name__ == "__main__":
    main()
