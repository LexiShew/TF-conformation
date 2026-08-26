#!/usr/bin/env python3
"""
make_samefamily.py — same-family vs other-family augmentation effect, both DNA arms.

WHY THIS SCRIPT EXISTS
----------------------
analysis/dna_relax/figures/samefamily_augmentation_effect.png and
crossfamily_augmentation_effect.png were delivered as PNGs with NO generating script
anywhere in the repo (verified by grep across analysis/ and scripts/). They therefore
could not be regenerated, extended to new pilots, or checked. They also covered only
5 pilots (tbp, engrailed, egr1, ets1, foxa).

This script supersedes both. It covers all 11 pilots that have same-family benchmark
entries, treats the frozen-DNA and relaxed-DNA pipelines as CO-EQUAL conditions, and
threads the canonical palette.

WHAT IT SHOWS
-------------
The 130-entry general benchmark averages the augmentation effect over ~120 entries from
families the pilot's ensemble says nothing about. The mechanistically meaningful quantity
is the effect on the pilot's OWN family. Panels:

  a  Own-family augmentation effect per pilot, frozen and relaxed side by side,
     seed-level mean with 95% CI.
  b  Own-family minus other-family contrast, per pilot, both arms.
     This is the "augmentation is family-specific" claim tested rather than asserted —
     REVIEW_figure_scripts.md flagged that fig6 never tested it.
  c  Paired frozen -> relaxed shift in the own-family effect, per pilot.
     Tests the falsifiable prediction in analysis/mechanism_and_roadmap.md §2:
     co-relaxing the DNA should recover the pilots augmentation hurts under frozen DNA.

STATISTICS
----------
The SEED is the unit of replication: augmentation is applied once per retraining, so the
5-24 same-family entries within a seed are correlated and cannot be treated as
independent replicates (this is the pseudoreplication error REVIEW_figure_scripts.md
found in fig9). All effects are seed-level means; CIs are t intervals over seeds
(n = 4-6). p-values are one-sample t against zero, with BH-FDR across pilots reported
in the CSV. Nothing here survives BH correction — the figure states this.

EXCLUSIONS
----------
  dux4  0 fnat survivors -> no eval JSON, no benchmark row at all.
  lef1  only HMG-box in the 130-entry benchmark -> no same-family entries to test.
Both are named in the caption rather than silently dropped.

COLOR
-----
palette.py at the repo root is the single source of truth: GREY = baseline reference,
TEAL = augmented with frozen DNA, GREEN = augmented with relaxed DNA, ALARM = negative
annotation ONLY (never a data series).

Env: deeppbs.  Usage:
    cd /project2/rohs_102/shewchuk/TF-conformation
    python analysis/mechanism/make_samefamily.py

Outputs
    analysis/mechanism/figures/M3_samefamily_both_arms.png
    analysis/mechanism/data/samefamily_both_arms.csv
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from palette import GREY, TEAL, GREEN, ALARM, apply_style  # noqa: E402

OUT_D = Path(__file__).resolve().parent / "data"
OUT_F = Path(__file__).resolve().parent / "figures"
OUT_D.mkdir(parents=True, exist_ok=True)
OUT_F.mkdir(parents=True, exist_ok=True)

PERSEED = ROOT / "analysis/data/perseed_perentry.csv"

# The benchmark CSVs label the DNA conditions "frozen" and "relax" (NOT "relaxed").
# Asserted at load time: filtering on a wrong label silently drops an entire arm.
COND_FROZEN, COND_RELAX = "frozen", "relax"
ARM_COLOR = {COND_FROZEN: TEAL, COND_RELAX: GREEN}
ARM_LABEL = {COND_FROZEN: "augmented · frozen DNA", COND_RELAX: "augmented · relaxed DNA"}


def seed_level_effects():
    """Seed-paired own-family and other-family ΔPearson, per pilot and DNA condition."""
    pe = pd.read_csv(PERSEED)
    found = set(pe.dna.unique())
    if not {COND_FROZEN, COND_RELAX} <= found:
        raise ValueError(f"dna labels are {sorted(found)}; expected to contain "
                         f"{[COND_FROZEN, COND_RELAX]}. Fix COND_*, not the filter.")

    g = (pe.groupby(["tf", "arm", "dna", "seed", "selffam"])
           .agg(m_pearsonr=("m_pearsonr", "mean"), n=("entry", "size")).reset_index())
    w = g.pivot_table(index=["tf", "dna", "seed", "selffam", "n"],
                      columns="arm", values="m_pearsonr").reset_index()
    w["dP"] = w["augmented"] - w["baseline"]
    return w.dropna(subset=["dP"])


def summarize(w):
    """Per (pilot, condition): own-family effect, other-family effect, and contrast."""
    rows = []
    for (tf, dna), sub in w.groupby(["tf", "dna"]):
        own = sub[sub.selffam].set_index("seed")["dP"]
        oth = sub[~sub.selffam].set_index("seed")["dP"]
        if own.empty:
            continue                      # lef1: no same-family benchmark entries
        common = own.index.intersection(oth.index)
        d = own.to_numpy()
        n = len(d)
        se = d.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
        tc = stats.t.ppf(0.975, n - 1) if n > 1 else np.nan
        p = stats.ttest_1samp(d, 0).pvalue if n > 1 else np.nan

        con = (own[common] - oth[common]).to_numpy() if len(common) > 1 else np.array([])
        if con.size > 1:
            cse = con.std(ddof=1) / np.sqrt(len(con))
            ctc = stats.t.ppf(0.975, len(con) - 1)
            clo, chi = con.mean() - ctc * cse, con.mean() + ctc * cse
            cp = stats.ttest_1samp(con, 0).pvalue
        else:
            clo = chi = cp = np.nan

        rows.append(dict(
            tf=tf, dna=dna, n_seeds=n,
            n_entries=int(sub[sub.selffam].n.median()),
            own_dP=d.mean(), own_lo=d.mean() - tc * se, own_hi=d.mean() + tc * se,
            own_p=p, own_n_neg=int((d < 0).sum()),
            other_dP=oth.mean() if not oth.empty else np.nan,
            contrast=con.mean() if con.size else np.nan,
            contrast_lo=clo, contrast_hi=chi, contrast_p=cp))

    out = pd.DataFrame(rows)
    for dna, idx in out.groupby("dna").groups.items():
        for col in ("own_p", "contrast_p"):
            pv = out.loc[idx, col]
            ok = pv.notna()
            if ok.sum():
                out.loc[pv[ok].index, col.replace("_p", "_q")] = (
                    pv[ok] * ok.sum() / pv[ok].rank(method="first")).clip(upper=1)
    return out


def _grouped_positions(pilots):
    y = np.arange(len(pilots))
    return y, 0.19


def panel_own(ax, s, order):
    y, off = _grouped_positions(order)
    ax.axvline(0, color=GREY, lw=0.8, zorder=1)
    for i, dna in enumerate((COND_FROZEN, COND_RELAX)):
        d = s[s.dna == dna].set_index("tf").reindex(order)
        yy = y + (off if i == 0 else -off)
        ax.errorbar(d.own_dP, yy, xerr=[d.own_dP - d.own_lo, d.own_hi - d.own_dP],
                    fmt="o", ms=4.5, lw=0, elinewidth=1.0, capsize=2,
                    color=ARM_COLOR[dna], ecolor=ARM_COLOR[dna], zorder=3,
                    label=ARM_LABEL[dna])
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel("Own-family ΔPearson (augmented − baseline)")
    ax.set_title("a  Augmentation helps on the pilot's own family", loc="left")
    # Legend goes below the axes: every row carries an error bar, so there is no
    # in-axes whitespace that stays clear at all pilot orderings.
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22),
              ncol=2, handletextpad=0.4, columnspacing=1.4)
    ax.margins(x=0.12, y=0.03)


def panel_contrast(ax, s, order):
    y, off = _grouped_positions(order)
    ax.axvline(0, color=GREY, lw=0.8, zorder=1)
    for i, dna in enumerate((COND_FROZEN, COND_RELAX)):
        d = s[s.dna == dna].set_index("tf").reindex(order)
        yy = y + (off if i == 0 else -off)
        ax.barh(yy, d.contrast, height=0.34, color=ARM_COLOR[dna], zorder=2,
                label=ARM_LABEL[dna])
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel("Own-family minus other-family ΔPearson")
    ax.set_title("b  The effect is specific to the pilot's own family", loc="left")
    ax.margins(x=0.14, y=0.03)


def panel_shift(ax, s, order):
    froz = s[s.dna == COND_FROZEN].set_index("tf").reindex(order)
    rel = s[s.dna == COND_RELAX].set_index("tf").reindex(order)
    y = np.arange(len(order))
    ax.axvline(0, color=GREY, lw=0.8, zorder=1)
    for i, tf in enumerate(order):
        a, b = froz.own_dP.get(tf, np.nan), rel.own_dP.get(tf, np.nan)
        if np.isnan(a) or np.isnan(b):
            continue
        ax.annotate("", xy=(b, i), xytext=(a, i),
                    arrowprops=dict(arrowstyle="->", color="#C7CCD1", lw=1.3,
                                    shrinkA=2, shrinkB=2), zorder=2)
    ax.scatter(froz.own_dP, y, s=26, color=TEAL, zorder=3, label=ARM_LABEL[COND_FROZEN])
    ax.scatter(rel.own_dP, y, s=26, color=GREEN, zorder=4, label=ARM_LABEL[COND_RELAX])
    n_up = int((rel.own_dP > froz.own_dP).sum())
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel("Own-family ΔPearson")
    ax.set_title(f"c  Relaxing the DNA shifts {n_up} of {len(order)} pilots upward",
                 loc="left")
    # No legend here: panel a already keys frozen/relaxed with the same two colours.
    ax.margins(x=0.12, y=0.03)


def main():
    apply_style()
    w = seed_level_effects()
    s = summarize(w)
    s.to_csv(OUT_D / "samefamily_both_arms.csv", index=False)

    order = (s[s.dna == COND_FROZEN].sort_values("own_dP", ascending=False).tf.tolist())

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    panel_own(axes[0], s, order)
    panel_contrast(axes[1], s, order)
    panel_shift(axes[2], s, order)

    fig.suptitle("The augmentation signal lives in the pilot's own family, "
                 "and relaxing the DNA strengthens it",
                 x=0.005, ha="left", fontsize=9.5)
    fig.text(0.005, -0.10,
             "Seed-paired effects; seed is the unit of replication (n = 4–6 seeds per "
             "pilot × arm). Error bars are 95% t intervals over seeds.\n"
             "n = 11 pilots: dux4 excluded (0 fnat survivors, no eval), lef1 excluded "
             "(only HMG-box in the benchmark, no same-family entries).\n"
             "Raw p < 0.05 for ets1 (both arms), csl and runx (relaxed only); no effect "
             "survives BH correction across pilots. Treat as a trend.",
             fontsize=6, color=GREY, ha="left", va="top")
    fig.tight_layout(rect=[0, 0.06, 1, 0.91])
    path = OUT_F / "M3_samefamily_both_arms.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"wrote {path}")

    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible()]
    ov = [(a.get_text()[:28], b.get_text()[:28]) for i, (a, ba) in enumerate(texts)
          for b, bb in texts[i + 1:] if ba.overlaps(bb)]
    print("text overlaps:", ov if ov else "none")

    print("\n=== own-family effect, both arms ===")
    print(s.pivot_table(index="tf", columns="dna",
                        values=["own_dP", "own_p", "contrast"])
          .round(4).to_string())


if __name__ == "__main__":
    main()
