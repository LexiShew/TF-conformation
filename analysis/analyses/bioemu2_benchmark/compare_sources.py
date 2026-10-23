#!/usr/bin/env python3
"""Does BioEmu2-sourced augmentation beat BioEmu1-sourced augmentation?

Both sources were run through the SAME Stages 4-7 (see lib/common.sh CONF_SRC),
so the only difference is where the augmenting conformations came from:

  bioemu1  Stage 1-2 BioEmu + HPacker + redock, then Stage 3 minimization
  bioemu2  BioEmu2 complex sampling, gated raw (no minimization - raw frames
           clear the fnat gate at 91-100% for these pilots)

Metric is the repo's own augmentation effect, matching
analyses/conformation/extract_metrics.py: per benchmark entry,
  dPearson = pearsonr(augmented) - pearsonr(baseline)
averaged over the 130-entry ID benchmark. Positive = augmentation helped.

Like-for-like seeds (important)
-------------------------------
id_benchmark_<tf>.json carries 12 checkpoints: the default-seed pair plus
_s1.._s5. id_benchmark_<tf>_be2.json carries only the default-seed pair, because
the BioEmu2 arm has so far been run at the single default seed 42. Comparing a
seed-AVERAGED BioEmu1 number against a single BioEmu2 run would conflate seed
variance with the source effect, so this script:

  * compares DEFAULT SEED against DEFAULT SEED as the headline, and
  * reports BioEmu1's across-seed spread (s1..s5) as the noise band, so a
    reader can see whether the source difference exceeds seed-to-seed scatter.

Any BioEmu2-vs-BioEmu1 gap smaller than that band is not yet interpretable and
needs the BioEmu2 arm re-run at 5 seeds.

Outputs (theme-local, per analysis/docs/LAYOUT.md - canonical tables untouched):
  data/bioemu2_vs_bioemu1_benchmark.csv
  data/bioemu2_vs_bioemu1_perentry.csv
  figures/B1_source_augmentation_effect.png
"""
import json
import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "common"))
from fig_common import TFCONF, EVAL_DIR, label_of  # noqa: E402

DATA = os.path.join(_HERE, "data")
FIGS = os.path.join(_HERE, "figures")
PILOTS = ["ets1", "tbp", "egr1", "engrailed", "foxa", "lef1", "csl", "err",
          "nfat", "runx", "hsf", "ndt80"]
METRIC = "pearsonr"


def arm_entries(res, arm, pilot, suffix):
    """Per-entry metric dict for exactly one checkpoint: <arm>_<pilot>_fold0<suffix>."""
    key = f"{arm}_{pilot}_fold0{suffix}"
    return res.get(key)


def seed_keys(res, arm, pilot):
    """The _sN checkpoints of one arm (excludes the unsuffixed default-seed run)."""
    pre = f"{arm}_{pilot}_fold0_s"
    return sorted(k for k in res if k.startswith(pre))


def delta(res, pilot, suffix, entries=None):
    """mean over entries of augmented - baseline for METRIC, plus the per-entry series."""
    base = arm_entries(res, "baseline", pilot, suffix)
    aug = arm_entries(res, "augmented", pilot, suffix)
    if base is None or aug is None:
        return None, None
    keys = sorted(set(base) & set(aug)) if entries is None else entries
    d = {k: aug[k][METRIC] - base[k][METRIC] for k in keys
         if METRIC in aug.get(k, {}) and METRIC in base.get(k, {})}
    if not d:
        return None, None
    return float(np.mean(list(d.values()))), d


def plot(df):
    """B1 - per-pilot augmentation effect by conformation source.

    BioEmu1's across-seed range is drawn as a grey span so the reader can judge
    each source difference against seed-to-seed scatter rather than against zero.
    Only a difference that clears that span is interpretable at one BioEmu2 seed.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    d = df.sort_values("diff_be2_minus_be1", ascending=False).reset_index(drop=True)
    x = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(1.15 * len(d) + 3.2, 5.0))

    for i, r in d.iterrows():
        if r["be1_seed_n"]:
            ax.add_patch(plt.Rectangle(
                (i - 0.34, r["be1_seed_min"]), 0.68,
                max(r["be1_seed_max"] - r["be1_seed_min"], 1e-6),
                facecolor="0.86", edgecolor="none", zorder=1))
    ax.axhline(0, color="0.35", lw=1.0, zorder=2)

    ax.scatter(x - 0.13, d["dP_bioemu1"], s=86, color="#2E7FA8",
               edgecolor="white", lw=1.1, zorder=4, label="BioEmu 1 (docked + minimized)")
    ax.scatter(x + 0.13, d["dP_bioemu2"], s=86, color="#D98BA5",
               edgecolor="white", lw=1.1, zorder=4, label="BioEmu2 (complex, unminimized)")
    for i, r in d.iterrows():
        ax.plot([i - 0.13, i + 0.13], [r["dP_bioemu1"], r["dP_bioemu2"]],
                color="0.55", lw=1.0, zorder=3)
        if r["exceeds_seed_band"] is True or str(r["exceeds_seed_band"]) == "True":
            ax.annotate("*", (i, max(r["dP_bioemu1"], r["dP_bioemu2"]) + 0.012),
                        ha="center", va="bottom", fontsize=15, color="#7B4B6B", zorder=5)

    ax.set_xticks(x)
    ax.set_xticklabels([p.upper() if p != "engrailed" else "engrailed"
                        for p in d["pilot"]], fontsize=10)
    ax.set_ylabel(r"augmentation effect  $\Delta$Pearson  (augmented $-$ baseline)",
                  fontsize=11)
    ax.set_title("Augmentation effect by conformation source, over the 130-entry ID benchmark",
                 fontsize=12.5, pad=11)
    handles = [Line2D([], [], marker="o", ls="", color="#2E7FA8", ms=9,
                      label="BioEmu 1 (docked + minimized)"),
               Line2D([], [], marker="o", ls="", color="#D98BA5", ms=9,
                      label="BioEmu2 (complex, unminimized)"),
               plt.Rectangle((0, 0), 1, 1, facecolor="0.86", edgecolor="none",
                             label="BioEmu 1 across-seed range (s1-s5)")]
    ax.legend(handles=handles, frameon=False, fontsize=9.5, loc="lower left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color="0.9", lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    fig.text(0.01, 0.034,
             "Both sources run through identical Stages 4-7 (lib/common.sh CONF_SRC); only the augmenting conformations differ. "
             "Positive = augmentation helped.",
             fontsize=7.8, color="0.3")
    fig.text(0.01, 0.010,
             "Default seed 42 for both arms. BioEmu2 has ONE seed, so * marks the only pilots whose source difference exceeds "
             "BioEmu 1's own across-seed range; the rest are within seed noise.",
             fontsize=7.8, color="0.3")
    fig.tight_layout(rect=[0, 0.075, 1, 1])
    out = os.path.join(FIGS, "B1_source_augmentation_effect.png")
    os.makedirs(FIGS, exist_ok=True)
    fig.savefig(out, dpi=200)
    print("wrote", out)


def main():
    rows, per_entry = [], []
    for tf in PILOTS:
        p1 = os.path.join(EVAL_DIR, f"id_benchmark_{tf}.json")
        p2 = os.path.join(EVAL_DIR, f"id_benchmark_{tf}_be2.json")
        if not (os.path.exists(p1) and os.path.exists(p2)):
            print(f"  {tf:11s} SKIP ("
                  f"{'no bioemu1 eval' if not os.path.exists(p1) else 'no bioemu2 eval'})")
            continue
        r1 = json.load(open(p1))["results"]
        r2 = json.load(open(p2))["results"]

        d1, s1 = delta(r1, tf, "")
        d2, s2 = delta(r2, tf, "_be2")
        if d1 is None or d2 is None:
            print(f"  {tf:11s} SKIP (missing a default-seed arm)")
            continue
        # common entry set, so the two sources are compared on identical entries
        common = sorted(set(s1) & set(s2))
        d1, s1 = delta(r1, tf, "", common)
        d2, s2 = delta(r2, tf, "_be2", common)

        # BioEmu1 seed band from the _sN pairs
        seed_deltas = []
        for k in seed_keys(r1, "augmented", tf):
            sfx = k.replace(f"augmented_{tf}_fold0", "")
            dk, _ = delta(r1, tf, sfx, common)
            if dk is not None:
                seed_deltas.append(dk)

        rows.append(dict(
            pilot=tf, n_entries=len(common),
            dP_bioemu1=round(d1, 5), dP_bioemu2=round(d2, 5),
            diff_be2_minus_be1=round(d2 - d1, 5),
            be1_seed_n=len(seed_deltas),
            be1_seed_mean=round(float(np.mean(seed_deltas)), 5) if seed_deltas else "",
            be1_seed_sd=round(float(np.std(seed_deltas, ddof=1)), 5) if len(seed_deltas) > 1 else "",
            be1_seed_min=round(float(np.min(seed_deltas)), 5) if seed_deltas else "",
            be1_seed_max=round(float(np.max(seed_deltas)), 5) if seed_deltas else "",
            exceeds_seed_band=(
                bool(abs(d2 - d1) > (np.max(seed_deltas) - np.min(seed_deltas)))
                if seed_deltas else ""),
        ))
        for k in common:
            per_entry.append(dict(pilot=tf, entry=k,
                                  dP_bioemu1=s1[k], dP_bioemu2=s2[k]))
        print(f"  {tf:11s} n={len(common):3d}  dP_be1={d1:+.4f}  dP_be2={d2:+.4f}  "
              f"diff={d2-d1:+.4f}  be1 seed range=["
              f"{min(seed_deltas):+.4f},{max(seed_deltas):+.4f}]" if seed_deltas else "")

    if not rows:
        print("nothing to compare")
        return
    df = pd.DataFrame(rows)
    os.makedirs(DATA, exist_ok=True)
    df.to_csv(os.path.join(DATA, "bioemu2_vs_bioemu1_benchmark.csv"), index=False)
    pd.DataFrame(per_entry).to_csv(
        os.path.join(DATA, "bioemu2_vs_bioemu1_perentry.csv"), index=False)
    plot(df)
    print()
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
