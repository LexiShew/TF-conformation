#!/usr/bin/env python
"""
RMSD of every BioEmu ensemble against its NMR bundle -- distributions, not just
the best-state counts.

plot_best_state.py answers "which deposited state is each frame nearest to".
This answers "how near is it at all", which the counts cannot show: an ensemble
that spreads its frames evenly over 20 states looks identical there whether it
sits 1 A or 6 A away.

Three figures, all from nmr_state_rmsd_{best,pairs}_<atoms>.csv:

  rmsd_overview   every ensemble on one axis, ordered by median d_min -- box +
                  strip per entry. The at-a-glance ranking of how well BioEmu
                  reproduces each deposited bundle.
  rmsd_panels     small multiples, one per ensemble: the distribution of d_min
                  (distance to the CLOSEST deposited state) drawn over the
                  distribution of all state-frame distances. The gap between the
                  two is how much the bundle's own spread is helping.
  rmsd_apo_holo   for pairs.csv entries only: d_min to the apo bundle and to the
                  holo bundle overlaid, one panel per pair. Uses the per-entry
                  CSVs, so it covers any pair whose two members were both run.

Usage:
    python plot_rmsd_distributions.py
    python plot_rmsd_distributions.py --atoms backbone
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from pair_core import load_pairs                            # noqa: E402

_TFCONF = _HERE.parents[2]
sys.path.insert(0, str(_TFCONF))
from palette import GREY, GREY_R, TEAL, TEAL_R, ALARM, apply_style   # noqa: E402


def load_best(path):
    """{(pdb, chain) -> array of per-frame d_min}"""
    out = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            out[(r["pdb_id"], r["chain"])].append(float(r["rmsd_A"]))
    return {k: np.array(v) for k, v in out.items()}


def load_all_pairs(path):
    """{(pdb, chain) -> array of every (state, frame) RMSD}"""
    out = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            out[(r["pdb_id"], r["chain"])].append(float(r["rmsd_A"]))
    return {k: np.array(v) for k, v in out.items()}


def _tag(args):
    """<atoms> or <atoms>_<label> - keeps trimmed and untrimmed outputs apart."""
    return f"{args.atoms}{('_' + args.label) if args.label else ''}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atoms", choices=("ca", "backbone"), default="ca")
    ap.add_argument("--label", default="", help="suffix of the CSVs to read / figures to write")
    ap.add_argument("--pairs", default=_HERE / "pairs.csv")
    ap.add_argument("--output-dir", default=_HERE / "plots")
    args = ap.parse_args()

    best_csv = _HERE / f"nmr_state_rmsd_best_{_tag(args)}.csv"
    pairs_csv = _HERE / f"nmr_state_rmsd_pairs_{_tag(args)}.csv"
    if not best_csv.is_file():
        sys.exit(f"missing {best_csv} -- run compute_nmr_state_rmsd.py first")
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    apply_style()

    best = load_best(best_csv)
    allr = load_all_pairs(pairs_csv) if pairs_csv.is_file() else {}
    keys = sorted(best, key=lambda k: np.median(best[k]))
    labels = [f"{p}/{c}" for p, c in keys]

    # ---- 1. overview: every ensemble, ordered by median d_min ------------
    fig, ax = plt.subplots(figsize=(7.2, 0.26 * len(keys) + 1.1))
    data = [best[k] for k in keys]
    bp = ax.boxplot(data, vert=False, widths=0.62, showfliers=False,
                    patch_artist=True, medianprops=dict(color="white", lw=1.1))
    for patch in bp["boxes"]:
        patch.set(facecolor=TEAL, edgecolor="none")
    for part in ("whiskers", "caps"):
        for ln in bp[part]:
            ln.set(color=GREY, lw=0.8)
    rng = np.random.default_rng(0)
    for i, d in enumerate(data, start=1):
        ax.scatter(d, i + rng.uniform(-0.17, 0.17, len(d)), s=1.6,
                   c=GREY_R[2], alpha=0.35, linewidths=0, zorder=3)
    ax.set_yticks(range(1, len(keys) + 1))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlabel(f"RMSD to the closest NMR state, {args.atoms} ({chr(197)})", fontsize=8)
    ax.set_title(f"BioEmu frame vs its deposited bundle - {len(keys)} ensembles, "
                 f"ordered by median", loc="left", fontsize=8)
    ax.axvline(2.0, color=ALARM, ls=":", lw=0.8)
    ax.annotate("2 " + chr(197), xy=(2.0, 1.005), xycoords=("data", "axes fraction"),
                fontsize=5.5, color=ALARM, ha="center")
    ax.set_xlim(left=0)
    ax.margins(y=0.008)
    fig.tight_layout()
    p = out_dir / f"rmsd_overview_{_tag(args)}.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print(f"wrote {p}")

    # ---- 2. per-ensemble panels -----------------------------------------
    ncols = 5
    nrows = -(-len(keys) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.05 * ncols, 1.5 * nrows),
                             squeeze=False)
    for ax, k in zip(axes.ravel(), keys):
        d = best[k]
        hi = max(np.percentile(allr[k], 99) if k in allr else d.max(), d.max()) * 1.05
        bins = np.linspace(0, hi, 46)
        if k in allr:
            ax.hist(allr[k], bins=bins, color=GREY_R[0], linewidth=0,
                    density=True, label="to every state")
        ax.hist(d, bins=bins, color=TEAL, linewidth=0, alpha=0.9,
                density=True, label="to closest state")
        ax.axvline(np.median(d), color=ALARM, ls="--", lw=0.8)
        ax.set_title(f"{k[0]}/{k[1]}", loc="left", fontsize=6.5)
        ax.annotate(f"med {np.median(d):.2f} " + chr(197), xy=(0.97, 0.9),
                    xycoords="axes fraction", fontsize=5.5, color="#666666", ha="right")
        ax.set_yticks([])
        ax.tick_params(labelsize=5.5)
    for ax in axes.ravel()[len(keys):]:
        ax.set_visible(False)
    axes[0][0].legend(fontsize=5, frameon=False, loc="upper left")
    fig.supxlabel(f"RMSD to NMR state, {args.atoms} ({chr(197)})", fontsize=8)
    fig.tight_layout(rect=(0, 0.012, 1, 1))
    p = out_dir / f"rmsd_panels_{_tag(args)}.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print(f"wrote {p}")

    # ---- 3. apo vs holo, for pairs whose members were both run -----------
    pairs = [r for r in load_pairs(args.pairs, statuses=None)
             if (r["apo_pdb"], r["apo_chain"]) in best and (r["holo_pdb"], r["holo_chain"]) in best]
    if not pairs:
        print("no pairs with both members computed; skipping rmsd_apo_holo")
        return
    ncols = 4
    nrows = -(-len(pairs) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.3 * ncols, 1.6 * nrows), squeeze=False)
    for ax, r in zip(axes.ravel(), pairs):
        a = best[(r["apo_pdb"], r["apo_chain"])]
        h = best[(r["holo_pdb"], r["holo_chain"])]
        hi = max(a.max(), h.max()) * 1.05
        bins = np.linspace(0, hi, 40)
        ax.hist(a, bins=bins, color=GREY_R[1], linewidth=0, alpha=0.8,
                density=True, label="apo ensemble vs apo bundle")
        ax.hist(h, bins=bins, color=TEAL_R[2], linewidth=0, alpha=0.75,
                density=True, label="holo ensemble vs holo bundle")
        ax.set_title(f"{r['pair_id']}  {r['label'][:22]}", loc="left", fontsize=6)
        ax.annotate(f"{np.median(a):.2f} / {np.median(h):.2f} " + chr(197),
                    xy=(0.97, 0.88), xycoords="axes fraction", fontsize=5.5,
                    color="#666666", ha="right")
        ax.set_yticks([])
        ax.tick_params(labelsize=5.5)
    for ax in axes.ravel()[len(pairs):]:
        ax.set_visible(False)
    axes[0][0].legend(fontsize=5, frameon=False, loc="upper right")
    fig.supxlabel(f"RMSD to the closest state of its OWN bundle, {args.atoms} ({chr(197)})",
                  fontsize=8)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    p = out_dir / f"rmsd_apo_holo_{_tag(args)}.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
