#!/usr/bin/env python
"""
Apo beside holo, for every pair: same axes, same scale, nothing overlaid.

plot_rmsd_distributions.py overlays the two members of a pair in one panel,
which hides the closest-state structure. Here each pair gets apo on the left and
holo on the right, on SHARED axis limits so the two panels are directly
comparable -- an apo ensemble that fits its bundle tightly next to a holo
ensemble that does not is then a visual fact, not an arithmetic one.

Per pair (pair_side_by_side_<pair>_<tag>.png):

  top     closest-state counts -- which deposited state each BioEmu frame is
          nearest to, with the uniform expectation marked. Shared y-limit.
  bottom  RMSD to that closest state. Shared x-limit, medians marked.

Across pairs (pair_side_by_side_all_<tag>.png): one row per pair, the apo and
holo RMSD distributions as a box pair, so the whole set ranks at a glance.

Both read the per-entry CSVs from compute_nmr_state_rmsd.py, so pass the same
--atoms/--label used to build them. Use --label core for the terminal-trimmed
numbers; untrimmed, entries with long tails measure their tails.

Usage:
    python plot_pair_side_by_side.py --label core
    python plot_pair_side_by_side.py --atoms backbone --label core
"""
import argparse
import csv
import sys
from collections import Counter, defaultdict
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
from palette import GREY, TEAL_R, ALARM, apply_style        # noqa: E402

APO, HOLO = TEAL_R[0], TEAL_R[2]      # BioEmu, light = apo construct, dark = holo


def _tag(args):
    return f"{args.atoms}{('_' + args.label) if args.label else ''}"


def load_best(path):
    """{(pdb, chain) -> dict(dmin=array, best=Counter, n_states=int)}"""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            k = (r["pdb_id"], r["chain"])
            e = out.setdefault(k, {"dmin": [], "best": Counter(), "n_states": 0})
            e["dmin"].append(float(r["rmsd_A"]))
            e["best"][int(r["best_state"])] += 1
            e["n_states"] = int(r["n_states"])
    for e in out.values():
        e["dmin"] = np.array(e["dmin"])
    return out


def draw_counts(ax, e, color, title):
    states = np.arange(1, e["n_states"] + 1)
    counts = np.array([e["best"].get(int(s), 0) for s in states])
    n = counts.sum()
    ax.bar(states, counts, width=0.8, color=color, linewidth=0)
    ax.axhline(n / e["n_states"], color=GREY, ls="--", lw=0.8)
    ax.set_title(title, loc="left", fontsize=6.5)
    ax.annotate(f"{n} frames · {int((counts > 0).sum())}/{e['n_states']} states used",
                xy=(0, 1.02), xycoords="axes fraction", fontsize=5, color="#666666")
    ax.set_xlim(0.4, e["n_states"] + 0.6)
    ax.tick_params(labelsize=5.5)
    return counts


def draw_dmin(ax, e, color, bins):
    d = e["dmin"]
    ax.hist(d, bins=bins, color=color, linewidth=0)
    ax.axvline(np.median(d), color=ALARM, ls="--", lw=0.9)
    ax.annotate(f"med {np.median(d):.2f} " + chr(197), xy=(0.97, 0.88),
                xycoords="axes fraction", fontsize=5.5, color="#666666", ha="right")
    ax.tick_params(labelsize=5.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atoms", choices=("ca", "backbone"), default="ca")
    ap.add_argument("--label", default="", help="suffix of the CSVs to read")
    ap.add_argument("--pairs", default=_HERE / "pairs.csv")
    ap.add_argument("--output-dir", default=_HERE / "plots")
    args = ap.parse_args()

    best_csv = _HERE / f"nmr_state_rmsd_best_{_tag(args)}.csv"
    if not best_csv.is_file():
        sys.exit(f"missing {best_csv} -- run compute_nmr_state_rmsd.py "
                 f"--atoms {args.atoms}"
                 + (f" --trim-termini 3.0 --label {args.label}" if args.label else ""))
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    apply_style()

    best = load_best(best_csv)
    pairs = [r for r in load_pairs(args.pairs, statuses=None)
             if (r["apo_pdb"], r["apo_chain"]) in best
             and (r["holo_pdb"], r["holo_chain"]) in best]
    if not pairs:
        sys.exit("no pair has both members in the CSV")
    print(f"{len(pairs)} pairs with both members computed")

    for r in pairs:
        a = best[(r["apo_pdb"], r["apo_chain"])]
        h = best[(r["holo_pdb"], r["holo_chain"])]
        fig, axes = plt.subplots(2, 2, figsize=(6.6, 3.4))
        ca = draw_counts(axes[0][0], a, APO, f"apo  {r['apo_pdb']}/{r['apo_chain']}")
        ch = draw_counts(axes[0][1], h, HOLO, f"holo  {r['holo_pdb']}/{r['holo_chain']}")
        ytop = max(ca.max(), ch.max()) * 1.15                 # shared, so bars compare
        for ax in axes[0]:
            ax.set_ylim(0, ytop)
            ax.set_xlabel("NMR state", fontsize=6.5)
        axes[0][0].set_ylabel("frames closest\nto this state", fontsize=6.5)

        hi = max(a["dmin"].max(), h["dmin"].max()) * 1.04     # shared x, so widths compare
        bins = np.linspace(0, hi, 42)
        draw_dmin(axes[1][0], a, APO, bins)
        draw_dmin(axes[1][1], h, HOLO, bins)
        for ax in axes[1]:
            ax.set_xlim(0, hi)
            ax.set_xlabel(f"RMSD to closest state, {args.atoms} ({chr(197)})", fontsize=6.5)
        axes[1][0].set_ylabel("frames", fontsize=6.5)
        ymax = max(ax.get_ylim()[1] for ax in axes[1])
        for ax in axes[1]:
            ax.set_ylim(0, ymax)

        fig.suptitle(f"{r['pair_id']}  ·  {r['label']}  ·  {r['family']}",
                     x=0.012, ha="left", fontsize=8)
        fig.tight_layout(rect=(0, 0, 1, 0.955))
        p = out_dir / f"pair_side_by_side_{r['pair_id']}_{_tag(args)}.png"
        fig.savefig(p, dpi=200)
        plt.close(fig)
        print(f"wrote {p}")

    # ---- all pairs on one axis, apo/holo as a box pair per row ----------
    order = sorted(pairs, key=lambda r: np.median(best[(r["holo_pdb"], r["holo_chain"])]["dmin"]))
    fig, ax = plt.subplots(figsize=(7.0, 0.42 * len(order) + 1.1))
    for i, r in enumerate(order):
        for off, key, color in ((+0.19, ("apo_pdb", "apo_chain"), APO),
                                (-0.19, ("holo_pdb", "holo_chain"), HOLO)):
            d = best[(r[key[0]], r[key[1]])]["dmin"]
            bp = ax.boxplot(d, positions=[i + off], vert=False, widths=0.32,
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="white", lw=1.0))
            bp["boxes"][0].set(facecolor=color, edgecolor="none")
            for part in ("whiskers", "caps"):
                for ln in bp[part]:
                    ln.set(color=GREY, lw=0.7)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{r['pair_id']}  {r['label'][:24]}" for r in order], fontsize=6)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel(f"RMSD to the closest state of its own bundle, {args.atoms} ({chr(197)})",
                  fontsize=8)
    ax.set_xlim(left=0)
    handles = [plt.Line2D([], [], marker="s", ls="", ms=5, color=APO, label="apo"),
               plt.Line2D([], [], marker="s", ls="", ms=5, color=HOLO, label="holo")]
    ax.legend(handles=handles, fontsize=6, frameon=False, loc="lower right")
    ax.set_title("Each BioEmu ensemble against its own NMR bundle, apo vs holo",
                 loc="left", fontsize=8)
    fig.tight_layout()
    p = out_dir / f"pair_side_by_side_all_{_tag(args)}.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
