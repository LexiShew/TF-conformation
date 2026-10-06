#!/usr/bin/env python
"""
How often is each NMR state the closest one to a BioEmu conformation?

Reads nmr_state_rmsd_best_<atoms>.csv (from compute_nmr_state_rmsd.py) and
draws, per protein, a bar chart of

    x = NMR state (1..n_states)
    y = number of BioEmu conformations whose best_state is that state

A dashed grey line marks the uniform expectation (n_conformations / n_states):
bars above it are states the ensemble preferentially lands near.  States never
chosen still get a tick, so gaps in the ensemble's coverage are visible.

Outputs (in --output-dir, default ./plots):
  <pdb>_chain<X>_best_state_<atoms>.png    one per protein
  all_best_state_<atoms>.png               small-multiple grid of the same panels
  best_state_counts_<atoms>.csv            the counts behind the plots

Usage:
    python plot_best_state.py                  # CA-based best states
    python plot_best_state.py --atoms backbone
"""
import argparse
import csv
import sys
from collections import Counter, OrderedDict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
_TFCONF = _HERE.parents[2]
sys.path.insert(0, str(_TFCONF))
from palette import TEAL, GREY, apply_style   # project-canonical colors


def load_best(path):
    """OrderedDict {(pdb_id, chain): {'counts': Counter, 'n_states', 'n_conf'}}."""
    entries = OrderedDict()
    with open(path) as fh:
        for row in csv.DictReader(fh):
            key = (row["pdb_id"], row["chain"])
            e = entries.setdefault(key, {"counts": Counter(), "n_states": 0, "n_conf": 0})
            e["counts"][int(row["best_state"])] += 1
            e["n_states"] = int(row["n_states"])
            e["n_conf"] += 1
    return entries


def draw_panel(ax, pdb_id, chain, entry, label_axes=True):
    n_states, n_conf = entry["n_states"], entry["n_conf"]
    states = np.arange(1, n_states + 1)
    counts = np.array([entry["counts"].get(int(s), 0) for s in states])
    expected = n_conf / n_states
    used = int((counts > 0).sum())

    ax.bar(states, counts, width=0.78, color=TEAL, linewidth=0)
    ax.axhline(expected, color=GREY, linestyle="--", linewidth=0.8, zorder=3)
    ax.annotate(f"uniform = {expected:.1f}", xy=(1.01, expected),
                xycoords=("axes fraction", "data"), va="center", ha="left",
                fontsize=6, color=GREY, annotation_clip=False)

    ax.set_title(f"{pdb_id} chain {chain}", loc="left", pad=13)
    ax.annotate(f"{n_conf} conformations · {used}/{n_states} states used · "
                f"top state {counts.max() / n_conf:.0%}",
                xy=(0, 1.012), xycoords="axes fraction", fontsize=6,
                color="#666666", va="bottom")
    if label_axes:
        ax.set_xlabel("NMR state")
        ax.set_ylabel("conformations closest to this state")
    # Thin the ticks for the 50-state entries; always show first and last.
    step = 1 if n_states <= 20 else (2 if n_states <= 30 else 5)
    ticks = list(range(1, n_states + 1, step))
    if ticks[-1] != n_states:
        ticks.append(n_states)
    ax.set_xticks(ticks)
    ax.set_xlim(0.4, n_states + 0.6)
    ax.set_ylim(0, max(counts.max(), expected) * 1.18)
    ax.margins(x=0)
    return counts


def _tag(args):
    """<atoms> or <atoms>_<label> - keeps trimmed and untrimmed outputs apart."""
    return f"{args.atoms}{('_' + args.label) if args.label else ''}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atoms", choices=("ca", "backbone"), default="ca")
    ap.add_argument("--label", default="", help="suffix of the CSVs to read / figures to write")
    ap.add_argument("--best-csv", default=None,
                    help="default: ./nmr_state_rmsd_best_<atoms>.csv")
    ap.add_argument("--output-dir", default=_HERE / "plots")
    args = ap.parse_args()

    best_csv = Path(args.best_csv) if args.best_csv else _HERE / f"nmr_state_rmsd_best_{_tag(args)}.csv"
    if not best_csv.is_file():
        sys.exit(f"missing {best_csv} — run compute_nmr_state_rmsd.py first")
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    apply_style()
    entries = load_best(best_csv)
    rows = []

    for (pdb_id, chain), entry in entries.items():
        fig, ax = plt.subplots(figsize=(max(3.4, entry["n_states"] * 0.16), 2.3))
        counts = draw_panel(ax, pdb_id, chain, entry)
        path = out_dir / f"{pdb_id}_chain{chain}_best_state_{_tag(args)}.png"
        fig.savefig(path)
        plt.close(fig)
        print(f"wrote {path}")
        for state, count in enumerate(counts, start=1):
            rows.append({"pdb_id": pdb_id, "chain": chain, "state": state,
                         "n_best": int(count), "n_conformations": entry["n_conf"],
                         "n_states": entry["n_states"],
                         "frac_best": f"{count / entry['n_conf']:.4f}"})

    # 35 ensembles at two columns is an unreadable 18-row strip; scale the grid.
    ncols = 2 if len(entries) <= 6 else (3 if len(entries) <= 12 else 4)
    nrows = -(-len(entries) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 1.85 * nrows),
                             squeeze=False)
    for ax, ((pdb_id, chain), entry) in zip(axes.ravel(), entries.items()):
        draw_panel(ax, pdb_id, chain, entry, label_axes=False)
    for ax in axes.ravel()[len(entries):]:
        ax.set_visible(False)
    fig.supxlabel("NMR state", fontsize=8)
    fig.supylabel("conformations closest to this state", fontsize=8)
    fig.tight_layout(rect=(0.015, 0.015, 1, 1))
    grid_path = out_dir / f"all_best_state_{_tag(args)}.png"
    fig.savefig(grid_path)
    plt.close(fig)
    print(f"wrote {grid_path}")

    counts_csv = out_dir / f"best_state_counts_{_tag(args)}.csv"
    with open(counts_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {counts_csv}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
