#!/usr/bin/env python
"""
Trim-sensitivity sweep: is the apo/holo verdict an artefact of where we cut?

compute_pair_rmsd.py trims disordered termini at a fixed --spread-cutoff (3.0 A).
That choice is a judgement call, and vnd -- the one pair whose BioEmu frames sit
closer to the BOUND bundle -- survives at only 36 core positions, so its sign may
be a property of the cut rather than of the ensemble.

This re-runs the pair analysis across a range of cutoffs and reports how the
core size, the (1) separation gate and the (3) delta move with it. A result that
is real should be flat in the middle of the range; one that flips with the cut is
a trimming artefact. The largest cutoff is the no-trim control.

Outputs (in --output-dir):
  trim_sweep_<atoms>.csv        one row per (pair, cutoff)
  plots/trim_sweep_<atoms>.png  core size / separation ratio / delta vs cutoff

Usage:
    python sweep_trim.py
    python sweep_trim.py --pair-ids vnd nhp6a rok --cutoffs 2 2.5 3 3.5 4 5
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from pair_core import load_pairs, load_pair, rmsd_matrix   # noqa: E402

_TFCONF = _HERE.parents[2]
sys.path.insert(0, str(_TFCONF))
from palette import TEAL, GREY, ALARM, apply_style          # noqa: E402

NO_TRIM = 1e3          # cutoff large enough that nothing is trimmed


def upper(m):
    i, j = np.triu_indices(m.shape[0], k=1)
    return m[i, j]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default=_HERE / "pairs.csv")
    ap.add_argument("--pair-ids", nargs="*", default=None)
    ap.add_argument("--atoms", choices=("ca", "backbone"), default="ca")
    ap.add_argument("--cutoffs", type=float, nargs="*",
                    default=(2.0, 2.5, 3.0, 3.5, 4.0, 5.0, NO_TRIM))
    ap.add_argument("--output-dir", default=_HERE)
    ap.add_argument("--no-fetch", action="store_true")
    args = ap.parse_args()

    pairs = load_pairs(args.pairs, args.pair_ids)
    out_dir = Path(args.output_dir)
    (out_dir / "plots").mkdir(parents=True, exist_ok=True)
    rows = []

    for pair in pairs:
        pid = pair["pair_id"]
        print(f"{pid} ({pair['label']})")
        for cutoff in args.cutoffs:
            try:
                P = load_pair(pair, atoms=args.atoms, spread_cutoff=cutoff,
                              fetch=not args.no_fetch, verbose=False)
            except Exception as exc:                    # noqa: BLE001
                print(f"  cutoff {cutoff:g}: SKIP ({exc})", file=sys.stderr)
                continue
            traj, idx = P["traj"], P["atom_idx"]

            apo_apo = rmsd_matrix(traj["apo_ref"], traj["apo_ref"],
                                  idx["apo_ref"], idx["apo_ref"])
            holo_holo = rmsd_matrix(traj["holo_ref"], traj["holo_ref"],
                                    idx["holo_ref"], idx["holo_ref"])
            apo_holo = rmsd_matrix(traj["holo_ref"], traj["apo_ref"],
                                   idx["holo_ref"], idx["apo_ref"])
            within = np.concatenate([upper(apo_apo), upper(holo_holo)])
            between = apo_holo.ravel()
            ratio = float(np.median(between) / np.median(within))

            row = {
                "pair_id": pid, "label": pair["label"],
                "spread_cutoff_A": "none" if cutoff >= NO_TRIM else f"{cutoff:g}",
                "n_intersect": P["n_intersect"], "n_core": P["n_core"],
                "n_atoms": P["n_atoms"],
                "median_within_A": f"{np.median(within):.3f}",
                "median_between_A": f"{np.median(between):.3f}",
                "separation_A": f"{np.median(between) - np.median(within):.3f}",
                "ratio_between_within": f"{ratio:.3f}",
            }
            for ens_role in ("apo", "holo"):
                ens = traj[f"{ens_role}_ens"]
                d = {}
                for bundle in ("apo", "holo"):
                    m = rmsd_matrix(ens, traj[f"{bundle}_ref"],
                                    idx[f"{ens_role}_ens"], idx[f"{bundle}_ref"])
                    d[bundle] = m.min(axis=0)
                delta = d["holo"] - d["apo"]
                row[f"{ens_role}_seq_dmin_apo_A"] = f"{np.median(d['apo']):.3f}"
                row[f"{ens_role}_seq_dmin_holo_A"] = f"{np.median(d['holo']):.3f}"
                row[f"{ens_role}_seq_delta_A"] = f"{np.median(delta):.3f}"
                row[f"{ens_role}_seq_frac_closer_apo"] = f"{(delta > 0).mean():.4f}"
            rows.append(row)
            print(f"  cutoff {row['spread_cutoff_A']:>4s}: core {P['n_core']:3d} pos, "
                  f"ratio {ratio:.2f}, delta apo-seq {row['apo_seq_delta_A']:>6s} A, "
                  f"holo-seq {row['holo_seq_delta_A']:>6s} A")

    if not rows:
        sys.exit("nothing computed")
    csv_path = out_dir / f"trim_sweep_{args.atoms}.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {csv_path}  ({len(rows)} rows)")

    # ---- plot -----------------------------------------------------------
    apply_style()
    pids = list(dict.fromkeys(r["pair_id"] for r in rows))
    fig, axes = plt.subplots(3, len(pids), figsize=(2.2 * len(pids), 5.6),
                             sharex=True, squeeze=False)

    def xs(sub):
        return [len(sub) - 1 if r["spread_cutoff_A"] == "none" else i
                for i, r in enumerate(sub)]

    for col, pid in enumerate(pids):
        sub = [r for r in rows if r["pair_id"] == pid]
        x = np.arange(len(sub))
        labels = [r["spread_cutoff_A"] for r in sub]

        ax = axes[0][col]
        ax.plot(x, [r["n_core"] for r in sub], "o-", color=GREY, ms=3)
        ax.set_title(sub[0]["label"], loc="left", fontsize=7)
        if col == 0:
            ax.set_ylabel("core positions", fontsize=7)

        ax = axes[1][col]
        ax.plot(x, [float(r["ratio_between_within"]) for r in sub], "o-", color=GREY, ms=3)
        ax.axhline(1.25, color=ALARM, ls="--", lw=0.8)
        if col == 0:
            ax.set_ylabel("between / within", fontsize=7)

        ax = axes[2][col]
        for role, ls in (("apo", "-"), ("holo", "--")):
            ax.plot(x, [float(r[f"{role}_seq_delta_A"]) for r in sub], "o" + ls,
                    color=TEAL, ms=3, label=f"BioEmu({role} seq)")
        ax.axhline(0, color=ALARM, ls="--", lw=0.8)
        if col == 0:
            ax.set_ylabel("median delta (A)\n>0 = closer to apo", fontsize=7)
            ax.legend(fontsize=5, frameon=False)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=6)
        ax.set_xlabel("trim cutoff (A)", fontsize=7)

    fig.tight_layout()
    path = out_dir / "plots" / f"trim_sweep_{args.atoms}.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
