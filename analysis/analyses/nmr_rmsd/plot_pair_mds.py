#!/usr/bin/env python
"""
One map per apo/holo pair: where do the BioEmu ensembles sit relative to the two
NMR bundles?

Scalar RMSDs cannot tell "covers the bound basin" from "sits between the two and
is equidistant from both". So pool everything -- apo states, holo states,
BioEmu(apo seq), BioEmu(holo seq) -- on the pair's shared core, build the full
pairwise RMSD matrix, and embed it.

Two views, both from the same matrix:

  left   classical MDS (Torgerson) of the pooled matrix. RMSD after optimal
         superposition is not exactly Euclidean, so the panel reports the share
         of positive eigenvalue mass the 2 axes capture -- read the geometry
         only as far as that number allows.

  right  the apo->holo reaction coordinate, from distances alone, against the
         displacement off that axis. With a = RMSD to the apo medoid, h = RMSD
         to the holo medoid and L their separation,

             t = (a^2 - h^2 + L^2) / (2 L^2)

         is the projection onto the apo->holo axis: t=0 at the apo medoid, t=1
         at the holo medoid. Medoids (the bundle's most central deposited state)
         are used rather than abstract centroids so both anchors are real
         structures.

         t alone is not enough: t~0.5 is also what a structure displaced
         PERPENDICULAR to the axis returns, however far away it is. So the panel
         plots t against the off-axis distance

             p = sqrt(a^2 - (t*L)^2)

         p/L << 1 means the structure genuinely lies on the apo->holo line and t
         can be read as "how far along"; p/L ~ 1 means the structure is somewhere
         else entirely and t is an artefact of being equidistant.

Outputs (in --output-dir):
  plots/pair_mds_<pair>_<atoms>.png   per-pair panel
  plots/pair_mds_all_<atoms>.png      all pairs stacked
  pair_projection_<atoms>.csv         per structure: group, a, h, t, p, p/L

Usage:
    python plot_pair_mds.py
    python plot_pair_mds.py --pair-ids vnd --atoms ca
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
from pair_core import load_pairs, load_pair, rmsd_matrix    # noqa: E402

_TFCONF = _HERE.parents[2]
sys.path.insert(0, str(_TFCONF))
from palette import GREY_R, TEAL_R, ALARM, apply_style      # noqa: E402

# grey = deposited reference structures, teal = BioEmu (palette.py SS 4.1);
# light/dark within each family = apo/holo.
GROUPS = (
    ("apo_ref",  "NMR apo",          GREY_R[0], "o", 26),
    ("holo_ref", "NMR holo",         GREY_R[2], "o", 26),
    ("apo_ens",  "BioEmu (apo seq)",  TEAL_R[0], ".", 9),
    ("holo_ens", "BioEmu (holo seq)", TEAL_R[2], ".", 9),
)


def pooled_matrix(P):
    """Full pairwise RMSD matrix over all four members, plus group slices."""
    traj, idx = P["traj"], P["atom_idx"]
    names = [g[0] for g in GROUPS]
    sizes = [traj[n].n_frames for n in names]
    offs = np.cumsum([0] + sizes)
    n = offs[-1]
    D = np.zeros((n, n))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j < i:
                continue
            block = rmsd_matrix(traj[b], traj[a], idx[b], idx[a])  # (frames of a, frames of b)
            D[offs[i]:offs[i + 1], offs[j]:offs[j + 1]] = block
            if i != j:
                D[offs[j]:offs[j + 1], offs[i]:offs[i + 1]] = block.T
    np.fill_diagonal(D, 0.0)
    D = 0.5 * (D + D.T)                       # symmetrise the tiny asymmetry
    slices = {n_: slice(offs[k], offs[k + 1]) for k, n_ in enumerate(names)}
    return D, slices


def classical_mds(D, n_dim=2):
    """Torgerson MDS. Returns (coords, fraction of positive eigen-mass kept)."""
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    pos = vals[vals > 0]
    keep = vals[:n_dim].clip(min=0)
    coords = vecs[:, :n_dim] * np.sqrt(keep)
    return coords, float(keep.sum() / pos.sum()) if pos.size else float("nan")


def medoid(D, sl):
    """Index (global) of the most central member of a group."""
    sub = D[sl, sl]
    return int(np.arange(D.shape[0])[sl][sub.sum(axis=1).argmin()])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default=_HERE / "pairs.csv")
    ap.add_argument("--pair-ids", nargs="*", default=None)
    ap.add_argument("--atoms", choices=("ca", "backbone"), default="ca")
    ap.add_argument("--spread-cutoff", type=float, default=3.0)
    ap.add_argument("--output-dir", default=_HERE)
    ap.add_argument("--no-fetch", action="store_true")
    args = ap.parse_args()

    pairs = load_pairs(args.pairs, args.pair_ids)
    out_dir = Path(args.output_dir)
    (out_dir / "plots").mkdir(parents=True, exist_ok=True)
    apply_style()
    proj_rows, panels = [], []

    for pair in pairs:
        pid = pair["pair_id"]
        print(f"{pid} ({pair['label']})")
        try:
            P = load_pair(pair, atoms=args.atoms, spread_cutoff=args.spread_cutoff,
                          fetch=not args.no_fetch, verbose=False)
        except Exception as exc:                        # noqa: BLE001 - report, continue
            print(f"  SKIP: {exc}", file=sys.stderr)
            continue
        D, sl = pooled_matrix(P)
        coords, kept = classical_mds(D)

        m_apo, m_holo = medoid(D, sl["apo_ref"]), medoid(D, sl["holo_ref"])
        L = D[m_apo, m_holo]
        a_all, h_all = D[m_apo], D[m_holo]
        t_all = (a_all ** 2 - h_all ** 2 + L ** 2) / (2 * L ** 2)
        # distance from the apo->holo line: what is left of `a` once the
        # on-axis component t*L is removed.
        p_all = np.sqrt(np.clip(a_all ** 2 - (t_all * L) ** 2, 0.0, None))

        for name, label, _, _, _ in GROUPS:
            for k, gi in enumerate(range(sl[name].start, sl[name].stop)):
                proj_rows.append({
                    "pair_id": pid, "group": name, "index": k,
                    "d_apo_medoid_A": f"{a_all[gi]:.4f}",
                    "d_holo_medoid_A": f"{h_all[gi]:.4f}",
                    "t": f"{t_all[gi]:.4f}",
                    "p_offaxis_A": f"{p_all[gi]:.4f}",
                    "p_over_L": f"{p_all[gi] / L:.4f}",
                })

        stats = {}
        for name, _, _, _, _ in GROUPS:
            t, pv = t_all[sl[name]], p_all[sl[name]]
            stats[name] = (float(np.median(t)), float(np.percentile(t, 25)),
                           float(np.percentile(t, 75)), float(np.median(pv)),
                           float(np.median(pv)) / L)
        print(f"  L(apo medoid -> holo medoid) = {L:.2f} A, MDS keeps {kept:.0%} "
              f"of positive eigen-mass")
        for name, _, _, _, _ in GROUPS:
            med, q1, q3, pmed, pol = stats[name]
            print(f"  {name:9s} t {med:+.2f} (IQR {q1:+.2f}..{q3:+.2f})  "
                  f"off-axis p {pmed:.2f} A  p/L {pol:.2f}")
        panels.append((pair, P, D, sl, coords, kept, t_all, p_all, L, stats))

    # ---- draw -----------------------------------------------------------
    def draw(axL, axR, panel, label_axes=True):
        pair, P, D, sl, coords, kept, t_all, p_all, L, stats = panel
        for name, label, color, marker, size in GROUPS:
            c = coords[sl[name]]
            axL.scatter(c[:, 0], c[:, 1], s=size, c=color, marker=marker,
                        linewidths=0, alpha=0.85, label=label, zorder=3 if "ref" in name else 2)
        axL.set_title(f"{pair['label']}", loc="left", fontsize=7)
        axL.annotate(f"{P['n_core']} core pos · 2D keeps {kept:.0%}",
                     xy=(0, 1.01), xycoords="axes fraction", fontsize=5.5, color="#666666")
        # A handful of far-flung frames otherwise squash the whole cloud into a
        # dot; clip the view to the bulk and say how many points fell outside.
        q1, q3 = np.percentile(coords, [25, 75], axis=0)
        iqr = np.maximum(q3 - q1, 1e-6)
        lo = np.maximum(q1 - 3 * iqr, coords.min(axis=0))
        hi = np.minimum(q3 + 3 * iqr, coords.max(axis=0))
        pad = 0.12 * np.maximum(hi - lo, 1e-6)
        lo, hi = lo - pad, hi + pad
        outside = int(((coords < lo) | (coords > hi)).any(axis=1).sum())
        axL.set_xlim(lo[0], hi[0])
        axL.set_ylim(lo[1], hi[1])
        axL.set_aspect("equal", adjustable="box")
        if outside:
            axL.annotate(f"{outside} points outside view", xy=(0.98, 0.02),
                         xycoords="axes fraction", fontsize=5, color="#666666", ha="right")
        if label_axes:
            axL.set_xlabel("MDS 1 (A)", fontsize=7)
            axL.set_ylabel("MDS 2 (A)", fontsize=7)

        for name, label, color, marker, size in GROUPS:
            axR.scatter(t_all[sl[name]], p_all[sl[name]] / L, s=size, c=color,
                        marker=marker, linewidths=0, alpha=0.85, label=label,
                        zorder=3 if "ref" in name else 2)
        for x, lab in ((0.0, "apo"), (1.0, "holo")):
            axR.axvline(x, color=ALARM, ls="--", lw=0.8)
            axR.annotate(lab, xy=(x, 1.005), xycoords=("data", "axes fraction"),
                         fontsize=5.5, color=ALARM, ha="center")
        axR.axhline(1.0, color=ALARM, ls=":", lw=0.8)
        axR.annotate("off-axis by a whole apo->holo step",
                     xy=(0.02, 1.02), xycoords=("axes fraction", "data"),
                     fontsize=5, color=ALARM, va="bottom")
        axR.set_xlim(-0.6, 1.6)
        # a few wild frames otherwise flatten every group onto y=0
        pl = p_all / L
        q1r, q3r = np.percentile(pl, [25, 75])
        ytop = min(float(pl.max()), float(q3r + 4 * (q3r - q1r))) * 1.12 + 1e-6
        n_out = int((pl > ytop).sum())
        axR.set_ylim(0, max(ytop, 1.25))
        if n_out:
            axR.annotate(f"{n_out} points above view", xy=(0.98, 0.02),
                         xycoords="axes fraction", fontsize=5, color="#666666", ha="right")
        axR.annotate(f"L = {L:.2f} A", xy=(0.02, 0.95), xycoords="axes fraction",
                     fontsize=5.5, color="#666666", ha="left")
        if label_axes:
            axR.set_xlabel("t  (apo -> holo projection)", fontsize=7)
            axR.set_ylabel("off-axis distance  p / L", fontsize=7)

    for panel in panels:
        pair = panel[0]
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(6.6, 2.9))
        draw(axL, axR, panel)
        axR.legend(fontsize=5.5, frameon=False, loc="upper right")
        fig.tight_layout()
        path = out_dir / "plots" / f"pair_mds_{pair['pair_id']}_{args.atoms}.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        print(f"wrote {path}")

    fig, axes = plt.subplots(len(panels), 2, figsize=(6.6, 2.7 * len(panels)), squeeze=False)
    for row, panel in enumerate(panels):
        draw(axes[row][0], axes[row][1], panel, label_axes=(row == len(panels) - 1))
    axes[0][1].legend(fontsize=5.5, frameon=False, loc="upper right")
    fig.tight_layout()
    path = out_dir / "plots" / f"pair_mds_all_{args.atoms}.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"wrote {path}")

    csv_path = out_dir / f"pair_projection_{args.atoms}.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(proj_rows[0]))
        w.writeheader()
        w.writerows(proj_rows)
    print(f"wrote {csv_path}  ({len(proj_rows)} rows)")


if __name__ == "__main__":
    main()
