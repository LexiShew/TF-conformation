#!/usr/bin/env python
"""
Cartesian PCA fitted on the DEPOSITED bundles; BioEmu projected into it.

Replaces the classical-MDS view. MDS embeds an RMSD matrix, and RMSD after
optimal superposition is not a Euclidean metric, so the 2D picture kept only
35-67% of the positive eigenvalue mass and had to be read as a sketch. Here the
space is built from superposed Cartesian coordinates instead: it is metric, the
variance explained is an honest number, and the axes have a defined meaning.

The subspace is fitted on the apo + holo NMR states ONLY. Those are the
conformations nature (and the depositors) actually put on the table, so they
define the coordinate system; the BioEmu ensembles are projected into it as test
data and never influence the axes. PC1 comes out as the apo->holo direction
whenever the pair is informative.

Two numbers per structure, both in Angstrom and directly comparable to an RMSD:

  PC1, PC2   position inside the deposited subspace -- where along apo->holo
  residual   the part of the structure the subspace CANNOT reproduce, i.e.
             || x - mean - P_k(x - mean) || / sqrt(n_atoms), with k the number
             of PCs holding --var-target of the bundle variance

A small residual means "this structure is a combination of motions the deposited
states already show". A residual on the scale of the apo-holo separation itself
means the structure is somewhere the deposited states never go, and its PC1
coordinate is not worth interpreting.

Bundle residuals are computed LEAVE-ONE-OUT (refit the subspace without that
state), otherwise the states that defined the space would trivially score zero
and the comparison against BioEmu would be rigged.

Outputs (in --output-dir):
  plots/pair_pca_<pair>_<atoms>.png   per-pair panel
  plots/pair_pca_all_<atoms>.png      all pairs stacked
  pair_pca_<atoms>.csv                per structure: PC1, PC2, residual
  pair_pca_summary_<atoms>.csv        per (pair, group): medians + apo/holo scale

Usage:
    python plot_pair_pca.py
    python plot_pair_pca.py --pair-ids vnd nhp6a rok
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
from pair_core import load_pairs, load_pair                 # noqa: E402

_TFCONF = _HERE.parents[2]
sys.path.insert(0, str(_TFCONF))
from palette import GREY_R, TEAL_R, ALARM, apply_style      # noqa: E402

# grey = deposited reference structures, teal = BioEmu (palette.py SS 4.1);
# light/dark within each family = apo/holo.
GROUPS = (
    ("apo_ref",  "NMR apo",           GREY_R[0], "o", 26),
    ("holo_ref", "NMR holo",          GREY_R[2], "o", 26),
    ("apo_ens",  "BioEmu (apo seq)",  TEAL_R[0], ".", 9),
    ("holo_ens", "BioEmu (holo seq)", TEAL_R[2], ".", 9),
)


def superposed_coords(P, ref_role="holo_ref", ref_frame=0):
    """{group -> (n_struct, 3*n_core_atoms)} superposed on one common reference."""
    traj, idx = P["traj"], P["atom_idx"]
    ref = traj[ref_role]
    out = {}
    for name, _, _, _, _ in GROUPS:
        t = traj[name].slice(range(traj[name].n_frames), copy=True)
        t.superpose(ref, frame=ref_frame,
                    atom_indices=idx[name], ref_atom_indices=idx[ref_role])
        out[name] = (t.xyz[:, idx[name]] * 10.0).reshape(t.n_frames, -1)
    return out


def fit_pca(X):
    """Mean + right singular vectors + explained-variance fractions."""
    mean = X.mean(axis=0)
    U, S, Vt = np.linalg.svd(X - mean, full_matrices=False)
    var = S ** 2
    return mean, Vt, var / var.sum() if var.sum() > 0 else var


def n_components(frac, target):
    """Smallest k whose cumulative explained variance reaches target."""
    return int(np.searchsorted(np.cumsum(frac), target) + 1)


def residual(X, mean, Vt, k, n_atoms):
    """Per-structure out-of-subspace deviation, in Angstrom per atom."""
    C = X - mean
    recon = (C @ Vt[:k].T) @ Vt[:k]
    return np.linalg.norm(C - recon, axis=1) / np.sqrt(n_atoms)


def loo_residual(X, k, n_atoms):
    """Leave-one-out residual: each row scored by a subspace fitted without it."""
    out = np.empty(len(X))
    for i in range(len(X)):
        keep = np.delete(np.arange(len(X)), i)
        mean, Vt, _ = fit_pca(X[keep])
        kk = min(k, Vt.shape[0])
        out[i] = residual(X[i:i + 1], mean, Vt, kk, n_atoms)[0]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default=_HERE / "pairs.csv")
    ap.add_argument("--pair-ids", nargs="*", default=None)
    ap.add_argument("--atoms", choices=("ca", "backbone"), default="ca")
    ap.add_argument("--spread-cutoff", type=float, default=3.0)
    ap.add_argument("--var-target", type=float, default=0.90,
                    help="bundle variance the subspace must hold (default: %(default)s)")
    ap.add_argument("--output-dir", default=_HERE)
    ap.add_argument("--no-fetch", action="store_true")
    args = ap.parse_args()

    pairs = load_pairs(args.pairs, args.pair_ids)
    out_dir = Path(args.output_dir)
    (out_dir / "plots").mkdir(parents=True, exist_ok=True)
    apply_style()
    rows, summary, panels = [], [], []

    for pair in pairs:
        pid = pair["pair_id"]
        print(f"{pid} ({pair['label']})")
        try:
            P = load_pair(pair, atoms=args.atoms, spread_cutoff=args.spread_cutoff,
                          fetch=not args.no_fetch, verbose=False)
        except Exception as exc:                        # noqa: BLE001 - report, continue
            print(f"  SKIP: {exc}", file=sys.stderr)
            continue
        n_atoms = P["n_atoms"]
        X = superposed_coords(P)

        bundle = np.vstack([X["apo_ref"], X["holo_ref"]])
        n_apo = len(X["apo_ref"])
        mean, Vt, frac = fit_pca(bundle)
        k = n_components(frac, args.var_target)
        print(f"  subspace from {len(bundle)} deposited states: "
              f"PC1 {frac[0]:.0%}, PC2 {frac[1]:.0%}, k={k} for {args.var_target:.0%}")

        # apo/holo separation along PC1 -- the scale everything else is read against
        pc_apo = (X["apo_ref"] - mean) @ Vt[0]
        pc_holo = (X["holo_ref"] - mean) @ Vt[0]
        sep_pc1 = float(abs(np.median(pc_holo) - np.median(pc_apo)))

        scores, resid = {}, {}
        for name, _, _, _, _ in GROUPS:
            scores[name] = (X[name] - mean) @ Vt[:2].T
            if name in ("apo_ref", "holo_ref"):
                sl = slice(0, n_apo) if name == "apo_ref" else slice(n_apo, len(bundle))
                resid[name] = loo_residual(bundle, k, n_atoms)[sl]
            else:
                resid[name] = residual(X[name], mean, Vt, k, n_atoms)

        for name, label, _, _, _ in GROUPS:
            for i in range(len(scores[name])):
                rows.append({"pair_id": pid, "group": name, "index": i,
                             "pc1_A": f"{scores[name][i, 0]:.4f}",
                             "pc2_A": f"{scores[name][i, 1]:.4f}",
                             "residual_A": f"{resid[name][i]:.4f}"})
            r = resid[name]
            summary.append({
                "pair_id": pid, "label": pair["label"], "group": name,
                "n": len(r), "k_components": k,
                "pc1_var_frac": f"{frac[0]:.4f}", "pc2_var_frac": f"{frac[1]:.4f}",
                "median_pc1_A": f"{np.median(scores[name][:, 0]):.3f}",
                "median_residual_A": f"{np.median(r):.3f}",
                "q3_residual_A": f"{np.percentile(r, 75):.3f}",
                "apo_holo_sep_pc1_A": f"{sep_pc1:.3f}",
                "residual_over_sep": f"{np.median(r) / sep_pc1:.3f}" if sep_pc1 else "",
            })
            print(f"  {name:9s} PC1 {np.median(scores[name][:, 0]):+6.2f} A   "
                  f"residual {np.median(r):.2f} A   "
                  f"residual/sep {np.median(r) / sep_pc1:.2f}" if sep_pc1 else "")
        panels.append((pair, P, scores, resid, frac, k, sep_pc1,
                       float(np.median(pc_apo)), float(np.median(pc_holo))))

    # ---- draw -----------------------------------------------------------
    def draw(axL, axR, panel, label_axes=True):
        pair, P, scores, resid, frac, k, sep, m_apo, m_holo = panel
        for name, label, color, marker, size in GROUPS:
            s = scores[name]
            axL.scatter(s[:, 0], s[:, 1], s=size, c=color, marker=marker,
                        linewidths=0, alpha=0.85, label=label,
                        zorder=3 if "ref" in name else 2)
            axR.scatter(s[:, 0], resid[name], s=size, c=color, marker=marker,
                        linewidths=0, alpha=0.85, zorder=3 if "ref" in name else 2)
        # A few wild frames otherwise squash every cluster into a dot. Clip to the
        # bulk of the data plus both bundles, and say how many fell outside.
        def clip(vals, anchors=()):
            q1, q3 = np.percentile(vals, [25, 75])
            iqr = max(q3 - q1, 1e-6)
            lo = max(q1 - 3 * iqr, float(np.min(vals)))
            hi = min(q3 + 3 * iqr, float(np.max(vals)))
            for a in anchors:                      # never clip away a bundle
                lo, hi = min(lo, a), max(hi, a)
            pad = 0.12 * max(hi - lo, 1e-6)
            return lo - pad, hi + pad

        all_pc1 = np.concatenate([scores[n][:, 0] for n, *_ in GROUPS])
        all_pc2 = np.concatenate([scores[n][:, 1] for n, *_ in GROUPS])
        all_res = np.concatenate([resid[n] for n, *_ in GROUPS])
        xlim = clip(all_pc1, (m_apo, m_holo))
        axL.set_xlim(*xlim)
        axL.set_ylim(*clip(all_pc2))
        axR.set_xlim(*xlim)
        axR.set_ylim(0, max(clip(all_res)[1], sep * 1.15))
        n_out = int(((all_pc1 < xlim[0]) | (all_pc1 > xlim[1])).sum())

        for ax in (axL, axR):
            for x, lab in ((m_apo, "apo"), (m_holo, "holo")):
                ax.axvline(x, color=ALARM, ls="--", lw=0.8)
                ax.annotate(lab, xy=(x, 0.015), xycoords=("data", "axes fraction"),
                            fontsize=5.5, color=ALARM, ha="center", va="bottom")
            if n_out:
                ax.annotate(f"{n_out} points outside view", xy=(0.99, 0.985),
                            xycoords="axes fraction", fontsize=5, color="#666666",
                            ha="right", va="top")
        axL.set_title(f"{pair['label']}  ·  {P['n_core']} core pos  ·  "
                      f"subspace from deposited states only", loc="left", fontsize=7)
        axR.axhline(sep, color=ALARM, ls=":", lw=0.8)
        axR.annotate("unexplained by a whole apo->holo step",
                     xy=(0.02, sep), xycoords=("axes fraction", "data"),
                     fontsize=5, color=ALARM, va="bottom")
        if label_axes:
            axL.set_xlabel(f"PC1 ({frac[0]:.0%} of deposited variance), A", fontsize=7)
            axL.set_ylabel(f"PC2 ({frac[1]:.0%}), A", fontsize=7)
            axR.set_xlabel("PC1, A", fontsize=7)
            axR.set_ylabel(f"residual outside the deposited\nsubspace (k={k}), A", fontsize=7)

    for panel in panels:
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(6.8, 2.9))
        draw(axL, axR, panel)
        axL.legend(fontsize=5.5, frameon=False, loc="best")
        fig.tight_layout()
        path = out_dir / "plots" / f"pair_pca_{panel[0]['pair_id']}_{args.atoms}.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        print(f"wrote {path}")

    fig, axes = plt.subplots(len(panels), 2, figsize=(6.8, 2.7 * len(panels)), squeeze=False)
    for row, panel in enumerate(panels):
        draw(axes[row][0], axes[row][1], panel, label_axes=(row == len(panels) - 1))
    axes[0][0].legend(fontsize=5.5, frameon=False, loc="best")
    fig.tight_layout()
    path = out_dir / "plots" / f"pair_pca_all_{args.atoms}.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"wrote {path}")

    for name, data in ((f"pair_pca_{args.atoms}.csv", rows),
                       (f"pair_pca_summary_{args.atoms}.csv", summary)):
        with open(out_dir / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)
        print(f"wrote {out_dir / name}  ({len(data)} rows)")


if __name__ == "__main__":
    main()
