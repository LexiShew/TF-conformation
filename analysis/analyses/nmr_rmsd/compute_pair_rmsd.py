#!/usr/bin/env python
"""
Apo/holo NMR pair analysis on a single shared comparison core.

Every distance below is computed on the SAME atom set for all four members of a
pair (apo bundle, holo bundle, BioEmu(apo seq), BioEmu(holo seq)); see
pair_core.py for how that core is built and trimmed.

  (1) Is the pair informative?
      Within-bundle (apo-apo, holo-holo) vs between-bundle (apo-holo) RMSD.
      If the apo/holo separation is not clearly above NMR bundle heterogeneity,
      the pair cannot support any claim about BioEmu and is gated out.

  (2) The 2x2.
      Both BioEmu ensembles scored against both bundles. The off-diagonal cells
      ask whether the ensemble reaches the other state; comparing the two
      ensembles isolates construct/sequence effects (BioEmu never sees DNA).

  (3) Delta per frame.
      delta = d_min(holo) - d_min(apo) for every BioEmu frame. Positive means
      the frame sits closer to the free state.

Outputs (in --output-dir, all suffixed _<atoms>):
  pair_core_summary       one row per pair: core size, trims, states, frames
  pair_core_residues      one row per core position: residue ids + bundle spread
  pair_bundle_rmsd        long: pair, comparison in {apo_apo, holo_holo, apo_holo}
  pair_bundle_summary     (1) medians, effect size, bootstrap CI, informative gate
  pair_frame_state_rmsd   long: pair, ensemble, bundle, state, frame, rmsd
  pair_frame_summary      (3) one row per BioEmu frame: d_min to each bundle, delta
  pair_2x2_summary        (2) one row per (pair, ensemble, bundle) cell

Usage:
    python compute_pair_rmsd.py                       # all pairs, CA atoms
    python compute_pair_rmsd.py --atoms backbone
    python compute_pair_rmsd.py --pair-ids trf1 rok
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from pair_core import load_pairs, load_pair, rmsd_matrix   # noqa: E402

BOOTSTRAP = 2000
SEED = 20260923


def write_csv(path, rows):
    if not rows:
        print(f"  (nothing to write for {path.name})", file=sys.stderr)
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path}  ({len(rows)} rows)")


def upper(m):
    """Off-diagonal upper triangle of a square matrix, flattened."""
    i, j = np.triu_indices(m.shape[0], k=1)
    return m[i, j]


def cles(a, b):
    """Common-language effect size P(a > b), ties counted as half."""
    from scipy.stats import mannwhitneyu
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    u = mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(u / (len(a) * len(b)))


def bootstrap_separation(apo_apo, holo_holo, apo_holo, rng, n_boot=BOOTSTRAP):
    """Percentile CI for median(between) - median(within), resampling STATES.

    States, not pairs, are the sampling unit: the pairs are not independent.
    Resampling with replacement makes some within-bundle pairs compare a state
    with itself; those are dropped rather than counted as exact zeros.
    """
    n_a, n_h = apo_apo.shape[0], holo_holo.shape[0]
    seps = np.empty(n_boot)
    for b in range(n_boot):
        ai = rng.integers(0, n_a, n_a)
        hi = rng.integers(0, n_h, n_h)
        wa = apo_apo[np.ix_(ai, ai)][np.not_equal.outer(ai, ai)]
        wh = holo_holo[np.ix_(hi, hi)][np.not_equal.outer(hi, hi)]
        within = np.concatenate([wa, wh])
        between = apo_holo[np.ix_(ai, hi)].ravel()
        seps[b] = np.median(between) - np.median(within) if within.size else np.nan
    seps = seps[~np.isnan(seps)]
    if seps.size == 0:
        return float("nan"), float("nan")
    return float(np.percentile(seps, 2.5)), float(np.percentile(seps, 97.5))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default=_HERE / "pairs.csv")
    ap.add_argument("--pair-ids", nargs="*", default=None, help="restrict to these pair ids")
    ap.add_argument("--atoms", choices=("ca", "backbone"), default="ca",
                    help="atom subset for superposition + RMSD (default: %(default)s); "
                         "BioEmu output is backbone+CB only, so sidechains never compare")
    ap.add_argument("--output-dir", default=_HERE)
    ap.add_argument("--spread-cutoff", type=float, default=3.0,
                    help="trim terminal core positions whose intra-bundle spread exceeds "
                         "this many Angstrom in either NMR bundle (default: %(default)s)")
    ap.add_argument("--min-ratio", type=float, default=1.25,
                    help="informative gate: median(between)/median(within) must reach this "
                         "and the bootstrap CI must exclude 0 (default: %(default)s)")
    ap.add_argument("--thresholds", type=float, nargs="*", default=(2.0, 3.0),
                    help="coverage/recall thresholds in Angstrom (default: 2 3)")
    ap.add_argument("--no-fetch", action="store_true", help="never download from RCSB")
    args = ap.parse_args()

    pairs = load_pairs(args.pairs, args.pair_ids)
    if not pairs:
        sys.exit("no pairs selected")
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    taus = tuple(args.thresholds)

    core_rows, residue_rows = [], []
    bundle_rows, bundle_summary = [], []
    fs_rows, frame_rows, cell_rows = [], [], []

    for pair in pairs:
        pid = pair["pair_id"]
        print(f"{pid}  ({pair['label']}): apo {pair['apo_pdb']}/{pair['apo_chain']} vs "
              f"holo {pair['holo_pdb']}/{pair['holo_chain']}")
        try:
            P = load_pair(pair, atoms=args.atoms, spread_cutoff=args.spread_cutoff,
                          fetch=not args.no_fetch)
        except Exception as exc:                        # noqa: BLE001 - report, continue
            print(f"  SKIP: {exc}", file=sys.stderr)
            continue

        traj, idx = P["traj"], P["atom_idx"]
        core_rows.append({
            "pair_id": pid, "label": pair["label"], "family": pair["family"],
            "apo_pdb": pair["apo_pdb"], "apo_chain": pair["apo_chain"],
            "apo_ref_chain": P["apo_ref_chain"],
            "holo_pdb": pair["holo_pdb"], "holo_chain": pair["holo_chain"],
            "holo_ref_chain": P["holo_ref_chain"],
            "n_apo_states": P["apo_n_states"], "n_holo_states": P["holo_n_states"],
            "n_apo_frames": P["apo_n_frames"], "n_holo_frames": P["holo_n_frames"],
            "n_intersect": P["n_intersect"], "n_trimmed": P["n_trimmed"],
            "n_core_positions": P["n_core"], "n_atoms": P["n_atoms"],
            "atoms": args.atoms, "spread_cutoff_A": args.spread_cutoff,
        })
        for r in P["core_residues"]:
            residue_rows.append({"pair_id": pid, **r})

        # ---- (1) bundle heterogeneity vs apo/holo separation -------------
        apo_apo = rmsd_matrix(traj["apo_ref"], traj["apo_ref"], idx["apo_ref"], idx["apo_ref"])
        holo_holo = rmsd_matrix(traj["holo_ref"], traj["holo_ref"], idx["holo_ref"], idx["holo_ref"])
        # rows = apo states, cols = holo states
        apo_holo = rmsd_matrix(traj["holo_ref"], traj["apo_ref"], idx["holo_ref"], idx["apo_ref"])

        for name, mat, tri in (("apo_apo", apo_apo, True),
                               ("holo_holo", holo_holo, True),
                               ("apo_holo", apo_holo, False)):
            if tri:
                i, j = np.triu_indices(mat.shape[0], k=1)
            else:
                i, j = (np.repeat(np.arange(mat.shape[0]), mat.shape[1]),
                        np.tile(np.arange(mat.shape[1]), mat.shape[0]))
            for a, b in zip(i, j):
                bundle_rows.append({"pair_id": pid, "comparison": name,
                                    "i": int(a) + 1, "j": int(b) + 1,
                                    "rmsd_A": f"{mat[a, b]:.4f}"})

        w_apo, w_holo, between = upper(apo_apo), upper(holo_holo), apo_holo.ravel()
        within = np.concatenate([w_apo, w_holo])
        sep = float(np.median(between) - np.median(within))
        ratio = float(np.median(between) / np.median(within)) if np.median(within) > 0 else float("inf")
        lo, hi = bootstrap_separation(apo_apo, holo_holo, apo_holo, rng)
        informative = bool(lo > 0 and ratio >= args.min_ratio)
        bundle_summary.append({
            "pair_id": pid, "label": pair["label"],
            "n_apo_states": P["apo_n_states"], "n_holo_states": P["holo_n_states"],
            "median_apo_apo_A": f"{np.median(w_apo):.3f}",
            "median_holo_holo_A": f"{np.median(w_holo):.3f}",
            "median_within_A": f"{np.median(within):.3f}",
            "median_between_A": f"{np.median(between):.3f}",
            "min_between_A": f"{between.min():.3f}",
            "separation_A": f"{sep:.3f}",
            "sep_ci_lo_A": f"{lo:.3f}", "sep_ci_hi_A": f"{hi:.3f}",
            "ratio_between_within": f"{ratio:.3f}",
            "p_between_gt_within": f"{cles(between, within):.3f}",
            "informative": int(informative),
        })
        print(f"  (1) within {np.median(within):.2f} A vs between {np.median(between):.2f} A "
              f"-> separation {sep:+.2f} A [{lo:+.2f}, {hi:+.2f}], "
              f"ratio {ratio:.2f} -> {'INFORMATIVE' if informative else 'gated out'}")

        # ---- (2)+(3) both ensembles vs both bundles ----------------------
        dmin = {}
        for ens_role in ("apo", "holo"):
            ens_name = f"{ens_role}_seq"
            ens = traj[f"{ens_role}_ens"]
            for bundle in ("apo", "holo"):
                ref = traj[f"{bundle}_ref"]
                # (n_states, n_frames)
                m = rmsd_matrix(ens, ref, idx[f"{ens_role}_ens"], idx[f"{bundle}_ref"])
                for s in range(m.shape[0]):
                    for f in range(m.shape[1]):
                        fs_rows.append({"pair_id": pid, "ensemble": ens_name, "bundle": bundle,
                                        "state": s + 1, "frame": f,
                                        "rmsd_A": f"{m[s, f]:.4f}"})
                best = m.argmin(axis=0)
                per_frame = m[best, np.arange(m.shape[1])]
                dmin[(ens_name, bundle)] = (per_frame, best, m)

                cell = {
                    "pair_id": pid, "ensemble": ens_name, "bundle": bundle,
                    "n_frames": int(m.shape[1]), "n_states": int(m.shape[0]),
                    "median_dmin_A": f"{np.median(per_frame):.3f}",
                    "q1_dmin_A": f"{np.percentile(per_frame, 25):.3f}",
                    "q3_dmin_A": f"{np.percentile(per_frame, 75):.3f}",
                    "min_dmin_A": f"{per_frame.min():.3f}",
                    "n_states_hit": int(len(np.unique(best))),
                }
                for tau in taus:
                    key = f"{tau:g}".replace(".", "p")
                    # frame coverage: frames that reach the bundle at all
                    cell[f"frac_frames_within_{key}A"] = f"{(per_frame <= tau).mean():.4f}"
                    # state recall: states reached by at least one frame
                    cell[f"state_recall_{key}A"] = f"{(m.min(axis=1) <= tau).mean():.4f}"
                cell_rows.append(cell)

        for ens_name in ("apo_seq", "holo_seq"):
            d_apo, best_apo, _ = dmin[(ens_name, "apo")]
            d_holo, best_holo, _ = dmin[(ens_name, "holo")]
            for f in range(len(d_apo)):
                delta = d_holo[f] - d_apo[f]
                frame_rows.append({
                    "pair_id": pid, "ensemble": ens_name, "frame": f,
                    "dmin_apo_A": f"{d_apo[f]:.4f}", "best_apo_state": int(best_apo[f]) + 1,
                    "dmin_holo_A": f"{d_holo[f]:.4f}", "best_holo_state": int(best_holo[f]) + 1,
                    "delta_A": f"{delta:.4f}",
                    "closer_to": "apo" if delta > 0 else "holo",
                })
            print(f"  (2/3) BioEmu({ens_name}): median d_min apo {np.median(d_apo):.2f} A, "
                  f"holo {np.median(d_holo):.2f} A, median delta {np.median(d_holo - d_apo):+.2f} A "
                  f"({(d_holo > d_apo).mean():.0%} of frames closer to apo)")

    suffix = f"_{args.atoms}.csv"
    write_csv(out_dir / f"pair_core_summary{suffix}", core_rows)
    write_csv(out_dir / f"pair_core_residues{suffix}", residue_rows)
    write_csv(out_dir / f"pair_bundle_rmsd{suffix}", bundle_rows)
    write_csv(out_dir / f"pair_bundle_summary{suffix}", bundle_summary)
    write_csv(out_dir / f"pair_frame_state_rmsd{suffix}", fs_rows)
    write_csv(out_dir / f"pair_frame_summary{suffix}", frame_rows)
    write_csv(out_dir / f"pair_2x2_summary{suffix}", cell_rows)


if __name__ == "__main__":
    main()
