#!/usr/bin/env python
"""
compute_reachability.py — the two coordinate passes that fill M1's reachability
axis and R2's per-residue profile for pilots the tabular extractor can't reach.

It REUSES compute_rmsds.py's loaders and atom-matching (import, not reimplement),
so the Cα correspondence between the BioEmu ensemble (resSeq 0..N) and the crystal
(real PDB numbering) is identical to the rest of rmsd_analysis.

Two outputs (both keyed on the same whole-protein Cα best-fit superposition the
original analysis used — "Cα-RMSD is whole-protein, superposed on all protein Cα
vs the crystal reference"):

  reachability.csv         one row per pilot:
      pilot,n,d_min,d_med,spread,reach_ratio,rmsf_mean
    d_min / d_med = min / median over the ensemble of the whole-protein Cα-RMSD
        to the crystal bound pose (how close the free ensemble REACHES the bound
        geometry). spread = mean pairwise Cα-RMSD within the ensemble. reach_ratio
        = d_min / spread. rmsf_mean = mean over residues of the per-residue Cα
        std across the ensemble (after best-fit) — the ensemble flexibility.
    -> feeds mechanism_apo_holo.csv's d_min/reach_ratio/rmsf_mean columns.

  ca_rmsd_perresidue.csv   one row per (pilot,resid_idx):
      pilot,resid_idx,per_res_rmsd,n_states
    per_res_rmsd = mean over states of the per-residue Cα deviation from crystal
        after whole-protein Cα best-fit. -> feeds R2.

ENSEMBLE SOURCE: which stage is the "free-state ensemble" is a --source flag.
Default 'stage3' reproduces the existing ca_rmsd_perresidue.csv (Stage-3 minimized
states, whole-protein Cα). For the apo/holo reachability the original memo used the
apo ensemble; --source stage1 gives that. The script prints a validation line
against known values (ets1 d_min≈0.87, tbp d_min≈0.59) so you can confirm the
recipe before trusting new pilots.

Run in the `pycurves` env (mdtraj), on a compute node for big ensembles:
  conda activate pycurves
  python compute_reachability.py --tfs <all> --source stage3 --out ../analysis/data
"""
import argparse, os, sys, glob
import numpy as np

# import compute_rmsds as a library (it guards __main__), reusing its loaders
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "rmsd_analysis"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compute_rmsds as cr


def _ca_pairs(traj_frame, ref, chain_b=None):
    """Cα index pairs (ensemble, crystal) via compute_rmsds' residue-position match."""
    idx_a, idx_b, info = cr.common_atom_indices(traj_frame, ref, chain_a=None, chain_b=chain_b)
    ca = [(a, b, k[0]) for a, b, k in zip(idx_a, idx_b, info) if k[2] == "CA"]
    if not ca:
        return None
    ia = np.array([c[0] for c in ca]); ib = np.array([c[1] for c in ca])
    rpos = np.array([c[2] for c in ca])
    return ia, ib, rpos


def _superpose_rmsd(coords, ref_coords):
    """Kabsch best-fit RMSD (Å) of one frame's Cα coords onto ref (both nm→Å inside)."""
    return cr.kabsch_rmsd(coords, ref_coords)  # returns Å


def _superpose_perres(coords, ref_coords):
    """Best-fit `coords` onto `ref_coords` (Kabsch) and return per-atom deviations (Å)."""
    A = coords - coords.mean(0); B = ref_coords - ref_coords.mean(0)
    H = A.T @ B
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d]); R = Vt.T @ D @ U.T
    A_rot = (R @ A.T).T
    return np.linalg.norm(A_rot - B, axis=1) * 10.0  # nm → Å


def load_frames(tf, pdb, source, n_states):
    """Return a LIST of single-frame trajs (never md.join — per-state minimized
    PDBs can have different atom counts, which join rejects)."""
    if source == "stage1":
        t = cr.load_stage1(pdb)
        return [t[i] for i in range(t.n_frames)] if t is not None else []
    stage_dir = {"stage2": "stage2_docked", "stage3": "stage3_min"}[source]
    frames = cr.load_stage_per_state(pdb, tf, stage_dir, n_states)
    return [f for f in frames if f is not None]


def process_pilot(tf, source):
    cfg = cr.load_pilot_config(tf)
    pdb = cfg.get("PDB_ID", "").lower()
    if not pdb:
        return None, None
    try:
        ref = cr.load_reference_full(pdb)
    except FileNotFoundError:
        print(f"  {tf}: no reference CIF", file=sys.stderr); return None, None
    n_states = 0
    sdir = f"{cr.CONF_ROOT}/stage3_min/{tf}"
    if os.path.isdir(sdir):
        n_states = len(glob.glob(f"{sdir}/{pdb}_state_*.pdb"))
    frames = load_frames(tf, pdb, source, n_states or 200)
    if not frames:
        print(f"  {tf}: no {source} ensemble", file=sys.stderr); return None, None

    # Build each frame's Cα coords in crystal-matched order, PER FRAME (topologies
    # differ across minimized states). Keep only residues common to ALL frames so
    # the per-residue profile and pairwise spread use a consistent residue set.
    ref_ca_by_rpos = {}   # rpos -> ref Cα xyz (nm)
    frame_ca = []         # list of {rpos: xyz}
    for fr in frames:
        pr = _ca_pairs(fr, ref)
        if pr is None:
            continue
        ia, ib, rpos = pr
        fca = {int(r): fr.xyz[0][a] for a, r in zip(ia, rpos)}
        for a_ref, r in zip(ib, rpos):
            ref_ca_by_rpos.setdefault(int(r), ref.xyz[0][a_ref])
        frame_ca.append(fca)
    if not frame_ca:
        print(f"  {tf}: no Cα correspondence", file=sys.stderr); return None, None
    common = sorted(set.intersection(*[set(d.keys()) for d in frame_ca]))
    if not common:
        print(f"  {tf}: no residues common to all frames", file=sys.stderr); return None, None
    rpos_arr = np.array(common)
    ref_ca = np.array([ref_ca_by_rpos[r] for r in common])          # (nres,3) nm
    ens_ca = np.array([[d[r] for r in common] for d in frame_ca])   # (nf,nres,3) nm
    nf = ens_ca.shape[0]

    dists = np.array([_superpose_rmsd(ens_ca[f], ref_ca) for f in range(nf)])
    d_min = float(np.nanmin(dists)); d_med = float(np.nanmedian(dists))
    perres = np.array([_superpose_perres(ens_ca[f], ref_ca) for f in range(nf)])
    per_res_mean = np.nanmean(perres, axis=0)

    if nf > 1:
        import itertools, random
        pairs = list(itertools.combinations(range(nf), 2))
        if len(pairs) > 2000:
            random.Random(0).shuffle(pairs); pairs = pairs[:2000]
        pw = np.array([_superpose_rmsd(ens_ca[i], ens_ca[j]) for i, j in pairs])
        spread = float(np.nanmean(pw))
    else:
        spread = float("nan")
    mean_frame = ens_ca.mean(0)
    devs = np.array([_superpose_perres(ens_ca[f], mean_frame) for f in range(nf)])
    rmsf_mean = float(np.nanmean(np.sqrt((devs ** 2).mean(0))))

    reach = {"pilot": tf, "n": int(nf), "d_min": round(d_min, 4),
             "d_med": round(d_med, 4), "spread": round(spread, 4),
             "reach_ratio": round(d_min / spread, 6) if spread and not np.isnan(spread) else float("nan"),
             "rmsf_mean": round(rmsf_mean, 4)}
    perres_rows = [{"pilot": tf, "resid_idx": int(r), "per_res_rmsd": round(float(v), 4),
                    "n_states": int(nf)} for r, v in zip(rpos_arr, per_res_mean)]
    return reach, perres_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tfs", nargs="+", required=True)
    ap.add_argument("--source", choices=["stage1", "stage2", "stage3"], default="stage3")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()
    cr.import_mdtraj()

    import pandas as pd
    reach_rows, perres_all = [], []
    for tf in args.tfs:
        print(f"[{tf}] source={args.source}")
        reach, perres = process_pilot(tf, args.source)
        if reach:
            reach_rows.append(reach)
            print(f"  d_min={reach['d_min']} d_med={reach['d_med']} "
                  f"spread={reach['spread']} reach_ratio={reach['reach_ratio']} n={reach['n']}")
        if perres:
            perres_all.extend(perres)

    os.makedirs(args.out, exist_ok=True)
    if reach_rows:
        rp = os.path.join(args.out, "reachability.csv")
        pd.DataFrame(reach_rows).to_csv(rp, index=False); print(f"wrote {rp}")
    if perres_all:
        pp = os.path.join(args.out, "ca_rmsd_perresidue.csv")
        pd.DataFrame(perres_all).to_csv(pp, index=False); print(f"wrote {pp} ({len(perres_all)} rows)")

    # validation hint against known values
    known = {"ets1": 0.87, "tbp": 0.59}
    for r in reach_rows:
        if r["pilot"] in known:
            exp = known[r["pilot"]]
            ok = abs(r["d_min"] - exp) < 0.15
            print(f"VALIDATION {r['pilot']}: d_min={r['d_min']} vs known≈{exp} "
                  f"[{'OK' if ok else 'MISMATCH — check --source'}]")


if __name__ == "__main__":
    main()
