#!/usr/bin/env python3
"""AF3 vs crystal RMSD, and BioEmu-vs-AF3 ensemble diversity.

Reuses the repo's own RMSD machinery (rmsd_analysis/compute_rmsds.py):
mdtraj loading, `common_atom_indices` (matches atoms by residue-ordinal +
name, robust to crystal/BioEmu/AF3 numbering differences) and `kabsch_rmsd`
(nm -> A). We restrict to Cα for a whole-protein backbone RMSD consistent
with the repo's ca_rmsd (= Cα-RMSD from the crystal bound pose).

Writes two CSVs into analysis/data/:
  - af3_rmsd_to_crystal.csv    per-AF3-model Cα-RMSD to the crystal chain
  - ensemble_diversity.csv     within-ensemble pairwise Cα-RMSD (BioEmu vs AF3)

Run in the `bioemu` conda env (has mdtraj 1.11). Cap BLAS threads on the
login node (see compute_details runbook); a compute node needs no cap.
"""
import os, sys, glob, csv, itertools, warnings
warnings.filterwarnings("ignore")
import numpy as np

BASE = "/project2/rohs_102/shewchuk/TF-conformation"
sys.path.insert(0, os.path.join(BASE, "rmsd_analysis"))
import compute_rmsds as cr      # repo module
cr.import_mdtraj()              # binds cr.md
md = cr.md

OUTD = os.path.join(BASE, "analysis", "data")

# pilot -> pdb id, af3 output subdir. Ordered by fnat pass-rate (repo convention).
PILOTS = [
    ("ets1",      "1k79", "ets1_1k79"),
    ("tbp",       "1tgh", "tbp_1tgh"),
    ("egr1",      "1aay", "egr1_1aay"),
    ("engrailed", "3hdd", "engrailed_3hdd"),
    ("foxa",      "1vtn", "foxa_1vtn"),
    ("lef1",      "2lef", "lef1_2lef"),
]


def ca_rmsd(traj_a, traj_b, chain_a=None, chain_b=None):
    """Cα-only Kabsch RMSD (Å) between single frames of two mdtraj trajs,
    matching atoms via the repo's residue-ordinal correspondence. chain_a/
    chain_b restrict each traj to one protein chain index (needed when the
    crystal reference carries a second protein copy; AF3/docked are single-
    chain so None is fine — DNA is skipped by residue name regardless)."""
    idx_a, idx_b, info = cr.common_atom_indices(traj_a, traj_b,
                                                chain_a=chain_a, chain_b=chain_b)
    # keep only Cα among the matched heavy atoms
    ca = np.array([traj_a.topology.atom(i).name == "CA" for i in idx_a])
    ca_a, ca_b = idx_a[ca], idx_b[ca]
    A = traj_a.xyz[0, ca_a, :]
    B = traj_b.xyz[0, ca_b, :]
    return cr.kabsch_rmsd(A, B), int(ca.sum())


def pairwise_diversity(trajs):
    """All unique-pair Cα-RMSDs plus summary stats."""
    vals = [ca_rmsd(trajs[i], trajs[j])[0]
            for i, j in itertools.combinations(range(len(trajs)), 2)]
    v = np.asarray(vals, dtype=float)
    stats = dict(n=len(trajs), n_pairs=len(v),
                 median=float(np.median(v)), mean=float(v.mean()),
                 p90=float(np.percentile(v, 90)), maxv=float(v.max()),
                 std=float(v.std()))
    return v, stats


def main():
    os.makedirs(OUTD, exist_ok=True)
    acc_rows, div_rows, pair_rows = [], [], []

    for tf, pdb, af3dir in PILOTS:
        ref = cr.load_reference_full(pdb)   # crystal .cif via repo loader
        cfg = cr.load_pilot_config(tf)
        prot_chain = int(cfg["PROTEIN_CHAIN"])  # crystal binding-chain index

        # --- AF3 models: 10 diffusion samples (seed-*/sample-*) ---
        af3_paths = sorted(glob.glob(f"{BASE}/af3/output/{af3dir}/seed-*/*_model.cif"))
        af3_trajs = []
        for p in af3_paths:
            base = os.path.basename(os.path.dirname(p))     # seed-1_sample-0
            seed = base.split("_")[0].replace("seed-", "")
            samp = base.split("_")[1].replace("sample-", "")
            t = md.load(p)
            r, n = ca_rmsd(t, ref, chain_a=None, chain_b=prot_chain)
            acc_rows.append(dict(pilot=tf, pdb=pdb, source="af3", model=base,
                                 seed=seed, sample=samp,
                                 ca_rmsd=round(r, 4), n_ca=n))
            af3_trajs.append(t)

        # --- BioEmu docked frames ---
        bio_paths = sorted(glob.glob(f"{BASE}/output/stage2_docked/{tf}/*.pdb"))
        bio_trajs = [md.load(p) for p in bio_paths]

        # --- diversity per ensemble ---
        for src, trajs in [("bioemu", bio_trajs), ("af3", af3_trajs)]:
            vals, s = pairwise_diversity(trajs)
            div_rows.append(dict(pilot=tf, source=src, n=s["n"],
                                 n_pairs=s["n_pairs"],
                                 median_pairwise=round(s["median"], 4),
                                 mean_pairwise=round(s["mean"], 4),
                                 p90_pairwise=round(s["p90"], 4),
                                 max_pairwise=round(s["maxv"], 4),
                                 std_pairwise=round(s["std"], 4)))
            for val in vals:
                pair_rows.append(dict(pilot=tf, source=src,
                                      pairwise_ca_rmsd=round(float(val), 4)))
        print(f"{tf:10s} AF3 n={len(af3_trajs)} BioEmu n={len(bio_trajs)} | "
              f"BioEmu medPair={div_rows[-2]['median_pairwise']:.3f}  "
              f"AF3 medPair={div_rows[-1]['median_pairwise']:.3f}")

    acc_csv = os.path.join(OUTD, "af3_rmsd_to_crystal.csv")
    with open(acc_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pilot", "pdb", "source", "model",
                                          "seed", "sample", "ca_rmsd", "n_ca"])
        w.writeheader(); w.writerows(acc_rows)

    div_csv = os.path.join(OUTD, "ensemble_diversity.csv")
    with open(div_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pilot", "source", "n", "n_pairs",
                                          "median_pairwise", "mean_pairwise",
                                          "p90_pairwise", "max_pairwise", "std_pairwise"])
        w.writeheader(); w.writerows(div_rows)

    pair_csv = os.path.join(OUTD, "ensemble_diversity_pairwise.csv")
    with open(pair_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pilot", "source", "pairwise_ca_rmsd"])
        w.writeheader(); w.writerows(pair_rows)

    print(f"\nwrote {acc_csv} ({len(acc_rows)} rows)")
    print(f"wrote {div_csv} ({len(div_rows)} rows)")
    print(f"wrote {pair_csv} ({len(pair_rows)} rows)")


if __name__ == "__main__":
    main()
