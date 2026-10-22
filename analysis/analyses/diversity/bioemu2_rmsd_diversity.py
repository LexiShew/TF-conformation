#!/usr/bin/env python3
"""BioEmu2 arm of the ensemble-diversity comparison (AF3 vs BioEmu vs BioEmu2).

The BioEmu2 counterpart to af3_rmsd_diversity.py, which produced
analysis/data/ensemble_diversity{,_pairwise}.csv for the BioEmu-1-vs-AF3 figure
(D1). It reuses the SAME machinery and the SAME quantity so the new arm plots on
the existing axis:

  * rmsd_analysis/compute_rmsds.py for mdtraj loading, `common_atom_indices`
    (atom correspondence by residue-ordinal + name) and `kabsch_rmsd` (nm -> A)
  * Calpha only, whole-protein, consistent with the repo's ca_rmsd
  * all unique within-ensemble pairs, same summary columns

Separate script rather than an edit to af3_rmsd_diversity.py because that
script's outputs are committed canonical tables in analysis/data/. This arm
writes to the theme's own data/ dir (analysis/docs/LAYOUT.md) so the canonical
tables are untouched and this arm regenerates independently as ensembles land.

THE ONE REAL DIFFERENCE, and it matters
---------------------------------------
BioEmu 1 frames are a single docked protein chain. BioEmu2 samples the COMPLEX
from the deposited assembly, so a frame can hold several protein chains - 1K79
and 3HDD both carry two protein copies. Pooling both protomers' Calpha would not
be comparable to AF3 (one chain) or BioEmu 1 (one chain), and would be inflated
by relative motion BETWEEN copies.

Each frame is therefore restricted to its FIRST PROTEIN CHAIN, found by residue
name via compute_rmsds.PROTEIN_RESNAMES rather than a fixed chain index - chain
order varies per structure (1SKN imports dna,dna,protein), so a hardcoded index
would silently select DNA. DNA needs no handling: common_atom_indices keeps only
protein residues.

Why the crystal column is gated
-------------------------------
common_atom_indices matches on (residue_ordinal, residue_name, atom_name), the
ordinal counting protein residues within the selected chain. Correct for BioEmu
1, whose frames came from the crystal's OBSERVED residues, so ordinals line up.
BioEmu2 builds from the deposited _entity_poly SEQRES and so models residues the
crystal never resolved (egr1 90 vs 85, engrailed 60 vs 55, tbp 185 vs 180). The
ordinals shift, residue names stop agreeing, and correspondence collapses to a
few coincidental atoms - the first run matched 1 Calpha for engrailed and
reported a median RMSD to crystal of exactly 0.000, which is the tell.

Rather than emit a plausible-looking wrong number, a frame whose correspondence
covers less than MIN_COVERAGE of the shorter chain gets an EMPTY ca_rmsd plus its
coverage, so the gap is explicit. Fixing it properly needs a SEQUENCE alignment
between model and crystal chain, not ordinal matching. Until then only the
within-ensemble diversity columns are trustworthy - and those are unaffected,
because every frame in an ensemble shares one topology (verified: matched Calpha
== full chain length for all 7 pilots measured).

Env: the `bioemu` conda env (mdtraj 1.11). On a login node cap BLAS threads
first; a compute node needs no cap. Prefer the sbatch wrapper:
    sbatch --partition=qcbr --account=rohs_102 \
      analysis/analyses/diversity/bioemu2_rmsd_diversity.sbatch

Outputs (analysis/analyses/diversity/data/)
    bioemu2_diversity.csv           same columns as ensemble_diversity.csv
    bioemu2_diversity_pairwise.csv  same columns as ensemble_diversity_pairwise.csv
    bioemu2_rmsd_to_crystal.csv     per-frame Calpha-RMSD to crystal, GATED
    bioemu2_diversity_skipped.csv   pilots with no usable ensemble, and why
"""
import argparse
import csv
import glob
import itertools
import json
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(BASE, "analysis", "analyses", "rmsd_analysis"))
import compute_rmsds as cr                                        # noqa: E402
cr.import_mdtraj()
md = cr.md

DATA = os.path.join(HERE, "data")

# Minimum fraction of the shorter chain that must correspond for a model-vs-crystal
# RMSD to count. 0.9 is strict on purpose: partial matches are exactly how a wrong
# number gets reported as a real one.
MIN_COVERAGE = 0.90

# Same pilot order as af3_rmsd_diversity.py (fnat pass-rate, repo convention).
PILOTS = [
    ("ets1", "1k79"), ("tbp", "1tgh"), ("egr1", "1aay"), ("engrailed", "3hdd"),
    ("foxa", "1vtn"), ("lef1", "2lef"), ("csl", "3brg"), ("err", "1lo1"),
    ("nfat", "1a66"), ("runx", "1hjc"), ("hsf", "5d5u"),
]


def first_protein_chain(traj):
    """Index of the first chain whose residues are protein, or None."""
    for ch in traj.topology.chains:
        if {r.name.strip() for r in ch.residues} & cr.PROTEIN_RESNAMES:
            return ch.index
    return None


def n_protein_residues(traj, chain_index):
    return sum(1 for r in traj.topology.chain(chain_index).residues
               if r.name in cr.PROTEIN_RESNAMES)


def ca_rmsd(traj_a, traj_b, chain_a=None, chain_b=None):
    """Calpha-only Kabsch RMSD (A) between single frames - same as the AF3 arm."""
    idx_a, idx_b, _info = cr.common_atom_indices(traj_a, traj_b,
                                                 chain_a=chain_a, chain_b=chain_b)
    ca = np.array([traj_a.topology.atom(i).name == "CA" for i in idx_a])
    if not ca.any():
        return float("nan"), 0
    A = traj_a.xyz[0, idx_a[ca], :]
    B = traj_b.xyz[0, idx_b[ca], :]
    return cr.kabsch_rmsd(A, B), int(ca.sum())


def complete_ensemble(pdb_id, want=100):
    """Frame paths for a BioEmu2 ensemble, only if its manifest says it finished."""
    out = os.path.join(BASE, "output", "bioemu2", f"{pdb_id}_complex")
    man = os.path.join(out, "sample_manifest.json")
    if not os.path.exists(man):
        return None, "no manifest"
    try:
        n = json.load(open(man)).get("n_frames_written", 0)
    except Exception:                                             # noqa: BLE001
        return None, "unreadable manifest"
    paths = sorted(glob.glob(os.path.join(out, "states", f"{pdb_id}_state_*.pdb")))
    if n != want or len(paths) != want:
        return None, f"incomplete ({n} in manifest, {len(paths)} state files)"
    return paths, "ok"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pilots", nargs="*", default=None)
    ap.add_argument("--num-samples", type=int, default=100)
    a = ap.parse_args()

    os.makedirs(DATA, exist_ok=True)
    pilots = [(t, p) for t, p in PILOTS if not a.pilots or t in a.pilots]

    div_rows, pair_rows, acc_rows, skipped = [], [], [], []
    for tf, pdb in pilots:
        paths, why = complete_ensemble(pdb, a.num_samples)
        if paths is None:
            print(f"{tf:10s} SKIP - {why}", flush=True)
            skipped.append(dict(pilot=tf, pdb=pdb, reason=why))
            continue

        trajs = [md.load(p) for p in paths]
        chains = [first_protein_chain(t) for t in trajs]
        if any(c is None for c in chains):
            print(f"{tf:10s} SKIP - no protein chain in some frames", flush=True)
            skipped.append(dict(pilot=tf, pdb=pdb, reason="no protein chain"))
            continue

        # --- model vs crystal, GATED on correspondence coverage ---
        ref = cr.load_reference_full(pdb)
        prot_chain = int(cr.load_pilot_config(tf)["PROTEIN_CHAIN"])
        expect = min(n_protein_residues(trajs[0], chains[0]),
                     n_protein_residues(ref, prot_chain))
        for pth, t, ch in zip(paths, trajs, chains):
            r, n_ca = ca_rmsd(t, ref, chain_a=ch, chain_b=prot_chain)
            cov = (n_ca / expect) if expect else 0.0
            ok = cov >= MIN_COVERAGE
            acc_rows.append(dict(pilot=tf, pdb=pdb, source="bioemu2",
                                 model=os.path.basename(pth),
                                 ca_rmsd=round(float(r), 4) if ok else "",
                                 n_ca=n_ca, n_expected=expect,
                                 coverage=round(cov, 3),
                                 gated="" if ok else "low_correspondence"))

        # --- within-ensemble diversity: one protein chain per frame ---
        vals = [ca_rmsd(trajs[i], trajs[j], chain_a=chains[i], chain_b=chains[j])[0]
                for i, j in itertools.combinations(range(len(trajs)), 2)]
        v = np.asarray(vals, dtype=float)
        v = v[~np.isnan(v)]
        div_rows.append(dict(pilot=tf, source="bioemu2", n=len(trajs), n_pairs=len(v),
                             median_pairwise=round(float(np.median(v)), 4),
                             mean_pairwise=round(float(v.mean()), 4),
                             p90_pairwise=round(float(np.percentile(v, 90)), 4),
                             max_pairwise=round(float(v.max()), 4),
                             std_pairwise=round(float(v.std()), 4)))
        pair_rows.extend(dict(pilot=tf, source="bioemu2",
                              pairwise_ca_rmsd=round(float(x), 4)) for x in v)
        g = acc_rows[-1]
        print(f"{tf:10s} n={len(trajs)} pairs={len(v)} "
              f"medPair={div_rows[-1]['median_pairwise']:.3f} "
              f"crystal_cov={g['coverage']:.2f}"
              f"{' GATED' if g['gated'] else ''}", flush=True)

    def dump(rows, name, cols):
        path = os.path.join(DATA, name)
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {path} ({len(rows)} rows)")

    dump(div_rows, "bioemu2_diversity.csv",
         ["pilot", "source", "n", "n_pairs", "median_pairwise", "mean_pairwise",
          "p90_pairwise", "max_pairwise", "std_pairwise"])
    dump(pair_rows, "bioemu2_diversity_pairwise.csv",
         ["pilot", "source", "pairwise_ca_rmsd"])
    dump(acc_rows, "bioemu2_rmsd_to_crystal.csv",
         ["pilot", "pdb", "source", "model", "ca_rmsd", "n_ca", "n_expected",
          "coverage", "gated"])
    if skipped:
        dump(skipped, "bioemu2_diversity_skipped.csv", ["pilot", "pdb", "reason"])


if __name__ == "__main__":
    main()
