# Compute Cα RMSD of each BioEmu-sampled conformation against the crystal
# structure and write a single CSV. No plotting, no structure rendering —
# see plot_rmsd.py and generate_structure_images.py for those.
#
# For each <MONOMER_DIR>/<PDB>_chains/<PDB>_conformations/ directory that
# contains both topology.pdb and samples.xtc, this script:
#   1. Loads the reference crystal protein and the conformational trajectory
#      into PyMOL.
#   2. For each sampled state, aligns it to the reference (sequence-based,
#      Cα atoms only, no outlier rejection) and records the RMSD.
#   3. Writes <MONOMER_DIR>/rmsd/rmsds.csv with one row per (pdb_id, state,
#      rmsd, pfam_family). The pfam_family column is populated from
#      deeppbs_tf_pfam_metadata.csv when --pfam-metadata is supplied; left
#      blank otherwise.
#
# Usage:
#   python compute_rmsds.py <MONOMER_DIR>
#   python compute_rmsds.py <MONOMER_DIR> 1abc 2def                 # filter to these PDB IDs
#   python compute_rmsds.py <MONOMER_DIR> --pfam-metadata path/to/deeppbs_tf_pfam_metadata.csv

import os
import sys
import glob
import csv
import argparse

import numpy as np
from pymol import cmd


def compute_rmsds_for_monomer(chains_dir, protein_pdb):
    """Compute Cα RMSD of each sampled state vs. the reference crystal.

    Returns a list of RMSD floats (one per state), or None if the trajectory
    is missing.
    """
    pdb_id = os.path.basename(chains_dir).replace("_chains", "")
    conf_dir = os.path.join(chains_dir, f"{pdb_id}_conformations")
    topo = os.path.join(conf_dir, "topology.pdb")
    xtc = os.path.join(conf_dir, "samples.xtc")
    if not (os.path.exists(topo) and os.path.exists(xtc)):
        return None

    cmd.delete("all")
    cmd.load(protein_pdb, "ref")
    cmd.load(topo, "confs")
    cmd.load_traj(xtc, "confs")
    n_states = cmd.count_states("confs")

    rmsds = []
    for state in range(1, n_states + 1):
        # cycles=0 disables outlier rejection so RMSD is reported over all
        # sequence-matched Cα atoms rather than a refined core subset.
        result = cmd.align(
            f"confs and name CA and state {state}",
            "ref and name CA",
            mobile_state=state,
            cycles=0,
        )
        rmsds.append(float(result[0]))

    cmd.delete("all")
    return rmsds


def compute_all(chains_dirs):
    """Compute RMSDs for each chains_dir; return dict pdb_id -> list[float]."""
    results = {}
    for chains_dir in chains_dirs:
        pdb_id = os.path.basename(chains_dir).replace("_chains", "")
        proteins = sorted(glob.glob(os.path.join(chains_dir, "*_chain*_protein.pdb")))
        if not proteins:
            print(f"[skip] {pdb_id}: no reference protein PDB")
            continue
        rmsds = compute_rmsds_for_monomer(chains_dir, proteins[0])
        if rmsds is None:
            print(f"[skip] {pdb_id}: missing topology.pdb or samples.xtc")
            continue
        results[pdb_id] = rmsds
        print(
            f"[{pdb_id}] n={len(rmsds)} "
            f"mean={np.mean(rmsds):.2f}Å median={np.median(rmsds):.2f}Å "
            f"max={np.max(rmsds):.2f}Å"
        )
    return results


def load_pfam_lookup(metadata_csv):
    """Return dict pdb_id (lowercase) -> primary Pfam family name.

    metadata_csv columns: PDB_ID, Chain, PWM_ID, UniProt_ID, TF_Class,
    TF_Family, Pfam_ID, Pfam_Name. Pfam_Name may be "|"-delimited when a
    chain has multiple domains; the first entry is taken as the primary
    family. Where multiple chains exist for a PDB, the alphabetically-first
    chain is used (matches the prior behavior).
    """
    chains_by_pdb = {}
    with open(metadata_csv) as f:
        for row in csv.DictReader(f):
            pid = (row.get("PDB_ID") or "").strip().lower()
            if not pid:
                continue
            chain = (row.get("Chain") or "").strip()
            names = (row.get("Pfam_Name") or "").strip()
            primary = names.split("|", 1)[0] if names else ""
            chains_by_pdb.setdefault(pid, []).append((chain, primary))

    return {
        pid: sorted(entries, key=lambda x: x[0])[0][1] or "Unclassified"
        for pid, entries in chains_by_pdb.items()
    }


def write_rmsds_csv(results, pfam_by_pdb, csv_path, preserve_others=False):
    """Write results to csv_path. Columns: pdb_id, state, rmsd_angstrom, pfam_family.

    If preserve_others, retain rows for pdb_ids not in `results`.
    """
    existing_rows = []
    if preserve_others and os.path.exists(csv_path):
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                if row["pdb_id"] not in results:
                    existing_rows.append(row)

    with open(csv_path, "w", newline="") as csvf:
        w = csv.writer(csvf)
        w.writerow(["pdb_id", "state", "rmsd_angstrom", "pfam_family"])
        for row in existing_rows:
            w.writerow([
                row["pdb_id"], row["state"], row["rmsd_angstrom"],
                row.get("pfam_family", ""),
            ])
        for pdb_id in sorted(results):
            family = pfam_by_pdb.get(pdb_id.lower(), "")
            for i, r in enumerate(results[pdb_id], start=1):
                w.writerow([pdb_id, i, f"{r:.4f}", family])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute Cα RMSDs vs. crystal structure and write a CSV.",
    )
    parser.add_argument(
        "monomer_dir",
        help="Directory containing <PDB>_chains/ subdirectories.",
    )
    parser.add_argument(
        "pdb_ids", nargs="*",
        help="Restrict to these PDB IDs (default: all monomers with a trajectory).",
    )
    parser.add_argument(
        "--pfam-metadata",
        help="Path to deeppbs_tf_pfam_metadata.csv. When provided, populates the "
             "pfam_family column in the output CSV.",
    )
    args = parser.parse_args()

    monomers_dir = os.path.abspath(args.monomer_dir)
    out_dir = os.path.join(monomers_dir, "rmsd")
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "rmsds.csv")

    pfam_by_pdb = {}
    if args.pfam_metadata:
        if not os.path.exists(args.pfam_metadata):
            sys.exit(f"Pfam metadata CSV not found: {args.pfam_metadata}")
        pfam_by_pdb = load_pfam_lookup(args.pfam_metadata)
        print(f"Loaded Pfam metadata for {len(pfam_by_pdb)} PDB IDs from {args.pfam_metadata}")

    requested = {x.upper() for x in args.pdb_ids}
    chains_dirs = sorted(glob.glob(os.path.join(monomers_dir, "*_chains")))
    if requested:
        chains_dirs = [
            d for d in chains_dirs
            if os.path.basename(d).replace("_chains", "").upper() in requested
        ]
    print(f"Processing {len(chains_dirs)} monomer directory(ies) from {monomers_dir}/")

    all_rmsds = compute_all(chains_dirs)
    write_rmsds_csv(all_rmsds, pfam_by_pdb, csv_path, preserve_others=bool(requested))
    print(f"Wrote {csv_path}")
