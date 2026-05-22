# Render one PNG per monomer showing the crystal structure (blue) overlaid
# with the best (green), 50th-percentile (yellow), and worst (red) sampled
# conformations. DNA, when present alongside the reference, is shown in
# salmon.
#
# The best/median/worst states are chosen from <MONOMER_DIR>/rmsd/rmsds.csv
# (written by compute_rmsds.py). PNGs are written to
# <MONOMER_DIR>/rmsd/structures/<PDB>.png.
#
# Usage:
#   python generate_structure_images.py <MONOMER_DIR>
#   python generate_structure_images.py <MONOMER_DIR> 1abc 2def    # filter to these PDB IDs
#   python generate_structure_images.py <MONOMER_DIR> --overwrite  # re-render existing images

import os
import sys
import csv
import glob
import argparse

import numpy as np
from pymol import cmd


COLOR_REF = "skyblue"
COLOR_BEST = "forest"
COLOR_MEDIAN = "yellow"
COLOR_WORST = "red"
COLOR_DNA = "salmon"


def load_rmsds_csv(csv_path, requested):
    """Return dict pdb_id -> list[(state:int, rmsd:float)] in CSV order."""
    results = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            pid = row["pdb_id"]
            if requested and pid.upper() not in requested:
                continue
            results.setdefault(pid, []).append(
                (int(row["state"]), float(row["rmsd_angstrom"]))
            )
    return results


def pick_states(state_rmsds):
    """Return (best_state, median_state, worst_state) by RMSD."""
    states = np.array([s for s, _ in state_rmsds])
    rmsds = np.array([r for _, r in state_rmsds])
    order = np.argsort(rmsds)
    best = int(states[order[0]])
    worst = int(states[order[-1]])
    median = int(states[order[len(order) // 2]])
    return best, median, worst


def render_overlay(protein_pdb, dna_pdb, topo, xtc, states, out_path, size=600):
    """Render ref + 3 picked sampled states (best/median/worst) as a single PNG."""
    best_state, median_state, worst_state = states

    cmd.delete("all")
    cmd.load(protein_pdb, "ref")
    cmd.load(topo, "_traj")
    cmd.load_traj(xtc, "_traj")

    cmd.create("best", "_traj", source_state=best_state, target_state=1)
    cmd.create("median", "_traj", source_state=median_state, target_state=1)
    cmd.create("worst", "_traj", source_state=worst_state, target_state=1)
    cmd.delete("_traj")

    for obj in ("best", "median", "worst"):
        cmd.align(f"{obj} and name CA", "ref and name CA", cycles=0)

    cmd.hide("everything")
    cmd.show("cartoon", "ref or best or median or worst")
    cmd.color(COLOR_REF, "ref")
    cmd.color(COLOR_BEST, "best")
    cmd.color(COLOR_MEDIAN, "median")
    cmd.color(COLOR_WORST, "worst")

    if dna_pdb and os.path.exists(dna_pdb):
        cmd.load(dna_pdb, "dna")
        cmd.show("cartoon", "dna")
        cmd.color(COLOR_DNA, "dna")

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_transparency", 0.0, "ref")
    cmd.orient()
    cmd.zoom("all", buffer=2)
    cmd.ray(size, size)
    cmd.png(out_path, dpi=150)
    cmd.delete("all")


def render_for_monomer(monomers_dir, structures_dir, pdb_id, state_rmsds, overwrite):
    out_path = os.path.join(structures_dir, f"{pdb_id}.png")
    if os.path.exists(out_path) and not overwrite:
        print(f"[skip] {pdb_id}: {out_path} already exists (use --overwrite to redraw)")
        return

    chains_dir = os.path.join(monomers_dir, f"{pdb_id}_chains")
    proteins = sorted(glob.glob(os.path.join(chains_dir, "*_chain*_protein.pdb")))
    if not proteins:
        print(f"[skip] {pdb_id}: no reference protein PDB in {chains_dir}")
        return
    conf_dir = os.path.join(chains_dir, f"{pdb_id}_conformations")
    topo = os.path.join(conf_dir, "topology.pdb")
    xtc = os.path.join(conf_dir, "samples.xtc")
    if not (os.path.exists(topo) and os.path.exists(xtc)):
        print(f"[skip] {pdb_id}: missing topology.pdb or samples.xtc in {conf_dir}")
        return

    dna_files = glob.glob(os.path.join(chains_dir, "*_dna.pdb"))
    dna_pdb = dna_files[0] if dna_files else None

    states = pick_states(state_rmsds)
    print(
        f"[{pdb_id}] best=state {states[0]} (RMSD {min(r for _, r in state_rmsds):.2f}Å) "
        f"median=state {states[1]} worst=state {states[2]} "
        f"(RMSD {max(r for _, r in state_rmsds):.2f}Å)"
    )
    try:
        render_overlay(proteins[0], dna_pdb, topo, xtc, states, out_path)
        print(f"  wrote {out_path}")
    except Exception as e:
        print(f"[warn] {pdb_id}: render failed ({e})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Render reference + best/median/worst-RMSD sampled structures per monomer.",
    )
    parser.add_argument(
        "monomer_dir",
        help="Directory containing <PDB>_chains/ subdirectories and rmsd/rmsds.csv.",
    )
    parser.add_argument(
        "pdb_ids", nargs="*",
        help="Restrict to these PDB IDs (default: all monomers in rmsds.csv).",
    )
    parser.add_argument(
        "--rmsds-csv",
        help="Path to rmsds.csv (default: <MONOMER_DIR>/rmsd/rmsds.csv).",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Re-render structure PNGs even if they already exist.",
    )
    args = parser.parse_args()

    monomers_dir = os.path.abspath(args.monomer_dir)
    out_dir = os.path.join(monomers_dir, "rmsd")
    structures_dir = os.path.join(out_dir, "structures")
    os.makedirs(structures_dir, exist_ok=True)
    csv_path = args.rmsds_csv or os.path.join(out_dir, "rmsds.csv")
    if not os.path.exists(csv_path):
        sys.exit(f"{csv_path} not found — run compute_rmsds.py first.")

    requested = {x.upper() for x in args.pdb_ids}
    by_pdb = load_rmsds_csv(csv_path, requested)
    if not by_pdb:
        sys.exit("No matching monomers in CSV.")

    print(f"Rendering {len(by_pdb)} monomer(s) → {structures_dir}/")
    for pdb_id, state_rmsds in by_pdb.items():
        render_for_monomer(monomers_dir, structures_dir, pdb_id, state_rmsds, args.overwrite)
