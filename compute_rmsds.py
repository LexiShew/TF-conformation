# Compute Cα RMSD of each BioEmu-sampled conformation against the crystal
# structure, then plot per-monomer and aggregated RMSD distributions.
#
# For each monomers/<PDB>_chains/ directory that contains both topology.pdb
# and samples.xtc, this script:
#   1. Renders a ray-traced PyMOL cartoon of the crystal protein (+ DNA if
#      present) to rmsd/structures/<PDB>.png.
#   2. Loads the reference crystal protein and the conformational trajectory
#      into PyMOL.
#   3. For each sampled state, aligns it to the reference (sequence-based,
#      Cα atoms only, no outlier rejection) and records the RMSD.
#   4. Writes rmsd/rmsds.csv with one row per (pdb_id, state, rmsd).
#   5. Writes one plot per monomer to rmsd/plots/<PDB>.png: structure image
#      on the left, RMSD histogram on the right. Plus two summary plots:
#      rmsd/summary_per_monomer.png (box plot per monomer, sorted by median)
#      and rmsd/summary_aggregated.png (pooled histogram across all monomers).
#
# Usage:
#   python compute_rmsds.py                 # all complete monomers
#   python compute_rmsds.py 1CIT 6PAX       # restrict to specific PDB IDs

import os
import sys
import glob
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
from pymol import cmd


def render_structure_image(protein_pdb, dna_pdb, out_path, size=500):
    """Ray-trace a PyMOL cartoon of the reference protein (+ DNA) to PNG."""
    cmd.delete("all")
    cmd.load(protein_pdb, "prot")
    cmd.hide("everything")
    cmd.show("cartoon", "prot")
    cmd.color("skyblue", "prot")
    if dna_pdb and os.path.exists(dna_pdb):
        cmd.load(dna_pdb, "dna")
        cmd.show("cartoon", "dna")
        cmd.color("salmon", "dna")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.orient()
    cmd.zoom("all", buffer=2)
    cmd.ray(size, size)
    cmd.png(out_path, dpi=150)
    cmd.delete("all")


def compute_rmsds_for_monomer(chains_dir, protein_pdb):
    """Compute Cα RMSD of each sampled state vs. the reference crystal.

    Returns:
        List of RMSD floats, one per state, or None if the trajectory is missing.
    """
    topo = os.path.join(chains_dir, "topology.pdb")
    xtc = os.path.join(chains_dir, "samples.xtc")
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


def plot_per_monomer(pdb_id, rmsds, structure_img_path, out_path):
    fig = plt.figure(figsize=(10, 4.2))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.5], wspace=0.15)

    ax_img = fig.add_subplot(gs[0])
    if structure_img_path and os.path.exists(structure_img_path):
        ax_img.imshow(mpimg.imread(structure_img_path))
    ax_img.axis("off")
    ax_img.set_title(pdb_id, fontsize=13, pad=6)

    ax_hist = fig.add_subplot(gs[1])
    ax_hist.hist(rmsds, bins=30, edgecolor="black", alpha=0.75, color="steelblue")
    mean = float(np.mean(rmsds))
    median = float(np.median(rmsds))
    ax_hist.axvline(mean, color="red", linestyle="--", linewidth=1.3,
                    label=f"mean = {mean:.2f} Å")
    ax_hist.axvline(median, color="orange", linestyle="--", linewidth=1.3,
                    label=f"median = {median:.2f} Å")
    ax_hist.set_xlabel("Cα RMSD to crystal structure (Å)")
    ax_hist.set_ylabel("Number of conformations")
    ax_hist.set_title(f"RMSD distribution (n={len(rmsds)})", fontsize=11)
    ax_hist.legend()

    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def plot_summary(all_rmsds, out_dir):
    """Write box-per-monomer and aggregated-histogram summary plots."""
    sorted_ids = sorted(all_rmsds.keys(), key=lambda k: np.median(all_rmsds[k]))
    values = [all_rmsds[k] for k in sorted_ids]

    width = max(8, len(sorted_ids) * 0.22)
    fig, ax = plt.subplots(figsize=(width, 5))
    ax.boxplot(values, showfliers=False, widths=0.6)
    ax.set_xticks(range(1, len(sorted_ids) + 1))
    ax.set_xticklabels(sorted_ids, rotation=90, fontsize=7)
    ax.set_ylabel("Cα RMSD to crystal structure (Å)")
    ax.set_title(
        f"RMSD distribution per monomer "
        f"(n={len(sorted_ids)} monomers, sorted by median)"
    )
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "summary_per_monomer.png"), dpi=120)
    plt.close(fig)

    flat = np.concatenate([np.asarray(v) for v in values])
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(flat, bins=60, edgecolor="black", alpha=0.75, color="steelblue")
    ax.axvline(float(np.mean(flat)), color="red", linestyle="--",
               label=f"mean = {np.mean(flat):.2f} Å")
    ax.axvline(float(np.median(flat)), color="orange", linestyle="--",
               label=f"median = {np.median(flat):.2f} Å")
    ax.set_xlabel("Cα RMSD to crystal structure (Å)")
    ax.set_ylabel("Number of conformations")
    ax.set_title(
        f"Aggregated RMSD distribution "
        f"(n={len(flat)} conformations across {len(sorted_ids)} monomers)"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "summary_aggregated.png"), dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    monomers_dir = os.path.join(script_dir, "monomers")
    out_dir = os.path.join(script_dir, "rmsd")
    plots_dir = os.path.join(out_dir, "plots")
    structures_dir = os.path.join(out_dir, "structures")
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(structures_dir, exist_ok=True)

    requested = {x.upper() for x in sys.argv[1:]}
    chains_dirs = sorted(glob.glob(os.path.join(monomers_dir, "*_chains")))
    if requested:
        chains_dirs = [
            d for d in chains_dirs
            if os.path.basename(d).replace("_chains", "") in requested
        ]
    print(f"Processing {len(chains_dirs)} monomer directory(ies) from {monomers_dir}/")

    all_rmsds = {}
    csv_path = os.path.join(out_dir, "rmsds.csv")
    with open(csv_path, "w", newline="") as csvf:
        writer = csv.writer(csvf)
        writer.writerow(["pdb_id", "state", "rmsd_angstrom"])

        for chains_dir in chains_dirs:
            pdb_id = os.path.basename(chains_dir).replace("_chains", "")
            proteins = sorted(glob.glob(os.path.join(chains_dir, "*_chain*_protein.pdb")))
            if not proteins:
                print(f"[skip] {pdb_id}: no reference protein PDB")
                continue

            dna_files = glob.glob(os.path.join(chains_dir, "*_dna.pdb"))
            dna_pdb = dna_files[0] if dna_files else None

            structure_img = os.path.join(structures_dir, f"{pdb_id}.png")
            if not os.path.exists(structure_img):
                try:
                    render_structure_image(proteins[0], dna_pdb, structure_img)
                except Exception as e:
                    print(f"[warn] {pdb_id}: structure render failed ({e})")
                    structure_img = None

            rmsds = compute_rmsds_for_monomer(chains_dir, proteins[0])
            if rmsds is None:
                print(f"[skip] {pdb_id}: missing topology.pdb or samples.xtc")
                continue

            all_rmsds[pdb_id] = rmsds
            for i, r in enumerate(rmsds, start=1):
                writer.writerow([pdb_id, i, f"{r:.4f}"])

            plot_per_monomer(
                pdb_id, rmsds, structure_img,
                os.path.join(plots_dir, f"{pdb_id}.png"),
            )
            print(
                f"[{pdb_id}] n={len(rmsds)} "
                f"mean={np.mean(rmsds):.2f}Å median={np.median(rmsds):.2f}Å "
                f"max={np.max(rmsds):.2f}Å"
            )

    if all_rmsds:
        plot_summary(all_rmsds, out_dir)
        print(f"\nWrote {csv_path}")
        print(f"Wrote {len(all_rmsds)} per-monomer plots to {plots_dir}/")
        print(f"Wrote summary plots to {out_dir}/")
    else:
        print("No complete monomer runs found.")
