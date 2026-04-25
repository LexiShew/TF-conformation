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
#   python compute_rmsds.py                             # compute + plot, all monomers
#   python compute_rmsds.py 1CIT 6PAX                   # compute + plot, these IDs only
#   python compute_rmsds.py --compute-only              # compute RMSDs + write CSV, no plots
#   python compute_rmsds.py --plot-only                 # skip compute; plot from rmsd/rmsds.csv
#   python compute_rmsds.py --plot-only 1CIT 6PAX       # plot from existing CSV, filtered

import os
import sys
import glob
import csv
import argparse

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
    ax_hist.set_xlim(left=0)
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
    ax.set_xlim(left=0)
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


def ensure_structure_image(chains_dir, structures_dir, pdb_id):
    """Return the path to a cached structure PNG, rendering it if missing."""
    out = os.path.join(structures_dir, f"{pdb_id}.png")
    if os.path.exists(out):
        return out
    proteins = sorted(glob.glob(os.path.join(chains_dir, "*_chain*_protein.pdb")))
    if not proteins:
        return None
    dna_files = glob.glob(os.path.join(chains_dir, "*_dna.pdb"))
    dna_pdb = dna_files[0] if dna_files else None
    try:
        render_structure_image(proteins[0], dna_pdb, out)
        return out
    except Exception as e:
        print(f"[warn] {pdb_id}: structure render failed ({e})")
        return None


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


def write_rmsds_csv(results, csv_path, preserve_others=False):
    """Write results to csv_path. If preserve_others, keep rows for pdb_ids not in results."""
    existing_rows = []
    if preserve_others and os.path.exists(csv_path):
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                if row["pdb_id"] not in results:
                    existing_rows.append(row)
    with open(csv_path, "w", newline="") as csvf:
        w = csv.writer(csvf)
        w.writerow(["pdb_id", "state", "rmsd_angstrom"])
        for row in existing_rows:
            w.writerow([row["pdb_id"], row["state"], row["rmsd_angstrom"]])
        for pdb_id in sorted(results):
            for i, r in enumerate(results[pdb_id], start=1):
                w.writerow([pdb_id, i, f"{r:.4f}"])


def load_rmsds_csv(csv_path, requested):
    """Load RMSDs from CSV, filtered to requested pdb_ids if non-empty."""
    results = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            pid = row["pdb_id"]
            if requested and pid not in requested:
                continue
            results.setdefault(pid, []).append(float(row["rmsd_angstrom"]))
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute and/or plot Cα RMSDs vs. crystal structure.",
    )
    parser.add_argument(
        "pdb_ids", nargs="*",
        help="Restrict to these PDB IDs (default: all monomers with a trajectory)",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--compute-only", action="store_true",
        help="Compute RMSDs and write rmsd/rmsds.csv; skip all plots.",
    )
    mode.add_argument(
        "--plot-only", action="store_true",
        help="Skip computation; re-draw plots from existing rmsd/rmsds.csv.",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    monomers_dir = os.path.join(script_dir, "monomers")
    out_dir = os.path.join(script_dir, "rmsd")
    plots_dir = os.path.join(out_dir, "plots")
    structures_dir = os.path.join(out_dir, "structures")
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(structures_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "rmsds.csv")

    requested = {x.upper() for x in args.pdb_ids}

    if args.plot_only:
        if not os.path.exists(csv_path):
            sys.exit(f"{csv_path} not found — run without --plot-only first.")
        all_rmsds = load_rmsds_csv(csv_path, requested)
        print(
            f"Loaded {sum(len(v) for v in all_rmsds.values())} RMSDs "
            f"for {len(all_rmsds)} monomer(s) from {csv_path}"
        )
    else:
        chains_dirs = sorted(glob.glob(os.path.join(monomers_dir, "*_chains")))
        if requested:
            chains_dirs = [
                d for d in chains_dirs
                if os.path.basename(d).replace("_chains", "") in requested
            ]
        print(f"Processing {len(chains_dirs)} monomer directory(ies) from {monomers_dir}/")
        all_rmsds = compute_all(chains_dirs)
        # When computing a subset, preserve rows for other monomers already in the CSV.
        write_rmsds_csv(all_rmsds, csv_path, preserve_others=bool(requested))
        print(f"Wrote {csv_path}")

    if args.compute_only:
        sys.exit(0)
    if not all_rmsds:
        print("No data to plot.")
        sys.exit(0)

    for pdb_id, rmsds in all_rmsds.items():
        chains_dir = os.path.join(monomers_dir, f"{pdb_id}_chains")
        structure_img = ensure_structure_image(chains_dir, structures_dir, pdb_id)
        plot_per_monomer(
            pdb_id, rmsds, structure_img,
            os.path.join(plots_dir, f"{pdb_id}.png"),
        )

    plot_summary(all_rmsds, out_dir)
    print(f"Wrote {len(all_rmsds)} per-monomer plots to {plots_dir}/")
    print(f"Wrote summary plots to {out_dir}/")
