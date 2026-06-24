# Render all RMSD plots from <MONOMER_DIR>/rmsd/rmsds.csv (written by
# compute_rmsds.py) and the cached structure PNGs in
# <MONOMER_DIR>/rmsd/structures/ (written by generate_structure_images.py).
#
# Outputs (under <MONOMER_DIR>/rmsd/):
#   plots/<PDB>.png            — structure overlay (left) + RMSD histogram (right)
#   summary_per_monomer.png    — box-per-monomer, sorted by median RMSD
#   summary_aggregated.png     — pooled RMSD histogram across all monomers
#   summary_by_family.png      — box-per-Pfam-family (only when the CSV's
#                                pfam_family column is populated)
#
# Usage:
#   python plot_rmsd.py <MONOMER_DIR>
#   python plot_rmsd.py <MONOMER_DIR> 1abc 2def              # filter to these PDB IDs
#   python plot_rmsd.py <MONOMER_DIR> --min-monomers 3       # require ≥3 monomers per family

import os
import sys
import csv
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec


def load_rmsds_csv(csv_path, requested):
    """Return (rmsds_by_pdb, family_by_pdb) loaded from the CSV.

    rmsds_by_pdb: pdb_id -> list[float]
    family_by_pdb: pdb_id -> str (may be "" if column blank)
    """
    rmsds_by_pdb = {}
    family_by_pdb = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            pid = row["pdb_id"]
            if requested and pid.upper() not in requested:
                continue
            rmsds_by_pdb.setdefault(pid, []).append(float(row["rmsd_angstrom"]))
            family_by_pdb.setdefault(pid, (row.get("pfam_family") or "").strip())
    return rmsds_by_pdb, family_by_pdb


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


def plot_by_family(all_rmsds, family_by_pdb, out_path, min_monomers=2):
    rmsds_by_family = {}
    monomers_by_family = {}
    for pdb_id, rmsds in all_rmsds.items():
        family = family_by_pdb.get(pdb_id) or "Unclassified"
        rmsds_by_family.setdefault(family, []).extend(rmsds)
        monomers_by_family.setdefault(family, set()).add(pdb_id)

    families = [
        f for f, mons in monomers_by_family.items() if len(mons) >= min_monomers
    ]
    if not families:
        print(
            f"[warn] no Pfam families with >= {min_monomers} monomers; "
            f"skipping by-family plot"
        )
        return False

    families.sort(key=lambda f: np.median(rmsds_by_family[f]))
    values = [rmsds_by_family[f] for f in families]
    labels = [f"{f}\n(n={len(monomers_by_family[f])})" for f in families]

    width = max(8, len(families) * 0.6)
    fig, ax = plt.subplots(figsize=(width, 5))
    ax.boxplot(values, showfliers=False, widths=0.6)
    ax.set_xticks(range(1, len(families) + 1))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Cα RMSD to crystal structure (Å)")
    ax.set_title(
        f"RMSD distribution per Pfam family "
        f"(families with >={min_monomers} monomers, sorted by median)"
    )
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot RMSD distributions from rmsds.csv + cached structure PNGs.",
    )
    parser.add_argument(
        "monomer_dir",
        help="Directory containing rmsd/rmsds.csv and rmsd/structures/.",
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
        "--min-monomers", type=int, default=2,
        help="Only include Pfam families with ≥N monomers in the by-family plot (default: 2).",
    )
    args = parser.parse_args()

    monomers_dir = os.path.abspath(args.monomer_dir)
    out_dir = os.path.join(monomers_dir, "rmsd")
    plots_dir = os.path.join(out_dir, "plots")
    structures_dir = os.path.join(out_dir, "structures")
    os.makedirs(plots_dir, exist_ok=True)
    csv_path = args.rmsds_csv or os.path.join(out_dir, "rmsds.csv")
    if not os.path.exists(csv_path):
        sys.exit(f"{csv_path} not found — run compute_rmsds.py first.")

    requested = {x.upper() for x in args.pdb_ids}
    all_rmsds, family_by_pdb = load_rmsds_csv(csv_path, requested)
    if not all_rmsds:
        sys.exit("No matching monomers in CSV.")
    print(
        f"Loaded {sum(len(v) for v in all_rmsds.values())} RMSDs "
        f"for {len(all_rmsds)} monomer(s) from {csv_path}"
    )

    for pdb_id, rmsds in all_rmsds.items():
        structure_img = os.path.join(structures_dir, f"{pdb_id}.png")
        if not os.path.exists(structure_img):
            print(
                f"[warn] {pdb_id}: no structure image at {structure_img} — "
                f"run generate_structure_images.py to create it"
            )
            structure_img = None
        plot_per_monomer(
            pdb_id, rmsds, structure_img,
            os.path.join(plots_dir, f"{pdb_id}.png"),
        )

    plot_summary(all_rmsds, out_dir)
    print(f"Wrote {len(all_rmsds)} per-monomer plots to {plots_dir}/")
    print(f"Wrote summary plots to {out_dir}/")

    if any(family_by_pdb.values()):
        plot_path = os.path.join(out_dir, "summary_by_family.png")
        if plot_by_family(all_rmsds, family_by_pdb, plot_path, min_monomers=args.min_monomers):
            print(f"Wrote {plot_path}")
    else:
        print(
            "pfam_family column is empty in CSV — skipping by-family plot. "
            "Re-run compute_rmsds.py with --pfam-metadata to populate it."
        )
