# Render all RMSD plots from rmsds_with_pfam.csv (written by
# merge_pfam_rmsd_info.py). This CSV is the only input; everything the
# plots need lives in its columns:
#   pdb_id, state, rmsd_angstrom, pfam_manual, TF Class_Jaspar, TF Family_Jaspar
#
# Outputs (under rmsd_analysis/plots/):
#   indiv_rmsd_plots/<PDB>.png — per-monomer RMSD histogram
#   summary_per_monomer.png    — box-per-monomer, sorted by median RMSD
#   summary_aggregated.png     — pooled RMSD histogram across all monomers
#   summary_by_family.png      — box per pfam_manual family
#   summary_by_tf_class.png    — box per TF Class_Jaspar
#   summary_by_tf_family.png   — box per TF Family_Jaspar
#
# Usage:
#   python plot_rmsd.py
#   python plot_rmsd.py 1abc 2def              # filter to these PDB IDs
#   python plot_rmsd.py --min-monomers 3       # require ≥3 monomers per group

import os
import sys
import csv
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# CSV column -> (output filename, human-readable noun) for the grouped plots.
GROUP_COLUMNS = [
    ("pfam_manual", "summary_by_family.png", "pfam_manual family"),
    ("TF Class_Jaspar", "summary_by_tf_class.png", "TF class (Jaspar)"),
    ("TF Family_Jaspar", "summary_by_tf_family.png", "TF family (Jaspar)"),
]


def load_rmsds_csv(csv_path, requested):
    """Return (rmsds_by_pdb, meta_by_pdb) loaded from rmsds_with_pfam.csv.

    rmsds_by_pdb: pdb_id -> list[float]
    meta_by_pdb:  pdb_id -> {column_name: label_str} for every GROUP_COLUMNS
                  column (blank cells become "" and are treated as
                  "Unclassified" downstream).
    """
    group_cols = [c for c, _, _ in GROUP_COLUMNS]
    rmsds_by_pdb = {}
    meta_by_pdb = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            pid = row["pdb_id"]
            if requested and pid.upper() not in requested:
                continue
            rmsds_by_pdb.setdefault(pid, []).append(float(row["rmsd_angstrom"]))
            if pid not in meta_by_pdb:
                meta_by_pdb[pid] = {
                    c: (row.get(c) or "").strip() for c in group_cols
                }
    return rmsds_by_pdb, meta_by_pdb


def plot_per_monomer(pdb_id, rmsds, meta, out_path):
    fig, ax_hist = plt.subplots(figsize=(7, 4.5))

    subtitle_bits = [meta.get("pfam_manual") or "Unclassified"]
    tf_class = meta.get("TF Class_Jaspar")
    if tf_class:
        subtitle_bits.append(tf_class)

    ax_hist.hist(rmsds, bins=30, edgecolor="black", alpha=0.75, color="steelblue")
    mean = float(np.mean(rmsds))
    median = float(np.median(rmsds))
    best = float(np.min(rmsds))
    ax_hist.axvline(mean, color="red", linestyle="--", linewidth=1.3,
                    label=f"mean = {mean:.2f} Å")
    ax_hist.axvline(median, color="orange", linestyle="--", linewidth=1.3,
                    label=f"median = {median:.2f} Å")
    ax_hist.axvline(best, color="green", linestyle="--", linewidth=1.3,
                    label=f"best = {best:.2f} Å")
    ax_hist.set_xlabel("Cα RMSD to crystal structure (Å)")
    ax_hist.set_xlim(left=0)
    ax_hist.set_ylabel("Number of conformations")
    ax_hist.set_title(
        f"{pdb_id}  ·  {' · '.join(subtitle_bits)}\n"
        f"RMSD distribution (n={len(rmsds)})",
        fontsize=11,
    )
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


def plot_by_group(all_rmsds, label_by_pdb, group_noun, out_path,
                  min_monomers=2):
    """Box-per-group plot. label_by_pdb maps pdb_id -> group label."""
    rmsds_by_group = {}
    monomers_by_group = {}
    for pdb_id, rmsds in all_rmsds.items():
        group = label_by_pdb.get(pdb_id) or "Unclassified"
        rmsds_by_group.setdefault(group, []).extend(rmsds)
        monomers_by_group.setdefault(group, set()).add(pdb_id)

    groups = [
        g for g, mons in monomers_by_group.items() if len(mons) >= min_monomers
    ]
    if not groups:
        print(
            f"[warn] no {group_noun} groups with >= {min_monomers} monomers; "
            f"skipping {os.path.basename(out_path)}"
        )
        return False

    groups.sort(key=lambda g: np.median(rmsds_by_group[g]))
    values = [rmsds_by_group[g] for g in groups]
    labels = [f"{g}\n(n={len(monomers_by_group[g])})" for g in groups]

    width = max(8, len(groups) * 0.6)
    fig, ax = plt.subplots(figsize=(width, 5.5))
    ax.boxplot(values, showfliers=False, widths=0.6)
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Cα RMSD to crystal structure (Å)")
    ax.set_title(
        f"RMSD distribution per {group_noun} "
        f"(groups with >={min_monomers} monomers, sorted by median)"
    )
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot RMSD distributions from rmsds_with_pfam.csv.",
    )
    parser.add_argument(
        "pdb_ids", nargs="*",
        help="Restrict to these PDB IDs (default: all monomers in the CSV).",
    )
    parser.add_argument(
        "--rmsds-csv",
        default=os.path.join(SCRIPT_DIR, "rmsds_with_pfam.csv"),
        help="Path to the merged CSV "
             "(default: rmsds_with_pfam.csv next to this script).",
    )
    parser.add_argument(
        "--output-dir",
        default=os.path.join(SCRIPT_DIR, "plots"),
        help="Where to write PNGs (default: rmsd_analysis/plots).",
    )
    parser.add_argument(
        "--min-monomers", type=int, default=2,
        help="Only include groups with ≥N monomers in the grouped plots "
             "(default: 2).",
    )
    args = parser.parse_args()

    csv_path = args.rmsds_csv
    if not os.path.exists(csv_path):
        sys.exit(f"{csv_path} not found — run merge_pfam_rmsd_info.py first.")

    out_dir = os.path.abspath(args.output_dir)
    indiv_dir = os.path.join(out_dir, "indiv_rmsd_plots")
    os.makedirs(indiv_dir, exist_ok=True)

    requested = {x.upper() for x in args.pdb_ids}
    all_rmsds, meta_by_pdb = load_rmsds_csv(csv_path, requested)
    if not all_rmsds:
        sys.exit("No matching monomers in CSV.")
    print(
        f"Loaded {sum(len(v) for v in all_rmsds.values())} RMSDs "
        f"for {len(all_rmsds)} monomer(s) from {csv_path}"
    )

    for pdb_id, rmsds in all_rmsds.items():
        plot_per_monomer(
            pdb_id, rmsds, meta_by_pdb.get(pdb_id, {}),
            os.path.join(indiv_dir, f"{pdb_id}.png"),
        )

    plot_summary(all_rmsds, out_dir)
    print(f"Wrote {len(all_rmsds)} per-monomer plots to {indiv_dir}/")
    print(f"Wrote summary plots to {out_dir}/")

    for col, fname, noun in GROUP_COLUMNS:
        label_by_pdb = {p: m.get(col, "") for p, m in meta_by_pdb.items()}
        if any(label_by_pdb.values()):
            out_path = os.path.join(out_dir, fname)
            if plot_by_group(all_rmsds, label_by_pdb, noun, out_path,
                             min_monomers=args.min_monomers):
                print(f"Wrote {out_path}")
        else:
            print(f"{col} column is empty in CSV — skipping {fname}.")
