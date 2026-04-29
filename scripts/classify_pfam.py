# Classify each protein chain in monomers/ against Pfam using HMMER (hmmscan),
# write one row per chain to rmsd/pfam_classifications.csv, then (if
# rmsd/rmsds.csv exists) plot the Cα RMSD distribution grouped by Pfam family.
#
# Requirements:
#   - HMMER 3 with `hmmscan` on PATH.
#   - Pfam-A.hmm, already pressed with `hmmpress` so .h3f/.h3i/.h3m/.h3p files
#     sit next to it. Pass its path via --pfam-db or $PFAM_DB.
#   - Optional: rmsd/rmsds.csv from compute_rmsds.py (needed only for the
#     by-family box plot).
#
# Usage:
#   export PFAM_DB=/path/to/Pfam-A.hmm
#   python classify_pfam.py                       # classify + plot
#   python classify_pfam.py --skip-plot           # classify only
#   python classify_pfam.py --min-monomers 3      # require ≥3 monomers/family
#   python classify_pfam.py --cpu 8

import os
import sys
import glob
import csv
import shutil
import subprocess
import argparse
import tempfile
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymol import cmd


# HMMER accepts the 20 standard AAs plus a few ambiguity codes; PyMOL emits '?'
# for residues it can't map to a single-letter code (MSE, modified residues,
# UNK, etc.). Anything outside this set gets rewritten to 'X' ("any residue")
# so hmmscan doesn't reject the input with "illegal character".
_HMMER_OK_AA = set("ACDEFGHIKLMNPQRSTVWYBZJUOX")


def extract_fasta_record(pdb_file, name):
    """Load a PDB in PyMOL and return a single-record FASTA string (or '')."""
    cmd.delete("all")
    cmd.load(pdb_file, "prot")
    fasta = cmd.get_fastastr("prot and polymer.protein")
    seq = "".join(ln for ln in fasta.splitlines() if not ln.startswith(">"))
    seq = "".join(c if c.upper() in _HMMER_OK_AA else "X" for c in seq)
    cmd.delete("all")
    return f">{name}\n{seq}\n" if seq else ""


def build_combined_fasta(monomers_dir, fasta_path):
    """Write one multi-FASTA covering every chain PDB under monomers_dir.

    Returns a dict mapping FASTA record name ``<PDB>_<CHAIN>`` to the tuple
    ``(pdb_id, chain_letter, pdb_path)``.
    """
    chain_map = {}
    with open(fasta_path, "w") as f:
        for chains_dir in sorted(glob.glob(os.path.join(monomers_dir, "*_chains"))):
            pdb_id = os.path.basename(chains_dir).replace("_chains", "")
            for pdb_file in sorted(
                glob.glob(os.path.join(chains_dir, "*_chain*_protein.pdb"))
            ):
                base = os.path.basename(pdb_file)
                # <PDB>_chain<X>_protein.pdb — chain letter between "_chain" and "_"
                chain_letter = base.split("_chain", 1)[1].split("_", 1)[0]
                name = f"{pdb_id}_{chain_letter}"
                record = extract_fasta_record(pdb_file, name)
                if record:
                    f.write(record)
                    chain_map[name] = (pdb_id, chain_letter, pdb_file)
    return chain_map


def run_hmmscan(pfam_db, fasta_path, out_tbl, cpu=4):
    """Run hmmscan against Pfam using the curated gathering threshold."""
    subprocess.run(
        [
            "hmmscan",
            "--cut_ga",
            "--cpu", str(cpu),
            "--domtblout", out_tbl,
            pfam_db,
            fasta_path,
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )


def parse_domtblout(path):
    """Yield dict per domain hit from hmmscan --domtblout output."""
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split(maxsplit=22)
            if len(parts) < 23:
                continue
            yield {
                "pfam_name": parts[0],
                "pfam_acc": parts[1],
                "query": parts[3],
                "full_evalue": float(parts[6]),
                "full_score": float(parts[7]),
                "i_evalue": float(parts[12]),
                "description": parts[22].strip(),
            }


def best_hit_per_chain(tbl_path):
    """Return dict query_name -> top hit (lowest full-sequence E-value)."""
    best = {}
    for hit in parse_domtblout(tbl_path):
        key = hit["query"]
        if key not in best or hit["full_evalue"] < best[key]["full_evalue"]:
            best[key] = hit
    return best


def write_classifications_csv(chain_map, best_hits, out_path):
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "pdb_id", "chain",
            "pfam_name", "pfam_acc", "pfam_description",
            "full_evalue", "bit_score",
        ])
        for name in sorted(chain_map):
            pdb_id, chain, _ = chain_map[name]
            hit = best_hits.get(name)
            if hit:
                w.writerow([
                    pdb_id, chain,
                    hit["pfam_name"], hit["pfam_acc"], hit["description"],
                    f"{hit['full_evalue']:.2e}", f"{hit['full_score']:.1f}",
                ])
            else:
                w.writerow([pdb_id, chain, "Unclassified", "", "", "", ""])


def plot_rmsd_by_family(classifications_csv, rmsds_csv, out_path, min_monomers=2):
    """Box plot of Cα RMSD grouped by the primary (first) chain's Pfam family."""
    primary_family = {}
    with open(classifications_csv) as f:
        for row in csv.DictReader(f):
            # rows are sorted by chain letter; keep the first one per pdb_id
            primary_family.setdefault(row["pdb_id"], row["pfam_name"])

    family_rmsds = defaultdict(list)
    family_monomers = defaultdict(set)
    with open(rmsds_csv) as f:
        for row in csv.DictReader(f):
            fam = primary_family.get(row["pdb_id"], "Unclassified")
            family_rmsds[fam].append(float(row["rmsd_angstrom"]))
            family_monomers[fam].add(row["pdb_id"])

    keep = {
        fam: vals for fam, vals in family_rmsds.items()
        if len(family_monomers[fam]) >= min_monomers
    }
    if not keep:
        print(f"No Pfam families with ≥{min_monomers} monomers; skipping family plot.")
        return

    sorted_fams = sorted(keep, key=lambda k: np.median(keep[k]))
    values = [keep[k] for k in sorted_fams]
    labels = [f"{k}\n(n={len(family_monomers[k])})" for k in sorted_fams]

    width = max(8, len(sorted_fams) * 0.55)
    fig, ax = plt.subplots(figsize=(width, 6))
    ax.boxplot(values, showfliers=False, widths=0.6)
    ax.set_xticks(range(1, len(sorted_fams) + 1))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Cα RMSD to crystal structure (Å)")
    ax.set_title(
        f"RMSD distribution by Pfam family "
        f"(≥{min_monomers} monomers/family, sorted by median)"
    )
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)

    # counts summary
    total_monomers = sum(len(family_monomers[f]) for f in sorted_fams)
    unclassified = len(family_monomers.get("Unclassified", set()))
    print(
        f"  {len(sorted_fams)} families plotted, "
        f"{total_monomers} monomers covered "
        f"({unclassified} unclassified)."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pfam-db", default=os.environ.get("PFAM_DB"),
                        help="Path to pressed Pfam-A.hmm (default: $PFAM_DB)")
    parser.add_argument("--skip-plot", action="store_true",
                        help="Skip the RMSD-by-family box plot")
    parser.add_argument("--min-monomers", type=int, default=2,
                        help="Only plot families with ≥N monomers (default: 2)")
    parser.add_argument("--cpu", type=int, default=4,
                        help="Threads for hmmscan (default: 4)")
    args = parser.parse_args()

    if not shutil.which("hmmscan"):
        sys.exit("hmmscan not found on PATH — install HMMER 3.")
    if not args.pfam_db:
        sys.exit("Pfam database required: --pfam-db or $PFAM_DB.")
    if not os.path.exists(args.pfam_db):
        sys.exit(f"Pfam database not found: {args.pfam_db}")
    if not os.path.exists(args.pfam_db + ".h3f"):
        sys.exit(
            f"{args.pfam_db} has not been pressed. "
            f"Run: hmmpress {args.pfam_db}"
        )

    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    monomers_dir = os.path.join(project_root, "monomers")
    out_dir = os.path.join(project_root, "rmsd")
    os.makedirs(out_dir, exist_ok=True)
    classifications_csv = os.path.join(out_dir, "pfam_classifications.csv")

    with tempfile.TemporaryDirectory() as wd:
        fasta_path = os.path.join(wd, "chains.fasta")
        tbl_path = os.path.join(wd, "pfam.domtblout")

        print("Extracting protein sequences from chain PDBs...")
        chain_map = build_combined_fasta(monomers_dir, fasta_path)
        print(f"  {len(chain_map)} chains → {fasta_path}")

        print(f"Running hmmscan (cpu={args.cpu}) against {args.pfam_db}...")
        run_hmmscan(args.pfam_db, fasta_path, tbl_path, cpu=args.cpu)

        best = best_hit_per_chain(tbl_path)
        print(f"  {len(best)}/{len(chain_map)} chains matched ≥1 Pfam family")

    write_classifications_csv(chain_map, best, classifications_csv)
    print(f"Wrote {classifications_csv}")

    if args.skip_plot:
        sys.exit(0)

    rmsds_csv = os.path.join(out_dir, "rmsds.csv")
    if not os.path.exists(rmsds_csv):
        print(f"No {rmsds_csv} — run compute_rmsds.py first to enable the family plot.")
        sys.exit(0)

    plot_path = os.path.join(out_dir, "summary_by_family.png")
    plot_rmsd_by_family(
        classifications_csv, rmsds_csv, plot_path,
        min_monomers=args.min_monomers,
    )
    print(f"Wrote {plot_path}")
