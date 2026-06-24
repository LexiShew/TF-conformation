#!/usr/bin/env python3
"""Link each PDB id to its JASPAR / HOCOMOCO PWM label(s).

Source of truth is the DeepPBS cluster-wise dataset
(jaspar_h11mo_cluster_wise_dna_containing_dataset.npy): an array of clusters,
each a list of [<pdb>_<chain>, [<pwm_label>, ...]] entries. PWM labels are
either HOCOMOCO ("*.H11MO.*") or JASPAR ("MA*.jaspar").

We aggregate per 4-char PDB id over all of its chains, split the labels into
jaspar/hocomoco columns, and write one row per id in the input list.

Usage:
  python map_pdb_to_pwm.py [--ids structures/deeppbs_pdb_ids.txt]
                           [--npy stage4_preprocess/jaspar_h11mo_cluster_wise_dna_containing_dataset.npy]
                           [--out config/deeppbs_pdb_pwms.csv]
"""
import argparse
import csv
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))  # scripts/pdb_prep -> repo root


def is_jaspar(label):
    return label.lower().endswith(".jaspar")


def is_hocomoco(label):
    return ".H11MO." in label.upper()


def build_pdb_to_pwms(npy_path):
    """pdb_id (lowercased 4-char) -> sorted list of distinct PWM labels."""
    data = np.load(npy_path, allow_pickle=True)
    mapping = {}
    for cluster in data:
        for entry in cluster:
            key, labels = entry[0], entry[1]
            pdb = key.split("_")[0].lower()
            mapping.setdefault(pdb, set()).update(labels)
    return {k: sorted(v) for k, v in mapping.items()}


def build_pdb_to_idtxt_pwms(idtxt_path):
    """pdb_id -> sorted PWM labels actually used by the benchmark id.txt entries.

    id.txt rows are <pdb>_<chain>_<PWM>.npz; the PWM itself may contain
    underscores (HOCOMOCO names), so we take everything after the 2nd '_'. This
    is the label a pilot config's PWM_LABEL must match — and it can differ from
    the cluster-dataset label (e.g. 1tgh is TBP_HUMAN.H11MO.0.A in the .npy but
    MA0343.1.jaspar in id.txt)."""
    mapping = {}
    if not idtxt_path or not os.path.exists(idtxt_path):
        return mapping
    with open(idtxt_path) as f:
        for line in f:
            name = line.strip()
            if not name:
                continue
            if name.endswith(".npz"):
                name = name[:-4]
            parts = name.split("_")
            if len(parts) < 3:
                continue
            pdb = parts[0].lower()
            pwm = "_".join(parts[2:])
            mapping.setdefault(pdb, set()).add(pwm)
    return {k: sorted(v) for k, v in mapping.items()}


def build_cluster_benchmark(npy_path, pdb_to_idtxt):
    """pdb_id -> sorted benchmark (id.txt) PWM labels found among the OTHER PDBs
    in the same npy TF-cluster.

    A pilot's reference PDB is usually not itself a benchmark entry; the held-out
    test entries are other structures of the same TF, which share its npy cluster.
    This surfaces the label(s) those cluster-mates use in id.txt — the candidate
    PWM_LABEL for a pilot built on this PDB. CAVEAT: clusters can be broad and span
    PWM variants, so treat this as candidates, not a single authoritative answer."""
    data = np.load(npy_path, allow_pickle=True)
    out = {}
    for cluster in data:
        pdbs = sorted({e[0].split("_")[0].lower() for e in cluster})
        labels = sorted({lbl for p in pdbs for lbl in pdb_to_idtxt.get(p, [])})
        for p in pdbs:
            out.setdefault(p, set()).update(labels)
    return {k: sorted(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", default=os.path.join(REPO, "structures", "deeppbs_pdb_ids.txt"))
    ap.add_argument("--npy", default=os.path.join(
        REPO, "stage4_preprocess", "jaspar_h11mo_cluster_wise_dna_containing_dataset.npy"))
    ap.add_argument("--out", default=os.path.join(REPO, "config", "deeppbs_pdb_pwms.csv"))
    ap.add_argument("--idtxt", default=os.path.join(
        REPO, "stage5_build_aug", "folds", "id.txt"),
        help="Benchmark id.txt; adds an id_txt_pwm column (the label a pilot "
             "config's PWM_LABEL must match). Pass '' to skip.")
    args = ap.parse_args()

    pdb_to_pwms = build_pdb_to_pwms(args.npy)
    pdb_to_idtxt = build_pdb_to_idtxt_pwms(args.idtxt)
    pdb_to_cluster_bench = build_cluster_benchmark(args.npy, pdb_to_idtxt)

    with open(args.ids) as f:
        ids = [line.strip().lower() for line in f if line.strip()]

    n_with, n_jaspar, n_hoco, n_idtxt = 0, 0, 0, 0
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pdb_id", "jaspar", "hocomoco", "id_txt_pwm", "in_benchmark",
                    "cluster_benchmark_pwms", "n_pwms", "all_pwms"])
        for pdb in ids:
            pwms = pdb_to_pwms.get(pdb, [])
            jaspar = [p for p in pwms if is_jaspar(p)]
            hoco = [p for p in pwms if is_hocomoco(p)]
            idtxt = pdb_to_idtxt.get(pdb, [])
            cluster_bench = pdb_to_cluster_bench.get(pdb, [])
            if pwms:
                n_with += 1
            n_jaspar += bool(jaspar)
            n_hoco += bool(hoco)
            n_idtxt += bool(idtxt)
            w.writerow([pdb, ";".join(jaspar), ";".join(hoco), ";".join(idtxt),
                        "yes" if idtxt else "no", ";".join(cluster_bench),
                        len(pwms), ";".join(pwms)])

    print(f"Wrote {args.out}")
    print(f"  {len(ids)} PDB ids; {n_with} mapped to >=1 PWM "
          f"({len(ids) - n_with} unmapped)")
    print(f"  {n_jaspar} have a JASPAR label, {n_hoco} have a HOCOMOCO label")
    print(f"  {n_idtxt} are in the benchmark id.txt (have an id_txt_pwm)")


if __name__ == "__main__":
    main()
