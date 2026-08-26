#!/usr/bin/env python3
"""Merge rmsds.csv with pivoted_pfams.csv into rmsds_with_pfam.csv.

Each rmsds.csv row (one per PDB_ID/state, ~100 states per PDB_ID) is kept.
The original `pfam_family` column is replaced by `pfam_manual` from
pivoted_pfams.csv, and the `TF Class_Jaspar` and `TF Family_Jaspar`
columns are appended.
"""

import argparse
import os

import pandas as pd

PFAM_COLS = ["pfam_manual", "TF Class_Jaspar", "TF Family_Jaspar"]


def merge(rmsds: pd.DataFrame, pivoted: pd.DataFrame) -> pd.DataFrame:
    info = pivoted[["PDB_ID"] + PFAM_COLS].rename(columns={"PDB_ID": "pdb_id"})

    missing = sorted(set(rmsds["pdb_id"]) - set(info["pdb_id"]))
    if missing:
        print(f"Warning: {len(missing)} PDB_ID(s) in rmsds.csv have no "
              f"pivoted_pfams entry (left blank): {missing}")

    merged = rmsds.merge(info, on="pdb_id", how="left")

    # Replace pfam_family in place with pfam_manual, keep column order.
    out_cols = []
    for col in rmsds.columns:
        if col == "pfam_family":
            out_cols.append("pfam_manual")
        else:
            out_cols.append(col)
    out_cols += [c for c in PFAM_COLS if c != "pfam_manual"]
    return merged[out_cols]


def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rmsds", default=os.path.join(here, "rmsds.csv"),
        help="Input rmsds.csv (default: next to this script)",
    )
    parser.add_argument(
        "--pivoted", default=os.path.join(here, "pivoted_pfams.csv"),
        help="Input pivoted_pfams.csv (default: next to this script)",
    )
    parser.add_argument(
        "-o", "--output", default=os.path.join(here, "rmsds_with_pfam.csv"),
        help="Output CSV (default: rmsds_with_pfam.csv next to this script)",
    )
    args = parser.parse_args()

    rmsds = pd.read_csv(args.rmsds)
    pivoted = pd.read_csv(args.pivoted)
    merged = merge(rmsds, pivoted)
    merged.to_csv(args.output, index=False)
    print(f"Wrote {len(merged)} rows to {args.output}")


if __name__ == "__main__":
    main()
