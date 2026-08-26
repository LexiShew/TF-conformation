#!/usr/bin/env python3
"""Pivot pfams_maybe_final.csv so each PDB_ID has a single row.

PWM_ID, UniProt, TF Class and TF Family are split into source-specific
columns: a "_Jaspar" suffix when the PWM_ID contains "jaspar" and a
"_Hocomoco" suffix when the PWM_ID contains "H11MO". When a PDB_ID has
more than one row for the same source, the values are joined with "|".
"""

import argparse
import os

import pandas as pd

# Columns that get split per PWM source -> output column base name.
SPLIT_COLS = {
    "PWM_ID": "PWM_ID",
    "UniProt_ID": "Uniprot",
    "TF_Class": "TF Class",
    "TF_Family": "TF Family",
}

# Columns kept once per PDB_ID (assumed consistent across rows of a PDB).
SHARED_COLS = ["Chain", "Pfam_ID", "Pfam_Name", "pfam_manual"]


def normalize_pfam_manual(value: str) -> str:
    """Split on commas, trim each token, sort lexicographically, rejoin.

    e.g. "Pou, Homeodomain" -> "Homeodomain | Pou".
    """
    if not isinstance(value, str):
        return value
    tokens = [t.strip() for t in value.split(",")]
    tokens = [t for t in tokens if t]
    return " | ".join(sorted(tokens))


def source_of(pwm_id: str) -> str | None:
    """Classify a PWM_ID as 'Jaspar', 'Hocomoco', or None."""
    if not isinstance(pwm_id, str):
        return None
    if "jaspar" in pwm_id:
        return "Jaspar"
    if "H11MO" in pwm_id:
        return "Hocomoco"
    return None


def pivot(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["_source"] = df["PWM_ID"].map(source_of)

    unclassified = df[df["_source"].isna()]
    if not unclassified.empty:
        ids = unclassified["PWM_ID"].unique().tolist()
        print(f"Warning: {len(unclassified)} row(s) with unclassified "
              f"PWM_ID dropped: {ids}")
        df = df[df["_source"].notna()]

    rows = []
    for pdb_id, group in df.groupby("PDB_ID", sort=False):
        row = {"PDB_ID": pdb_id}

        for col in SHARED_COLS:
            series = group[col].dropna().astype(str)
            if col == "pfam_manual":
                series = series.map(normalize_pfam_manual)
            row[col] = "|".join(series.unique().tolist())

        for source in ("Jaspar", "Hocomoco"):
            sub = group[group["_source"] == source]
            for src_col, out_base in SPLIT_COLS.items():
                out_col = f"{out_base}_{source}"
                vals = sub[src_col].dropna().astype(str).unique().tolist()
                row[out_col] = "|".join(vals)

        rows.append(row)

    ordered_cols = ["PDB_ID"] + SHARED_COLS
    for out_base in SPLIT_COLS.values():
        ordered_cols += [f"{out_base}_Jaspar", f"{out_base}_Hocomoco"]

    return pd.DataFrame(rows)[ordered_cols]


def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i", "--input",
        default=os.path.join(here, "pfams_maybe_final.csv"),
        help="Input CSV (default: pfams_maybe_final.csv next to this script)",
    )
    parser.add_argument(
        "-o", "--output",
        default=os.path.join(here, "pivoted_pfams.csv"),
        help="Output CSV (default: pivoted_pfams.csv next to this script)",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    pivoted = pivot(df)
    pivoted.to_csv(args.output, index=False)
    print(f"Wrote {len(pivoted)} PDB_ID rows to {args.output}")


if __name__ == "__main__":
    main()
