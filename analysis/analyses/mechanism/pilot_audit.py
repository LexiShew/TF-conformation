#!/usr/bin/env python
"""
pilot_audit.py -- which pilots are present at each stage of the mechanism analysis, and why any are missing.

The pilot roster has grown over the project (6 -> 10 -> 12 -> 13 configs). Different
analyses legitimately cover different subsets, and a silent drop is easy to mistake for
a result. This script prints the pilot set at every stage the mechanism analysis touches
and names the reason for each exclusion, so any n in a downstream figure can be traced.

Usage
-----
    source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
    cd /project2/rohs_102/shewchuk/TF-conformation
    python analysis/mechanism/pilot_audit.py

Output
------
    analysis/mechanism/data/pilot_coverage.csv   pilot x stage presence matrix
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / "data"
OUT.mkdir(parents=True, exist_ok=True)


def show(label, pilots):
    p = sorted(set(pilots))
    print("{:44s} n={:2d}  {}".format(label, len(p), " ".join(p)))
    return p


def main():
    stages = {}

    stages["config: base pilot files"] = show(
        "config: base pilot files",
        [f.stem for f in (ROOT / "config/pilots").glob("*.sh")
         if not f.stem.endswith("_dnarelax")])

    stages["stage7: standard eval JSON"] = show(
        "stage7: standard eval JSON",
        [f.stem.replace("id_benchmark_", "")
         for f in (ROOT / "output/stage7_eval").glob("id_benchmark_*.json")
         if "_dnarelax" not in f.stem])

    stages["stage7: dnarelax eval JSON"] = show(
        "stage7: dnarelax eval JSON",
        [f.stem.replace("id_benchmark_", "").replace("_dnarelax", "")
         for f in (ROOT / "output/stage7_eval").glob("id_benchmark_*_dnarelax.json")])

    stages["reachability.csv (protein axis)"] = show(
        "reachability.csv (protein axis)",
        pd.read_csv(ROOT / "analysis/data/reachability.csv").pilot)

    ps = pd.read_csv(ROOT / "analysis/data/perseed_summary.csv")
    stages["perseed_summary frozen"] = show(
        "perseed_summary frozen", ps.query("dna=='frozen'").tf)
    stages["perseed_summary relaxed"] = show(
        "perseed_summary relaxed", ps.query("dna=='relaxed'").tf)

    pe = pd.read_csv(ROOT / "analysis/data/perseed_perentry.csv",
                     usecols=["tf", "selffam"])
    stages["perseed_perentry: has own-family entries"] = show(
        "perseed_perentry: has own-family entries", pe.query("selffam").tf)

    pc = pd.read_csv(ROOT / "analysis/analyses/dna_relax/data/pycurves_all_perstructure.csv")
    stages["pycurves: crystal DNA geometry"] = show(
        "pycurves: crystal DNA geometry", pc.query("cond=='crystal'").tf)
    stages["pycurves: frozen ensemble"] = show(
        "pycurves: frozen ensemble", pc.query("cond=='frozen'").tf)

    stages["mgw_fl_summary (whole molecule)"] = show(
        "mgw_fl_summary (whole molecule)",
        pd.read_csv(ROOT / "analysis/analyses/dna_relax/data/mgw_fl_summary.csv").pilot)
    stages["iface_mgwfl (interface restricted)"] = show(
        "iface_mgwfl (interface restricted)",
        pd.read_csv(ROOT / "analysis/analyses/dna_relax/data/iface_mgwfl_vs_accuracy.csv").pilot)

    m = pd.read_csv(OUT / "mechanism_table.csv")
    stages["mechanism_table: all rows"] = show("mechanism_table: all rows", m.tf)
    stages["mechanism_table: has dPearson"] = show(
        "mechanism_table: has dPearson", m.dropna(subset=["froz_dP_mean"]).tf)

    sf = pd.read_csv(OUT / "selffam_effects.csv")
    stages["own-family effect (frozen)"] = show(
        "own-family effect (frozen)", sf.query("dna=='frozen' and selffam").tf)

    universe = sorted(set().union(*stages.values()))
    mat = pd.DataFrame({s: [p in v for p in universe] for s, v in stages.items()},
                       index=universe)
    mat.index.name = "pilot"
    mat.to_csv(OUT / "pilot_coverage.csv")

    print("\n=== pilots missing from any stage ===")
    for p in universe:
        missing = [s for s, v in stages.items() if p not in v]
        if missing:
            print("  {:10s} absent from: {}".format(p, "; ".join(missing)))


if __name__ == "__main__":
    main()
