#!/usr/bin/env python
"""Regenerate pilot_mechanism.csv from fnat_pilots_summary.csv.

This table previously had no writer in the tree: fnat_on_pilots.py emits only
fnat_pilots_summary.csv, and the column `nmr_pass` appears in no script. The
rule below was reconstructed from the summary and verified to reproduce all
three derived columns exactly for all five pilots.

Columns
-------
pilot     pilot id
n_apo     number of deposited apo ENTRIES (not models)
apo_med   median over the per-entry median fnats, all apo entries
apo_pass  mean pass rate over the CRYSTAL apo entries, as a percentage.
          Falls back to the NMR entries when a pilot has no crystal apo
          (nfat is the only such case).
nmr_pass  mean pass rate over the NMR apo bundles, as a percentage.
          NaN when a pilot has no NMR apo entry.
bm        BioEmu median fnat
bm_pass   BioEmu pass rate, as a percentage
bm_max    BioEmu maximum fnat over all frames

Why the apo_pass/nmr_pass split matters: it is the whole content of the runx
reclassification. RUNX has six crystal apo structures at fnat 0.667-0.875
(100% pass) and two NMR bundles, 1CMO (43 models) and 1CO1 (10), at 0.250 and
0.375 (0% pass). Its free state reads selection-like by crystallography and
induced-fit-like by NMR -- a 0.459 fnat gap by experimental method. runx is the
only pilot where the two disagree.

Usage:
    python regenerate_pilot_mechanism.py [--summary PATH] [--out PATH] [--check]
"""
import argparse
import csv
import statistics
import sys
from pathlib import Path


def is_nmr(method):
    return "NMR" in (method or "").upper()


def build(summary_path):
    rows = list(csv.DictReader(open(summary_path)))
    pilots = sorted({r["pilot"] for r in rows})
    out = []
    for p in pilots:
        g = [r for r in rows if r["pilot"] == p]
        apo = [r for r in g if r["source"] != "bioemu"]
        bm = [r for r in g if r["source"] == "bioemu"]
        if not apo or not bm:
            print(f"  skip {p}: apo={len(apo)} bioemu={len(bm)}", file=sys.stderr)
            continue
        bm = bm[0]
        cry = [r for r in apo if not is_nmr(r["method"])]
        nmr = [r for r in apo if is_nmr(r["method"])]
        pass_src = cry if cry else nmr          # crystal preferred; NMR fallback
        out.append({
            "pilot":    p,
            "n_apo":    len(apo),
            "apo_med":  round(statistics.median([float(r["fnat_median"]) for r in apo]), 4),
            "apo_pass": round(100 * statistics.fmean([float(r["pass_rate"]) for r in pass_src]), 2),
            "nmr_pass": (round(100 * statistics.fmean([float(r["pass_rate"]) for r in nmr]), 2)
                         if nmr else ""),
            "bm":       round(float(bm["fnat_median"]), 4),
            "bm_pass":  round(100 * float(bm["pass_rate"]), 2),
            "bm_max":   round(float(bm["fnat_max"]), 4),
        })
    return out


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=str(here / "fnat_pilots" / "fnat_pilots_summary.csv"))
    ap.add_argument("--out", default=str(here / "fnat_pilots" / "pilot_mechanism.csv"))
    ap.add_argument("--check", action="store_true",
                    help="compare against the existing file instead of overwriting")
    args = ap.parse_args()

    rows = build(args.summary)
    if not rows:
        sys.exit("nothing built")

    if args.check:
        old = {r["pilot"]: r for r in csv.DictReader(open(args.out))}
        bad = 0
        for r in rows:
            o = old.get(r["pilot"])
            if o is None:
                print(f"MISSING in existing: {r['pilot']}"); bad += 1; continue
            for k, v in r.items():
                ov = o.get(k, "")
                same = (str(v) == str(ov)) or (
                    v not in ("", None) and ov not in ("", None)
                    and abs(float(v) - float(ov)) < 0.011)
                if not same:
                    print(f"DIFF {r['pilot']}.{k}: regenerated={v} existing={ov}"); bad += 1
        print("MATCH" if not bad else f"{bad} differences")
        return

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out} ({len(rows)} pilots)")


if __name__ == "__main__":
    main()
