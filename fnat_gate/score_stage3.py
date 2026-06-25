#!/usr/bin/env python3
"""
score_stage3.py — score a directory of single-model state PDBs against a crystal
reference, with the state label parsed straight from the filename (so it can
never drift from the file the way per-file --tag did).

Globs <PDB>_state_<NNN>.pdb, scores each with interface_rmsd.score (fnat vs the
model's OWN DNA — correct frame for docked/minimized states), writes one CSV row
per state (schema matches batch_redock_score.py), prints a summary, and
optionally writes a pass-list of states clearing an fnat floor.

With --compare-dir, also scores a second directory (e.g. stage2_docked) and
prints a PAIRED comparison by state index — use it to check whether Stage 3
minimization is improving or loosening the interface relative to Stage 2.

Usage:
  python score_stage3.py --ref 1vtn.cif --dir .../foxa/stage3_min --floor 0.5
  python score_stage3.py --ref 1vtn.cif --dir .../foxa/stage3_min \
      --compare-dir .../foxa/stage2_docked

Store next to interface_rmsd.py (it imports it). Requires: biopython, numpy.
"""
import argparse, os, re, sys, glob
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import interface_rmsd as ir

STATE_RE = re.compile(r"_state_(\d+)", re.IGNORECASE)


def score_dir(R, dirpath, contact_cutoff, pdb_id=None):
    """Return (pdb_id, {state:int -> metrics dict}) for *_state_*.pdb in dirpath."""
    files = glob.glob(os.path.join(dirpath, "*_state_*.pdb"))
    items = []
    for f in files:
        m = STATE_RE.search(os.path.basename(f))
        if m:
            items.append((int(m.group(1)), f))
        else:
            print(f"skip (no _state_NNN): {os.path.basename(f)}", file=sys.stderr)
    items.sort(key=lambda x: x[0])
    if not items:
        return None, {}
    pdb = pdb_id or os.path.basename(items[0][1]).split("_state_")[0]
    out = {}
    for state, f in items:
        try:
            out[state] = ir.score(R, ir.load_models(f)[0], contact_cutoff, use_model_dna=True)
        except Exception as e:
            print(f"{dirpath} state {state}: SCORE_ERROR {e}", file=sys.stderr)
    return pdb, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="crystal reference (.cif/.pdb)")
    ap.add_argument("--dir", required=True, help="primary dir of <PDB>_state_*.pdb (e.g. stage3_min)")
    ap.add_argument("--compare-dir", help="second dir to compare against (e.g. stage2_docked)")
    ap.add_argument("--pdb-id", help="defaults to the prefix before _state_ in the first file")
    ap.add_argument("--out", help="CSV path (default <dir>/<pdb>_fnat.csv)")
    ap.add_argument("--floor", type=float, default=None,
                    help="if set, also write a pass-list of states with fnat >= floor")
    ap.add_argument("--pass-out", help="pass-list path (default <dir>/<pdb>_pass.txt)")
    ap.add_argument("--iface-cutoff", type=float, default=5.0)
    ap.add_argument("--contact-cutoff", type=float, default=4.5)
    ap.add_argument("--gap", type=int, default=4)
    a = ap.parse_args()

    R = ir.ref_side(ir.load_models(a.ref)[0], a.iface_cutoff, a.contact_cutoff, a.gap)

    pdb, primary = score_dir(R, a.dir, a.contact_cutoff, a.pdb_id)
    if not primary:
        sys.exit(f"no *_state_*.pdb found in {a.dir}")

    out = a.out or os.path.join(a.dir, f"{pdb}_fnat.csv")
    with open(out, "w") as fh:
        fh.write(ir.COLS + "\n")
        for s in sorted(primary):
            fh.write(ir.row(pdb, s, primary[s]) + "\n")
    print(f"wrote {out}  ({len(primary)} states)")

    pf = np.array([primary[s]["fnat"] for s in sorted(primary) if primary[s]["fnat"] == primary[s]["fnat"]])
    ps = np.array([primary[s]["smax"] for s in sorted(primary) if primary[s]["smax"] == primary[s]["smax"]])
    floor = a.floor if a.floor is not None else 0.5
    if len(pf):
        print(f"[{os.path.basename(a.dir.rstrip('/'))}] fnat median={np.median(pf):.3f} "
              f"min={pf.min():.3f} max={pf.max():.3f} | seg_max median={np.median(ps):.3f} "
              f"| pct fnat>={floor:g}: {(pf>=floor).mean()*100:.1f}%")

    if a.floor is not None:
        pass_out = a.pass_out or os.path.join(a.dir, f"{pdb}_pass.txt")
        passing = [s for s in sorted(primary) if primary[s]["fnat"] >= a.floor]
        with open(pass_out, "w") as fh:
            for s in passing:
                fh.write(f"{pdb}_state_{s:03d}\n")
        print(f"wrote {pass_out}  ({len(passing)}/{len(primary)} states >= {a.floor})")

    # ---------- paired comparison ----------
    if a.compare_dir:
        _, comp = score_dir(R, a.compare_dir, a.contact_cutoff, pdb)
        common = sorted(set(primary) & set(comp))
        if not common:
            print("compare: no common state indices", file=sys.stderr); return
        pfat = np.array([primary[s]["fnat"] for s in common])
        cfat = np.array([comp[s]["fnat"] for s in common])
        psm  = np.array([primary[s]["smax"] for s in common])
        csm  = np.array([comp[s]["smax"]  for s in common])
        d = pfat - cfat
        A = os.path.basename(a.dir.rstrip("/"))
        B = os.path.basename(a.compare_dir.rstrip("/"))
        print(f"\n--- paired comparison over {len(common)} common states ---")
        print(f"                  {A:>16}  {B:>16}")
        print(f"fnat median       {np.median(pfat):16.3f}  {np.median(cfat):16.3f}")
        print(f"seg_max median    {np.median(psm):16.3f}  {np.median(csm):16.3f}")
        print(f"mean Δfnat ({A} - {B}): {d.mean():+.3f}   median Δ: {np.median(d):+.3f}")
        print(f"states where {A} > {B}: {(d>0).sum()}/{len(common)}")
        print(f"legend: Δ<0 ⇒ {A} has LOWER fnat (e.g. minimization loosened the interface)")


if __name__ == "__main__":
    main()