#!/usr/bin/env python
"""
compare_to_pipeline.py — is a BioEmu2 complex ensemble better than stages 1-3?

Scores a BioEmu2 protein-DNA ensemble and the existing pipeline's minimized
states against the SAME crystal reference with the SAME interface metrics, so the
two are comparable rather than merely both reported.

Everything measured here is measured by existing code:
  * fnat / interface RMSD   fnat_gate/interface_rmsd.py  (ref_side, score)
  * effect size (CLES)      analysis/analyses/nmr_rmsd/compute_pair_rmsd.py:cles
  * paired stats, BH-FDR    analysis/common/stats.py
Nothing statistical is reimplemented here; this script only arranges the inputs
and reports.

The question it answers
----------------------
The pipeline's frozen arm pairs every protein frame with crystal-pose DNA, so its
fnat is high partly by construction. BioEmu2 samples the DNA too, so its fnat is
an honest measurement of whether the model reproduces the native interface. A
LOWER BioEmu2 fnat is therefore not automatically worse - read it together with
the DNA-geometry columns and the caveat printed at the end.

Usage
-----
    python bioemu2/compare_to_pipeline.py --tf-name ets1
    python bioemu2/compare_to_pipeline.py --tf-name ets1 --arms frozen dnarelax
    python bioemu2/compare_to_pipeline.py --tf-name ets1 --out-csv /tmp/cmp.csv
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_TFCONF = _HERE.parent

# Reuse, don't reimplement.
sys.path.insert(0, str(_TFCONF / "fnat_gate"))
sys.path.insert(0, str(_TFCONF / "analysis" / "common"))
import interface_rmsd as ir                                    # noqa: E402
import score_stage3 as s3                                      # noqa: E402
from stats import bh_fdr                                       # noqa: E402

try:
    sys.path.insert(0, str(_TFCONF / "analysis" / "analyses" / "nmr_rmsd"))
    from compute_pair_rmsd import cles                          # noqa: E402
    _HAVE_CLES = True
except Exception:                                               # noqa: BLE001
    _HAVE_CLES = False

try:
    from scipy.stats import mannwhitneyu
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False


def pilot_env(tf_name):
    cfg = _TFCONF / "config" / "pilots" / f"{tf_name}.sh"
    if not cfg.exists():
        sys.exit(f"ERROR: no pilot config at {cfg}")
    r = subprocess.run(["bash", "-c", f'set -a; source "{cfg}" >/dev/null 2>&1; env'],
                       capture_output=True, text=True, cwd=_TFCONF)
    return dict(l.split("=", 1) for l in r.stdout.splitlines() if "=" in l)


def dna_count_mismatches(R, state_dir):
    """How many states have a DNA residue count != the reference?

    This is a pre-check, not a rescoring. It matters because
    interface_rmsd.score() silently falls back to the REFERENCE DNA when the
    counts disagree (finding F10 of docs/codebase_review_2026-10-15.md), which
    would quietly turn a model-frame fnat into a cross-frame one. The existing
    pipeline carries crystal DNA so it never trips this; a BioEmu2 ensemble
    samples the DNA and can, so it is counted and reported rather than silent.
    """
    import glob
    n = 0
    files = sorted(glob.glob(os.path.join(state_dir, "*_state_*.pdb")))
    for f in files:
        try:
            _, _, md_dna = ir.ordered(ir.load_models(f)[0])
            if len(md_dna) != len(R["rd"]):
                n += 1
        except Exception:                                       # noqa: BLE001
            continue
    return n, len(files)


def score_states(R, state_dir, contact_cut):
    """Score a directory of states using the EXISTING gate, unchanged.

    fnat_gate/score_stage3.py:score_dir already globs <pdb>_state_NNN.pdb,
    parses the state index from the filename (deliberately, so the label cannot
    drift from the file), and scores each with use_model_dna=True. That is
    exactly what is wanted here, so it is called rather than reimplemented -
    which also guarantees these numbers are produced the same way the pipeline's
    own fnat gate produces them.
    """
    _pdb, by_state = s3.score_dir(R, state_dir, contact_cut)
    return sorted(by_state.items())


def summarize(name, rows, metric):
    vals = np.array([r[1][metric] for r in rows
                     if r[1].get(metric) is not None
                     and not np.isnan(r[1][metric])], dtype=float)
    if vals.size == 0:
        return None
    return dict(name=name, metric=metric, n=int(vals.size),
                median=float(np.median(vals)), mean=float(vals.mean()),
                p25=float(np.percentile(vals, 25)),
                p75=float(np.percentile(vals, 75)),
                vals=vals)


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tf-name", required=True)
    p.add_argument("--bioemu2-dir", default=None,
                   help="default: output/bioemu2/<pdb>_complex/states")
    p.add_argument("--arms", nargs="*", default=["frozen"],
                   choices=["frozen", "dnarelax"],
                   help="which pipeline arms to compare against")
    p.add_argument("--iface-cutoff", type=float, default=5.0)
    p.add_argument("--contact-cutoff", type=float, default=4.5)
    p.add_argument("--gap", type=int, default=4)
    p.add_argument("--out-csv", default=None)
    a = p.parse_args()

    env = pilot_env(a.tf_name)
    pdb_id = (env.get("PDB_ID") or "").lower()
    ref_cif = _TFCONF / "structures" / "source_chains" / f"{pdb_id}_chains" / f"{pdb_id}.cif"
    if not ref_cif.exists():
        sys.exit(f"ERROR: reference cif not found: {ref_cif}")

    prot = env.get("PROTEIN_CHAIN")
    dna = env.get("DNA_CHAINS")
    protein_chain = int(str(prot).split("#")[0].strip()) if prot else None
    dna_chains = [int(x) for x in dna.split(",") if x.strip()] if dna else None

    print(f"=== {a.tf_name} ({pdb_id}) ===")
    print(f"reference : {ref_cif}")
    print(f"chains    : protein={protein_chain} dna={dna_chains}")

    # Reference interface, defined once by the existing gate's own logic.
    R = ir.ref_side(ir.load_models(str(ref_cif))[0], a.iface_cutoff,
                    a.contact_cutoff, a.gap,
                    protein_chain=protein_chain, dna_chains=dna_chains)
    print(f"interface : {len(R['iface'])} residues, {len(R['segs'])} segments, "
          f"{len(R['native'])} native contacts, {len(R['rd'])} DNA residues")

    sets = []

    b2 = Path(a.bioemu2_dir) if a.bioemu2_dir else \
        _TFCONF / "output" / "bioemu2" / f"{pdb_id}_complex" / "states"
    if b2.is_dir():
        print(f"\nscoring BioEmu2 complex states: {b2}")
        mism, n_files = dna_count_mismatches(R, str(b2))
        rows = score_states(R, str(b2), a.contact_cutoff)
        print(f"  {len(rows)} states scored (of {n_files} files)")
        if mism:
            print(f"  WARNING: {mism}/{len(rows)} states have a DNA residue count "
                  f"!= the reference ({len(R['rd'])}). interface_rmsd.score() "
                  f"silently falls back to the REFERENCE DNA for those, so their "
                  f"fnat is NOT a model-frame measurement. See finding F10.",
                  file=sys.stderr)
        sets.append(("bioemu2_complex", rows))
    else:
        print(f"\nNOTE: no BioEmu2 states at {b2}")
        print(f"      run: python bioemu2/sample_complex.py --tf-name {a.tf_name}")

    for arm in a.arms:
        sfx = "" if arm == "frozen" else "_dnarelax"
        d = _TFCONF / "output" / f"stage3_min{sfx}" / a.tf_name
        if not d.is_dir():
            print(f"\nNOTE: no pipeline states at {d}")
            continue
        print(f"\nscoring pipeline states ({arm}): {d}")
        # The frozen arm's DNA *is* the crystal DNA, so model-frame and
        # reference-frame fnat coincide there; use_model_dna=True is still the
        # right call because it is what fnat_gate itself uses for stage 3.
        mism, n_files = dna_count_mismatches(R, str(d))
        rows = score_states(R, str(d), a.contact_cutoff)
        print(f"  {len(rows)} states scored (of {n_files} files)")
        if mism:
            print(f"  WARNING: {mism}/{len(rows)} DNA-count mismatches (F10)",
                  file=sys.stderr)
        sets.append((f"pipeline_{arm}", rows))

    if len(sets) < 1:
        sys.exit("nothing to compare")

    # ---- report ---------------------------------------------------------
    METRICS = [("fnat", "higher is better"), ("g", "iRMSD global, lower better"),
               ("smax", "iRMSD worst segment, lower better"),
               ("ident", "seq identity to reference")]
    print("\n" + "=" * 78)
    print(f"{'set':<20}{'metric':<8}{'n':>5}{'median':>10}{'IQR':>18}  note")
    print("-" * 78)
    summaries = {}
    for name, rows in sets:
        for metric, note in METRICS:
            s = summarize(name, rows, metric)
            if not s:
                continue
            summaries[(name, metric)] = s
            iqr = f"[{s['p25']:.3f}, {s['p75']:.3f}]"
            print(f"{name:<20}{metric:<8}{s['n']:>5}{s['median']:>10.3f}"
                  f"{iqr:>18}  {note}")

    # ---- head-to-head ---------------------------------------------------
    if ("bioemu2_complex", "fnat") in summaries:
        print("\n" + "=" * 78)
        print("BioEmu2 complex vs pipeline, per metric")
        print("-" * 78)
        pvals, labels = [], []
        for arm in a.arms:
            other = f"pipeline_{arm}"
            for metric, _ in METRICS:
                k1, k2 = ("bioemu2_complex", metric), (other, metric)
                if k1 not in summaries or k2 not in summaries:
                    continue
                x, y = summaries[k1]["vals"], summaries[k2]["vals"]
                line = (f"  {metric:<8} vs {arm:<9} "
                        f"median {np.median(x):+.3f} vs {np.median(y):+.3f}")
                if _HAVE_CLES:
                    line += f"  CLES={cles(x, y):.2f}"
                if _HAVE_SCIPY:
                    u, pv = mannwhitneyu(x, y, alternative="two-sided")
                    pvals.append(pv); labels.append(f"{metric} vs {arm}")
                    line += f"  p={pv:.3g}"
                print(line)
        if pvals:
            qs = bh_fdr(np.array(pvals))
            print("\n  BH-FDR across these comparisons "
                  "(the metrics are correlated, so this is a floor):")
            for lab, pv, q in zip(labels, pvals, qs):
                print(f"    {lab:<24} p={pv:.3g}  q={q:.3g}")

    # Effective n: states of one ensemble are not independent structures.
    for name, rows in sets:
        print(f"\n  {name}: {len(rows)} states from 1 structure - these are "
              f"correlated samples of one system, not {len(rows)} independent "
              f"observations. Any p-value above is a within-system comparison.")
        break

    if a.out_csv:
        import csv
        with open(a.out_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["set", "state"] + ir.COLS.split(",")[2:])
            for name, rows in sets:
                for state, m in rows:
                    w.writerow([name, state, m["n_iface"], m["n_seg"],
                                ir.fmt(m["ident"]), ir.fmt(m["g"]),
                                ir.fmt(m["smax"]), ir.fmt(m["smean"]),
                                ir.fmt(m["fnat"])])
        print(f"\nwrote {a.out_csv}")

    print("\n" + "=" * 78)
    print("READING THIS HONESTLY")
    print("-" * 78)
    print("The frozen pipeline arm pairs every frame with CRYSTAL-POSE DNA, so its")
    print("fnat is high partly by construction - the DNA was never free to move.")
    print("BioEmu2 samples the DNA as well, so its fnat is an actual measurement.")
    print("A lower BioEmu2 fnat is therefore not automatically a worse model; it")
    print("is a different, harder question being asked. What would be decisive is")
    print("BioEmu2 reproducing the native interface while ALSO reproducing the")
    print("crystal's DNA deformation - that is what stages 1-3 cannot test at all.")
    print("")
    print("Also: BioEmu2's model card states sampling frequency 'must not be")
    print("interpreted as a calibrated Boltzmann probability'. Do not read these")
    print("state counts as populations.")


if __name__ == "__main__":
    main()
