#!/usr/bin/env python
"""
The apo-fnat mechanism read-out, applied to the pipeline's own pilots.

On the NMR pairs, docking a DEPOSITED APO structure onto the holo DNA and
scoring fnat split the set cleanly: proteins whose real free state reproduces
the bound interface (selection-like) from those whose does not (induced-fit-
like). That read-out needs no BioEmu at all, which is what makes it useful --
Ch. 17 currently infers mechanism from the SIGN of the augmentation effect, and
this measures it independently.

Here the same thing is done for the pilots. The reference is the pilot's own
crystal (config/pilots/<tf>.sh: PDB_ID, PROTEIN_CHAIN, DNA_CHAINS), and the
docked source is an apo structure of the same protein found by 95% sequence
cluster. Three sources, all through the pipeline's own Stage-2 dock:

  apo      the deposited apo structure(s). NMR bundles give a pass RATE over
           models; a crystal apo gives a single fnat.
  bioemu   the pilot's own Stage-1 ensemble, on the same residues, so the
           BioEmu number is comparable to the apo one for that pilot.
  (the holo ceiling is trivially 1.0 here -- unlike the NMR case there is only
  ONE deposited holo structure and it IS the reference, so it is not scored.)

Everything for a pilot is trimmed to one residue set -- the pilot's reference
protein intersected with every apo candidate and with the BioEmu topology -- so
the native contact set is identical across sources. Stage 2 also hard-exits on
a Ca-count mismatch, so this trimming is required, not optional.

Outputs (in --output-dir, default ./fnat_pilots):
  refs/<pilot>_ref.pdb          trimmed docking target (protein + DNA)
  <pilot>_<source>_fnat.csv     per-state fnat
  fnat_pilots_summary.csv       per (pilot, source): median fnat, pass rate

Usage:
    python fnat_on_pilots.py --apo-json <pilot_apo.json>
    python fnat_on_pilots.py --pilots ets1 runx
"""
import argparse
import csv
import difflib
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from pair_core import residue_sequence                                  # noqa: E402
from compute_nmr_state_rmsd import load_reference, select_reference_chain  # noqa: E402
from fnat_on_pairs import _clean_serials, _match, _seq, DNA_RES, run     # noqa: E402

_TFCONF = _HERE.parents[2]
STAGE2 = _TFCONF / "stage2_redock" / "stage2_redock.py"
SCORE = _TFCONF / "fnat_gate" / "score_stage3.py"
PILOTS = _TFCONF / "config" / "pilots"
SRC_CHAINS = _TFCONF / "structures" / "source_chains"
STAGE1 = _TFCONF / "output" / "stage1_bioemu"


def read_pilot(name):
    """PDB_ID / BINDING_CHAIN / PROTEIN_CHAIN / DNA_CHAINS from the pilot config."""
    # Only real assignment lines: several configs discuss alternative values in
    # comments (ets1's warns 'if A actually binds E/F, set DNA_CHAINS="2,3"'),
    # and a whole-file regex happily picks the commented one up.
    lines = [l.split("#", 1)[0] for l in (PILOTS / f"{name}.sh").read_text().splitlines()]
    text = "\n".join(l for l in lines if l.strip())

    def g(k, quoted=True):
        m = re.search(rf'^\s*(?:export\s+)?{k}="([^"]+)"' if quoted
                      else rf"^\s*(?:export\s+){k}=(\d+)", text, re.M)
        return m.group(1) if m else None
    return dict(name=name, pdb=g("PDB_ID"), binding_chain=g("BINDING_CHAIN"),
                protein_chain=g("PROTEIN_CHAIN", False), dna_chains=g("DNA_CHAINS"))


def fetch_cif(pdb_id, cache):
    import urllib.request
    cache.mkdir(parents=True, exist_ok=True)
    dest = cache / f"{pdb_id.lower()}.cif"
    if not dest.is_file() or not dest.stat().st_size:
        urllib.request.urlretrieve(
            f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", dest)
    return dest


class BrokenReference(RuntimeError):
    """The scoring reference has no interface, so nothing can be scored against it."""


def _fnat_values(path):
    """fnat column as floats, plus whether EVERY state was unscorable.

    The gate writes "NA" when fnat is nan, and interface_rmsd.score returns nan
    only when the REFERENCE has no native protein-DNA contacts at all
    (`fnat = rec/len(nat) if nat else nan`). So NA is a property of the
    reference, not of the docked state: if one state is NA they all are, and it
    means the reference protein/DNA chain selection is wrong -- not that the
    ensemble scored zero. Silently dropping those rows hides a broken reference,
    so the caller is told.
    """
    out, total = [], 0
    for r in csv.DictReader(open(path)):
        total += 1
        try:
            v = float(r["fnat"])
        except (TypeError, ValueError):
            continue
        if v == v:
            out.append(v)
    if total and not out:
        raise BrokenReference(
            f"{Path(path).name}: all {total} states scored NA -- the reference has no "
            f"protein-DNA contacts, so its chain selection is wrong")
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apo-json", required=True,
                    help="JSON list of {pilot, apo, chain, method, models} candidates")
    ap.add_argument("--pilots", nargs="*", default=None)
    ap.add_argument("--output-dir", default=_HERE / "fnat_pilots")
    ap.add_argument("--floor", type=float, default=0.5)
    ap.add_argument("--min-frac-length", type=float, default=0.85,
                    help="drop apo candidates shorter than this fraction of the pilot "
                         "protein: a partial construct (2P81 is a 44-residue fragment of "
                         "a 60-residue homeodomain) shrinks the shared core for every "
                         "other source and scores badly for reasons unrelated to apo/holo")
    ap.add_argument("--python", default=sys.executable)
    args = ap.parse_args()

    import mdtraj as md
    cand = json.load(open(args.apo_json))
    by_pilot = {}
    for c in cand:
        by_pilot.setdefault(c["pilot"], []).append(c)
    if args.pilots:
        by_pilot = {k: v for k, v in by_pilot.items() if k in set(args.pilots)}

    out = Path(args.output_dir)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    cache = out / "cif"
    rows = []

    for pilot, apos in sorted(by_pilot.items()):
        cfg = read_pilot(pilot)
        print(f"\n{pilot}  (crystal {cfg['pdb']}, protein chainid {cfg['protein_chain']}, "
              f"DNA {cfg['dna_chains']})")
        ref_cif = SRC_CHAINS / f"{cfg['pdb']}_chains" / f"{cfg['pdb']}.cif"
        if not ref_cif.is_file():
            print(f"  SKIP: no {ref_cif}", file=sys.stderr)
            continue
        ref = load_reference(ref_cif, verbose=False)[0]

        chains = list(ref.top.chains)
        pc = int(cfg["protein_chain"])
        ref_res = [(r, aa) for r, aa in residue_sequence(ref.top)
                   if r.chain.index == pc]
        if not ref_res:
            print(f"  SKIP: chainid {pc} has no protein residues", file=sys.stderr)
            continue
        dna_idx = [a.index for c in cfg["dna_chains"].split(",")
                   for a in chains[int(c)].atoms
                   if a.residue.name.strip().upper() in DNA_RES]
        if not dna_idx:
            print(f"  SKIP: chainids {cfg['dna_chains']} hold no DNA", file=sys.stderr)
            continue

        # BioEmu ensemble for the pilot's binding chain
        stem = f"{cfg['pdb']}_chain{cfg['binding_chain']}"
        d = STAGE1 / f"{stem}_conformations"
        top = next((p for p in (d / f"{stem}_topology.pdb", d / "topology.pdb")
                    if p.is_file()), None)
        xtc = next((p for p in (d / f"{stem}_samples.xtc", d / "samples.xtc")
                    if p.is_file()), None)
        ens = md.load(str(xtc), top=str(top)) if (top and xtc) else None

        # one shared residue set: reference INTERSECT every apo INTERSECT BioEmu
        keep_ref = ref_res
        loaded = {}
        n_ref = len(ref_res)
        short = [c for c in apos if (c.get("L") or 0) < args.min_frac_length * n_ref]
        if short:
            print(f"  dropping {len(short)} partial construct(s): "
                  + ", ".join(f"{c['apo']}({c['L']}aa)" for c in short))
        apos = [c for c in apos if c not in short]
        for c in apos:
            try:
                traj = load_reference(fetch_cif(c["apo"], cache), verbose=False)
                ares, _, _ = select_reference_chain(traj.top, c["chain"], _seq(ref_res))
            except Exception as exc:                    # noqa: BLE001
                print(f"  {c['apo']}: SKIP load ({str(exc)[:60]})", file=sys.stderr)
                continue
            loaded[c["apo"]] = (traj, ares, c)
            keep_ref, _ = _match(_seq(keep_ref), keep_ref, _seq(ares), ares)
        if ens is not None:
            eres = residue_sequence(ens.top)
            keep_ref, _ = _match(_seq(keep_ref), keep_ref, _seq(eres), eres)
        if not loaded or len(keep_ref) < 4:
            print(f"  SKIP: {len(keep_ref)} residues shared across sources", file=sys.stderr)
            continue
        print(f"  shared core: {len(keep_ref)}/{len(ref_res)} reference residues "
              f"across {len(loaded)} apo + BioEmu")

        keep = {r.index for r, _ in keep_ref}
        prot_idx = [a.index for a in ref.top.atoms if a.residue.index in keep]
        combined = _clean_serials(ref.atom_slice(prot_idx).stack(ref.atom_slice(dna_idx)))
        for n, ch in enumerate(combined.top.chains):
            ch.chain_id = chr(ord("A") + n)
        ref_pdb = out / "refs" / f"{pilot}_ref.pdb"
        ref_pdb.parent.mkdir(parents=True, exist_ok=True)
        combined.save_pdb(str(ref_pdb))
        dna_ids = ",".join(str(i) for i in range(1, combined.top.n_chains))

        sources = []
        for pdb_id, (traj, ares, c) in loaded.items():
            _, keep_a = _match(_seq(keep_ref), keep_ref, _seq(ares), ares)
            idx = {r.index for r, _ in keep_a}
            sl = _clean_serials(traj.atom_slice(
                [a.index for a in traj.top.atoms if a.residue.index in idx]))
            sources.append((f"apo_{pdb_id}", sl, c))
        if ens is not None:
            _, keep_e = _match(_seq(keep_ref), keep_ref, _seq(eres), eres)
            idx = {r.index for r, _ in keep_e}
            sl = _clean_serials(ens.atom_slice(
                [a.index for a in ens.top.atoms if a.residue.index in idx]))
            sources.append(("bioemu", sl, None))

        for tag, sl, c in sources:
            t_pdb = out / "refs" / f"{pilot}_{tag}.pdb"
            t_xtc = out / "refs" / f"{pilot}_{tag}.xtc"
            sl[0].save_pdb(str(t_pdb))
            sl.save_xtc(str(t_xtc))
            docked = out / "docked" / f"{pilot}_{tag}"
            rc = run([args.python, STAGE2, "--ref", ref_pdb, "--traj", t_xtc,
                      "--top", t_pdb, "--pdb-id", pilot, "--out-dir", docked,
                      "--protein-chain", 0, "--dna-chains", dna_ids,
                      "--mismatch-action", "trim", "--allow-multimer"],
                     out / "logs" / f"{pilot}_{tag}_stage2.log")
            n = len(list(docked.glob(f"{pilot}_state_*.pdb"))) if docked.is_dir() else 0
            if rc != 0 or not n:
                print(f"  {tag}: SKIP dock (rc={rc}) -- logs/{pilot}_{tag}_stage2.log",
                      file=sys.stderr)
                continue
            csv_path = out / f"{pilot}_{tag}_fnat.csv"
            rc = run([args.python, SCORE, "--ref", ref_pdb, "--dir", docked,
                      "--pdb-id", pilot, "--out", csv_path,
                      "--protein-chain", 0, "--dna-chains", dna_ids],
                     out / "logs" / f"{pilot}_{tag}_fnat.log")
            if rc != 0 or not csv_path.is_file():
                print(f"  {tag}: SKIP score -- logs/{pilot}_{tag}_fnat.log", file=sys.stderr)
                continue
            try:
                f = _fnat_values(csv_path)
            except BrokenReference as exc:
                print(f"  BROKEN REFERENCE: {exc}", file=sys.stderr)
                continue
            if not f.size:
                continue
            rows.append({"pilot": pilot, "source": tag,
                         "method": (c or {}).get("method", "BioEmu"),
                         "n": len(f), "fnat_median": f"{np.median(f):.4f}",
                         "fnat_max": f"{f.max():.4f}",
                         "pass_rate": f"{(f >= args.floor).mean():.4f}"})
            print(f"  {tag:14s} n={len(f):3d}  fnat median {np.median(f):.3f}  "
                  f"pass {(f >= args.floor).mean():.0%}")

    if not rows:
        sys.exit("nothing scored")
    p = out / "fnat_pilots_summary.csv"
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {p}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
