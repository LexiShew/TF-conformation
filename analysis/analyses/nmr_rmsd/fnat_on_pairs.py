#!/usr/bin/env python
"""
Run the pipeline's Stage-2 dock + fnat gate on the NMR pairs.

The holo member of every pair is a protein-DNA complex, so each pair carries its
own docking target. That lets the pipeline's own interface-fidelity machinery be
pointed at structures where the apo/holo answer is known independently -- which
is the bridge between this analysis and the rest of the project:

    does where a BioEmu ensemble lands on the apo->holo axis (%across, from
    plot_pair_pca.py) predict how many of its frames survive the fnat gate?

If it does, fnat pass rate is a pipeline-side proxy for %across on the TFs that
have no NMR pair. If it does not, the gate is filtering on something unrelated
to how close the ensemble got to the bound state.

Per pair:
  1. write a single-model reference PDB from the holo NMR entry -- model 1,
     the matched protein chain, plus every DNA chain (protein first, so the
     reference chainids are 0 = protein, 1.. = DNA, as Stage 2 expects)

--source picks what gets docked, and the three together are the point:

  bioemu       the BioEmu ensemble (the thing being tested)
  holo_states  the DEPOSITED holo states, DNA stripped, put through the SAME
               dock. This is the reachable ceiling: a structure that IS the
               bound state, paying only the cost of rigid-body docking. Scoring
               deposited states WITH their own DNA is near-circular -- fnat's
               native contact set is defined from model 1 of that same bundle --
               so it measures deposition self-consistency, not what a docked
               structure can reach.
  apo_states   the DEPOSITED apo states docked onto the holo DNA. What a genuine
               experimental free state scores. This is the comparison BioEmu
               should actually be held to, since it is meant to be an apo
               ensemble: if real apo structures score no better, the gate is
               rejecting free-state conformations as such, not BioEmu's errors.

All three are trimmed to ONE residue set shared by the reference, the BioEmu
ensemble and the apo bundle, so the native contact set is identical across
sources and the numbers are comparable.
  2. stage2_redock.py: interface-aligned Kabsch of each BioEmu frame onto that
     reference's DNA, exactly as the pipeline docks onto crystal DNA
  3. fnat_gate/score_stage3.py: fnat per docked state vs the reference

Outputs (in --output-dir, default ./fnat):
  refs/<pair>_ref.pdb          the single-model docking target
  docked/<pair>/               per-state docked PDBs
  <pair>_fnat.csv              per-state fnat from the gate
  fnat_summary.csv             per pair: median fnat, pass rate at --floor

Usage:
    python fnat_on_pairs.py                     # every ready pair
    python fnat_on_pairs.py --pair-ids nhp6a
"""
import argparse
import csv
import difflib
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from pair_core import (CIF_SEARCH_DIRS, NMR_DIR, ensemble_paths, load_pairs,   # noqa: E402
                       residue_sequence)
from compute_nmr_state_rmsd import (load_reference, reference_cif,             # noqa: E402
                                    select_reference_chain)

_TFCONF = _HERE.parents[2]
STAGE2 = _TFCONF / "stage2_redock" / "stage2_redock.py"
SCORE = _TFCONF / "fnat_gate" / "score_stage3.py"
DNA_RES = {"DA", "DC", "DG", "DT", "DI", "DU", "A", "C", "G", "T", "U"}


def _clean_serials(traj):
    """mdtraj's cif reader can leave atom.serial as a str, which its own PDB
    writer then tries to use in an integer format ("not all arguments converted
    during string formatting"). Clearing it makes the writer fall back to the
    atom index."""
    for a in traj.top.atoms:
        a.serial = None
    return traj


def _match(seq_a, res_a, seq_b, res_b):
    """Residues of a and b that align, as two parallel lists."""
    ka, kb = [], []
    for i, j, n in difflib.SequenceMatcher(a=seq_a, b=seq_b, autojunk=False).get_matching_blocks():
        for k in range(n):
            ka.append(res_a[i + k])
            kb.append(res_b[j + k])
    return ka, kb


def _seq(res):
    return "".join(aa for _, aa in res)


def build_reference(pair, out_pdb, verbose=True):
    """Reference PDB + the residue set shared by reference, BioEmu and apo bundle."""
    import mdtraj as md
    pdb_id, chain = pair["holo_pdb"], pair["holo_chain"]
    top_path, xtc_path = ensemble_paths(pdb_id, chain)
    ens = md.load(str(xtc_path), top=str(top_path))
    ens_res = residue_sequence(ens.top)

    cif = reference_cif(pdb_id, CIF_SEARCH_DIRS, NMR_DIR / "cif", fetch=True)
    ref_all = load_reference(cif, verbose=False)
    ref = ref_all[0]
    prot_res, ref_chain, identity = select_reference_chain(ref.top, chain, _seq(ens_res))

    # apo bundle, so its residues constrain the shared set too
    apo_top, apo_xtc = ensemble_paths(pair["apo_pdb"], pair["apo_chain"])
    apo_ens = md.load(str(apo_xtc), top=str(apo_top))
    apo_cif = reference_cif(pair["apo_pdb"], CIF_SEARCH_DIRS, NMR_DIR / "cif", fetch=True)
    apo_all = load_reference(apo_cif, verbose=False)
    apo_res, _, _ = select_reference_chain(apo_all.top, pair["apo_chain"],
                                           _seq(residue_sequence(apo_ens.top)))

    # Stage 2 hard-exits on a Ca-count mismatch and checks before any trimming,
    # so every source must carry exactly the reference's residues. Intersect
    # reference / BioEmu / apo bundle once and trim all of them to it.
    r1, e1 = _match(_seq(prot_res), prot_res, _seq(ens_res), ens_res)
    r2, a1 = _match(_seq(r1), r1, _seq(apo_res), apo_res)
    keep_ref = r2
    ens_map = {id(x): y for x, y in zip(r1, e1)}
    keep_ens = [ens_map[id(x)] for x in keep_ref]
    keep_apo = a1
    if len(keep_ref) < 4:
        raise ValueError(f"{pdb_id}: only {len(keep_ref)} residues shared across sources")

    keep = {r.index for r, _ in keep_ref}
    prot_idx = [a.index for a in ref.top.atoms if a.residue.index in keep]
    dna_idx = [a.index for a in ref.top.atoms if a.residue.name.strip().upper() in DNA_RES]
    if not dna_idx:
        raise ValueError(f"{pdb_id}: no DNA found in the reference")

    combined = _clean_serials(ref.atom_slice(prot_idx).stack(ref.atom_slice(dna_idx)))
    # mdtraj can emit the same chain_id twice; biopython then MERGES those chains
    # and the gate is handed fewer than it was told ("chain index 2 out of range").
    for n, ch in enumerate(combined.top.chains):
        ch.chain_id = chr(ord("A") + n)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    combined.save_pdb(str(out_pdb))
    if verbose:
        print(f"  ref: chain {ref_chain} (id {identity:.2f}), {len(keep_ref)} shared residues, "
              f"{combined.top.n_chains - 1} DNA chain(s)")
    return dict(dna_ids=list(range(1, combined.top.n_chains)),
                ens=ens, keep_ens=keep_ens,
                ref_all=ref_all, keep_ref=keep_ref,
                apo_all=apo_all, keep_apo=keep_apo, out_dir=out_pdb.parent)


def source_ensemble(pair, info, source):
    """Write (topology, xtc) for the requested source, trimmed to the shared residues."""
    if source == "bioemu":
        traj, res = info["ens"], info["keep_ens"]
    elif source == "holo_states":
        traj, res = info["ref_all"], info["keep_ref"]
    elif source == "apo_states":
        traj, res = info["apo_all"], info["keep_apo"]
    else:
        raise ValueError(source)
    idx = {r.index for r, _ in res} if isinstance(res[0], tuple) else {r.index for r in res}
    sl = _clean_serials(traj.atom_slice([a.index for a in traj.top.atoms
                                        if a.residue.index in idx]))
    top = info["out_dir"] / f"{pair['pair_id']}_{source}.pdb"
    xtc = info["out_dir"] / f"{pair['pair_id']}_{source}.xtc"
    sl[0].save_pdb(str(top))
    sl.save_xtc(str(xtc))
    return top, xtc, sl.n_frames


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


def run(cmd, log):
    with open(log, "w") as fh:
        return subprocess.run([str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT).returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default=_HERE / "pairs.csv")
    ap.add_argument("--pair-ids", nargs="*", default=None)
    ap.add_argument("--output-dir", default=_HERE / "fnat")
    ap.add_argument("--floor", type=float, default=0.5,
                    help="fnat pass floor, matching the pipeline default (FNAT_FLOOR)")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--source", choices=("bioemu", "holo_states", "apo_states"),
                    default="bioemu", help="what to dock and score (see module docstring)")
    args = ap.parse_args()

    out = Path(args.output_dir)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    rows = []
    for pair in load_pairs(args.pairs, args.pair_ids):
        pid = pair["pair_id"]
        print(f"{pid}  ({pair['label']}) holo {pair['holo_pdb']}/{pair['holo_chain']}")
        try:
            info = build_reference(pair, out / "refs" / f"{pid}_ref.pdb")
        except Exception as exc:                        # noqa: BLE001
            print(f"  SKIP ref: {exc}", file=sys.stderr)
            continue

        docked = out / args.source / pid
        try:
            top_path, xtc_path, n_src = source_ensemble(pair, info, args.source)
        except Exception as exc:                        # noqa: BLE001
            print(f"  SKIP source: {exc}", file=sys.stderr)
            continue
        rc = run([args.python, STAGE2, "--ref", out / "refs" / f"{pid}_ref.pdb",
                  "--traj", xtc_path, "--top", top_path, "--pdb-id", pid,
                  "--out-dir", docked, "--protein-chain", 0,
                  "--dna-chains", ",".join(map(str, info["dna_ids"])),
                  "--mismatch-action", "trim", "--allow-multimer"],
                 out / "logs" / f"{pid}_{args.source}_stage2.log")
        n = len(list(docked.glob(f"{pid}_state_*.pdb"))) if docked.is_dir() else 0
        if rc != 0 or not n:
            print(f"  SKIP dock (rc={rc}, {n}/{n_src} states) -- "
                  f"logs/{pid}_{args.source}_stage2.log", file=sys.stderr)
            continue
        print(f"  docked {n} {args.source} states")

        csv_path = out / f"{pid}_fnat_{args.source}.csv"
        rc = run([args.python, SCORE, "--ref", out / "refs" / f"{pid}_ref.pdb",
                  "--dir", docked, "--pdb-id", pid, "--out", csv_path,
                  "--protein-chain", 0,
                  "--dna-chains", ",".join(map(str, info["dna_ids"]))],
                 out / "logs" / f"{pid}_{args.source}_fnat.log")
        if rc != 0 or not csv_path.is_file():
            print(f"  SKIP score (rc={rc}) -- see logs/{pid}_fnat.log", file=sys.stderr)
            continue

        try:
            f = _fnat_values(csv_path)
        except BrokenReference as exc:
            print(f"  BROKEN REFERENCE: {exc}", file=sys.stderr)
            continue
        if not f.size:
            print("  SKIP: no fnat values", file=sys.stderr)
            continue
        rows.append({"pair_id": pid, "label": pair["label"], "family": pair["family"],
                     "holo_pdb": pair["holo_pdb"], "n_states": len(f),
                     "fnat_median": f"{np.median(f):.4f}",
                     "fnat_q1": f"{np.percentile(f, 25):.4f}",
                     "fnat_q3": f"{np.percentile(f, 75):.4f}",
                     "fnat_max": f"{f.max():.4f}",
                     "pass_floor": args.floor,
                     "pass_rate": f"{(f >= args.floor).mean():.4f}",
                     "source": args.source})
        print(f"  fnat median {np.median(f):.3f}, pass rate at {args.floor:g}: "
              f"{(f >= args.floor).mean():.0%}")

    if not rows:
        sys.exit("nothing scored")
    p = out / f"fnat_summary_{args.source}.csv"
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {p}  ({len(rows)} pairs)")


if __name__ == "__main__":
    main()
