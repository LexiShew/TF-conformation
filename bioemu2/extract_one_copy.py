#!/usr/bin/env python
"""
extract_one_copy.py — reduce each BioEmu2 frame to ONE protein chain + its cognate duplex.

Why this exists
---------------
`System.from_cif()` imports the whole deposited assembly. For 1K79 that is TWO
copies of the complex, so every frame BioEmu2 writes holds 4 DNA chains
(4 x 15 nt) and 2 protein chains (2 x 110 aa). The repo's fnat gate assumes one
protein chain and one duplex per frame: `interface_rmsd.ordered()` takes no chain
selection, pools all 60 DNA residues, and `score()` then fails its
`len(md) == len(R['rd'])` test and silently falls back to the REFERENCE DNA
(finding F10 of docs/codebase_review_2026-10-15.md). Because BioEmu2 frames are
never superposed onto the crystal, that fallback measures protein-DNA contacts
ACROSS two coordinate frames, which is why the first ets1 run reported fnat
median 0.000 (87/100 frames exactly zero). That number is arithmetic, not a
measurement of the model.

This script writes a parallel `states_1copy/` tree holding one protein chain and
its two cognate strands per frame. Contacts are internal distances, so a
one-copy frame is scored in the model's own coordinate frame with
`use_model_dna=True` and no superposition is needed — the fnat it yields IS a
measurement of the sampled interface.

What it reuses
--------------
  * chain classification, heavy-atom filter   fnat_gate/interface_rmsd.py
  * reference interface + strand order        interface_rmsd.ref_side (same
                                              config chain selection the gate
                                              itself uses)
Only the copy selection is new, because nothing in the repo does it: the
pipeline's frames never had a second copy to choose between.

Strand ORDER is load-bearing, not cosmetic
------------------------------------------
`score()` indexes native contacts as (i, j) with j into `R['rd']`, the reference
DNA residue list. The model's DNA list is built the same way — chain iteration
order, residues ascending within a chain — so if the two strands are emitted in
the opposite order to the reference every j is off by one strand length and fnat
is garbage again, silently. The two strands are therefore matched to the
reference strands by nucleotide sequence before writing, and the per-frame
report records the match so the choice is auditable rather than assumed.

Usage
-----
    python bioemu2/extract_one_copy.py --tf-name ets1
    python bioemu2/compare_to_pipeline.py --tf-name ets1 --arms frozen dnarelax \
        --bioemu2-dir output/bioemu2/1k79_complex/states_1copy

Outputs
-------
    <states-dir>/../states_1copy/<pdb>_state_NNN.pdb   one protein + one duplex
    <states-dir>/../extract_one_copy_report.csv        per-frame selection audit
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_TFCONF = _HERE.parent

sys.path.insert(0, str(_TFCONF / "fnat_gate"))
import interface_rmsd as ir                                       # noqa: E402

from Bio.PDB import PDBIO                                         # noqa: E402
from Bio.PDB.Model import Model as BioModel                       # noqa: E402
from Bio.PDB.Structure import Structure as BioStructure           # noqa: E402

STATE_RE = re.compile(r"_state_(\d+)\.pdb$")


def pilot_env(tf_name):
    cfg = _TFCONF / "config" / "pilots" / f"{tf_name}.sh"
    if not cfg.exists():
        sys.exit(f"ERROR: no pilot config at {cfg}")
    r = subprocess.run(["bash", "-c", f'set -a; source "{cfg}" >/dev/null 2>&1; env'],
                       capture_output=True, text=True, cwd=_TFCONF)
    return dict(l.split("=", 1) for l in r.stdout.splitlines() if "=" in l)


def classify_chains(model):
    """Split a frame's chains into protein and DNA, by residue content.

    Same residue vocabulary interface_rmsd.ordered() uses, so a chain this calls
    DNA is a chain the gate will also call DNA.
    """
    prot, dna = [], []
    for ch in model:
        names = [r.get_resname().strip() for r in ch]
        n_dna = sum(1 for n in names if n in ir.DNA_RES)
        n_aa = sum(1 for n in names if n in ir.AA3TO1)
        if n_dna and n_dna >= n_aa:
            dna.append(ch)
        elif n_aa:
            prot.append(ch)
    return prot, dna


def heavy_coords(chain):
    return np.array([a.coord for r in chain for a in ir.heavy(r)], dtype=float)


def contact_count(a_xyz, b_xyz, cutoff):
    """Number of heavy-atom pairs within cutoff. Chunked to bound memory."""
    if not len(a_xyz) or not len(b_xyz):
        return 0
    n = 0
    for i in range(0, len(a_xyz), 2000):
        blk = a_xyz[i:i + 2000]
        d2 = ((blk[:, None, :] - b_xyz[None, :, :]) ** 2).sum(-1)
        n += int((d2 <= cutoff * cutoff).sum())
    return n


def strand_seq(chain):
    return [r.get_resname().strip()
            for r in sorted(chain, key=lambda r: r.id[1])
            if r.get_resname().strip() in ir.DNA_RES]


def seq_match(a, b):
    """Fraction of identical positions over the overlap, 0 if either is empty."""
    m = min(len(a), len(b))
    if not m:
        return 0.0
    return sum(1 for i in range(m) if a[i] == b[i]) / m


def reference_strands(R):
    """Reference DNA strands as resname lists, in R['rd'] order.

    R['rd'] is already the config-selected, chain-ordered reference DNA, so the
    strand boundaries are just the points where the parent chain changes.
    """
    out = []
    last = object()
    for r in R["rd"]:
        cid = r.get_parent().id
        if cid != last:
            out.append((cid, []))
            last = cid
        out[-1][1].append(r.get_resname().strip())
    return out


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tf-name", required=True)
    p.add_argument("--states-dir", default=None,
                   help="default: output/bioemu2/<pdb>_complex/states")
    p.add_argument("--out-dir", default=None,
                   help="default: <states-dir>/../states_1copy")
    p.add_argument("--report", default=None,
                   help="default: <states-dir>/../extract_one_copy_report.csv")
    p.add_argument("--iface-cutoff", type=float, default=5.0)
    p.add_argument("--contact-cutoff", type=float, default=4.5)
    p.add_argument("--gap", type=int, default=4)
    a = p.parse_args()

    env = pilot_env(a.tf_name)
    pdb_id = (env.get("PDB_ID") or "").lower()
    ref_cif = _TFCONF / "structures" / "source_chains" / f"{pdb_id}_chains" / f"{pdb_id}.cif"
    if not ref_cif.exists():
        sys.exit(f"ERROR: reference cif not found: {ref_cif}")

    prot_cfg = env.get("PROTEIN_CHAIN")
    dna_cfg = env.get("DNA_CHAINS")
    protein_chain = int(str(prot_cfg).split("#")[0].strip()) if prot_cfg else None
    dna_chains = [int(x) for x in dna_cfg.split(",") if x.strip()] if dna_cfg else None

    states = Path(a.states_dir) if a.states_dir else \
        _TFCONF / "output" / "bioemu2" / f"{pdb_id}_complex" / "states"
    out = Path(a.out_dir) if a.out_dir else states.parent / "states_1copy"
    report = Path(a.report) if a.report else states.parent / "extract_one_copy_report.csv"
    out.mkdir(parents=True, exist_ok=True)

    # The reference side, defined by the gate's own logic and the pilot config —
    # the same R the comparison will score against.
    R = ir.ref_side(ir.load_models(str(ref_cif))[0], a.iface_cutoff,
                    a.contact_cutoff, a.gap,
                    protein_chain=protein_chain, dna_chains=dna_chains)
    ref_str = reference_strands(R)
    n_ref_dna = len(R["rd"])
    print(f"=== {a.tf_name} ({pdb_id}) ===")
    print(f"reference  : {ref_cif}")
    print(f"ref strands: " + ", ".join(f"{cid}={''.join(s[-1] for s in seq)}"
                                       for cid, seq in ref_str))
    print(f"ref DNA    : {n_ref_dna} residues in {len(ref_str)} strands")
    if len(ref_str) != 2:
        sys.exit(f"ERROR: the reference interface has {len(ref_str)} DNA strands, not 2. "
                 f"Strand matching below assumes a single duplex; a tandem or "
                 f"multi-duplex pilot needs the selection generalized first.")
    print(f"states in  : {states}")
    print(f"states out : {out}")

    files = sorted(glob.glob(os.path.join(str(states), "*_state_*.pdb")))
    if not files:
        sys.exit(f"ERROR: no *_state_*.pdb in {states}")

    io = PDBIO()
    rows, n_written, n_flipped, n_nonduplex, n_badcount = [], 0, 0, 0, 0
    for f in files:
        m = STATE_RE.search(os.path.basename(f))
        state = int(m.group(1)) if m else -1
        model = ir.load_models(f)[0]
        prots, dnas = classify_chains(model)
        if not prots or len(dnas) < 2:
            print(f"  state {state}: SKIP ({len(prots)} protein, {len(dnas)} DNA chains)",
                  file=sys.stderr)
            continue

        pxyz = {ch.id: heavy_coords(ch) for ch in prots}
        dxyz = {ch.id: heavy_coords(ch) for ch in dnas}

        # 1. the protein copy that is actually bound: most protein-DNA contacts.
        pc = {ch.id: sum(contact_count(pxyz[ch.id], dxyz[d.id], a.contact_cutoff)
                         for d in dnas) for ch in prots}
        prot_id = max(pc, key=pc.get)

        # 2. its two best-contacting strands.
        dc = {d.id: contact_count(pxyz[prot_id], dxyz[d.id], a.contact_cutoff)
              for d in dnas}
        ranked = sorted(dc, key=dc.get, reverse=True)
        s1, s2 = ranked[0], ranked[1]

        # 3. are those two strands actually base-paired to each other? Reported,
        #    not enforced: a melted duplex is a real BioEmu2 outcome and dropping
        #    it would bias the ensemble.
        inter = contact_count(dxyz[s1], dxyz[s2], a.contact_cutoff)
        if inter == 0:
            n_nonduplex += 1

        # 4. strand order must match the reference, or every native-contact index
        #    is off by one strand length (see module docstring).
        chains = {ch.id: ch for ch in dnas}
        q1, q2 = strand_seq(chains[s1]), strand_seq(chains[s2])
        keep = (seq_match(ref_str[0][1], q1) + seq_match(ref_str[1][1], q2)) / 2
        swap = (seq_match(ref_str[0][1], q2) + seq_match(ref_str[1][1], q1)) / 2
        flipped = swap > keep
        if flipped:
            s1, s2 = s2, s1
            q1, q2 = q2, q1
            n_flipped += 1
        order_match = max(keep, swap)

        n_dna_res = len(q1) + len(q2)
        if n_dna_res != n_ref_dna:
            n_badcount += 1
            print(f"  state {state}: WARNING {n_dna_res} DNA residues != reference "
                  f"{n_ref_dna}; the gate will still fall back to reference DNA",
                  file=sys.stderr)

        # 5. emit protein, strand 1, strand 2 — in that order.
        new_s = BioStructure("one_copy")
        new_m = BioModel(0)
        new_s.add(new_m)
        for cid in (prot_id, s1, s2):
            src = chains[cid] if cid in chains else next(c for c in prots if c.id == cid)
            cp = src.copy()
            cp.detach_parent()
            new_m.add(cp)
        io.set_structure(new_s)
        io.save(str(out / os.path.basename(f)))
        n_written += 1

        rows.append(dict(state=state, protein_chain=prot_id,
                         dna_chain_1=s1, dna_chain_2=s2,
                         prot_dna_contacts=pc[prot_id],
                         contacts_strand_1=dc[s1], contacts_strand_2=dc[s2],
                         interstrand_contacts=inter, n_dna_res=n_dna_res,
                         strand_order_flipped=int(flipped),
                         strand_seq_match=round(order_match, 3)))

    if not rows:
        sys.exit(f"ERROR: no frame yielded a protein chain plus two DNA chains "
                 f"({len(files)} files read from {states}).")

    with open(report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- summary ---------------------------------------------------------
    picks = {}
    for r in rows:
        key = (r["protein_chain"], r["dna_chain_1"], r["dna_chain_2"])
        picks[key] = picks.get(key, 0) + 1
    print(f"\nwrote {n_written}/{len(files)} one-copy frames -> {out}")
    print(f"report -> {report}")
    print(f"copy selection (protein, strand1, strand2) x frames:")
    for key, n in sorted(picks.items(), key=lambda kv: -kv[1]):
        print(f"    {key}  x{n}")
    print(f"strand order flipped to match reference : {n_flipped}/{n_written}")
    print(f"chosen strands NOT in contact (melted)  : {n_nonduplex}/{n_written}")
    print(f"DNA residue count != reference          : {n_badcount}/{n_written}")
    if n_badcount:
        print("  Those frames will STILL trip the F10 fallback and their fnat is not a\n"
              "  model-frame measurement. Exclude them before reading fnat.",
              file=sys.stderr)
    inter_vals = np.array([r["interstrand_contacts"] for r in rows], dtype=float)
    pd_vals = np.array([r["prot_dna_contacts"] for r in rows], dtype=float)
    print(f"interstrand contacts  median {np.median(inter_vals):.0f} "
          f"[{np.percentile(inter_vals, 25):.0f}, {np.percentile(inter_vals, 75):.0f}]")
    print(f"protein-DNA contacts  median {np.median(pd_vals):.0f} "
          f"[{np.percentile(pd_vals, 25):.0f}, {np.percentile(pd_vals, 75):.0f}]")
    print("\nNow score them with the existing comparison, which will take the\n"
          "use_model_dna path because the DNA residue count now matches:\n"
          f"    python bioemu2/compare_to_pipeline.py --tf-name {a.tf_name} "
          f"--arms frozen dnarelax \\\n        --bioemu2-dir {out}")


if __name__ == "__main__":
    main()
