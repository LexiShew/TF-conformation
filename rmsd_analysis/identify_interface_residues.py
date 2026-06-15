#!/usr/bin/env python
"""
Identify protein–DNA interface residues from a reference structure and write
the result as JSON keyed by 0-based residue ordinal among protein residues
(the same numbering convention compute_rmsds.py uses internally).

Three methods:

  --method distance   : a protein residue is interface iff any of its heavy
                        atoms is within --cutoff Å of any DNA heavy atom
                        (default 4.5 Å). Fast, deterministic.
  --method sasa       : a protein residue is interface iff its
                        ΔSASA = SASA(apo) − SASA(complex) ≥ --cutoff Å²
                        (default 5 Å²). Catches 2nd-shell residues that
                        get buried without making direct vdW contact.
  --method dssr       : stub — invokes x3dna-dssr and parses its
                        nt_aa_interactions block to get base/backbone
                        split. Implemented if dssr is on PATH; falls back
                        to distance with a warning otherwise.

Usage (pilot-aware, batch — matches compute_rmsds.py UX):
    python identify_interface_residues.py \\
        --pilots egr1 dux4 tbp \\
        --output-dir ./interface_residues \\
        --method distance --cutoff 4.5

Single-PDB mode (no pilot config needed):
    python identify_interface_residues.py \\
        --reference-cif path/to/1aay.cif --protein-chain 2 \\
        --output ./interface_residues/egr1.json --method distance
"""
import argparse
import json
import os
import sys
from pathlib import Path
import numpy as np


def import_mdtraj():
    global md
    import mdtraj as md_mod
    md = md_mod


PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
}
DNA_RESNAMES = {"DA", "DC", "DG", "DT", "DU", "A", "C", "G", "T", "U",
                "DA3", "DA5", "DC3", "DC5", "DG3", "DG5", "DT3", "DT5"}

# Defaults match compute_rmsds.py (env-overridable).
PROJECT_ROOT = os.environ.get("DEEPPBS_PROJECT_ROOT",
                              "/project2/rohs_102/shewchuk")
DEEPPBS_DIR  = os.environ.get("DEEPPBS_DIR", f"{PROJECT_ROOT}/DeepPBS")
BIOEMU_ROOT  = os.environ.get(
    "DEEPPBS_BIOEMU_ROOT",
    f"{PROJECT_ROOT}/TF-conformation/deeppbs_pdbs/monomer_chains")
PILOTS_DIR   = os.environ.get("DEEPPBS_PILOTS_DIR",
                              f"{DEEPPBS_DIR}/run/jobs/config/pilots")


def load_pilot_config(tf_name):
    path = os.path.join(PILOTS_DIR, f"{tf_name}.sh")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Pilot config not found: {path}")
    cfg = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("export ") and "=" in line:
                key_val = line[len("export "):].split("=", 1)
                key = key_val[0].strip()
                val = key_val[1].strip().strip('"').split('#')[0].strip().strip('"')
                cfg[key] = val
    return cfg


def reference_cif_path(pdb_id):
    return f"{BIOEMU_ROOT}/{pdb_id}_chains/{pdb_id}.cif"


# ---------------------------------------------------------------------------
# Method: distance ----------------------------------------------------------
def find_interface_distance(ref, protein_chain, cutoff_A):
    """Per protein residue, min distance to any DNA heavy atom (Å)."""
    cutoff_nm = cutoff_A / 10.0

    # DNA heavy atoms (any chain)
    dna_idx = [a.index for a in ref.topology.atoms
               if a.residue.name in DNA_RESNAMES
               and a.element is not None and a.element.symbol != "H"]
    if not dna_idx:
        raise ValueError("No DNA heavy atoms found in reference structure.")
    dna_coords = ref.xyz[0, np.array(dna_idx, dtype=int)]   # (n_dna, 3) nm

    interface = []
    res_pos = 0
    for chain in ref.topology.chains:
        if chain.index != protein_chain:
            continue
        for res in chain.residues:
            if res.name not in PROTEIN_RESNAMES:
                continue
            heavy = [a.index for a in res.atoms
                     if a.element is not None and a.element.symbol != "H"]
            if heavy:
                prot_coords = ref.xyz[0, np.array(heavy, dtype=int)]
                diffs = prot_coords[:, None, :] - dna_coords[None, :, :]
                dists = np.sqrt((diffs ** 2).sum(axis=-1))
                min_d_nm = float(dists.min())
                if min_d_nm < cutoff_nm:
                    interface.append({
                        "position":       res_pos,
                        "resname":        res.name,
                        "resseq":         int(res.resSeq),
                        "min_distance_A": round(min_d_nm * 10.0, 3),
                    })
            res_pos += 1
    return interface


# ---------------------------------------------------------------------------
# Method: ΔSASA -------------------------------------------------------------
def find_interface_sasa(ref, protein_chain, threshold_A2):
    """ΔSASA = SASA(apo) − SASA(complex), per protein residue, in Å²."""
    threshold_nm2 = threshold_A2 / 100.0   # 1 nm² = 100 Å²

    sasa_complex = md.shrake_rupley(ref, mode="residue")[0]  # nm²

    # Build apo (drop DNA atoms)
    keep = [a.index for a in ref.topology.atoms
            if a.residue.name not in DNA_RESNAMES]
    apo = ref.atom_slice(np.array(keep, dtype=int))
    sasa_apo = md.shrake_rupley(apo, mode="residue")[0]

    # Match per-residue between complex and apo by (chain_index, resSeq, resname).
    apo_by_key = {}
    for i, res in enumerate(apo.topology.residues):
        apo_by_key[(res.chain.index, int(res.resSeq), res.name)] = float(sasa_apo[i])

    interface = []
    res_pos = 0
    for chain in ref.topology.chains:
        if chain.index != protein_chain:
            continue
        for res in chain.residues:
            if res.name not in PROTEIN_RESNAMES:
                continue
            key = (chain.index, int(res.resSeq), res.name)
            s_complex = float(sasa_complex[res.index])
            s_apo = apo_by_key.get(key)
            if s_apo is None:
                res_pos += 1
                continue
            dsasa_nm2 = s_apo - s_complex
            if dsasa_nm2 >= threshold_nm2:
                interface.append({
                    "position": res_pos,
                    "resname":  res.name,
                    "resseq":   int(res.resSeq),
                    "dSASA_A2": round(dsasa_nm2 * 100.0, 3),
                })
            res_pos += 1
    return interface


# ---------------------------------------------------------------------------
# Method: DSSR (optional) ---------------------------------------------------
def find_interface_dssr(ref_path, protein_chain, cutoff_A):
    """Use x3dna-dssr if available to get base/backbone split. Falls back
    to the distance method if dssr isn't on PATH."""
    import shutil, subprocess, tempfile, json as _json
    if shutil.which("x3dna-dssr") is None:
        print("  WARN: x3dna-dssr not found on PATH; falling back to "
              "--method distance", file=sys.stderr)
        return None  # caller will retry with distance

    with tempfile.TemporaryDirectory() as tmp:
        json_path = os.path.join(tmp, "out.json")
        try:
            subprocess.run(
                ["x3dna-dssr", f"-i={ref_path}", f"-o={json_path}", "--json"],
                check=True, capture_output=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"  WARN: x3dna-dssr failed ({e}); falling back to "
                  f"distance method", file=sys.stderr)
            return None
        with open(json_path) as f:
            data = _json.load(f)

    # nt_aa_interactions: list of {nt: ..., aa: ..., type: ..., level: ...}
    # type in {"base", "backbone", "sugar"}; level in {"H-bond", "vdW", ...}
    interactions = data.get("nt_aa_interactions") or data.get("ntaa_interactions") or []
    # Reduce to protein chain residue resSeq -> (set of contact_types)
    by_resseq = {}
    for it in interactions:
        aa = it.get("aa") or it.get("aa_id") or ""
        # Format examples: "C.ARG18" (chain.resname+resseq) — parse defensively
        if "." not in aa:
            continue
        try:
            resseq = int("".join(c for c in aa.split(".")[1] if c.isdigit()))
        except ValueError:
            continue
        by_resseq.setdefault(resseq, set()).add(it.get("type", "unknown"))

    # Walk the reference topology to map resseq -> 0-based protein-residue pos.
    interface = []
    res_pos = 0
    for chain in ref_ref.topology.chains if False else []:
        # placeholder — actual walk happens in the caller after we re-load
        pass
    return by_resseq  # caller assembles


# ---------------------------------------------------------------------------
# Top-level per-pilot driver
# ---------------------------------------------------------------------------
def identify_for_pilot(pdb_id, protein_chain, method, cutoff,
                       reference_cif=None):
    cif = reference_cif or reference_cif_path(pdb_id)
    if not os.path.isfile(cif):
        raise FileNotFoundError(f"Reference not found: {cif}")
    ref = md.load(cif)

    n_protein_residues = sum(
        1 for chain in ref.topology.chains
        if chain.index == protein_chain
        for res in chain.residues
        if res.name in PROTEIN_RESNAMES
    )

    if method == "distance":
        interface = find_interface_distance(ref, protein_chain, cutoff)
        method_args = {"cutoff_A": cutoff}
    elif method == "sasa":
        interface = find_interface_sasa(ref, protein_chain, cutoff)
        method_args = {"threshold_A2": cutoff}
    elif method == "dssr":
        dssr_resseq = find_interface_dssr(cif, protein_chain, cutoff)
        if dssr_resseq is None:
            print("  Falling back to distance method.", file=sys.stderr)
            return identify_for_pilot(pdb_id, protein_chain, "distance",
                                       cutoff or 4.5, reference_cif)
        # Map resseq -> position
        interface = []
        res_pos = 0
        for chain in ref.topology.chains:
            if chain.index != protein_chain:
                continue
            for res in chain.residues:
                if res.name not in PROTEIN_RESNAMES:
                    continue
                types = dssr_resseq.get(int(res.resSeq))
                if types:
                    interface.append({
                        "position":      res_pos,
                        "resname":       res.name,
                        "resseq":        int(res.resSeq),
                        "contact_types": sorted(types),
                    })
                res_pos += 1
        method_args = {"cutoff_A": cutoff}
    else:
        raise ValueError(f"Unknown method: {method}")

    return {
        "pdb_id":                     pdb_id,
        "protein_chain":              protein_chain,
        "method":                     method,
        "method_args":                method_args,
        "n_total_protein_residues":   n_protein_residues,
        "n_interface_residues":       len(interface),
        "interface_residue_positions": [r["position"] for r in interface],
        "interface_residue_details":   interface,
    }


def main():
    global BIOEMU_ROOT, PILOTS_DIR
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--pilots", nargs="+",
                   help="Pilot names (egr1, dux4, tbp); each must have a "
                        "<tf>.sh in --pilots-dir.")
    g.add_argument("--reference-cif",
                   help="Path to a single reference structure (with DNA). "
                        "Requires --protein-chain and --output.")

    p.add_argument("--protein-chain", type=int,
                   help="0-based chain index of the protein chain "
                        "(only when using --reference-cif).")
    p.add_argument("--output",
                   help="Output JSON path (only when using --reference-cif).")
    p.add_argument("--output-dir", default="./interface_residues",
                   help="Where per-pilot JSONs go (default: "
                        "./interface_residues).")

    p.add_argument("--method", choices=("distance", "sasa", "dssr"),
                   default="distance",
                   help="Identification method (default: distance).")
    p.add_argument("--cutoff", type=float, default=None,
                   help="Distance: cutoff in Å (default 4.5). "
                        "SASA: ΔSASA threshold in Å² (default 5).")

    p.add_argument("--bioemu-root", default=BIOEMU_ROOT,
                   help=f"(default: {BIOEMU_ROOT}; env DEEPPBS_BIOEMU_ROOT)")
    p.add_argument("--pilots-dir", default=PILOTS_DIR,
                   help=f"(default: {PILOTS_DIR}; env DEEPPBS_PILOTS_DIR)")
    args = p.parse_args()

    BIOEMU_ROOT = args.bioemu_root
    PILOTS_DIR  = args.pilots_dir

    if args.cutoff is None:
        args.cutoff = 4.5 if args.method != "sasa" else 5.0

    import_mdtraj()

    if args.reference_cif:
        if args.protein_chain is None or args.output is None:
            sys.exit("--reference-cif requires both --protein-chain and --output")
        pdb_id = Path(args.reference_cif).stem.lower()
        out = identify_for_pilot(pdb_id, args.protein_chain, args.method,
                                  args.cutoff, reference_cif=args.reference_cif)
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        with open(args.output, "w") as f:
            json.dump(out, f, indent=2)
        print(f"  wrote {args.output}  ({out['n_interface_residues']}/"
              f"{out['n_total_protein_residues']} residues)")
        return

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for tf in args.pilots:
        try:
            cfg = load_pilot_config(tf)
            pdb_id = cfg["PDB_ID"]
            chain  = int(cfg["PROTEIN_CHAIN"])
            print(f"=== {tf}  (PDB={pdb_id}, chain={chain}) ===")
            out = identify_for_pilot(pdb_id, chain, args.method, args.cutoff)
            out["tf"] = tf
            json_path = output_dir / f"{tf}.json"
            with open(json_path, "w") as f:
                json.dump(out, f, indent=2)
            print(f"  wrote {json_path}  "
                  f"({out['n_interface_residues']}/"
                  f"{out['n_total_protein_residues']} residues)")
        except Exception as e:
            print(f"  ERROR for {tf}: {e}", file=sys.stderr)
            import traceback; traceback.print_exc()


if __name__ == "__main__":
    main()
