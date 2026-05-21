"""
Stage 2: Cα-align each Stage 1 frame onto a reference structure, carry the
reference DNA and any structural metal ions across, write per-state PDBs.

Carrying metals: BioEmu doesn't sample metals and HPACKER doesn't restore
them. To preserve metal coordination geometry in Stage 3, we extract metal
positions from the reference crystal and apply the same Kabsch transform
that we apply to the reference DNA. The metals are written as HETATM
records in each docked PDB so Stage 3 sees them.

This version also validates that BioEmu's protein and the reference crystal's
protein are the same residue sequence — not just the same count of Cα atoms.
A mere count-match misses cases where the crystal has a disordered gap that
BioEmu modeled across, producing positional misalignment (this caused the
DUX4 5z6z pilot to dock DNA against the wrong protein region).

Auto-detects protein vs DNA chains in the reference unless --protein-chain
and --dna-chains are explicitly provided. Use --inspect-only to see the
detected layout without writing any output.
"""
import argparse
import os
import sys
import numpy as np
import mdtraj as md


PROTEIN_RESNAMES = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
    "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","HID","HIE","HIP","HSD","HSE","HSP",
}
DNA_RESNAMES = {
    "DA","DG","DC","DT",
    "DA5","DG5","DC5","DT5","DA3","DG3","DC3","DT3",
}
HISTIDINE_VARIANTS = {"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"}
# Structural metal ions worth carrying across to docked frames so Stage 3
# can set up coordination-cage restraints.
STRUCTURAL_METAL_RESNAMES = {"ZN", "MG", "MN", "FE", "CA", "CO", "NI", "CU"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb-id", required=True, help="Used in output filenames, e.g. '1aay'")
    p.add_argument("--ref", required=True, help="Reference structure (cif or pdb)")
    p.add_argument("--traj", help="Trajectory of relaxed protein frames (xtc)")
    p.add_argument("--top", help="Topology for --traj (pdb)")
    p.add_argument("--out-dir", help="Where per-state PDBs are written")
    p.add_argument("--protein-chain", type=int, default=None,
                   help="0-based chainid of the protein in --ref")
    p.add_argument("--dna-chains", default=None,
                   help="Comma-separated 0-based chainids of DNA strands")
    p.add_argument("--inspect-only", action="store_true",
                   help="Print detected chain layout and exit")
    p.add_argument("--require-sequence-match", action="store_true", default=True,
                   help="Verify residue sequences match between traj and ref (default: True)")
    p.add_argument("--mismatch-action", choices=["fail", "warn", "trim"], default="fail",
                   help="What to do if residue sequences differ: "
                        "fail = exit (default), "
                        "warn = continue with misalignment risk, "
                        "trim = slice both to the longest common-sequence run")
    p.add_argument("--max-mismatches", type=int, default=0,
                   help="Maximum allowed sequence mismatches before action triggers. "
                        "Set higher to tolerate occasional discrepancies "
                        "(e.g. HIS protonation variant labels).")
    return p.parse_args()


def detect_chains(ref):
    classification = {"protein": [], "dna": [], "other": []}
    for chain in ref.topology.chains:
        residues = list(chain.residues)
        n_prot = sum(1 for r in residues if r.name in PROTEIN_RESNAMES)
        n_dna  = sum(1 for r in residues if r.name in DNA_RESNAMES)
        n_total = len(residues)
        if n_prot >= 10 and n_prot >= n_dna:
            classification["protein"].append((chain.index, n_total, f"protein (n_res={n_total}, n_prot={n_prot})"))
        elif n_dna >= 3 and n_dna >= n_prot:
            classification["dna"].append((chain.index, n_total, f"DNA (n_res={n_total}, n_dna={n_dna})"))
        else:
            sample = ",".join(r.name for r in residues[:3])
            classification["other"].append((chain.index, n_total, f"other (n_res={n_total}, sample={sample})"))
    return classification


def normalize_resname(name):
    """Map His-protonation variants to one canonical name."""
    if name in HISTIDINE_VARIANTS: return "HIS"
    return name


def check_sequence_match(traj, ref, protein_chain, max_mismatches):
    """Compare residue sequences between traj (all chains protein) and ref protein chain.

    Returns (n_mismatches, mismatch_positions, traj_seq, ref_seq).
    """
    traj_residues = [normalize_resname(r.name) for r in traj.topology.residues]
    ref_chain = list(ref.topology.chains)[protein_chain]
    ref_residues = [normalize_resname(r.name) for r in ref_chain.residues]

    if len(traj_residues) != len(ref_residues):
        return None, None, traj_residues, ref_residues  # length mismatch, handled separately

    mismatches = []
    for i, (t, r) in enumerate(zip(traj_residues, ref_residues)):
        if t != r:
            mismatches.append((i, t, r))
    return len(mismatches), mismatches, traj_residues, ref_residues


def find_reference_metals(ref):
    """Find structural metal ions in the reference structure.

    Returns:
        ref_metal_atom_idx (np.ndarray): atom indices of every metal atom
        metal_info (list): list of (resname, atom_index) tuples for logging
    """
    ref_metal_atom_idx = []
    metal_info = []
    for chain in ref.topology.chains:
        for res in chain.residues:
            if res.name.strip().upper() in STRUCTURAL_METAL_RESNAMES:
                for atom in res.atoms:
                    ref_metal_atom_idx.append(atom.index)
                    metal_info.append((res.name.strip().upper(), atom.index))
    return np.array(ref_metal_atom_idx, dtype=int), metal_info


def report_reference_metals(ref_metal_atom_idx, metal_info):
    print(f"Reference structural metal atoms: {len(ref_metal_atom_idx)}")
    for resname, idx in metal_info:
        print(f"  {resname} (atom index {idx})")


def kabsch(P, Q):
    p_centroid = P.mean(axis=0)
    q_centroid = Q.mean(axis=0)
    H = (P - p_centroid).T @ (Q - q_centroid)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = q_centroid - R @ p_centroid
    return R, t


def find_common_subsequence(seq_a, seq_b):
    """Return (start_a, start_b, length) of the longest contiguous match."""
    n, m = len(seq_a), len(seq_b)
    best = (0, 0, 0)
    for i in range(n):
        for j in range(m):
            k = 0
            while i + k < n and j + k < m and seq_a[i + k] == seq_b[j + k]:
                k += 1
            if k > best[2]:
                best = (i, j, k)
    return best


def main():
    args = parse_args()
    ref = md.load(args.ref)
    detection = detect_chains(ref)

    print("=== Chain detection ===")
    for kind in ("protein", "dna", "other"):
        for chainid, n, label in detection[kind]:
            print(f"  chainid {chainid}: {label}")

    if args.protein_chain is not None:
        protein_chain = args.protein_chain
        print(f"\nUsing user-specified protein chain: {protein_chain}")
    else:
        if not detection["protein"]:
            print("ERROR: No protein chain detected. Specify --protein-chain manually.", file=sys.stderr)
            sys.exit(1)
        protein_chain = max(detection["protein"], key=lambda x: x[1])[0]
        print(f"\nAuto-selected protein chain: {protein_chain} (longest)")

    if args.dna_chains is not None:
        dna_chain_ids = [int(c) for c in args.dna_chains.split(",")]
        print(f"Using user-specified DNA chains: {dna_chain_ids}")
    else:
        if not detection["dna"]:
            print("ERROR: No DNA chain detected. Specify --dna-chains manually.", file=sys.stderr)
            sys.exit(1)
        dna_chain_ids = [c[0] for c in detection["dna"]]
        print(f"Auto-selected DNA chains: {dna_chain_ids}")

    # Identify structural metals (carried across to docked frames during the
    # main per-frame loop). Reported early so --inspect-only users can verify.
    ref_metal_atom_idx, metal_info = find_reference_metals(ref)
    report_reference_metals(ref_metal_atom_idx, metal_info)

    if args.inspect_only:
        print("\n--inspect-only: stopping here.")
        return

    for required in ("traj", "top", "out_dir"):
        if getattr(args, required) is None:
            print(f"ERROR: --{required.replace('_', '-')} is required unless --inspect-only", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)

    dna_selector = " or ".join(f"chainid {c}" for c in dna_chain_ids)
    ref_prot_idx = ref.topology.select(f"chainid {protein_chain} and name CA")
    ref_dna_idx  = ref.topology.select(f"({dna_selector}) and not element H")

    print(f"\nReference Cα atoms (chainid {protein_chain}): {len(ref_prot_idx)}")
    print(f"Reference DNA heavy atoms (chainids {dna_chain_ids}): {len(ref_dna_idx)}")

    traj = md.load(args.traj, top=args.top)
    print(f"Loaded trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")

    traj_prot_ca_idx = traj.topology.select("name CA")
    print(f"Trajectory Cα atoms: {len(traj_prot_ca_idx)}")

    if len(traj_prot_ca_idx) != len(ref_prot_idx):
        print(
            f"ERROR: Cα count mismatch — traj has {len(traj_prot_ca_idx)}, ref has {len(ref_prot_idx)}.",
            file=sys.stderr,
        )
        sys.exit(1)

    # =========================================================================
    # CRITICAL: sequence-identity check
    # =========================================================================
    # A matching Cα count is not sufficient. If the crystal has a disordered
    # loop that BioEmu modeled across, the residue *positions* will line up
    # but the residue *identities* won't. This is exactly what bit us in DUX4
    # (5z6z had a 9-residue crystal gap that BioEmu filled in).
    if args.require_sequence_match:
        result = check_sequence_match(traj, ref, protein_chain, args.max_mismatches)
        n_mis, mismatches, traj_seq, ref_seq = result

        if n_mis is None:
            print(f"\nERROR: residue count mismatch", file=sys.stderr)
            sys.exit(1)

        if n_mis > args.max_mismatches:
            print(f"\n⚠ Sequence mismatch detected: {n_mis} residues differ between "
                  f"BioEmu trajectory and reference (max allowed: {args.max_mismatches})")
            print(f"First 5 mismatches (idx, traj, ref):")
            for i, t, r in mismatches[:5]:
                print(f"  position {i}: traj={t}, ref={r}")

            # Check if it's a contiguous gap-style mismatch (suggests crystal gap)
            if mismatches and mismatches[0][0] > 0:
                # Look at where mismatches cluster
                positions = [m[0] for m in mismatches]
                clusters = []
                cluster_start = positions[0]
                cluster_end = positions[0]
                for p in positions[1:]:
                    if p == cluster_end + 1:
                        cluster_end = p
                    else:
                        clusters.append((cluster_start, cluster_end))
                        cluster_start = cluster_end = p
                clusters.append((cluster_start, cluster_end))
                print(f"\nMismatch clusters: {clusters}")
                if len(clusters) <= 3 and any(c[1]-c[0] >= 3 for c in clusters):
                    print("  Pattern suggests crystal-gap modeling by BioEmu.")
                    print("  Consider using --mismatch-action trim to drop the gap-modeled region.")

            if args.mismatch_action == "fail":
                print("\nExiting (use --mismatch-action warn or trim to override)", file=sys.stderr)
                sys.exit(1)
            elif args.mismatch_action == "warn":
                print("\n⚠ Continuing despite mismatch (--mismatch-action warn). "
                      "Output may be biologically meaningless in mismatch regions.")
            elif args.mismatch_action == "trim":
                start_a, start_b, length = find_common_subsequence(traj_seq, ref_seq)
                print(f"\nTrimming to longest common run: "
                      f"traj[{start_a}:{start_a+length}], ref[{start_b}:{start_b+length}], "
                      f"length={length}")
                if length < 0.5 * len(traj_seq):
                    print(f"WARNING: trimmed region is less than half the full sequence; "
                          f"likely lots is being dropped.", file=sys.stderr)
                # Restrict the Cα indices we use for alignment
                traj_prot_ca_idx = traj_prot_ca_idx[start_a:start_a+length]
                ref_prot_idx = ref_prot_idx[start_b:start_b+length]
                print(f"After trim: traj Cα={len(traj_prot_ca_idx)}, ref Cα={len(ref_prot_idx)}")
        else:
            print(f"✓ Sequence match verified ({n_mis} mismatches within tolerance)")

    if len(ref_prot_idx) == 0 or len(ref_dna_idx) == 0:
        print("ERROR: zero atoms selected for protein or DNA after filtering", file=sys.stderr)
        sys.exit(1)

    # Metal detection already done above (before --inspect-only exit). We
    # already have ref_metal_atom_idx and metal_info in scope.

    ref_dna_traj = ref.atom_slice(ref_dna_idx)
    ref_prot_xyz = ref.xyz[0, ref_prot_idx, :]
    ref_dna_xyz  = ref.xyz[0, ref_dna_idx, :]
    if len(ref_metal_atom_idx) > 0:
        ref_metal_traj = ref.atom_slice(ref_metal_atom_idx)
        ref_metal_xyz  = ref.xyz[0, ref_metal_atom_idx, :]
    else:
        ref_metal_traj = None
        ref_metal_xyz  = None

    n_written, n_failed = 0, 0
    for i in range(traj.n_frames):
        try:
            frame_ca_xyz = traj.xyz[i, traj_prot_ca_idx, :]
            R, t = kabsch(ref_prot_xyz, frame_ca_xyz)

            # Apply the same Kabsch transform to DNA and metals
            dna_xyz_aligned = (R @ ref_dna_xyz.T).T + t
            ref_dna_aligned = ref_dna_traj[0]
            ref_dna_aligned.xyz = dna_xyz_aligned[np.newaxis, :, :]

            combined = traj[i].stack(ref_dna_aligned)

            if ref_metal_traj is not None:
                metal_xyz_aligned = (R @ ref_metal_xyz.T).T + t
                ref_metal_aligned = ref_metal_traj[0]
                ref_metal_aligned.xyz = metal_xyz_aligned[np.newaxis, :, :]
                combined = combined.stack(ref_metal_aligned)

            # If we trimmed, we want to keep the full BioEmu protein but only use
            # the trimmed Cα subset for alignment. Write the full protein frame +
            # aligned DNA. (Stage 3 will handle any residues we kept that aren't
            # in the reference.)
            out_path = os.path.join(args.out_dir, f"{args.pdb_id}_state_{i+1:03d}.pdb")
            combined.save_pdb(out_path)
            n_written += 1
        except Exception as e:
            print(f"Frame {i}: FAILED — {e}")
            n_failed += 1

    print(f"\nDone. Written: {n_written}, Failed: {n_failed}")
    print(f"Output: {args.out_dir}/{args.pdb_id}_state_*.pdb")


if __name__ == "__main__":
    main()