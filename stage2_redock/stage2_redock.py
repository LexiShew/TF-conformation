"""
Stage 2: align each Stage 1 frame onto a reference structure, carry the
reference DNA and any structural metal ions across, write per-state PDBs.

ALIGNMENT MODE (new):
  --align-mode interface (default) : Kabsch-fit on the DNA-contacting protein
      Cα only. Global all-Cα fitting minimizes whole-chain RMSD, which spreads
      placement error onto the interface and systematically loses native
      protein-DNA contacts (observed: fnat capped ~0.47 even for sub-Å folds).
      Fitting on interface Cα places the recognition residues on their cognate
      DNA contacts instead of averaging error across the whole chain.
  --align-mode all : original behaviour (all-Cα fit), kept for A/B comparison.
  --align-mode per_domain : DIAGNOSTIC ONLY. Segments the interface Cα into
      domains by residue-number gap, fits each domain independently, and anchors
      the carried DNA to the largest interface domain. Valid only for multidomain
      folds whose subdomains contact independent DNA subsites (e.g. C2H2 arrays);
      for cooperative/tandem binders it discards inter-domain geometry and
      launders a real failure into deceptively high fnat. Never the batch default.

The BioEmu protein conformation is left rigid in all modes; only the
superposition target changes.

Monomer scope: by default Stage 2 refuses any assembly where more than one
protein chain contacts the DNA (pass --allow-multimer to override), so it never
silently docks one chain's ensemble against a multi-protein complex site.

Carrying metals: BioEmu doesn't sample metals and HPACKER doesn't restore
them. We extract metal positions from the reference crystal and apply the same
Kabsch transform used for the reference DNA, writing them as HETATM so Stage 3
sees them.

Validates that BioEmu's protein and the reference crystal's protein are the
same residue sequence (not just the same Cα count); a mere count-match misses
crystal gaps BioEmu modeled across (the DUX4 5z6z case).
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
STRUCTURAL_METAL_RESNAMES = {"ZN", "MG", "MN", "FE", "CA", "CO", "NI", "CU"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb-id", help="Used in output filenames, e.g. '1aay'. "
                                    "Required unless --inspect-only.")
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
                        "fail = exit (default), warn = continue, "
                        "trim = slice both to the longest common-sequence run")
    p.add_argument("--max-mismatches", type=int, default=0,
                   help="Maximum allowed sequence mismatches before action triggers.")
    # --- alignment controls ---
    p.add_argument("--align-mode", choices=["interface", "all", "per_domain"], default="interface",
                   help="Kabsch fit on DNA-contacting Cα ('interface', default), "
                        "all Cα ('all', original behaviour), or 'per_domain' "
                        "(DIAGNOSTIC ONLY — fit each interface domain independently, "
                        "anchor DNA to the largest; valid only for multidomain folds "
                        "whose subdomains contact independent DNA subsites, e.g. C2H2 "
                        "arrays — see module docstring).")
    p.add_argument("--iface-cutoff", type=float, default=5.0,
                   help="Interface defn for --align-mode interface/per_domain: protein "
                        "residue within this many Å of any DNA heavy atom (default 5.0).")
    p.add_argument("--domain-gap", type=int, default=10,
                   help="--align-mode per_domain: split interface Cα into domains "
                        "wherever consecutive interface residues' resSeq differ by more "
                        "than this (default 10).")
    p.add_argument("--allow-multimer", action="store_true",
                   help="Skip the monomer guard (B3). By default Stage 2 refuses an "
                        "assembly with >1 protein chain contacting the DNA, to avoid "
                        "silently docking one chain's ensemble against a complex site.")
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
    if name in HISTIDINE_VARIANTS: return "HIS"
    return name


def check_sequence_match(traj, ref, protein_chain, max_mismatches):
    traj_residues = [normalize_resname(r.name) for r in traj.topology.residues]
    ref_chain = list(ref.topology.chains)[protein_chain]
    ref_residues = [normalize_resname(r.name) for r in ref_chain.residues]
    if len(traj_residues) != len(ref_residues):
        return None, None, traj_residues, ref_residues
    mismatches = [(i, t, r) for i, (t, r) in enumerate(zip(traj_residues, ref_residues)) if t != r]
    return len(mismatches), mismatches, traj_residues, ref_residues


def find_reference_metals(ref):
    ref_metal_atom_idx, metal_info = [], []
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


def kabsch_rmsd(P, Q):
    """RMSD between P and Q after the optimal Kabsch superposition of P onto Q."""
    R, t = kabsch(P, Q)
    Pa = (R @ P.T).T + t
    return float(np.sqrt(((Pa - Q) ** 2).sum(axis=1).mean()))


def dna_contacting_protein_chains(ref, ref_dna_idx, cutoff_ang):
    """0-based chain indices of protein chains with any heavy atom within cutoff
    of the reference DNA. Used by the monomer guard (B3)."""
    cutoff_nm = cutoff_ang / 10.0
    xyz = ref.xyz[0]
    dna_xyz = xyz[ref_dna_idx]
    contacting = []
    for chain in ref.topology.chains:
        n_prot = sum(1 for r in chain.residues if r.name in PROTEIN_RESNAMES)
        if n_prot < 10:
            continue
        heavy = ref.topology.select(f"chainid {chain.index} and not element H")
        if len(heavy) == 0:
            continue
        p = xyz[heavy]
        dmin = np.sqrt(((p[:, None, :] - dna_xyz[None, :, :]) ** 2).sum(-1)).min()
        if dmin <= cutoff_nm:
            contacting.append(chain.index)
    return contacting


def segment_into_domains(ref, ref_prot_idx, iface_pos, gap):
    """Split interface positions into domains wherever consecutive interface
    residues' resSeq jump by more than `gap`. Returns a list of position-lists
    (positions index into ref_prot_idx), ordered by resSeq."""
    items = sorted(
        ((ref.topology.atom(int(ref_prot_idx[k])).residue.resSeq, k) for k in iface_pos),
        key=lambda x: x[0],
    )
    domains, cur, prev = [], [], None
    for resseq, k in items:
        if prev is not None and resseq - prev > gap:
            domains.append(cur)
            cur = []
        cur.append(k)
        prev = resseq
    if cur:
        domains.append(cur)
    return domains


def interface_positions(ref, protein_chain, ref_prot_idx, ref_dna_idx, cutoff_ang):
    """Positions within ref_prot_idx whose residue is within cutoff of DNA.

    ref_prot_idx is the (possibly trimmed) array of protein Cα atom indices, in
    residue order. Returns a list of positions k into that array. mdtraj xyz is
    in nanometres, so the Å cutoff is converted here.
    """
    cutoff_nm = cutoff_ang / 10.0
    xyz = ref.xyz[0]
    dna_xyz = xyz[ref_dna_idx]
    # residues that contact DNA (by residue topology index)
    iface_resids = set()
    for res in list(ref.topology.chains)[protein_chain].residues:
        heavy = [a.index for a in res.atoms if a.element.symbol != "H"]
        if not heavy:
            continue
        p = xyz[heavy]
        dmin = np.sqrt(((p[:, None, :] - dna_xyz[None, :, :]) ** 2).sum(-1)).min()
        if dmin <= cutoff_nm:
            iface_resids.add(res.index)
    # map back to positions within the aligned Cα array
    pos = [k for k in range(len(ref_prot_idx))
           if ref.topology.atom(int(ref_prot_idx[k])).residue.index in iface_resids]
    return pos


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

    ref_metal_atom_idx, metal_info = find_reference_metals(ref)
    report_reference_metals(ref_metal_atom_idx, metal_info)

    if args.inspect_only:
        print("\n--inspect-only: stopping here.")
        return

    for required in ("pdb_id", "traj", "top", "out_dir"):
        if getattr(args, required) is None:
            print(f"ERROR: --{required.replace('_', '-')} is required unless --inspect-only", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)

    dna_selector = " or ".join(f"chainid {c}" for c in dna_chain_ids)
    ref_prot_idx = ref.topology.select(f"chainid {protein_chain} and name CA")
    ref_dna_idx  = ref.topology.select(f"({dna_selector}) and not element H")

    print(f"\nReference Cα atoms (chainid {protein_chain}): {len(ref_prot_idx)}")
    print(f"Reference DNA heavy atoms (chainids {dna_chain_ids}): {len(ref_dna_idx)}")

    # --- Monomer guard (B3) ---
    contacting = dna_contacting_protein_chains(ref, ref_dna_idx, args.iface_cutoff)
    print(f"Protein chains contacting DNA (within {args.iface_cutoff} Å): {contacting}")
    if len(contacting) > 1 and not args.allow_multimer:
        print(f"\nERROR: not a functional monomer — {len(contacting)} protein chains "
              f"contact the DNA site (chainids {contacting}).", file=sys.stderr)
        print("This pipeline targets monomeric TF–DNA complexes; docking a single "
              "chain's ensemble against a multi-protein site would mislabel the data.",
              file=sys.stderr)
        print("Pass --allow-multimer to override (diagnostic only).", file=sys.stderr)
        sys.exit(1)

    traj = md.load(args.traj, top=args.top)
    print(f"Loaded trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")

    traj_prot_ca_idx = traj.topology.select("name CA")
    print(f"Trajectory Cα atoms: {len(traj_prot_ca_idx)}")

    if len(traj_prot_ca_idx) != len(ref_prot_idx):
        print(f"ERROR: Cα count mismatch — traj has {len(traj_prot_ca_idx)}, ref has {len(ref_prot_idx)}.",
              file=sys.stderr)
        sys.exit(1)

    # --- sequence-identity check (unchanged) ---
    if args.require_sequence_match:
        result = check_sequence_match(traj, ref, protein_chain, args.max_mismatches)
        n_mis, mismatches, traj_seq, ref_seq = result
        if n_mis is None:
            print(f"\nERROR: residue count mismatch", file=sys.stderr); sys.exit(1)
        if n_mis > args.max_mismatches:
            print(f"\n⚠ Sequence mismatch: {n_mis} residues differ (max allowed: {args.max_mismatches})")
            for i, t, r in mismatches[:5]:
                print(f"  position {i}: traj={t}, ref={r}")
            if args.mismatch_action == "fail":
                print("\nExiting (use --mismatch-action warn or trim to override)", file=sys.stderr)
                sys.exit(1)
            elif args.mismatch_action == "warn":
                print("\n⚠ Continuing despite mismatch (--mismatch-action warn).")
            elif args.mismatch_action == "trim":
                start_a, start_b, length = find_common_subsequence(traj_seq, ref_seq)
                print(f"\nTrimming to longest common run: traj[{start_a}:{start_a+length}], "
                      f"ref[{start_b}:{start_b+length}], length={length}")
                if length < 0.5 * len(traj_seq):
                    print("WARNING: trimmed region < half the sequence.", file=sys.stderr)
                traj_prot_ca_idx = traj_prot_ca_idx[start_a:start_a+length]
                ref_prot_idx = ref_prot_idx[start_b:start_b+length]
                print(f"After trim: traj Cα={len(traj_prot_ca_idx)}, ref Cα={len(ref_prot_idx)}")
        else:
            print(f"✓ Sequence match verified ({n_mis} mismatches within tolerance)")

    if len(ref_prot_idx) == 0 or len(ref_dna_idx) == 0:
        print("ERROR: zero atoms selected for protein or DNA after filtering", file=sys.stderr)
        sys.exit(1)

    ref_dna_traj = ref.atom_slice(ref_dna_idx)
    ref_prot_xyz = ref.xyz[0, ref_prot_idx, :]
    ref_dna_xyz  = ref.xyz[0, ref_dna_idx, :]
    if len(ref_metal_atom_idx) > 0:
        ref_metal_traj = ref.atom_slice(ref_metal_atom_idx)
        ref_metal_xyz  = ref.xyz[0, ref_metal_atom_idx, :]
    else:
        ref_metal_traj = None
        ref_metal_xyz  = None

    # =========================================================================
    # Choose alignment atom subset: interface Cα (default) or all Cα.
    # =========================================================================
    domains = None  # set only in per_domain mode, for per-frame diagnostics
    if args.align_mode == "interface":
        align_pos = interface_positions(ref, protein_chain, ref_prot_idx,
                                        ref_dna_idx, args.iface_cutoff)
        if len(align_pos) < 3:
            print(f"⚠ Only {len(align_pos)} interface Cα found at {args.iface_cutoff} Å; "
                  f"falling back to all-Cα alignment.", file=sys.stderr)
            align_pos = list(range(len(ref_prot_idx)))
            mode_used = "all (fallback)"
        else:
            mode_used = "interface"
    elif args.align_mode == "per_domain":
        iface = interface_positions(ref, protein_chain, ref_prot_idx,
                                    ref_dna_idx, args.iface_cutoff)
        if len(iface) < 3:
            print(f"⚠ Only {len(iface)} interface Cα; per_domain falling back to all-Cα.",
                  file=sys.stderr)
            align_pos = list(range(len(ref_prot_idx)))
            mode_used = "all (fallback)"
        else:
            domains = segment_into_domains(ref, ref_prot_idx, iface, args.domain_gap)
            print(f"\nper_domain: {len(domains)} interface domain(s) "
                  f"(split where consecutive interface resSeq gap > {args.domain_gap}):")
            for di, dom in enumerate(domains):
                rs = [ref.topology.atom(int(ref_prot_idx[k])).residue.resSeq for k in dom]
                print(f"  domain {di}: {len(dom)} Cα, resSeq {min(rs)}-{max(rs)}")
            largest = max(domains, key=len)
            align_pos = largest
            mode_used = (f"per_domain (DNA anchored to largest of {len(domains)} "
                         f"domains, {len(largest)} Cα)")
    else:
        align_pos = list(range(len(ref_prot_idx)))
        mode_used = "all"
    align_pos = np.asarray(align_pos, dtype=int)
    ref_align_xyz = ref_prot_xyz[align_pos]
    print(f"\nAlignment mode: {mode_used} — fitting on {len(align_pos)}/{len(ref_prot_idx)} Cα")
    per_domain_rmsd = ([[] for _ in domains] if domains is not None else None)

    n_written, n_failed = 0, 0
    for i in range(traj.n_frames):
        try:
            frame_ca_xyz = traj.xyz[i, traj_prot_ca_idx, :]
            # Kabsch on the chosen subset only; transform applies to everything.
            R, t = kabsch(ref_align_xyz, frame_ca_xyz[align_pos])

            # per_domain diagnostic: how well each interface domain fits on its own.
            if per_domain_rmsd is not None:
                for di, dom in enumerate(domains):
                    dpos = np.asarray(dom, dtype=int)
                    per_domain_rmsd[di].append(
                        kabsch_rmsd(ref_prot_xyz[dpos], frame_ca_xyz[dpos]) * 10.0)

            dna_xyz_aligned = (R @ ref_dna_xyz.T).T + t
            ref_dna_aligned = ref_dna_traj[0]
            ref_dna_aligned.xyz = dna_xyz_aligned[np.newaxis, :, :]
            combined = traj[i].stack(ref_dna_aligned)

            if ref_metal_traj is not None:
                metal_xyz_aligned = (R @ ref_metal_xyz.T).T + t
                ref_metal_aligned = ref_metal_traj[0]
                ref_metal_aligned.xyz = metal_xyz_aligned[np.newaxis, :, :]
                combined = combined.stack(ref_metal_aligned)

            out_path = os.path.join(args.out_dir, f"{args.pdb_id}_state_{i+1:03d}.pdb")
            combined.save_pdb(out_path)
            n_written += 1
        except Exception as e:
            print(f"Frame {i}: FAILED — {e}")
            n_failed += 1

    if per_domain_rmsd is not None:
        print("\nper_domain independent-fit Cα RMSD (Å), mean over frames:")
        for di, vals in enumerate(per_domain_rmsd):
            if vals:
                print(f"  domain {di}: {np.mean(vals):.2f} (min {np.min(vals):.2f}, "
                      f"max {np.max(vals):.2f})")

    print(f"\nDone. Written: {n_written}, Failed: {n_failed}")
    print(f"Output: {args.out_dir}/{args.pdb_id}_state_*.pdb")


if __name__ == "__main__":
    main()