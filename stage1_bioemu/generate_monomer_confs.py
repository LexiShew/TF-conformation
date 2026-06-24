# Generate protein monomer conformations using BioEmu and (optionally) dock DNA.
#
# Processes ONE protein per invocation. Two mutually exclusive input modes:
#
#   1. --chains-dir <DIR>
#        <DIR> is a <PDB>_chains/ directory containing:
#          <PDB>_chain<X>_protein.pdb   -- protein chain from crystal structure
#          <PDB>_dna.pdb                -- DNA from crystal structure (optional)
#        Default output root: <DIR> itself (output nests inside <PDB>_chains/).
#        DNA docking runs when *_dna.pdb exists.
#
#   2. --pdb-file <FILE>
#        Path to a single protein PDB file. The PDB ID is parsed from the
#        filename (the part before the first underscore). Default output root:
#        current dir. DNA docking is skipped.
#
# In both modes, output is written to <output-root>/<PDB>_conformations/,
# containing:
#   topology.pdb                    -- reference backbone-only topology
#   samples.xtc                     -- backbone-only trajectory
#   samples_sidechain_rec.pdb/.xtc  -- (optional) reconstructed full-atom
#                                      structures (with --reconstruct-sidechains)
#   samples_md_equil.pdb/.xtc       -- (optional) MD-equilibrated full-atom
#                                      structures (with --md-equil, implies
#                                      sidechain reconstruction)
#   <PDB>_docked.pdb                -- single multi-state PDB with protein+DNA
#                                      per state (mode 1, when DNA exists)
#   <PDB>_docked.pml                -- PyMOL loader script that opens
#                                      <PDB>_docked.pdb and assigns per-state
#                                      secondary structure. Open with
#                                      `pymol <PDB>_docked.pml` to view with
#                                      cartoons rendered.
#
# BioEmu only generates protein backbones (N, Cα, C, O); side-chain atoms in
# samples.xtc are placeholders, not physically meaningful positions. Use
# --reconstruct-sidechains to run HPacker for proper side-chain placement
# before docking. This requires conda + CUDA12; see BioEmu's docs.
#
# By default, BioEmu's intermediate batch_*.npz files are removed after
# sampling and docked output is written as a single multi-state PDB. Pass
# --keep-intermediates to retain the .npz files (useful for resume) and write
# per-state PDB files under docked/ instead.
#
# --output-dir / -o overrides the default output root in either mode.

import argparse
import glob
import os
import re
import sys
import tempfile

from pymol import cmd
from bioemu.sample import main as sample


# ---------------------------------------------------------------------------
# Sequence + sampling
# ---------------------------------------------------------------------------

def get_protein_sequence(pdb_file):
    """Extract the single-letter amino acid sequence from a PDB file."""
    cmd.delete("all")
    cmd.load(pdb_file, "protein")
    fasta = cmd.get_fastastr("protein and polymer.protein")
    sequence = "".join(
        line for line in fasta.splitlines() if not line.startswith(">")
    )
    cmd.delete("all")
    return sequence


def generate_conformations(protein_pdb, num_conformations, out_dir,
                           cache_dir=None, keep_intermediates=False):
    """Run BioEmu on a protein PDB and write outputs to out_dir.

    Removes BioEmu's intermediate batch_*.npz files after sampling unless
    keep_intermediates is True.
    """
    os.makedirs(out_dir, exist_ok=True)

    print(f"Extracting sequence from: {protein_pdb}")
    seq = get_protein_sequence(protein_pdb)
    print(f"Sequence: {seq[:60]}... (length: {len(seq)})")

    sample(
        sequence=seq,
        num_samples=num_conformations,
        output_dir=out_dir,
        cache_embeds_dir=cache_dir,
    )

    if not keep_intermediates:
        npz_files = glob.glob(os.path.join(out_dir, "*.npz"))
        for npz in npz_files:
            os.remove(npz)
        if npz_files:
            print(f"Removed {len(npz_files)} intermediate .npz files")

    print(f"Generated {num_conformations} conformations in {out_dir}")


def reconstruct_sidechains(out_dir, md_equil=False):
    """Run HPacker side-chain reconstruction on the BioEmu output, optionally
    followed by short MD equilibration.

    Reads topology.pdb and samples.xtc from out_dir and writes:
        samples_sidechain_rec.pdb / .xtc   -- HPacker output (always)
        samples_md_equil.pdb / .xtc        -- MD-relaxed output (if md_equil)

    Returns (topology_path, trajectory_path) for the latest output (MD-relaxed
    if md_equil else sidechain-reconstructed), so the caller can use those for
    downstream docking instead of the backbone-only originals.

    Requires the bioemu[md] extras and a conda-installed HPacker environment.
    See BioEmu's setup_sidechain_relax.sh for installation.
    """
    # Imported lazily so the script still runs without bioemu[md] installed
    # for users who don't pass --reconstruct-sidechains.
    from bioemu.sidechain_relax import main as sidechain_relax

    topo_in = os.path.join(out_dir, "topology.pdb")
    xtc_in = os.path.join(out_dir, "samples.xtc")
    if not (os.path.exists(topo_in) and os.path.exists(xtc_in)):
        sys.exit(f"Cannot reconstruct side chains: missing topology/samples in {out_dir}")

    print(f"Reconstructing side chains via HPacker (md_equil={md_equil})...")
    sidechain_relax(
        pdb_path=topo_in,
        xtc_path=xtc_in,
        outpath=out_dir,
        md_equil=md_equil,
    )

    if md_equil:
        topo_out = os.path.join(out_dir, "samples_md_equil.pdb")
        xtc_out = os.path.join(out_dir, "samples_md_equil.xtc")
    else:
        topo_out = os.path.join(out_dir, "samples_sidechain_rec.pdb")
        xtc_out = os.path.join(out_dir, "samples_sidechain_rec.xtc")

    if not (os.path.exists(topo_out) and os.path.exists(xtc_out)):
        sys.exit(
            f"Sidechain reconstruction completed but expected outputs not found: "
            f"{topo_out}, {xtc_out}"
        )
    print(f"Reconstructed structures: {topo_out}, {xtc_out}")
    return topo_out, xtc_out


# ---------------------------------------------------------------------------
# DNA docking
# ---------------------------------------------------------------------------

def _reassign_chains(protein_obj, dna_obj):
    """Force protein onto chain A and DNA onto a non-clashing chain.

    Preserves DNA's existing chain IDs unless they clash with the protein
    (e.g. dsDNA already on chains B/C is left alone; DNA on chain A or with
    no chain ID is moved to B).
    """
    cmd.alter(protein_obj, "chain='A'")

    dna_chains = set()
    cmd.iterate(
        dna_obj, "dna_chains.add(chain)", space={"dna_chains": dna_chains}
    )
    if "A" in dna_chains or not dna_chains or dna_chains == {""}:
        cmd.alter(dna_obj, "chain='B'")

    cmd.sort()


def _atom_lines(pdb_path):
    """Yield ATOM/HETATM/TER lines from a single-state PDB, skipping headers
    and MODEL/ENDMDL/END records."""
    with open(pdb_path) as f:
        for line in f:
            tag = line[:6].rstrip()
            if tag in ("ATOM", "HETATM", "TER"):
                yield line


def _ss_header_lines(pdb_path):
    """Yield HELIX/SHEET records from a PDB file's header. PyMOL emits these
    when saving an object with ss labels assigned (via cmd.dss); we hoist
    them to the top of the multi-state file so secondary structure is
    preserved on reload."""
    with open(pdb_path) as f:
        for line in f:
            tag = line[:6].rstrip()
            if tag in ("HELIX", "SHEET"):
                yield line
            elif tag in ("ATOM", "HETATM", "MODEL"):
                # Headers come before atoms; we can stop scanning.
                break


def _write_multistate_pdb(num_states, confs_obj, dna_obj, out_path):
    """Write a multi-state PDB with protein conformation + DNA in every MODEL.

    For each state, runs cmd.dss on that state, saves protein and DNA to
    separate temp PDBs, then emits a MODEL block with per-state HELIX/SHEET
    records followed by protein and DNA ATOM lines. Per-state SS lets you
    track helix/sheet changes across conformations (e.g. unfolding events).
    """
    with tempfile.TemporaryDirectory() as tmpdir, open(out_path, "w") as out:
        # DNA is the same in every state; save it once.
        dna_path = os.path.join(tmpdir, "dna.pdb")
        cmd.save(dna_path, dna_obj, state=1)
        dna_lines = list(_atom_lines(dna_path))

        for state in range(1, num_states + 1):
            # Re-run dss on this specific state so SS labels reflect this
            # conformation's geometry, not whatever was assigned last.
            cmd.dss(confs_obj, state=state)

            prot_path = os.path.join(tmpdir, f"prot_{state}.pdb")
            cmd.save(prot_path, confs_obj, state=state)

            out.write(f"MODEL     {state}\n")
            for line in _ss_header_lines(prot_path):
                out.write(line)
            for line in _atom_lines(prot_path):
                out.write(line)
            for line in dna_lines:
                out.write(line)
            out.write("ENDMDL\n")
        out.write("END\n")


def _write_per_state_pdbs(num_states, confs_obj, dna_obj, docked_dir, pdb_id):
    """Write one PDB per state, each containing protein + DNA with per-state
    HELIX/SHEET headers reflecting that conformation's geometry."""
    os.makedirs(docked_dir, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        dna_path = os.path.join(tmpdir, "dna.pdb")
        cmd.save(dna_path, dna_obj, state=1)
        dna_lines = list(_atom_lines(dna_path))

        for state in range(1, num_states + 1):
            cmd.dss(confs_obj, state=state)
            prot_path = os.path.join(tmpdir, f"prot_{state}.pdb")
            cmd.save(prot_path, confs_obj, state=state)

            out_path = os.path.join(
                docked_dir, f"{pdb_id}_docked_{state:04d}.pdb"
            )
            with open(out_path, "w") as out:
                for line in _ss_header_lines(prot_path):
                    out.write(line)
                for line in _atom_lines(prot_path):
                    out.write(line)
                for line in dna_lines:
                    out.write(line)
                out.write("END\n")


def _write_loader_pml(out_dir, pdb_filename, pml_filename, object_name):
    """Write a companion .pml script that loads the multi-state PDB and runs
    dss on every state, ensuring per-state secondary structure is assigned
    in PyMOL regardless of how the parser handles per-MODEL HELIX/SHEET
    records. Users open this with `pymol <pml_filename>` instead of the .pdb.

    The script resolves the .pdb path relative to the .pml file's own
    location, so it works regardless of the user's current working directory.
    """
    pml_path = os.path.join(out_dir, pml_filename)
    with open(pml_path, "w") as f:
        f.write("python\n")
        f.write("import os\n")
        # __script__ is set by PyMOL to the absolute path of the .pml when
        # run via `pymol path/to/script.pml`. Fall back to cwd if unset
        # (e.g. when the script is @-included from another script).
        f.write("_script_dir = os.path.dirname(os.path.abspath(__script__)) "
                "if '__script__' in dir() else os.getcwd()\n")
        f.write(f"_pdb_path = os.path.join(_script_dir, {pdb_filename!r})\n")
        f.write(f"cmd.load(_pdb_path, {object_name!r})\n")
        f.write("# Assign secondary structure per state (geometry-based)\n")
        f.write(f"for s in range(1, cmd.count_states({object_name!r}) + 1):\n")
        f.write(f"    cmd.dss({object_name!r}, state=s)\n")
        f.write("python end\n")
        f.write("show cartoon\n")
    return pml_path


def dock_dna(out_dir, protein_pdb, dna_pdb, pdb_id,
             topo_path=None, xtc_path=None, multistate=True):
    """Dock DNA onto each BioEmu conformation via template-based superposition.

    For each conformation in the trajectory, aligns it to the crystal-structure
    protein. If `multistate` is True, writes a single multi-state
    <pdb_id>_docked.pdb containing protein+DNA in each MODEL with protein and
    DNA on separate chains. Otherwise writes per-state PDBs under
    out_dir/<pdb_id>_docked/.

    By default uses out_dir/topology.pdb and out_dir/samples.xtc. Pass explicit
    topo_path/xtc_path to dock onto sidechain-reconstructed or MD-equilibrated
    structures instead.
    """
    cmd.delete("all")

    if topo_path is None:
        topo_path = os.path.join(out_dir, "topology.pdb")
    if xtc_path is None:
        xtc_path = os.path.join(out_dir, "samples.xtc")
    if not os.path.exists(topo_path) or not os.path.exists(xtc_path):
        print(f"Missing {topo_path} or {xtc_path}, skipping docking.")
        return

    cmd.load(protein_pdb, "ref_protein")
    cmd.load(dna_pdb, "dna")
    cmd.load(topo_path, "confs")
    # Replace topology (state 1) with first trajectory frame; subsequent frames
    # extend states 2..N. Without state=1, load_traj appends after state 1,
    # leaving the topology as a stray first state.
    cmd.load_traj(xtc_path, "confs", state=1)
    num_states = cmd.count_states("confs")

    # Align every conformation onto the crystal-structure protein so each
    # state lands in the DNA's frame of reference.
    for state in range(1, num_states + 1):
        cmd.align(f"confs and state {state}", "ref_protein", mobile_state=state)

    # Make protein and DNA show as distinct chains in the sequence viewer.
    _reassign_chains("confs", "dna")

    # Secondary structure is assigned per-state inside the writer functions
    # so each MODEL gets HELIX/SHEET records reflecting its own geometry.

    print(f"Docking DNA onto {num_states} conformations...")
    if multistate:
        pdb_filename = f"{pdb_id}_docked.pdb"
        pml_filename = f"{pdb_id}_docked.pml"
        object_name = f"{pdb_id}_docked"
        out_path = os.path.join(out_dir, pdb_filename)
        _write_multistate_pdb(num_states, "confs", "dna", out_path)
        pml_path = _write_loader_pml(out_dir, pdb_filename, pml_filename, object_name)
        cmd.delete("all")
        print(f"Saved {num_states}-state docked complex to {out_path}")
        print(f"Companion loader written to {pml_path} (open with `pymol {pml_filename}`)")
    else:
        docked_dir = os.path.join(out_dir, f"{pdb_id}_docked")
        _write_per_state_pdbs(num_states, "confs", "dna", docked_dir, pdb_id)
        cmd.delete("all")
        print(f"Saved {num_states} docked complexes to {docked_dir}/")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description="Generate monomer conformations with BioEmu and optionally dock DNA."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--chains-dir", "-c",
        help="Path to a <PDB>_chains/ directory (mode 1).",
    )
    group.add_argument(
        "--pdb-file", "-p",
        help="Path to a single protein PDB file (mode 2).",
    )
    parser.add_argument(
        "--output-dir", "-o",
        help="Output root for <PDB>_conformations/. "
             "Defaults to --chains-dir itself, or cwd for --pdb-file.",
    )
    parser.add_argument(
        "-n", "--num-conformations", type=int, default=100,
        help="Number of conformations to generate (default: 100).",
    )
    parser.add_argument(
        "--cache-dir",
        default=os.environ.get("BIOEMU_CACHE_DIR"),
        help="Directory for BioEmu's ColabFold embeddings cache. "
             "Defaults to $BIOEMU_CACHE_DIR, then BioEmu's built-in default "
             "(~/.bioemu_embeds_cache). Point this at a filesystem with quota "
             "headroom if home is full.",
    )
    parser.add_argument(
        "--keep-intermediates", action="store_true",
        help="Keep BioEmu's batch_*.npz files (useful for resuming partial "
             "runs) and write per-state docked PDBs under docked/ instead of "
             "a single multi-state docked.pdb.",
    )
    parser.add_argument(
        "--reconstruct-sidechains", action="store_true",
        help="Run HPacker side-chain reconstruction on BioEmu's backbone-only "
             "output before docking. BioEmu only generates backbone atoms, so "
             "side chains in the raw samples are placeholders. Required for "
             "any analysis that depends on side-chain positions (interface "
             "contacts, DeepPBS, etc.). Requires bioemu[md] + conda + CUDA12.",
    )
    parser.add_argument(
        "--md-equil", action="store_true",
        help="After side-chain reconstruction, run a short MD equilibration to "
             "relax the structure. Implies --reconstruct-sidechains. Slower "
             "but gives more physically realistic structures.",
    )
    parser.add_argument(
        "--redock-only", action="store_true",
        help="Skip sampling; re-run only the DNA docking step using existing "
             "topology.pdb and samples.xtc in the output dir. Mode 1 only.",
    )
    parser.add_argument(
        "--all-chains", action="store_true",
        help="Chains-dir mode: run BioEmu (+ HPacker with --reconstruct-"
             "sidechains) on EVERY *_chain*_protein.pdb in the dir, not just "
             "the first. Each chain writes its own <PDB>_chain<X>_"
             "conformations/. DNA docking is skipped (that's a later, "
             "per-complex stage). Chains whose output already exists are "
             "skipped, so the job is resumable. Used by stage1_bioemu.",
    )
    args = parser.parse_args()

    if args.redock_only and args.pdb_file:
        sys.exit("--redock-only requires --chains-dir (no DNA in --pdb-file mode).")

    if args.all_chains and not args.chains_dir:
        sys.exit("--all-chains requires --chains-dir.")
    if args.all_chains and args.redock_only:
        sys.exit("--all-chains does not dock DNA, so --redock-only is moot.")

    # --md-equil implies --reconstruct-sidechains
    if args.md_equil:
        args.reconstruct_sidechains = True

    return args


def _run_chains_mode(args):
    chains_dir = os.path.abspath(args.chains_dir.rstrip("/"))
    if not os.path.isdir(chains_dir):
        sys.exit(f"Not a directory: {chains_dir}")

    pdb_id = os.path.basename(chains_dir).split("_")[0]
    protein_files = glob.glob(os.path.join(chains_dir, "*_chain*_protein.pdb"))
    if not protein_files:
        sys.exit(f"No *_chain*_protein.pdb file found in {chains_dir}")
    protein_pdb = protein_files[0]

    output_root = args.output_dir or chains_dir
    out_dir = os.path.join(output_root, f"{pdb_id}_conformations")

    if not args.redock_only:
        generate_conformations(
            protein_pdb, args.num_conformations, out_dir,
            args.cache_dir, args.keep_intermediates,
        )

    # Optional side-chain reconstruction. If enabled, downstream docking uses
    # the reconstructed structures so docked.pdb has full-atom protein.
    # In --redock-only mode, skip the (slow) reconstruction step but pick up
    # existing reconstructed files if they're on disk.
    topo_path, xtc_path = None, None
    if args.reconstruct_sidechains:
        if args.redock_only:
            suffix = "md_equil" if args.md_equil else "sidechain_rec"
            existing_topo = os.path.join(out_dir, f"samples_{suffix}.pdb")
            existing_xtc = os.path.join(out_dir, f"samples_{suffix}.xtc")
            if os.path.exists(existing_topo) and os.path.exists(existing_xtc):
                topo_path, xtc_path = existing_topo, existing_xtc
                print(f"Using existing reconstructed structures: {topo_path}")
            else:
                topo_path, xtc_path = reconstruct_sidechains(
                    out_dir, md_equil=args.md_equil
                )
        else:
            topo_path, xtc_path = reconstruct_sidechains(
                out_dir, md_equil=args.md_equil
            )

    dna_files = glob.glob(os.path.join(chains_dir, "*_dna.pdb"))
    if dna_files:
        dock_dna(
            out_dir, protein_pdb, dna_files[0], pdb_id,
            topo_path=topo_path, xtc_path=xtc_path,
            multistate=not args.keep_intermediates,
        )
    else:
        print(f"No DNA file found in {chains_dir}, skipping docking.")


def _run_all_chains_mode(args):
    """Run BioEmu (+ optional HPacker) on EVERY protein chain in a <PDB>_chains/
    dir, each to its own <PDB>_chain<X>_conformations/ dir.

    This is the "all chains, not just monomers" path used by stage1_bioemu:
    the default --chains-dir mode only processes protein_files[0] (fine for
    single-chain monomers), whereas multi-chain crystal structures need every
    chain sampled. DNA docking is intentionally skipped here — stage1 produces
    the per-chain conformation library; docking onto DNA is a later stage.

    Chains whose expected output already exists are skipped, so a partially
    completed batch (e.g. a pre-empted SLURM array task) resumes cleanly.
    """
    chains_dir = os.path.abspath(args.chains_dir.rstrip("/"))
    if not os.path.isdir(chains_dir):
        sys.exit(f"Not a directory: {chains_dir}")

    pdb_id = os.path.basename(chains_dir).split("_")[0]
    protein_files = sorted(
        glob.glob(os.path.join(chains_dir, "*_chain*_protein.pdb"))
    )
    if not protein_files:
        sys.exit(f"No *_chain*_protein.pdb file found in {chains_dir}")

    output_root = args.output_dir or chains_dir
    # The file written last in each chain's pipeline; its presence means done.
    done_marker = (
        "samples_md_equil.xtc" if args.md_equil
        else "samples_sidechain_rec.xtc" if args.reconstruct_sidechains
        else "samples.xtc"
    )

    print(f"[all-chains] {pdb_id}: {len(protein_files)} protein chain(s)")
    for protein_pdb in protein_files:
        base = os.path.basename(protein_pdb)
        # <PDB>_chain<X>_protein.pdb -> chain tag "chain<X>"
        m = re.search(r"_(chain[A-Za-z0-9]+)_protein\.pdb$", base)
        chain_tag = m.group(1) if m else os.path.splitext(base)[0]
        out_dir = os.path.join(output_root, f"{pdb_id}_{chain_tag}_conformations")

        if os.path.exists(os.path.join(out_dir, done_marker)):
            print(f"[all-chains] SKIP {pdb_id} {chain_tag}: {done_marker} exists")
            continue

        print(f"[all-chains] {pdb_id} {chain_tag} -> {out_dir}")
        generate_conformations(
            protein_pdb, args.num_conformations, out_dir,
            args.cache_dir, args.keep_intermediates,
        )
        if args.reconstruct_sidechains:
            reconstruct_sidechains(out_dir, md_equil=args.md_equil)


def _run_pdb_file_mode(args):
    pdb_path = os.path.abspath(args.pdb_file)
    if not os.path.isfile(pdb_path):
        sys.exit(f"Not a file: {pdb_path}")

    pdb_id = os.path.basename(pdb_path).split("_")[0]
    output_root = args.output_dir or os.getcwd()
    out_dir = os.path.join(output_root, f"{pdb_id}_conformations")

    generate_conformations(
        pdb_path, args.num_conformations, out_dir,
        args.cache_dir, args.keep_intermediates,
    )

    if args.reconstruct_sidechains:
        reconstruct_sidechains(out_dir, md_equil=args.md_equil)


if __name__ == "__main__":
    args = _parse_args()
    if args.all_chains:
        _run_all_chains_mode(args)
    elif args.chains_dir:
        _run_chains_mode(args)
    else:
        _run_pdb_file_mode(args)