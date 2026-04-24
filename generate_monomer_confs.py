# Generate protein monomer conformations using BioEmu and dock DNA onto them.
#
# Given a PDB ID, this script:
#   1. Reads the protein chain PDB file from monomers/<PDB>_chains/
#   2. Extracts the amino acid sequence using PyMOL
#   3. Runs BioEmu to generate conformational samples (saved as topology.pdb + samples.xtc)
#   4. If a DNA PDB file exists, docks the DNA onto each conformation by aligning
#      the conformation to the crystal structure protein and saving combined PDB files
#
# Expected directory structure:
#   monomers/<PDB>_chains/
#     <PDB>_chain<X>_protein.pdb   -- protein chain extracted from crystal structure
#     <PDB>_dna.pdb                -- DNA extracted from crystal structure (optional)
#
# Output (written to the same <PDB>_chains/ directory):
#   topology.pdb    -- reference topology for the generated conformations
#   samples.xtc     -- trajectory of generated conformations
#   docked/         -- combined protein+DNA PDB files for each conformation (if DNA exists)
#
# Usage: python generate_monomer_confs.py <PDB_ID> [NUM_CONFORMATIONS]
#   PDB_ID             -- 4-letter PDB identifier (case-insensitive)
#   NUM_CONFORMATIONS  -- number of conformations to generate (default: 100)
#
# This script is called by run_conf.sh via SLURM. To submit all monomers at once,
# use submit_all.sh.

import os
import sys
import glob
from pymol import cmd
from bioemu.sample import main as sample


def get_protein_sequence(pdb_file):
    """Extract the amino acid sequence from a local PDB file.

    Loads the PDB into PyMOL, extracts the FASTA sequence for protein
    polymers, and returns it as a plain string.

    Args:
        pdb_file: Path to a PDB file containing a protein structure.

    Returns:
        Single-letter amino acid sequence string.
    """
    cmd.delete("all")
    cmd.load(pdb_file, "protein")
    fasta = cmd.get_fastastr("protein and polymer.protein")
    sequence = "".join(line for line in fasta.splitlines() if not line.startswith(">"))
    cmd.delete("all")
    return sequence


def dock_dna(chains_dir, protein_pdb, dna_pdb):
    """Dock DNA onto each BioEmu conformation via template-based superposition.

    For each conformation in the samples trajectory, aligns it to the crystal
    structure protein, then saves the conformation + DNA as a combined PDB file
    in a docked/ subdirectory.

    Args:
        chains_dir: Path to the <PDB>_chains/ directory containing
            topology.pdb and samples.xtc from BioEmu.
        protein_pdb: Path to the crystal structure protein PDB file
            (used as the alignment reference).
        dna_pdb: Path to the crystal structure DNA PDB file.
    """
    cmd.delete("all")

    topo_path = os.path.join(chains_dir, "topology.pdb")
    xtc_path = os.path.join(chains_dir, "samples.xtc")

    if not os.path.exists(topo_path) or not os.path.exists(xtc_path):
        print(f"Missing topology.pdb or samples.xtc in {chains_dir}, skipping docking.")
        return

    cmd.load(protein_pdb, "ref_protein")
    cmd.load(dna_pdb, "dna")

    cmd.load(topo_path, "confs")
    cmd.load_traj(xtc_path, "confs")
    num_states = cmd.count_states("confs")

    docked_dir = os.path.join(chains_dir, "docked")
    os.makedirs(docked_dir, exist_ok=True)

    print(f"Docking DNA onto {num_states} conformations...")
    for state in range(1, num_states + 1):
        cmd.align(f"confs and state {state}", "ref_protein", mobile_state=state)
        out_path = os.path.join(docked_dir, f"docked_{state:04d}.pdb")
        cmd.save(out_path, "confs or dna", state=state)

    cmd.delete("all")
    print(f"Saved {num_states} docked complexes to {docked_dir}/")


if __name__ == "__main__":
    pdb_id = sys.argv[1] if len(sys.argv) > 1 else sys.exit("Usage: python generate_monomer_confs.py <PDB_ID> [NUM_CONFORMATIONS]")
    num_conformations = int(sys.argv[2]) if len(sys.argv) > 2 else 100

    script_dir = os.path.dirname(os.path.abspath(__file__))
    chains_dir = os.path.join(script_dir, "monomers", f"{pdb_id.upper()}_chains")

    # find protein and DNA PDB files
    protein_files = glob.glob(os.path.join(chains_dir, "*_chain*_protein.pdb"))
    dna_files = glob.glob(os.path.join(chains_dir, "*_dna.pdb"))

    if not protein_files:
        sys.exit(f"No protein PDB file found in {chains_dir}")

    protein_pdb = protein_files[0]
    print(f"Extracting sequence from: {protein_pdb}")
    seq = get_protein_sequence(protein_pdb)
    print(f"Sequence for {pdb_id}: {seq[:60]}... (length: {len(seq)})")

    sample(sequence=seq, num_samples=num_conformations, output_dir=chains_dir)

    if dna_files:
        dock_dna(chains_dir, protein_pdb, dna_files[0])
    else:
        print(f"No DNA file found in {chains_dir}, skipping docking.")

    print(f"Generated {num_conformations} conformations for {pdb_id} in {chains_dir}")
