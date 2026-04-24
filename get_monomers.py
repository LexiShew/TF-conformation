# Fetch PDB structures from RCSB and split them into separate monomer files.
#
# For each PDB ID, this script:
#   1. Fetches the structure from RCSB PDB
#   2. Saves all DNA nucleotides as a single file (<PDB>_dna.pdb, if DNA is present)
#   3. Saves each protein chain as a separate file (<PDB>_chain<X>_protein.pdb)
#
# Output is written to monomers/<PDB>_chains/ (one directory per structure),
# resolved relative to this script so it works from any CWD.
#
# Usage:
#   python get_monomers.py                   # run on the default hardcoded list below
#   python get_monomers.py <pdb_ids.txt>     # run on a file with one PDB ID per line
#   python get_monomers.py 1CIT 6PAX 1TC3    # run on PDB IDs passed on the command line

from pymol import cmd
import os
import sys

# Default list of PDB IDs (DNA-binding proteins of interest, dedup-sorted).
DEFAULT_MONOMERS = [
    "1A02", "1A1F", "1A3Q", "1A66", "1A6Y", "1AM9", "1AN2", "1AN4", "1AU7", "1AWC",
    "1B72", "1BC7", "1BF5", "1BG1", "1BY4", "1C7U", "1CDW", "1CF7", "1CIT", "1DH3",
    "1DP7", "1DSZ", "1DUX", "1E3O", "1GAT", "1GJI", "1GLU", "1GT0", "1GTW", "1H6F",
    "1HCQ", "1HJC", "1HLV", "1HRY", "1IC8", "1IF1", "1IG7", "1IMH", "1IV6", "1JNM",
    "1K78", "1KB2", "1LE5", "1LFU", "1LO1", "1MDY", "1MHD", "1MSE", "1N6J", "1NWQ",
    "1ODH", "1PUE", "1PUF", "1R4I", "1SRS", "1T2K", "1TF3", "1TSR", "1UBD", "1VFC",
    "1VTN", "1XBR", "1YO5", "2A07", "2A66", "2C6Y", "2C7A", "2D5V", "2EZD", "2F8X",
    "2FF0", "2GLI", "2H1K", "2H8R", "2HDC", "2IRF", "2JP9", "2JX1", "2KMK", "2KO0",
    "2LD5", "2LEF", "2LKX", "2LT7", "2ME6", "2O6G", "2OEH", "2QL2", "2UZK", "2WBS",
    "2WTY", "2X6V", "2XSD", "3A5T", "3CBB", "3CMY", "3CO6", "3D1N", "3DFV", "3DZU",
    "3F27", "3G73", "3GNA", "3GUT", "3JTG", "3L1P", "3L2C", "3MLN", "3MVA", "3QMI",
    "3QRF", "3QSV", "3QYM", "3RKQ", "3U2B", "3UK3", "3VD0", "3ZP5", "4A04", "4A75",
    "4ATI", "4BNC", "4BQA", "4DA4", "4EOT", "4EUW", "4GZN", "4H10", "4HP1", "4IRI",
    "4J19", "4MHG", "4PZI", "4RBO", "4RDU", "4TNT", "4UNO", "4UUV", "4XRM", "4XRS",
    "4Y5W", "4Y60", "4YJ0", "4YO2", "4ZKG", "4ZPK", "4ZPR", "5D5U", "5D8K", "5E8I",
    "5EGB", "5FD3", "6PAX",
]


def get_monomers(pdb_ids, out_root):
    """Fetch PDB structures and split into separate DNA and protein chain files.

    For each PDB ID, creates <out_root>/<PDB>_chains/ and saves:
      - <PDB>_dna.pdb: all DNA nucleotides (if present)
      - <PDB>_chain<X>_protein.pdb: one file per protein chain

    Args:
        pdb_ids: Iterable of 4-letter PDB identifiers to fetch.
        out_root: Directory under which <PDB>_chains/ subdirectories are created.
    """
    os.makedirs(out_root, exist_ok=True)
    start_dir = os.getcwd()
    for pdb_id in pdb_ids:
        dir_name = os.path.join(out_root, f"{pdb_id}_chains")
        os.makedirs(dir_name, exist_ok=True)
        os.chdir(dir_name)

        cmd.fetch(pdb_id)
        obj_name = cmd.get_object_list()[0]

        if cmd.count_atoms(f"{obj_name} and resn DA+DT+DG+DC") > 0:
            cmd.save(f"{obj_name}_dna.pdb", f"{obj_name} and resn DA+DT+DG+DC")

        for chain in cmd.get_chains(obj_name):
            sel_name = f"chain_{chain}"
            cmd.select(sel_name, f"{obj_name} and chain {chain} and polymer.protein")
            if cmd.count_atoms(sel_name) > 0:
                cmd.save(f"{obj_name}_chain{chain}_protein.pdb", sel_name)
            cmd.delete(sel_name)

        cmd.delete(obj_name)
        os.chdir(start_dir)


def _parse_args(argv):
    """Resolve the list of PDB IDs from CLI args, a file, or the default list."""
    if len(argv) == 1:
        return DEFAULT_MONOMERS
    if len(argv) == 2 and os.path.isfile(argv[1]):
        with open(argv[1]) as f:
            return [line.strip() for line in f if line.strip() and not line.startswith("#")]
    return argv[1:]


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_root = os.path.join(script_dir, "monomers")
    pdb_ids = _parse_args(sys.argv)
    print(f"Fetching {len(pdb_ids)} structure(s) into {out_root}/")
    get_monomers(pdb_ids, out_root)
