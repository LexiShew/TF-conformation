#!/usr/bin/env python
"""Diagnose why common_atom_indices returns zero matches across stages."""
import mdtraj as md
import sys

PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
}

# Use EGR1 / 1aay as the test case
REF_CIF = "/project2/rohs_102/shewchuk/TF-conformation/deeppbs_pdbs/monomer_chains/1aay_chains/1aay.cif"
BIOEMU_TOPO = "/project2/rohs_102/shewchuk/TF-conformation/deeppbs_pdbs/monomer_chains/1aay_chains/1aay_conformations/topology.pdb"
S2_PDB = "/project2/rohs_102/shewchuk/DeepPBS/data/conformations/egr1/stage2_docked/1aay_state_001.pdb"
S3_PDB = "/project2/rohs_102/shewchuk/DeepPBS/data/conformations/egr1/stage3_min/1aay_state_001.pdb"

def report(label, traj, protein_chain=None):
    print(f"\n=== {label} ===")
    print(f"  n_atoms={traj.n_atoms}, n_residues={traj.topology.n_residues}, "
          f"n_chains={traj.topology.n_chains}")
    chains = list(traj.topology.chains)
    for chain in chains[:5]:
        res = list(chain.residues)
        if not res:
            continue
        # Filter to protein only for the reporting
        prot_res = [r for r in res if r.name in PROTEIN_RESNAMES]
        if not prot_res:
            kind = "non-protein"
            print(f"  chain {chain.index} (id='{chain.chain_id if hasattr(chain, 'chain_id') else '?'}'): "
                  f"{kind}, {len(res)} residues, first={res[0].name}")
            continue
        first = prot_res[0]
        last = prot_res[-1]
        sample_atoms = [a.name for a in first.atoms][:8]
        print(f"  chain {chain.index} (id='{chain.chain_id if hasattr(chain, 'chain_id') else '?'}'): "
              f"protein, {len(prot_res)} prot residues, "
              f"resSeq range {first.resSeq}–{last.resSeq}, "
              f"first residue: {first.name} resSeq={first.resSeq} sample_atoms={sample_atoms}")

ref = md.load(REF_CIF)
report("Reference 1aay.cif", ref)

s0 = md.load(BIOEMU_TOPO)
report("Stage 0 BioEmu topology.pdb", s0)

s2 = md.load(S2_PDB)
report("Stage 2 docked", s2)

s3 = md.load(S3_PDB)
report("Stage 3 minimized", s3)

# Specifically: do the residue keys (resSeq, resname, atomname) intersect?
def get_keys(traj, restrict_chain=None):
    keys = set()
    for atom in traj.topology.atoms:
        if atom.residue.name not in PROTEIN_RESNAMES:
            continue
        if atom.element is None or atom.element.symbol == "H":
            continue
        if restrict_chain is not None and atom.residue.chain.index != restrict_chain:
            continue
        keys.add((atom.residue.resSeq, atom.residue.name, atom.name))
    return keys

ref_chain2 = get_keys(ref, restrict_chain=2)  # 1aay protein is chain 2
s2_keys = get_keys(s2)
s3_keys = get_keys(s3)
print(f"\nReference protein-chain keys: {len(ref_chain2)}")
print(f"Stage 2 protein keys:         {len(s2_keys)}")
print(f"Stage 3 protein keys:         {len(s3_keys)}")
print(f"Intersection ref ∩ s2: {len(ref_chain2 & s2_keys)}")
print(f"Intersection ref ∩ s3: {len(ref_chain2 & s3_keys)}")
print(f"Intersection s2 ∩ s3:  {len(s2_keys & s3_keys)}")

# Show first 5 keys from each set for inspection
print("\nFirst 5 ref keys:")
for k in sorted(ref_chain2)[:5]:
    print(f"  {k}")
print("First 5 stage 2 keys:")
for k in sorted(s2_keys)[:5]:
    print(f"  {k}")
print("First 5 stage 3 keys:")
for k in sorted(s3_keys)[:5]:
    print(f"  {k}")