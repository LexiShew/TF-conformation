import csv
import numpy as np

# Load the dataset
# allow_pickle=True is required because the file contains nested lists/objects
data = np.load('jaspar_h11mo_cluster_wise_dna_containing_dataset.npy', allow_pickle=True)

# Each cluster is a list of entries shaped like:
#   ['7eg8_A', ['TAF1_HUMAN.H11MO.0.A', ...]]
# i.e. (PDB_chain, [motif_id, ...]). Motif IDs encode a UniProt-style name
# (e.g. TAF1_HUMAN) before the first dot.

rows = []
pdb_ids = set()

for cluster_idx, cluster in enumerate(data):
    for entry in cluster:
        pdb_chain = entry[0]
        motif_ids = list(entry[1]) if len(entry) > 1 and entry[1] is not None else []

        if '_' in pdb_chain:
            pdb_id, chain = pdb_chain.split('_', 1)
        else:
            pdb_id, chain = pdb_chain, ''
        pdb_id = pdb_id.lower()
        pdb_ids.add(pdb_id)

        if not motif_ids:
            rows.append([cluster_idx, pdb_id, chain, '', ''])
            continue

        for motif_id in motif_ids:
            uniprot_name = motif_id.split('.')[0] if isinstance(motif_id, str) else ''
            rows.append([cluster_idx, pdb_id, chain, motif_id, uniprot_name])

unique_pdbs = sorted(pdb_ids)

print(f"Successfully parsed {len(data)} clusters.")
print(f"Found {len(unique_pdbs)} unique PDB structures.")
print(f"Expanded to {len(rows)} (pdb_chain, motif) rows.")
print("-" * 30)
print("Preview of PDB IDs:", unique_pdbs[:10])

# Full export: one row per (pdb_chain, motif) pair.
# TODO: add TF class/family and Pfam family columns when that data is available.
with open("dataset_entries.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["cluster_index", "pdb_id", "chain", "motif_id", "uniprot_name"])
    writer.writerows(rows)
print(f"Saved {len(rows)} rows to 'dataset_entries.csv'.")

# Keep the simple unique-PDB list for batch downloads / lookups.
with open("pdb_ids.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["pdb_id"])
    for pid in unique_pdbs:
        writer.writerow([pid])
print(f"Saved {len(unique_pdbs)} unique PDB IDs to 'pdb_ids.csv'.")
