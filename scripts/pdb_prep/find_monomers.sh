
SRC_DIR="deeppbs_pdbs/pdb_chains"
DEST_DIR="deeppbs_pdbs/monomer_chains"

mkdir -p "$DEST_DIR"

for chain_dir in "$SRC_DIR"/*_chains; do
    # Skip if no matching directories exist
    [[ -d "$chain_dir" ]] || continue

    protein_count=$(find "$chain_dir" -maxdepth 1 -type f -name "*_protein.pdb" | wc -l)

    if [[ "$protein_count" -eq 1 ]]; then
        echo "Copying: $chain_dir"
        cp -R "$chain_dir" "$DEST_DIR/"
    else
        echo "Skipping: $chain_dir ($protein_count protein PDB files)"
    fi
done