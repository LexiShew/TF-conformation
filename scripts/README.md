# scripts/

Standalone utilities that are **not** part of the staged pilot pipeline
(`stage1_bioemu` … `stage7_eval`). Grouped by purpose:

| Subfolder | Contents | Purpose |
|---|---|---|
| `analysis/` | `compute_rmsds.py`, `plot_rmsd.py`, `generate_structure_images.py` | Per-conformation RMSD computation + plots + rendered structure images. |
| `pymol/` | `color_protein_states.py`, `gradient_protein_states.py`, `gradient_protein_split.py`, `spectrum_states.py` | PyMOL helpers that color multi-state objects by conformation index. |
| `classification/` | `classify_pfam.py`, `classify_DBDs.py` | Assign Pfam families / DNA-binding-domain classes to structures. |
| `pdb_prep/` | `process_pdb_chains.py`, `find_monomers.sh` | Split raw PDBs into protein/DNA chains; identify monomers. |
| `deprecated/` | `run_conf.sh`, `submit_all.sh` | **Superseded by `stage1_bioemu/`.** These called the old single-chain `generate_monomer_confs.py` (now moved into `stage1_bioemu/` and generalized to all chains via `--all-chains`). Kept for reference only; do not use. |

These are invoked by hand (`python scripts/analysis/compute_rmsds.py …`), not
by the SLURM pipeline. Cross-references between scripts are documentation only
(no imports), so the grouping above is purely organizational.
