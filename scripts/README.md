# scripts/

Scripts that live outside the per-stage directories. Grouped by purpose:

| Subfolder | Contents | Purpose |
|---|---|---|
| `pipeline/` | `run_pilot.sh`, `run_multiseed_pilot.sh`, `run_legacy_ab.sh` | **Pipeline launchers.** Submit the per-TF SLURM DAG (Stage 1 → 7) via `wrappers/`. Self-locate `TFCONF_DIR`; run from anywhere, e.g. `./scripts/pipeline/run_pilot.sh egr1`. |
| `analysis/` | `compute_rmsds.py`, `plot_rmsd.py`, `generate_structure_images.py` | Per-conformation RMSD computation + plots + rendered structure images. |
| `pymol/` | `color_protein_states.py`, `gradient_protein_states.py`, `gradient_protein_split.py`, `spectrum_states.py` | PyMOL helpers that color multi-state objects by conformation index. |
| `classification/` | `classify_pfam.py`, `classify_DBDs.py` | Assign Pfam families / DNA-binding-domain classes to structures. |
| `pdb_prep/` | `process_pdb_chains.py`, `find_monomers.sh` | Split raw PDBs into protein/DNA chains; identify monomers. |
| `maintenance/` | `finish_cleanup.sh` | One-off repo-maintenance utility (gitignore/LFS hygiene); hardcodes its repo path. |
| `deprecated/` | `run_conf.sh`, `submit_all.sh` | **Superseded by `stage1_bioemu/`.** These called the old single-chain `generate_monomer_confs.py` (now moved into `stage1_bioemu/` and generalized to all chains via `--all-chains`). Kept for reference only; do not use. |

The launchers in `pipeline/` ARE the pipeline entry points; everything else
here is invoked by hand (e.g. `python scripts/analysis/compute_rmsds.py …`),
not by the SLURM jobs. Cross-references between the utility scripts are
documentation only (no imports), so the grouping is purely organizational.
