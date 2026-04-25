# BioEmu DNA-Binding Protein Pipeline

Generate conformational ensembles for DNA-binding protein monomers using
[BioEmu](https://github.com/microsoft/bioemu), then dock the crystal-structure
DNA onto each sampled conformation by template-based superposition on the
reference protein.

## Pipeline overview

```
          RCSB PDB
              │
              │  get_monomers.py
              ▼
   monomers/<PDB>_chains/
     <PDB>_chain<X>_protein.pdb   ← reference crystal chain
     <PDB>_dna.pdb                ← crystal DNA (if present)
              │
              │  generate_monomer_confs.py  (submitted via run_conf.sh / submit_all.sh)
              ▼
   monomers/<PDB>_chains/
     topology.pdb, samples.xtc    ← BioEmu conformations
     batch_*.npz                  ← BioEmu intermediate batches
     docked/docked_0001.pdb, ...  ← conformation + DNA, aligned on reference
```

## Scripts

All scripts live at the top level. `monomers/` is data-only.

| Script | Purpose |
| --- | --- |
| `get_monomers.py` | Fetch PDB structures from RCSB and split each into one `<PDB>_dna.pdb` plus one `<PDB>_chain<X>_protein.pdb` per protein chain, under `monomers/<PDB>_chains/`. |
| `generate_monomer_confs.py` | For a single PDB ID: read the protein chain from `monomers/<PDB>_chains/`, extract its FASTA sequence via PyMOL, run BioEmu sampling, then dock the crystal DNA onto each conformation. |
| `run_conf.sh` | SLURM batch wrapper that runs `generate_monomer_confs.py` for one PDB ID on one GPU. |
| `submit_all.sh` | Submit one `run_conf.sh` job per `monomers/<PDB>_chains/` directory. |

### Requirements

- `bioemu` Python package (with its CUDA/JAX dependencies)
- `pymol` Python module (used for sequence extraction, chain splitting, DNA superposition)
- SLURM cluster with a GPU partition if you want to use `run_conf.sh` /
  `submit_all.sh` as-is. The scripts currently request:
  `-p rohs --account=rohs_102 --gres=gpu:rtx5000:1 --time=24:00:00 -n 8 -N 1`.
  Edit the `#SBATCH` lines in `run_conf.sh` for other clusters.

## How to run

### 1. Fetch and split structures

```bash
# Use the built-in default list of ~157 DNA-binding proteins:
python get_monomers.py

# ...or pass PDB IDs directly:
python get_monomers.py 1CIT 6PAX 1TC3

# ...or read PDB IDs from a file (one per line, # comments allowed):
python get_monomers.py my_targets.txt
```

Creates `monomers/<PDB>_chains/` for each ID, populated with the reference
protein and DNA PDB files.

### 2. Generate conformations + dock DNA

**One protein:**
```bash
python generate_monomer_confs.py <PDB_ID> [NUM_CONFORMATIONS]   # default: 100
```

**One protein on SLURM:**
```bash
sbatch --job-name="<PDB_ID>_<N>" run_conf.sh <PDB_ID> <N>
```

**All proteins in `monomers/` on SLURM:**
```bash
./submit_all.sh [NUM_CONFORMATIONS]   # default: 100
```

Each job writes into its own `monomers/<PDB>_chains/` directory, so the jobs
don't contend for output paths.

### 3. Re-running partial / failed jobs

BioEmu resumes from the `batch_*.npz` checkpoints already in the output
directory, so re-submitting the same `run_conf.sh` for a PDB ID will pick up
where it left off. If you want to re-run from scratch, delete the `batch_*.npz`
files (and `topology.pdb`, `samples.xtc`, `docked/`) first.

## What's currently in `monomers/`

- **65** `<PDB>_chains/` directories
- **55** complete (have `topology.pdb`, `samples.xtc`, and a populated
  `docked/`)
- **10** incomplete — only the reference protein/DNA files and partial
  `batch_*.npz` checkpoints, no `topology.pdb`/`samples.xtc` yet. Re-submit
  `run_conf.sh` for these to finish them. They are:
  `1A1F, 1A66, 1AN2, 1BC7, 1BF5, 1BG1, 1CDW, 1CIT, 1DP7, 1E3O`.

## Other artifacts at the top level

These are older standalone outputs kept for reference — the canonical pipeline
output now lives inside `monomers/<PDB>_chains/`.

| Path | What it is |
| --- | --- |
| `1cit_output/`, `1tc3_output/`, `1skh_output/`, `6pax_output/` | Early BioEmu runs (pre-refactor) with `topology.pdb`, `samples.xtc`, and batch NPZ files. `6pax_output/` additionally contains a `docked/` directory with 96 docked complexes. |
| `combined_docked.pse` | PyMOL session with multiple docked complexes loaded together (visualization only). |
| `1cit.cif`, `1tc3.cif`, `6pax.cif`, `pdb1cit.ent` | Raw crystal structures fetched earlier; equivalents now live in `monomers/<PDB>_chains/` when the corresponding PDB ID has been processed. |
| `topology.pdb`, `samples.xtc` | Stray BioEmu outputs from an early test run (kept for reproducibility; not referenced by the current pipeline). |
