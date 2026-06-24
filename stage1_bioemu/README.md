# stage1_bioemu

Stage 1 of the conformation pipeline: generate a backbone ensemble with
**BioEmu** and rebuild full-atom side chains with **HPacker**, for **every
protein chain** of every DeepPBS structure.

This replaces the old `stage1_hpacker/`, which (a) only ran HPacker on
BioEmu samples that some earlier step had already produced, and (b) only
handled the single monomer chain. `stage1_bioemu` does the sampling *and* the
side-chain reconstruction in one stage, over all chains of multi-chain
complexes too.

## Inputs

The source-chain library, one dir per structure:

```
structures/source_chains/<PDB>_chains/
    <PDB>_chain<X>_protein.pdb     # one or more protein chains
    <PDB>_dna.pdb                  # (used later, by stage2_redock)
    <PDB>.cif                      # reference crystal structure
```

## Outputs

All chains from all structures land flat in one output directory,
`structures/stage1_bioemu_output/` (override with `STAGE1_OUTPUT_DIR`). Each
chain's dir name is globally unique, so there are no collisions:

```
structures/stage1_bioemu_output/<PDB>_chain<X>_conformations/
    topology.pdb                   # BioEmu backbone-only reference
    samples.xtc                    # BioEmu backbone-only trajectory
    samples_sidechain_rec.pdb/.xtc # HPacker full-atom reconstruction
```

> **Note on naming.** The old monomer pipeline wrote a single
> `<PDB>_conformations/` nested inside each `<PDB>_chains/` dir. We now sample
> *every* chain into a shared output dir, one `<PDB>_chain<X>_conformations/`
> per chain. Downstream stages (stage2_redock onward, and `BIOEMU_DIR` in
> `lib/common.sh`) still assume the old layout and will need updating when you
> wire them to this per-chain library.

DNA docking is **not** done here — that's `stage2_redock`, per complex.

## Running

Single structure (one GPU job):

```bash
sbatch stage1_bioemu/stage1_bioemu.sh \
    structures/source_chains/1a1f_chains 100
```

Whole library (one job per structure):

```bash
./stage1_bioemu/submit_stage1_bioemu.sh           # defaults to source_chains, 100 confs
./stage1_bioemu/submit_stage1_bioemu.sh structures/source_chains 200
```

Jobs are **resumable**: a chain whose `samples_sidechain_rec.xtc` already
exists is skipped, so re-submitting only fills in what's missing.

## Engine

Both scripts call `scripts/generate_monomer_confs.py --chains-dir <dir>
--all-chains --reconstruct-sidechains`. The `--all-chains` flag is what
iterates over every `*_chain*_protein.pdb` (the default mode processes only
the first chain).
