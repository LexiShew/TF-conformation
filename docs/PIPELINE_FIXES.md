# TF-conformation pipeline — fix spec for Claude Code

## Context

Pipeline: per-TF DAG that (Stage 1) samples protein backbone ensembles with **BioEmu**
and rebuilds side chains with **HPacker**, (Stage 2) docks each conformation onto the
crystal DNA, (Stage 3+) minimizes, preprocesses, and builds augmented training data for
**DeepPBS** (structure → binding-specificity PWM). Goal of the project: augment DeepPBS
training with conformational ensembles of monomeric TF–DNA complexes.

**Root problem:** Stage 1 was refactored from a single-monomer layout to a per-chain
"sample every protein chain" library, but `lib/common.sh`, the Stage 2 wrapper, and the
pilot configs still encode the OLD single-monomer layout. The two halves are disconnected.
The three pilots that "work" (1aay, 1tgh, 5z6z) only work because their Stage 1 outputs were
pre-populated by the old pipeline and pointed at by hand. A genuinely new TF run through the
new Stage 1 will not connect to Stage 2.

Files referenced (paths relative to repo root `TF-conformation/`):
- `lib/common.sh` — path resolution + pilot config loader
- `config/pilots/<tf>.sh` — per-TF config (sets PDB_ID, FOLD, PWM_LABEL, PROTEIN_CHAIN, DNA_CHAINS, …)
- `stage1_bioemu/stage1_bioemu.sh` — Stage 1 SLURM logic
- `stage1_bioemu/generate_monomer_confs.py` — Stage 1 engine (BioEmu + HPacker; also contains a legacy PyMOL docker)
- `stage1_bioemu/submit_stage1_bioemu.sh` — fan-out submitter
- `stage2_redock/stage2_redock.sh` — Stage 2 SLURM logic (sourced by wrapper)
- `stage2_redock/stage2_redock.py` — Stage 2 engine (mdtraj Kabsch dock)
- `wrappers/stage{1,2}_*.sh` — SBATCH wrappers
- Scoring/plotting (separate tree): `TF_conf_init_outputs/scripts/interface_rmsd.py`, `plot_interface_metrics.py`

---

## INVARIANTS — do not regress these

1. **Frame preservation.** The dock must keep the BioEmu protein in (or aligned to) the
   crystal protein's coordinate frame and carry the crystal DNA in that same frame. The
   canonical Stage 2 engine (`stage2_redock.py`) does this via a per-frame Kabsch transform
   applied to the protein-fit AND to the carried DNA/metals. fnat is scored with
   `--use_model_dna` (model protein vs the DNA embedded in the same docked file).
2. **Interface alignment is canonical, not global.** Global all-Cα Kabsch caps fnat at
   ~0.47 even for sub-Å folds (it averages placement error across the whole chain).
   Interface-Cα alignment raises fnat to ~0.72. `--align-mode interface` is the default.
3. **`per_domain` mode is a DIAGNOSTIC, never the batch default.** It fits each domain
   independently and discards inter-domain geometry — only valid for multidomain folds whose
   subdomains contact independent DNA subsites (e.g. C2H2 arrays). For cooperative/tandem
   binders it launders a real failure into deceptively high fnat.
4. **Keep the existing Stage 2 safety logic:** sequence-match/trim guard, metal carry-across,
   `--inspect-only`.
5. **Monomer-only scope.** The project filters to functionally monomeric TFs. The pipeline
   must enforce this, not silently dock chain 0 of an assembly.

---

## BLOCKERS (fix before the next new-TF run)

### B1 — Stage 1 output location/naming does not match what Stage 2 reads
- **Stage 1 writes:** `output/stage1_bioemu/<PDB>_chain<X>_conformations/`
  containing `topology.pdb`, `samples.xtc`, `samples_sidechain_rec.{pdb,xtc}`
  (see `generate_monomer_confs.py::_run_all_chains_mode` and `reconstruct_sidechains`).
- **`common.sh` says:** `STAGE1_DIR="${WORK_DIR}/stage1_relax"` (under
  `DeepPBS/data/conformations/<TF>/`) and
  `BIOEMU_DIR="${BIOEMU_RAW_ROOT}/${PDB_ID}_chains/${PDB_ID}_conformations"` — both OLD paths.
- **`stage2_redock/stage2_redock.sh` reads:** `${STAGE1_DIR}/${PDB_ID}_sidechain_rec.xtc`
  and `..._sidechain_rec.pdb` — wrong prefix (`${PDB_ID}_` vs `samples_`) AND wrong dir.
- **Fix:** Pick ONE canonical Stage 1 output location and naming, then update `common.sh`
  (`STAGE1_DIR`, `BIOEMU_DIR`) and the Stage 2 wrapper's `--traj/--top` to derive from it,
  using the `samples_sidechain_rec.*` filenames. Per-chain dir must include the chain token
  (see B2).
- **Accept:** From a clean run, `stage2_redock.sh` resolves `--traj/--top` to existing files
  with no hand-editing.

### B2 — No "which chain binds DNA" selector (chain identity dropped)
- Stage 1 now produces `<PDB>_chain<X>_conformations/` for EVERY protein chain.
  `common.sh` keys everything on `${PDB_ID}` with no chain token, so Stage 2 cannot know which
  chain's ensemble to use. Separately, `stage2_redock.py --protein-chain` is a **0-based mdtraj
  chainid into the .cif**, whereas Stage 1 tags dirs by the **source-filename chain letter**
  (`chain<X>`). These two indexings are not guaranteed to agree.
- **Fix:** Add `BINDING_CHAIN` to each pilot config and thread it through:
  (a) select the correct `<PDB>_chain<BINDING_CHAIN>_conformations/` dir, and
  (b) pass the matching 0-based `--protein-chain` to `stage2_redock.py`.
  Establish and document one mapping from chain-letter tag → .cif 0-based chainid; assert they
  agree (fail loudly if not).
- **Accept:** A pilot config sets `BINDING_CHAIN`; Stage 2 docks that chain's ensemble against
  the DNA using the matching `.cif` chain, verified by a non-empty interface + sane fnat.

### B3 — Monomer assumption is unenforced
- `_run_all_chains_mode` samples all chains; Stage 2 docks a single `--protein-chain` against
  the full DNA. For a multi-protein-chain assembly this silently docks one chain's ensemble
  and labels it with the complex PWM (the 5z6z/dimer failure mode) — no error.
- **Fix:** Add an explicit assertion (in Stage 2 or a pre-flight check) that the complex is a
  functional monomer: exactly one protein chain in contact with the DNA site. Refuse/skip
  multi-protein-chain assemblies with a clear message rather than docking chain 0.
- **Accept:** Running a known assembly (e.g. a dimer) aborts with a clear "not a monomer"
  message; monomers proceed.

### B4 — Two dockers; the inferior one is still reachable
- `generate_monomer_confs.py::dock_dna` is a SECOND docker using PyMOL `cmd.align` (global
  superposition) writing `<PDB>_docked.pdb`. The interface/per-domain fix lives ONLY in
  `stage2_redock.py`. `--all-chains` mode skips `dock_dna` (good), but plain `--chains-dir`
  mode still calls it — a foot-gun that produces globally-docked (fnat ~0.47) states.
- **Fix:** Make `stage2_redock.py` the sole production docker. Either gate `dock_dna` behind an
  explicit `--legacy-pymol-dock` flag (off by default) or remove it. Ensure no production path
  calls it.
- **Accept:** No production code path produces `<PDB>_docked.pdb` via PyMOL; all docked states
  come from `stage2_redock.py`.

### B5 — `per_domain` mode is missing from the repo copy of `stage2_redock.py`
- The committed `stage2_redock.py` has only `interface` and `all` modes. The `per_domain` mode
  (segment interface Cα by residue-number gap, fit each domain independently, anchor DNA to the
  largest interface domain) is needed for the C2H2 / tandem diagnostics.
- **Fix:** Add `--align-mode per_domain` and `--domain-gap` (default 10). Keep `interface` the
  default. See INVARIANT 3 — diagnostic only.
- **Accept:** `--align-mode per_domain` runs on a multidomain case and prints the per-domain
  Cα split; `interface`/`all` outputs are byte-identical to before.

### B6 — `--align-mode` not exposed through the SLURM wrapper
- `stage2_redock/stage2_redock.sh` never passes `--align-mode`, so the batch path always uses
  the default and there's no way to run `all`/`per_domain` baselines.
- **Fix:** Add a `STAGE2_ALIGN_MODE` env passthrough (mirror the existing `MISMATCH_ACTION`
  pattern), default `interface`.
- **Accept:** `STAGE2_ALIGN_MODE=all sbatch …` runs the global baseline.

### B7 — fnat gate before training data is built — DONE
- **Implemented:** a single structural-quality gate at **Stage 3 → Stage 4** (`fnat_gate/`,
  vendored scorer `interface_rmsd.py` + `score_stage3.py`). It scores every **post-minimization**
  state (`fnat` vs the model's own DNA) and builds `${STAGE3_DIR}_pass/` — a symlink mirror of
  only the states with `fnat >= FNAT_FLOOR`. Stage 4 reads ONLY that pass dir.
- **Why Stage 3, not Stage 2:** minimization moves per-state fnat both ways, so the gate must
  score the minimized pose. Stage 2 now carries EVERY docked state forward — it never drops,
  moves, or quarantines (the earlier `rejected_fnat/` + `fnat_drops.log` Stage-2 rejection was
  removed; it shrank the denominator before minimization could rescue near-misses).
- **Floor:** `FNAT_FLOOR=0.5` default (`lib/common.sh`), overridable globally or per-pilot
  (`config/pilots/<tf>.sh`).
- **Fail-loud:** empty pass-list ⇒ gate exits non-zero; Stage 4 depends `afterok` on the gate,
  so the DAG halts (no training on an empty/ungated set).
- **Accept:** Stage 4 consumes only passing states; a per-state fnat CSV + pass-list are produced;
  `FNAT_FLOOR=0.99` halts the DAG. Verified locally on FOXA (65/72 minimized states pass @ 0.5;
  0/72 @ 0.99 → exit 1). Full-ensemble yield requires the cluster re-run (see below).

---

## SHOULD-FIX / CLEANUP

- **S1 — `--no-md-equil` flag does not exist.** `generate_monomer_confs.py` has only
  `--md-equil` (absent = no equilibration). Any wrapper still passing `--no-md-equil` (the old
  `stage1_hpacker.sh` did) will argparse-error. Grep for `--no-md-equil` and remove.
- **S2 — Guard required pilot vars.** In `common.sh::load_pilot_config`, add `require_var` for
  `PDB_ID`, `FOLD`, `PWM_LABEL` (and now `BINDING_CHAIN`) so a malformed pilot config fails with
  a clear message instead of a confusing downstream path error.
- **S3 — Verify `REF_CIF` exists.** `common.sh` sets
  `REF_CIF="${BIOEMU_RAW_ROOT}/${PDB_ID}_chains/${PDB_ID}.cif"`. Confirm the `source_chains`
  refactor actually copied each `.cif` into its `<PDB>_chains/` dir; add an existence check.
- **S4 — Single-chain input assert.** `get_protein_sequence` concatenates all-chain FASTA with
  no separator; harmless in all-chains mode (one chain per file) but would silently fuse
  sequences on a multi-chain PDB. Assert the input `*_chain*_protein.pdb` is single-chain.
- **S5 — Comment the intentional LEGACY non-suffix on `STAGE2_DIR`.** Docking is metal-
  independent, so the legacy A/B intentionally does NOT suffix `STAGE2_DIR`. Add a one-line
  comment so it isn't "fixed" later.
- **S6 — `STAGE2_INSPECT` uses `return 0`.** Correct because the script is `source`d, but add a
  guard so running it directly doesn't silently no-op.
- **S7 — Scoring `--tag` must not contain commas.** `interface_rmsd.py` emits CSV; a tag like
  `1tgh,1` injects a column. Use `pdb_id_state` style tags (the current code sanitizes commas to
  `_`; keep that).

---

## Suggested order
B1 + B2 together (same root cause), then B3, then B4/B5/B6 (Stage 2 surface), then B7 (gate).
Cleanups S1–S7 anytime. After B1–B3, run one NEW pilot (e.g. engrailed 1hdd, a true monomeric
homeodomain) end-to-end Stage 1 → Stage 2 → fnat with zero hand-editing as the acceptance test.

## Definition of done
A new monomeric pilot, added only via a `config/pilots/<tf>.sh` (with `BINDING_CHAIN`), runs
Stage 1 → Stage 2 → fnat scoring with no manual path edits; Stage 2 uses interface alignment;
fnat is scored with `--use_model_dna`; sub-floor states are filtered before Stage 5; and an
attempt to add a non-monomeric assembly aborts with a clear message.
