# Retired pilots

Pilots withdrawn from the study, with the reason and the scope of removal. Any number in an
older document that involves a retired pilot is **withdrawn** — the docs are kept as a record
of what was actually run, not corrected in place.

---

## dux4 (5z6z) — retired 2026-10-15

**Why.** Two independent disqualifiers, both already noted in the code before retirement:

1. **Not a monomer.** DUX4 binds DNA as a homodimer, violating the pipeline's core
   TRUE MONOMER assumption (Stage 2's monomer guard, `README.md` Stage 2). This was known —
   `analysis/common/fig_common.py` carried `NONMONOMER = {"dux4"}` and labelled the family
   "Homeodomain (dimer)".
2. **It never produced a result.** `0 fnat survivors`, so no training data, no
   `output/stage7_eval/id_benchmark_dux4.json`, no benchmark row at any point. Recorded in
   `make_dnaflex.py`, `make_samefamily.py`, and `config/pilots/dux4.sh` itself
   ("Kept only for reference").

So nothing downstream depended on it, and nothing is lost by removing it. Its only footprint
was in structural figures (RMSD progression, ensemble diversity), where it appeared as a
contrast case.

## irf (1if1) — retired 2026-10-15

**Why.** The source structure is not usable with this pipeline:

1. **Its duplex resolves as two helical segments.** 3DNA reports
   `helix count problem 2` for 1if1 — 139 occurrences across the stage-4 logs. The
   `assert len(dna_helices) == 1` in `process_co_crystal.py` then skips the state. The result
   was **91–94% silent attrition**: 78 gated states → 5 npz (frozen arm), 74 → 7 (relaxed).
   The pilot trained on 5 conformers. `config/pilots/ndt80.sh` already referred to this as
   the "IRF-style missing strand" when justifying its own source choice.
2. **The loss was biased, not random.** The failure triggers on states whose DNA deviates
   enough from a single ideal helix to be split in two, so it preferentially removed the
   most-deformed conformations — the diversity the project is about. The 5 survivors are the
   most crystal-like frames.
3. **No usable on-target test.** `TEST_PWM_FILTER="MA1419|IRF"` matched 4 entries that are
   all four chains of a *single* PDB (7jm4). Effective structural n = 1, not 4.

`analysis/README.md` reported IRF as one of only two net-positive families
("IRF (n=4, +0.039)"). That number was 5 training conformers scored on 4 chains of one
crystal, and is withdrawn.

---

## What was removed

**Configs** — `config/pilots/{dux4,dux4_dnarelax,irf,irf_dnarelax}.sh`

**Code** — dux4/irf dropped from every hardcoded pilot list (23 files): the seven
`conf_vs_crystal_pwm/` scripts, `dna_relax/{interface_mgwfl,pycurves_batch/*}`,
`diversity/{af3_rmsd_diversity,montage_ensembles,plot_diversity,render_ensembles}`,
`conformation/{render_prepost_min,make_M}`,
`rmsd_analysis/{compute_rmsds,plot_rmsd_progression,identify_interface_residues,rmsd_regen.sbatch}`,
`mechanism/{make_dnaflex,make_dnaflex_labeled,make_samefamily,make_mechanism_fig}`,
`af3/{run_af3_batch2,af3_pycurves}.sh`, `stage7_eval/predict_ab_pwms.py`,
`scripts/structure_viz/render_state_stages.sh`, `palette.py`, `analysis/common/fig_common.py`.

`fig_common.py` now carries `RETIRED = {"dux4", "irf"}` so that any old CSV still holding
those rows renders a label instead of raising `KeyError`. They must not appear in new figures.
`ndt80` took the palette slot dux4 vacated.

**Data** — 53 CSVs filtered on their pilot column; notable row drops:
`analysis/data/perseed_perentry.csv` −2,600, `ensemble_diversity_pairwise.csv` −8,926,
`rmsd_analysis/per_state_rmsds.csv` −1,311, `stage3/data/traj_sample.csv` −2,172,
`perentry_accuracy.csv` −130.

**Outputs** — ~112 MB across `output/stage{1,2,3,4,5,6,7}*`, `structures/source_chains/`,
`af3/output/`, including 42 `stage6_train` configs and all run dirs.

**Other** — `af3/input/{dux4_5z6z,irf_1if1}.json`,
`analysis/analyses/rmsd_analysis/interface_residues/dux4.json`, the `5z6z_*` /`dux4_*`
pycurves viz files, `stage7_eval/confcmp_idfiles/id_irf_confcmp*.txt`,
`stage5_build_aug/folds_aug/train0_aug{,_legacy}_dux4.txt`.

---

## What was deliberately KEPT, and why

Retiring a pilot removes it as an **augmentation source**. It does not touch the shared
130-entry evaluation benchmark or the base DeepPBS training set — doing so would change what
every *other* pilot is measured against and invalidate all 264 existing checkpoints.

| kept | what it is |
|---|---|
| 9 entries in `data/folds/id.txt` | DUX4 (`6u81`, `6a8r`, `5z2t`, `6dfy` ×2) and IRF (`7jm4` ×4) **benchmark** entries — different structures from the retired pilot sources |
| `1if1_A_IRF1_MOUSE…npz` in train0/1/3/4 + valid2 | the 1if1 **crystal** as a base training entry (1 entry, not 100 conformers) |
| `5z6z_A_DUX4_HUMAN…npz` in train0/2/3/4 + valid1 | likewise for 5z6z |
| `structures/source_chains/2irf_chains/` | **2IRF is a different PDB** (IRF2), benchmark data — not the retired `1if1` pilot |
| `*DUX4*` / `*IRF*` npz inside other pilots' `combined_assembly_*` | base DeepPBS training entries that every pilot's assembly includes |
| `analysis/common/fig_common.py: RETIRED` | label fallback for historical CSVs |

The distinction is checkable: pilot augmentation is always named `<pdb>_state_NNN_<pwm>.npz`.
Zero such files remain for `1if1` or `5z6z`. Base crystal entries are `<pdb>_<chain>_<pwm>.npz`
and are untouched.

---

## Remaining pilots (12)

`csl · egr1 · engrailed · err · ets1 · foxa · hsf · lef1 · ndt80 · nfat · runx · tbp`

Three of these have no usable on-target subset and now declare it explicitly via
`ON_TARGET_USABLE=no`: **engrailed** (own motif absent from the benchmark; a 10-entry
Homeodomain family-transfer subset is available if opted into), **lef1** (the only HMG-box in
the benchmark — no family fallback exists), and **tbp** (n=1, and that one entry is a histone
chain; see `docs/tbp_label_error_downstream.md`). **foxa**'s filter was corrected from a
0-match regex to a 5-entry Forkhead family-transfer subset.
