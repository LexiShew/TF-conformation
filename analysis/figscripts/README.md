# Figure regeneration — F / I / R / S / P / M series

Pilot-agnostic scripts that rebuild the book's structural/accuracy figure series
from pipeline output. **Nothing hardcodes the pilot list** — pilots are
discovered from disk (`fig_common.discover_pilots`), so a newly-run TF appears in
every figure automatically once its pipeline output and eval JSON exist.

Location: `analysis/figscripts/` in the repo. Palette: `palette.py` at the repo
root (one hue per entity — see `COLOR_CONSISTENCY.md`).

## Files

| File | What it does |
|---|---|
| `fig_common.py` | Shared foundation: pilot auto-discovery, pilot→PDB/family maps (`PILOT_META`), pass-rate ordering, labels, palette import, `savefig()`. **Edit `PILOT_META` when you add a pilot.** |
| `extract_metrics.py` | Pipeline output → metric CSVs in `analysis/data/`: `perstate_metrics.csv`, `perentry_accuracy.csv`, `mechanism_apo_holo.csv`. Auto-discovers pilots. |
| `compute_reachability.py` | The coordinate pass that fills M1's reachability axis (`d_min`/`reach_ratio`/`rmsf_mean` → `reachability.csv`) and R2's per-residue profile (`ca_rmsd_perresidue.csv`), for all pilots. Reuses `compute_rmsds.py`'s loaders. **`pycurves` env.** |
| `make_F.py` | F1 fnat distributions · F2 pass-rate bars · F3 fnat-vs-iRMSD · F4 interface-size-vs-fidelity |
| `make_I.py` | I1 iRMSD distributions · I2 seg max-vs-mean (localized distortion) · I4 interface geometry |
| `make_P.py` | P1 pilot→family map · P2 baseline accuracy by family · P3 augmentation effect by family (⭐ key figure) |
| `make_R.py` | R1 Cα-RMSD stage2-vs-stage3 · R2 per-residue profiles · R3 minimization delta |
| `make_M.py` | M1 apo/holo mechanism (reachability × spread, aug-sign colored) |
| `make_S.py` + `render_S.sbatch` | PyMOL renders: S1 best/worst-fnat frame, S2 crystal/docked/minimized. **Must run via sbatch on a compute node.** |

## Environment

All matplotlib scripts run in the **`deeppbs`** conda env; `compute_rmsds.py`
(R-data regen) needs **`pycurves`** (it imports `mdtraj`, absent from `deeppbs`);
`make_S.py` needs **`pymol`**. On the login node always cap BLAS threads first:

```bash
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate deeppbs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Regeneration workflow

```bash
cd analysis/figscripts

# 1. (re)build the metric CSVs for ALL discovered pilots
python extract_metrics.py
#    → analysis/data/{perstate_metrics,perentry_accuracy,mechanism_apo_holo}.csv
#    (family_annotation.csv is NOT rebuilt — the 130-entry benchmark is fixed and
#     its motif→family map is pilot-independent.)

# 2. the cheap matplotlib series (seconds each), all pilots:
python make_F.py
python make_I.py
python make_P.py
python make_M.py
#    → analysis/figures/{F1..F4,I1,I2,I4,P1,P2,P3,M1}.png
```

### R series — needs a Cα-RMSD regen first
`make_R.py` reads `rmsd_analysis/per_state_rmsds.csv`. That file is produced by
the already-pilot-agnostic `rmsd_analysis/compute_rmsds.py`. To cover all pilots:

```bash
conda activate pycurves          # mdtraj lives here, not in deeppbs
cd rmsd_analysis
python compute_rmsds.py --tfs csl dux4 egr1 engrailed err ets1 foxa hsf irf lef1 nfat runx tbp \
    --interface-dir interface_residues --output per_state_rmsds.csv
#    ~2 min on a rohs node; submit via srun if the login node is loaded:
#    srun -p rohs -A rohs_102 -c 4 -t 00:30:00 python compute_rmsds.py --tfs <all> ...
cd ../analysis/figscripts && conda activate deeppbs && python make_R.py
```

**R2 per-residue profile:** produced by `compute_reachability.py` (not
`compute_rmsds.py`, which is whole-protein only). It now covers **all 13 pilots**.
To (re)build it for the current pilot set:

```bash
conda activate pycurves          # mdtraj
cd analysis/figscripts
python compute_reachability.py --tfs csl dux4 egr1 engrailed err ets1 foxa hsf irf lef1 nfat runx tbp \
    --source stage3 --out ../data
#    writes ../data/ca_rmsd_perresidue.csv (R2) AND ../data/reachability.csv (M1).
#    --source stage3 reproduces the original 6-pilot d_min values to 4 decimals
#    (ets1 0.8694, tbp 0.587); it prints a VALIDATION line so you can confirm.
```

### S series — PyMOL, compute node only
`cmd.png()` segfaults on the login node (no GL context) — always render via sbatch
on the `rohs` partition with `ray=1`. The edu license auto-loads from
`~/.pymol/license.lic` (valid through 2027-02-01).

```bash
sbatch analysis/figscripts/render_S.sbatch          # renders the default pilots
#    → analysis/figures/pymol/{S1_bestworst_<tf>,S2_stage_progression_<tf>}.png
```
`make_S.py` takes a pilot list, so the successor can render the rest:
`pymol -cq make_S.py -- <tf1> <tf2> …` (inside the sbatch wrapper). Verified today
for ets1 (S1+S2) and lef1 (S1).

### M series — reachability axis
`M1` plots two axes. `spread` (BioEmu ensemble diversity) and the own-family
`aug_dP` are rebuilt for every pilot by `extract_metrics.py`. The **reachability**
axis (`d_min`, `reach_ratio`) comes from `compute_reachability.py` →
`reachability.csv`, which `extract_metrics.py` merges into
`mechanism_apo_holo.csv` (reachability.csv WINS over any carried-forward value).
This now covers **all 13 pilots** (`reach_status = computed`). `make_M.py` caps
panel A's x-axis at 5 Å and draws a far-unreachable pilot (dux4 ≈ 10 Å) as an
off-scale annotation so the resolved 0.5–4.5 Å cluster stays legible; panel B
shows spread for all pilots with the augmentation-effect sign dot. A new pilot
enters M1 automatically once it appears in `reachability.csv` — just re-run
`compute_reachability.py` including it, then `extract_metrics.py`, then `make_M.py`.
If a pilot ever lacks a reachability row it is flagged `needs_coord_pass` and
shown spread-only in panel B.

## Adding a new pilot — checklist

1. Run the pipeline for the TF (stage1→7) so `output/stage3_min/<tf>/` and
   `output/stage7_eval/id_benchmark_<tf>.json` exist, and `config/pilots/<tf>.sh`
   is present.
2. Add one row to `PILOT_META` in `fig_common.py` (family, short label, curated
   `dna_deform`).
3. `python extract_metrics.py` → the new pilot is auto-discovered into every CSV.
4. Re-run `make_F/I/P/M.py` (instant). For R, re-run `compute_rmsds.py --tfs …`
   including the new TF, then `make_R.py`. For S, add it to the `render_S.sbatch`
   pilot list.
5. For M1's reachability axis and R2's per-residue profile, run
   `compute_reachability.py --tfs <all incl. new> --source stage3 --out ../data`
   (pycurves env), then re-run `extract_metrics.py` and `make_M.py` / `make_R.py`.
   The new pilot enters M1 panel A and R2 automatically once it's in
   `reachability.csv` / `ca_rmsd_perresidue.csv`.

## Current coverage (13 pilots discovered)

`csl dux4 egr1 engrailed err ets1 foxa hsf irf lef1 nfat runx tbp`
— F/I/P + R1/R2/R3 + M1 (both axes): **all 13** (dux4 has 0 fnat survivors, shown
as the 0% reference; irf lacks a `spread` value because the diversity pass wasn't
run for it, but its reachability is computed). Benchmark accuracy (P2/P3): 12
pilots (all but dux4, which has no eval).
